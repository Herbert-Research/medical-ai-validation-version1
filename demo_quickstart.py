"""Headless entry point that reproduces the survival and calibration workflow.

Usage:
    python demo_quickstart.py
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Literal, Sequence, Tuple

import numpy as np
import pandas as pd

from config import (
    BASE_HAZARD_MONTHLY,
    CENSORING_RATE,
    COEF_AGE_PER_YEAR,
    COEF_AI_BENEFIT_BASE,
    COEF_ALBUMIN_PER_UNIT,
    COEF_LYMPH_POSITIVE,
    COEF_MARGIN_POSITIVE,
    COEF_NERVE_INVASION,
    COEF_TUMOR_STAGE,
    COEF_VASCULAR_INVASION,
    DEFAULT_METASTASIS_SAMPLES,
    DEFAULT_N_PATIENTS,
    DEFAULT_SEED,
    ECE_N_BINS,
    GB_LEARNING_RATE,
    GB_N_ESTIMATORS,
    LYMPH_POSITIVE_RATE,
    MARGIN_POSITIVE_RATE,
    NERVE_INVASION_RATE,
    VASCULAR_INVASION_RATE,
)

# Configure logging for pipeline execution
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

# Direct matplotlib cache to a writable workspace path to avoid warnings when $HOME is locked down.
MPL_CONFIG_DIR = Path.cwd() / "outputs" / ".mplconfig"
MPL_CONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(MPL_CONFIG_DIR)

import matplotlib
import seaborn as sns
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
from lifelines.utils import concordance_index

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from scipy.stats import chi2
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss
from sklearn.model_selection import train_test_split

CALIBRATION_BIN_STRATEGY: Literal["uniform", "quantile"] = "quantile"
STAGE_SUFFIX_OFFSETS = {"A": -1, "B": 0}
GRADE_MAPPING = {"G1": 1, "G2": 2, "G3": 3, "G4": 4}
DEFAULT_SURVIVAL_COVARIATES = [
    "tumor_stage",
    "lymph_positive",
    "nerve_invasion",
    "vascular_invasion",
    "margin_positive",
    "age_years",
    "albumin_g_dl",
]


@dataclass
class RunConfig:
    n_patients: int = DEFAULT_N_PATIENTS
    metastasis_samples: int = DEFAULT_METASTASIS_SAMPLES
    seed: int = DEFAULT_SEED
    output_dir: Path = field(default_factory=lambda: Path("outputs"))
    tcga_path: Path = field(
        default_factory=lambda: Path("data/tcga_2018_clinical_data.tsv")
    )
    calibration_bin_strategy: Literal["uniform", "quantile"] = CALIBRATION_BIN_STRATEGY
    skip_tcga: bool = False


@dataclass
class SurvivalResults:
    treatment_label: str
    control_label: str
    median_treatment: float
    median_control: float
    logrank_p_value: float
    hazard_ratio: float
    hazard_ratio_ci: Tuple[float, float]
    c_index: float
    ph_global_p_value: float


def survival_results_to_frame(
    results: SurvivalResults, cohort_label: str
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Cohort": [cohort_label],
            "Treatment label": [results.treatment_label],
            "Control label": [results.control_label],
            "Median DFS (treatment)": [results.median_treatment],
            "Median DFS (control)": [results.median_control],
            "Log-rank p-value": [results.logrank_p_value],
            "Hazard ratio": [results.hazard_ratio],
            "HR CI lower": [results.hazard_ratio_ci[0]],
            "HR CI upper": [results.hazard_ratio_ci[1]],
            "C-index": [results.c_index],
            "PH global p-value": [results.ph_global_p_value],
        }
    )


def write_dataframe(df: pd.DataFrame, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def parse_status_flag(value: str) -> float:
    if isinstance(value, str) and ":" in value:
        flag = value.split(":", 1)[0].strip()
        if flag in {"0", "1"}:
            return float(flag)
    return float(np.nan)


def stage_code_to_numeric(code: str) -> float:
    if not isinstance(code, str):
        return float(np.nan)
    cleaned = code.strip().upper()

    if cleaned == "TIS":
        return 0.0
    if cleaned == "T0":
        return 1.0
    if not cleaned.startswith("T"):
        return float(np.nan)

    numeric_part = "".join(ch for ch in cleaned[1:] if ch.isdigit())
    if not numeric_part:
        return float(np.nan)

    base_stage = int(numeric_part)
    base_value = 3 * (base_stage - 1) + 2
    suffix = "".join(ch for ch in cleaned[1:] if ch.isalpha())
    offset = STAGE_SUFFIX_OFFSETS.get(suffix[:1], 0)

    return float(base_value + offset)


def grade_to_numeric(grade: str) -> float:
    if not isinstance(grade, str):
        return float(np.nan)
    value = GRADE_MAPPING.get(grade.strip().upper())
    return float(value) if value is not None else float(np.nan)


def build_survival_dataset(
    n_patients: int = DEFAULT_N_PATIENTS, seed: int = DEFAULT_SEED
) -> pd.DataFrame:
    """Generate synthetic survival cohort with clinical covariates.

    Parameters
    ----------
    n_patients : int
        Total number of patients. Must be >= 2 to have at least one per group.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        Synthetic survival dataset with patient characteristics and outcomes.

    Raises
    ------
    ValueError
        If n_patients < 2.
    TypeError
        If seed is not an integer.
    """
    if n_patients < 2:
        raise ValueError(
            f"n_patients must be >= 2 to have at least one patient per group, got {n_patients}"
        )
    if not isinstance(seed, int):
        raise TypeError(f"seed must be an integer, got {type(seed).__name__}")

    rng = np.random.default_rng(seed)
    n_control = n_patients // 2
    n_treatment = n_patients - n_control
    group = np.array(["Standard"] * n_control + ["AI-guided"] * n_treatment)
    rng.shuffle(group)
    ai_indicator = (group == "AI-guided").astype(int)

    tumor_stage = rng.integers(1, 5, size=n_patients)  # I-IV
    lymph_positive = rng.binomial(1, LYMPH_POSITIVE_RATE, size=n_patients)
    nerve_invasion = rng.binomial(1, NERVE_INVASION_RATE, size=n_patients)
    vascular_invasion = rng.binomial(1, VASCULAR_INVASION_RATE, size=n_patients)
    margin_positive = rng.binomial(1, MARGIN_POSITIVE_RATE, size=n_patients)
    age_years = rng.normal(62, 9, size=n_patients)
    albumin_g_dl = rng.normal(3.9, 0.35, size=n_patients)

    ai_benefit = (
        COEF_AI_BENEFIT_BASE * ai_indicator * (0.65 + 0.35 * (tumor_stage >= 3))
    )

    linear_predictor = (
        COEF_TUMOR_STAGE * (tumor_stage - 1)
        + COEF_LYMPH_POSITIVE * lymph_positive
        + COEF_NERVE_INVASION * nerve_invasion
        + COEF_VASCULAR_INVASION * vascular_invasion
        + COEF_MARGIN_POSITIVE * margin_positive
        + COEF_AGE_PER_YEAR * (age_years - 60)
        + COEF_ALBUMIN_PER_UNIT * (albumin_g_dl - 4.0)
        + ai_benefit
        + rng.normal(0.0, 0.48, size=n_patients)
    )

    base_hazard = BASE_HAZARD_MONTHLY
    survival_time_months = -np.log(rng.uniform(size=n_patients)) / (
        base_hazard * np.exp(linear_predictor)
    )

    censoring_time = -np.log(rng.uniform(size=n_patients)) / CENSORING_RATE
    observed_time = np.minimum(survival_time_months, censoring_time)
    event_observed = (survival_time_months <= censoring_time).astype(int)

    return pd.DataFrame(
        {
            "patient_id": range(1, n_patients + 1),
            "group": group,
            "time_months": observed_time,
            "event": event_observed,
            "tumor_stage": tumor_stage,
            "lymph_positive": lymph_positive,
            "nerve_invasion": nerve_invasion,
            "vascular_invasion": vascular_invasion,
            "margin_positive": margin_positive,
            "age_years": age_years,
            "albumin_g_dl": albumin_g_dl,
        }
    )


def build_synthetic_gastric_node_dataset(
    n_patients: int = 220, seed: int = 7
) -> pd.DataFrame:
    """Simulated gastric cancer cohort focused on node status discrimination.

    Parameters
    ----------
    n_patients : int
        Total number of patients. Must be >= 2.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        Synthetic gastric cancer survival dataset.

    Raises
    ------
    ValueError
        If n_patients < 2.
    TypeError
        If seed is not an integer.
    """
    if n_patients < 2:
        raise ValueError(f"n_patients must be >= 2, got {n_patients}")
    if not isinstance(seed, int):
        raise TypeError(f"seed must be an integer, got {type(seed).__name__}")

    rng = np.random.default_rng(seed)
    node_positive = rng.binomial(1, 0.45, size=n_patients)
    group = np.where(node_positive == 1, "Node-positive", "Node-negative")

    tumor_stage = rng.integers(1, 5, size=n_patients)  # I-IV
    grade = rng.integers(1, 4, size=n_patients)  # 1-3
    lymph_ratio = np.clip(
        rng.normal(0.22 + 0.25 * node_positive, 0.12, size=n_patients), 0, 1
    )
    perineural_invasion = rng.binomial(1, 0.35, size=n_patients)
    vascular_invasion = rng.binomial(1, 0.28, size=n_patients)
    age_years = rng.normal(63, 8, size=n_patients)
    albumin_g_dl = rng.normal(3.8, 0.35, size=n_patients)

    node_effect = 0.38 * node_positive * (0.55 + 0.45 * (tumor_stage >= 3))

    linear_predictor = (
        0.55 * (tumor_stage - 1)
        + 0.30 * (grade - 1)
        + 1.0 * lymph_ratio
        + 0.35 * perineural_invasion
        + 0.32 * vascular_invasion
        + 0.012 * (age_years - 60)
        - 0.35 * (albumin_g_dl - 4.0)
        + node_effect
        + rng.normal(0.0, 0.60, size=n_patients)
    )

    base_hazard = 0.02
    survival_time_months = -np.log(rng.uniform(size=n_patients)) / (
        base_hazard * np.exp(linear_predictor)
    )
    censoring_time = -np.log(rng.uniform(size=n_patients)) / 0.014
    observed_time = np.minimum(survival_time_months, censoring_time)
    event_observed = (survival_time_months <= censoring_time).astype(int)

    return pd.DataFrame(
        {
            "patient_id": range(1, n_patients + 1),
            "group": group,
            "time_months": observed_time,
            "event": event_observed,
            "tumor_stage": tumor_stage,
            "grade": grade,
            "lymph_ratio": lymph_ratio,
            "perineural_invasion": perineural_invasion,
            "vascular_invasion": vascular_invasion,
            "age_years": age_years,
            "albumin_g_dl": albumin_g_dl,
        }
    )


def _plot_kaplan_meier(
    data: pd.DataFrame,
    treatment_label: str,
    control_label: str,
    output_path: Path,
    plot_title: str,
) -> Dict[str, float]:
    """Generate Kaplan-Meier survival curves with confidence intervals."""
    sns.set_context("talk")
    sns.set_style("whitegrid")
    plt.figure(figsize=(10, 6))

    kmf = KaplanMeierFitter()
    medians: Dict[str, float] = {}
    default_colors = {"Standard": "#2A9D8F", "AI-guided": "#E76F51"}
    palette_iter = iter(sns.color_palette("Set2", max(2, len(set(data["group"])))))
    color_lookup = {}
    for group in [control_label, treatment_label]:
        if group in default_colors:
            color_lookup[group] = default_colors[group]
        else:
            color_lookup[group] = next(palette_iter)

    for group in [control_label, treatment_label]:
        subset = data[data["group"] == group]
        kmf.fit(subset["time_months"], event_observed=subset["event"], label=group)
        kmf.plot_survival_function(
            ci_show=True, ci_alpha=0.15, lw=2, color=color_lookup[group]
        )
        medians[group] = kmf.median_survival_time_

    plt.title(plot_title, fontsize=16, fontweight="bold")
    plt.xlabel("Time since gastrectomy (months)")
    plt.ylabel("Disease-free survival probability")
    plt.ylim(0, 1.0)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    return medians


def _fit_cox_model(
    data: pd.DataFrame,
    treatment_label: str,
    covariate_cols: Iterable[str],
) -> Tuple[CoxPHFitter, float, Tuple[float, float], float, float]:
    """Fit Cox proportional hazards model and extract key metrics."""
    cox_df = data[["time_months", "event"] + list(covariate_cols)].copy()
    cox_df["treatment"] = (data["group"] == treatment_label).astype(int)

    cph = CoxPHFitter()
    cph.fit(cox_df, duration_col="time_months", event_col="event")

    hazard_ratio = float(np.exp(cph.params_["treatment"]))
    ci_lower, ci_upper = np.exp(cph.confidence_intervals_.loc["treatment"])

    risk_scores = cph.predict_partial_hazard(cox_df).values.ravel()
    c_index = float(
        concordance_index(
            event_times=cox_df["time_months"],
            # negate partial hazards so higher scores = longer survival
            predicted_scores=-risk_scores,
            event_observed=cox_df["event"],
        )
    )
    ph_test = proportional_hazard_test(cph, cox_df, time_transform="rank")
    global_stat = float(ph_test.test_statistic.sum())
    df_ph = len(ph_test.test_statistic)
    ph_p_value = float(chi2.sf(global_stat, df_ph))

    return cph, hazard_ratio, (float(ci_lower), float(ci_upper)), c_index, ph_p_value


def analyse_survival(
    data: pd.DataFrame,
    output_path: Path,
    *,
    treatment_label: str = "AI-guided",
    control_label: str = "Standard",
    covariate_cols: Iterable[str] | None = None,
    plot_title: str = "Kaplan-Meier Survival Curves with 95% CI",
) -> SurvivalResults:
    covariate_cols = list(covariate_cols or DEFAULT_SURVIVAL_COVARIATES)
    if not {treatment_label, control_label}.issubset(set(data["group"])):
        raise ValueError("Both treatment and control groups must be present.")

    medians = _plot_kaplan_meier(
        data, treatment_label, control_label, output_path, plot_title
    )
    control = data[data["group"] == control_label]
    treatment = data[data["group"] == treatment_label]
    logrank = logrank_test(
        control["time_months"],
        treatment["time_months"],
        control["event"],
        treatment["event"],
    )
    _, hazard_ratio, hazard_ratio_ci, c_index, ph_p_value = _fit_cox_model(
        data, treatment_label, covariate_cols
    )

    return SurvivalResults(
        treatment_label=treatment_label,
        control_label=control_label,
        median_treatment=medians[treatment_label],
        median_control=medians[control_label],
        logrank_p_value=float(logrank.p_value),
        hazard_ratio=hazard_ratio,
        hazard_ratio_ci=hazard_ratio_ci,
        c_index=c_index,
        ph_global_p_value=ph_p_value,
    )


def simulate_metastasis_dataset(
    n_samples: int = DEFAULT_METASTASIS_SAMPLES, seed: int = DEFAULT_SEED
) -> pd.DataFrame:
    """Generate synthetic metastasis prediction dataset.

    Parameters
    ----------
    n_samples : int
        Number of samples to generate. Must be >= 10 for meaningful splits.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        Synthetic dataset with features and metastasis outcome.

    Raises
    ------
    ValueError
        If n_samples < 10.
    TypeError
        If seed is not an integer.
    """
    if n_samples < 10:
        raise ValueError(
            f"n_samples must be >= 10 for meaningful train/test splits, got {n_samples}"
        )
    if not isinstance(seed, int):
        raise TypeError(f"seed must be an integer, got {type(seed).__name__}")

    rng = np.random.default_rng(seed)
    tumor_stage = rng.integers(1, 4, size=n_samples)
    lymph_nodes = rng.poisson(lam=6, size=n_samples)
    station_distance = rng.normal(loc=3.0, scale=1.0, size=n_samples)
    histology_score = rng.normal(loc=0.0, scale=1.0, size=n_samples)
    bmi = rng.normal(loc=24, scale=4, size=n_samples)

    linear_pred = (
        0.6 * tumor_stage
        + 0.12 * lymph_nodes
        - 0.35 * station_distance
        - 0.25 * histology_score
        + 0.02 * (bmi - 25)
    )

    base_risk = 1 / (1 + np.exp(-(linear_pred - 2.5)))
    metastasis = rng.binomial(1, base_risk)

    return pd.DataFrame(
        {
            "tumor_stage": tumor_stage,
            "lymph_nodes": lymph_nodes,
            "station_distance": station_distance,
            "histology_score": histology_score,
            "bmi": bmi,
            "metastasis": metastasis,
        }
    )


def expected_calibration_error(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    n_bins: int = ECE_N_BINS,
    *,
    strategy: Literal["uniform", "quantile"] = CALIBRATION_BIN_STRATEGY,
) -> float:
    if strategy not in {"uniform", "quantile"}:
        raise ValueError("strategy must be 'uniform' or 'quantile'")

    y_true = np.asarray(y_true)
    y_prob = np.asarray(y_prob)

    if strategy == "quantile":
        quantiles = np.linspace(0.0, 1.0, n_bins + 1)
        bins = np.quantile(y_prob, quantiles)
        bins = np.unique(bins)
        if len(bins) <= 1:
            bins = np.linspace(0.0, 1.0, n_bins + 1)
    else:
        bins = np.linspace(0.0, 1.0, n_bins + 1)

    bin_ids = np.digitize(y_prob, bins, right=False) - 1
    bin_ids = np.clip(bin_ids, 0, len(bins) - 2)

    ece = 0.0
    for idx in range(len(bins) - 1):
        mask = bin_ids == idx
        if np.any(mask):
            observed = y_true[mask].mean()
            predicted = y_prob[mask].mean()
            ece += abs(observed - predicted) * mask.mean()
    return float(ece)


def calibration_slope_intercept(
    y_true: np.ndarray, y_prob: np.ndarray
) -> Tuple[float, float]:
    eps = 1e-6
    clipped = np.clip(y_prob, eps, 1 - eps)
    logits = np.log(clipped / (1 - clipped)).reshape(-1, 1)
    lr = LogisticRegression(penalty=None, solver="lbfgs", max_iter=1000)
    lr.fit(logits, y_true)
    return float(lr.coef_[0][0]), float(lr.intercept_[0])


def evaluate_calibration_models(
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
    y_train: pd.Series,
    y_test: pd.Series,
    *,
    strategy: Literal["uniform", "quantile"] = CALIBRATION_BIN_STRATEGY,
    n_bins: int = ECE_N_BINS,
) -> Tuple[pd.DataFrame, Dict[str, np.ndarray]]:
    def make_estimator() -> GradientBoostingClassifier:
        return GradientBoostingClassifier(
            random_state=DEFAULT_SEED,
            n_estimators=GB_N_ESTIMATORS,
            learning_rate=GB_LEARNING_RATE,
        )

    base_model = make_estimator()
    base_model.fit(X_train, y_train)
    raw_probs = base_model.predict_proba(X_test)[:, 1]

    sigmoid_model = CalibratedClassifierCV(
        estimator=make_estimator(), method="sigmoid", cv=5
    )
    sigmoid_model.fit(X_train, y_train)
    sigmoid_probs = sigmoid_model.predict_proba(X_test)[:, 1]

    isotonic_model = CalibratedClassifierCV(
        estimator=make_estimator(), method="isotonic", cv=5
    )
    isotonic_model.fit(X_train, y_train)
    isotonic_probs = isotonic_model.predict_proba(X_test)[:, 1]

    labels = y_test.to_numpy()

    def summarise(label: str, probs: np.ndarray) -> Dict[str, float | str]:
        slope, intercept = calibration_slope_intercept(labels, probs)
        return {
            "Model": label,
            "Brier Score": float(brier_score_loss(labels, probs)),
            f"ECE ({n_bins}-bin)": expected_calibration_error(
                labels, probs, n_bins=n_bins, strategy=strategy
            ),
            "Calibration Slope": slope,
            "Calibration Intercept": intercept,
        }

    summary = pd.DataFrame(
        [
            summarise("Gradient Boosting (uncalibrated)", raw_probs),
            summarise("Platt Scaling (sigmoid)", sigmoid_probs),
            summarise("Isotonic Regression", isotonic_probs),
        ]
    )

    return summary, {
        "Uncalibrated": raw_probs,
        "Platt scaling": sigmoid_probs,
        "Isotonic": isotonic_probs,
    }


def plot_calibration_curves(
    probabilities: Dict[str, np.ndarray],
    y_test: pd.Series,
    output_path: Path,
    *,
    strategy: Literal["uniform", "quantile"] = CALIBRATION_BIN_STRATEGY,
) -> None:
    plt.figure(figsize=(10, 6))
    plt.plot([0, 1], [0, 1], "k--", label="Perfect calibration")

    for label, probs in probabilities.items():
        frac_pos, mean_pred = calibration_curve(
            y_test, probs, n_bins=ECE_N_BINS, strategy=strategy
        )
        plt.plot(
            mean_pred, frac_pos, marker="o", linewidth=2, markersize=7, label=label
        )

    plt.xlabel("Predicted probability")
    plt.ylabel("Observed frequency")
    plt.title("Calibration Curves for Station-Level Metastasis Risk")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def load_tcga_clinical_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", na_values=["NA", "na", "NaN", "null"])


def prepare_tcga_survival_cohort(clinical_df: pd.DataFrame) -> pd.DataFrame:
    df = clinical_df.copy()
    required_cols = [
        "Patient ID",
        "Progress Free Survival (Months)",
        "Disease Free (Months)",
        "Progression Free Status",
        "Disease Free Status",
        "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code",
        "American Joint Committee on Cancer Tumor Stage Code",
        "Neoplasm Histologic Grade",
        "Diagnosis Age",
    ]
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise KeyError(f"TCGA table is missing required columns: {missing}")
    pfs = pd.to_numeric(df["Progress Free Survival (Months)"], errors="coerce")
    dfs = pd.to_numeric(df["Disease Free (Months)"], errors="coerce")
    time_months = pfs.fillna(dfs)

    pfs_status = df["Progression Free Status"].apply(parse_status_flag)
    dfs_status = df["Disease Free Status"].apply(parse_status_flag)
    event_indicator = pfs_status.fillna(dfs_status)

    n_stage = (
        df["Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code"]
        .str.upper()
        .str.strip()
    )
    valid_n_stage = n_stage.isin({"N0", "N1", "N2", "N3", "N3A", "N3B"})
    group = np.where(n_stage == "N0", "Node-negative", "Node-positive")

    tumor_stage = df["American Joint Committee on Cancer Tumor Stage Code"].apply(
        stage_code_to_numeric
    )
    grade = df["Neoplasm Histologic Grade"].apply(grade_to_numeric)
    diagnosis_age = pd.to_numeric(df["Diagnosis Age"], errors="coerce")

    cohort = pd.DataFrame(
        {
            "patient_id": df["Patient ID"],
            "group": group,
            "time_months": time_months,
            "event": event_indicator,
            "tumor_stage_num": tumor_stage,
            "grade_num": grade,
            "diagnosis_age": diagnosis_age,
        }
    )

    cohort = cohort[valid_n_stage].dropna(subset=["time_months", "event"])
    cohort = cohort[cohort["time_months"] > 0]
    cohort["event"] = cohort["event"].astype(int)
    cohort = cohort.dropna(subset=["tumor_stage_num", "grade_num", "diagnosis_age"])

    return cohort


def parse_cli_args(argv: Sequence[str] | None = None) -> RunConfig:
    """Expose a lightweight CLI so the workflow can be parameterised."""

    parser = argparse.ArgumentParser(
        description="Quantitative survival validation and calibration quickstart.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--n-patients",
        type=int,
        default=RunConfig().n_patients,
        help="Number of synthetic patients to include in the survival simulation.",
    )
    parser.add_argument(
        "--metastasis-samples",
        type=int,
        default=RunConfig().metastasis_samples,
        help="Number of samples to draw for the metastasis calibration dataset.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=RunConfig().seed,
        help="Random seed applied to both synthetic data generators.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=RunConfig().output_dir,
        help="Directory where figures and CSV artefacts will be written.",
    )
    parser.add_argument(
        "--tcga-path",
        type=Path,
        default=RunConfig().tcga_path,
        help="Path to the TCGA clinical TSV used for the pilot external validation.",
    )
    parser.add_argument(
        "--calibration-strategy",
        choices=["uniform", "quantile"],
        default=RunConfig().calibration_bin_strategy,
        help="Binning strategy for calibration curves and ECE computation.",
    )
    parser.add_argument(
        "--skip-tcga",
        action="store_true",
        help="Disable the TCGA pilot analysis (useful when the TSV is not available).",
    )
    args = parser.parse_args(argv)
    return RunConfig(
        n_patients=args.n_patients,
        metastasis_samples=args.metastasis_samples,
        seed=args.seed,
        output_dir=args.output_dir,
        tcga_path=args.tcga_path,
        calibration_bin_strategy=args.calibration_strategy,
        skip_tcga=args.skip_tcga,
    )


def run(config: RunConfig | None = None) -> int:
    """Execute the validation pipeline.

    Parameters
    ----------
    config : RunConfig | None
        Pipeline configuration. If None, uses default settings.

    Returns
    -------
    int
        Exit code: 0 for success, 1 for failure.
    """
    cfg = config or RunConfig()

    # --- Output directory setup ---
    try:
        output_dir = Path(cfg.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Output directory: {output_dir.resolve()}")
    except PermissionError as e:
        logger.error(f"Cannot create output directory '{cfg.output_dir}': {e}")
        return 1
    except OSError as e:
        logger.error(f"Failed to create output directory '{cfg.output_dir}': {e}")
        return 1

    output_km = output_dir / "kaplan_meier_example.png"
    output_calibration = output_dir / "calibration_curve_example.png"
    output_gastric = output_dir / "synthetic_gastric_kaplan_meier.png"
    output_tcga = output_dir / "tcga_kaplan_meier.png"

    # --- Synthetic survival cohort generation ---
    try:
        logger.info(
            f"Generating synthetic survival cohort (n={cfg.n_patients}, seed={cfg.seed})..."
        )
        survival_df = build_survival_dataset(n_patients=cfg.n_patients, seed=cfg.seed)
        logger.info(
            f"Generated {len(survival_df)} patients with {survival_df['event'].sum()} events"
        )
    except (ValueError, TypeError) as e:
        logger.error(f"Invalid parameters for survival dataset: {e}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error generating survival dataset: {e}")
        return 1

    # --- Survival analysis ---
    try:
        logger.info("Running survival analysis...")
        survival_results = analyse_survival(survival_df, output_km)
        logger.info(
            f"Survival analysis complete. HR={survival_results.hazard_ratio:.2f}, p={survival_results.logrank_p_value:.4f}"
        )
    except Exception as e:
        logger.error(f"Survival analysis failed: {e}")
        return 1
    synthetic_summary_path = output_dir / "synthetic_survival_summary.csv"
    try:
        write_dataframe(
            survival_results_to_frame(
                survival_results, "Synthetic (Standard vs AI-guided)"
            ),
            synthetic_summary_path,
        )
    except (IOError, OSError) as e:
        logger.error(f"Failed to write synthetic survival summary: {e}")
        return 1

    # --- Synthetic gastric cancer cohort ---
    try:
        logger.info("Generating synthetic gastric cancer cohort...")
        gastric_df = build_synthetic_gastric_node_dataset(seed=cfg.seed + 101)
        logger.info(f"Generated {len(gastric_df)} gastric cancer patients")
    except (ValueError, TypeError) as e:
        logger.error(f"Invalid parameters for gastric dataset: {e}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error generating gastric dataset: {e}")
        return 1
    try:
        logger.info("Running gastric cancer survival analysis...")
        gastric_results = analyse_survival(
            gastric_df,
            output_gastric,
            treatment_label="Node-positive",
            control_label="Node-negative",
            covariate_cols=[
                "tumor_stage",
                "grade",
                "lymph_ratio",
                "perineural_invasion",
                "vascular_invasion",
                "age_years",
                "albumin_g_dl",
            ],
            plot_title="Synthetic Gastric Survival by Node Status",
        )
        logger.info(f"Gastric analysis complete. HR={gastric_results.hazard_ratio:.2f}")
    except Exception as e:
        logger.error(f"Gastric cancer survival analysis failed: {e}")
        return 1

    synthetic_gastric_summary_path = (
        output_dir / "synthetic_gastric_cancer_survival_summary.csv"
    )
    try:
        write_dataframe(
            survival_results_to_frame(
                gastric_results, "Synthetic gastric cancer (Node+ vs Node-)"
            ),
            synthetic_gastric_summary_path,
        )
    except (IOError, OSError) as e:
        logger.error(f"Failed to write gastric survival summary: {e}")
        return 1

    # --- Metastasis calibration analysis ---
    try:
        logger.info(f"Generating metastasis dataset (n={cfg.metastasis_samples})...")
        metastasis_df = simulate_metastasis_dataset(
            n_samples=cfg.metastasis_samples, seed=cfg.seed
        )
    except (ValueError, TypeError) as e:
        logger.error(f"Invalid parameters for metastasis dataset: {e}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error generating metastasis dataset: {e}")
        return 1

    try:
        X = metastasis_df.drop(columns="metastasis")
        y = metastasis_df["metastasis"]
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.25, stratify=y, random_state=DEFAULT_SEED
        )
        logger.info(
            f"Train/test split: {len(X_train)} train, {len(X_test)} test samples"
        )
    except Exception as e:
        logger.error(f"Failed to prepare train/test split: {e}")
        return 1

    try:
        logger.info("Evaluating calibration models...")
        calibration_summary, probs = evaluate_calibration_models(
            X_train,
            X_test,
            y_train,
            y_test,
            strategy=cfg.calibration_bin_strategy,
        )
        logger.info(f"Calibration evaluation complete ({len(probs)} models compared)")
    except Exception as e:
        logger.error(f"Calibration model evaluation failed: {e}")
        return 1

    try:
        plot_calibration_curves(
            probs, y_test, output_calibration, strategy=cfg.calibration_bin_strategy
        )
        logger.info(f"Calibration curves saved to {output_calibration}")
    except Exception as e:
        logger.error(f"Failed to generate calibration plots: {e}")
        return 1

    # --- TCGA pilot analysis ---
    tcga_path = Path(cfg.tcga_path)
    tcga_results: SurvivalResults | None = None
    tcga_group_summary: pd.DataFrame | None = None
    if cfg.skip_tcga:
        logger.info("Skipping TCGA pilot analysis (--skip-tcga flag enabled)")
    else:
        if tcga_path.exists():
            try:
                logger.info(f"Loading TCGA clinical data from {tcga_path}...")
                tcga_raw = load_tcga_clinical_table(tcga_path)
                tcga_cohort = prepare_tcga_survival_cohort(tcga_raw)
            except Exception as e:
                logger.error(f"Failed to load/prepare TCGA data: {e}")
                return 1

            if not tcga_cohort.empty:
                try:
                    logger.info(
                        f"Running TCGA survival analysis (n={len(tcga_cohort)})..."
                    )
                    tcga_results = analyse_survival(
                        tcga_cohort,
                        output_tcga,
                        treatment_label="Node-positive",
                        control_label="Node-negative",
                        covariate_cols=[
                            "tumor_stage_num",
                            "grade_num",
                            "diagnosis_age",
                        ],
                        plot_title="TCGA-STAD Progression-Free Survival by N-stage",
                    )
                    logger.info(
                        f"TCGA analysis complete. HR={tcga_results.hazard_ratio:.2f}"
                    )
                except Exception as e:
                    logger.error(f"TCGA survival analysis failed: {e}")
                    return 1

                try:
                    group_rows = []
                    for group_label, subset in tcga_cohort.groupby("group"):
                        kmf = KaplanMeierFitter()
                        kmf.fit(subset["time_months"], event_observed=subset["event"])
                        group_rows.append(
                            {
                                "group": group_label,
                                "patients": len(subset),
                                "median_months": kmf.median_survival_time_,
                                "event_rate": subset["event"].mean(),
                            }
                        )
                    tcga_group_summary = pd.DataFrame(group_rows).round(
                        {"median_months": 1, "event_rate": 3}
                    )
                except Exception as e:
                    logger.warning(f"Failed to compute TCGA group summary: {e}")
                    tcga_group_summary = None
            else:
                logger.warning(
                    "TCGA cohort preparation yielded 0 patients after filtering; "
                    "skipping pilot analysis."
                )
        else:
            logger.warning(
                f"TCGA clinical file not found at {tcga_path}. "
                "Place it under data/ or point --tcga-path to a valid TSV."
            )

    # --- Final summary output ---
    logger.info("Pipeline execution successful. Generating summary reports...")

    print("\n" + "=" * 60)
    print("=== Survival Analysis Summary ===")
    print("=" * 60)
    print(
        f"Median DFS ({survival_results.treatment_label}): "
        f"{survival_results.median_treatment:.1f} months"
    )
    print(
        f"Median DFS ({survival_results.control_label}):  "
        f"{survival_results.median_control:.1f} months"
    )
    print(f"Log-rank p-value:      {survival_results.logrank_p_value:.4f}")
    print(
        f"Hazard ratio ({survival_results.treatment_label} vs {survival_results.control_label}): "
        f"{survival_results.hazard_ratio:.2f} "
        f"(95% CI {survival_results.hazard_ratio_ci[0]:.2f}–{survival_results.hazard_ratio_ci[1]:.2f})"
    )
    print(f"C-index:                {survival_results.c_index:.3f}")
    print(f"PH global p-value:      {survival_results.ph_global_p_value:.3f}")

    calibration_summary_path = output_dir / "calibration_summary.csv"
    try:
        write_dataframe(calibration_summary, calibration_summary_path)
    except (IOError, OSError) as e:
        logger.error(f"Failed to write calibration summary: {e}")
        return 1

    print("\n=== Calibration Summary ===")
    print(calibration_summary.round(3).to_string(index=False))
    print(
        f"\nFigures saved to: {output_km.resolve()}, {output_calibration.resolve()}. "
        f"Tabular outputs written to: "
        f"{synthetic_summary_path.resolve()}, "
        f"{synthetic_gastric_summary_path.resolve()}, "
        f"{calibration_summary_path.resolve()}."
    )

    print("\n=== Synthetic Gastric Cancer (Node-positive vs Node-negative) ===")
    print(
        f"Median DFS (Node-positive): {gastric_results.median_treatment:.1f} months | "
        f"Median DFS (Node-negative): {gastric_results.median_control:.1f} months"
    )
    print(
        f"Log-rank p-value: {gastric_results.logrank_p_value:.4f} | "
        f"PH global p-value: {gastric_results.ph_global_p_value:.3f}"
    )
    print(
        f"Hazard ratio: {gastric_results.hazard_ratio:.2f} "
        f"(95% CI {gastric_results.hazard_ratio_ci[0]:.2f}–{gastric_results.hazard_ratio_ci[1]:.2f})"
    )
    print(f"C-index: {gastric_results.c_index:.3f}")
    print(f"Synthetic gastric figure saved to: {output_gastric.resolve()}.")
    print(
        f"Synthetic gastric summary saved to: "
        f"{synthetic_gastric_summary_path.resolve()}."
    )

    if tcga_results is not None and tcga_group_summary is not None:
        tcga_summary_path = output_dir / "tcga_survival_summary.csv"
        try:
            write_dataframe(
                survival_results_to_frame(
                    tcga_results, "TCGA-STAD (Node-positive vs Node-negative)"
                ),
                tcga_summary_path,
            )
            tcga_group_summary_path = output_dir / "tcga_group_summary.csv"
            write_dataframe(tcga_group_summary, tcga_group_summary_path)
        except (IOError, OSError) as e:
            logger.error(f"Failed to write TCGA summaries: {e}")
            return 1

        print("\n=== TCGA-STAD Pilot (Node-positive vs Node-negative) ===")
        print(
            f"Median PFS (Node-positive): {tcga_results.median_treatment:.1f} months | "
            f"Median PFS (Node-negative): {tcga_results.median_control:.1f} months"
        )
        print(
            f"Log-rank p-value: {tcga_results.logrank_p_value:.4f} | "
            f"PH global p-value: {tcga_results.ph_global_p_value:.3f}"
        )
        print(
            f"Hazard ratio: {tcga_results.hazard_ratio:.2f} "
            f"(95% CI {tcga_results.hazard_ratio_ci[0]:.2f}–{tcga_results.hazard_ratio_ci[1]:.2f})"
        )
        print("\nGroup-level summary (TCGA):")
        print(tcga_group_summary.to_string(index=False))
        print(
            f"\nTCGA figure saved to: {output_tcga.resolve()}. "
            f"Summaries saved to: {tcga_summary_path.resolve()}, {tcga_group_summary_path.resolve()}."
        )

    logger.info("Pipeline completed successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(run(parse_cli_args()))
