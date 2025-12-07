from pathlib import Path

import numpy as np
import pytest

import demo_quickstart as dq


def test_survival_dataset_structure():
    df = dq.build_survival_dataset(n_patients=10, seed=123)
    expected_cols = {
        "patient_id",
        "group",
        "time_months",
        "event",
        "tumor_stage",
        "lymph_positive",
        "nerve_invasion",
        "vascular_invasion",
        "margin_positive",
        "age_years",
        "albumin_g_dl",
    }
    assert set(df.columns) == expected_cols
    assert df["group"].value_counts().to_dict() == {"Standard": 5, "AI-guided": 5}
    assert df["event"].isin({0, 1}).all()


def test_survival_dataset_handles_odd_patient_counts():
    df = dq.build_survival_dataset(n_patients=11, seed=321)
    counts = df["group"].value_counts().to_dict()
    assert counts["Standard"] == 5
    assert counts["AI-guided"] == 6
    assert counts["Standard"] + counts["AI-guided"] == 11


def test_calibration_metrics_are_finite():
    df = dq.simulate_metastasis_dataset(n_samples=200, seed=99)
    X = df.drop(columns="metastasis")
    y = df["metastasis"]
    X_train, X_test, y_train, y_test = dq.train_test_split(
        X, y, test_size=0.5, stratify=y, random_state=99
    )

    summary, probs = dq.evaluate_calibration_models(X_train, X_test, y_train, y_test)
    assert {"Uncalibrated", "Platt scaling", "Isotonic"} <= set(probs.keys())
    assert np.isfinite(summary["Brier Score"]).all()
    assert np.isfinite(summary["ECE (10-bin)"]).all()

    raw_brier = summary.loc[
        summary["Model"] == "Gradient Boosting (uncalibrated)", "Brier Score"
    ].item()
    iso_brier = summary.loc[
        summary["Model"] == "Isotonic Regression", "Brier Score"
    ].item()
    assert iso_brier <= raw_brier + 1e-6


def test_survival_metrics_are_stable(tmp_path):
    df = dq.build_survival_dataset(n_patients=240, seed=123)
    results = dq.analyse_survival(df, tmp_path / "km.png")
    assert 0.4 < results.hazard_ratio < 0.8
    assert results.logrank_p_value < 0.01
    assert results.ph_global_p_value > 0.05


def test_tcga_cohort_preparation():
    data_path = Path("data/tcga_2018_clinical_data.tsv")
    if not data_path.exists():
        pytest.skip("TCGA pilot cohort not available.")
    clinical = dq.load_tcga_clinical_table(data_path)
    cohort = dq.prepare_tcga_survival_cohort(clinical)
    assert not cohort.empty
    assert {"Node-negative", "Node-positive"} <= set(cohort["group"])
    assert (cohort["time_months"] > 0).all()


def test_expected_calibration_error_perfect_calibration():
    """ECE should be zero for perfectly calibrated predictions."""
    y_true = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
    y_prob = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.9, 0.9, 0.9, 0.9, 0.9])
    ece = dq.expected_calibration_error(y_true, y_prob, n_bins=2, strategy="uniform")
    assert ece < 0.15  # Allow some binning noise


def test_expected_calibration_error_worst_case():
    """ECE should be high for completely miscalibrated predictions."""
    y_true = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0])
    y_prob = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.9, 0.9, 0.9, 0.9, 0.9])
    ece = dq.expected_calibration_error(y_true, y_prob, n_bins=2, strategy="uniform")
    assert ece > 0.7


def test_calibration_slope_intercept_returns_tuple():
    """Calibration slope/intercept should return two floats."""
    y_true = np.array([0, 0, 1, 1, 0, 1, 0, 1])
    y_prob = np.array([0.2, 0.3, 0.7, 0.8, 0.4, 0.6, 0.1, 0.9])
    slope, intercept = dq.calibration_slope_intercept(y_true, y_prob)
    assert isinstance(slope, float)
    assert isinstance(intercept, float)


def test_cli_args_parsing():
    """CLI parser should accept valid arguments."""
    config = dq.parse_cli_args(["--n-patients", "500", "--seed", "123"])
    assert config.n_patients == 500
    assert config.seed == 123


def test_cli_args_defaults():
    """CLI parser should use defaults when no arguments provided."""
    config = dq.parse_cli_args([])
    assert config.n_patients == 240
    assert config.seed == 42


def test_metastasis_dataset_class_balance():
    """Metastasis dataset should have reasonable class balance."""
    df = dq.simulate_metastasis_dataset(n_samples=1000, seed=42)
    metastasis_rate = df["metastasis"].mean()
    assert 0.15 < metastasis_rate < 0.40  # Expect ~25% metastasis rate


def test_synthetic_gastric_dataset_structure():
    """Synthetic gastric dataset should have expected columns."""
    df = dq.build_synthetic_gastric_node_dataset(n_patients=50, seed=42)
    expected_cols = {"patient_id", "group", "time_months", "event", "tumor_stage"}
    assert expected_cols.issubset(set(df.columns))


def test_survival_results_dataclass_fields():
    """SurvivalResults should contain all expected fields."""
    from dataclasses import fields

    field_names = {f.name for f in fields(dq.SurvivalResults)}
    expected = {"treatment_label", "control_label", "hazard_ratio", "c_index"}
    assert expected.issubset(field_names)
