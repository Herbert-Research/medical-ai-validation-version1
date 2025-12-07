#!/usr/bin/env python3
"""
Repository validation helper that checks for required artefacts and basic data consistency.

Usage:
    python3 validate_repository.py
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT / "outputs"


def find_file(name: str) -> Optional[Path]:
    """Return the first matching path for a file across common locations."""
    candidates = [
        ROOT / name,
        OUTPUTS / name,
        ROOT / "tests" / name,
        ROOT / "docs" / name,
    ]
    for path in candidates:
        if path.exists():
            return path
    return None


def check_required_files(required: Iterable[str]) -> Tuple[List[str], List[str]]:
    """Verify required files exist; return (found, missing) lists."""
    found: List[str] = []
    missing: List[str] = []
    print("Checking required files...")
    for name in required:
        path = find_file(name)
        if path:
            print(f"[OK] {name} -> {path}")
            found.append(name)
        else:
            print(f"[MISSING] {name}")
            missing.append(name)
    return found, missing


def as_number(value: object) -> float:
    """Convert values that may be strings, numpy scalars, or Python scalars to float."""
    if isinstance(value, (float, int, np.floating, np.integer)):
        return float(value)
    try:
        return float(str(value))
    except Exception:
        return float("nan")


def medians_consistent(group_value: object, survival_value: object) -> bool:
    """Check whether two median survival values agree within tolerance, handling inf/NaN."""
    g = as_number(group_value)
    s = as_number(survival_value)
    if (np.isnan(g) and np.isnan(s)) or (np.isinf(g) and np.isinf(s)):
        return True
    if np.isnan(g) != np.isnan(s):
        return False
    if np.isinf(g) != np.isinf(s):
        return False
    return abs(g - s) < 0.5


def check_survival_dataset(
    group_file: str,
    survival_file: str,
    treatment_label: str,
    control_label: str,
) -> List[str]:
    """Validate survival summaries against group-level medians."""
    errors: List[str] = []
    group_path = find_file(group_file)
    survival_path = find_file(survival_file)

    if not survival_path:
        return errors  # Nothing to check for this dataset

    print(f"\nChecking survival dataset: {survival_file}")
    try:
        survival_df = pd.read_csv(survival_path)
    except Exception as exc:
        errors.append(f"{survival_file}: load failed ({exc})")
        print(f"[ERROR] Unable to load {survival_file}: {exc}")
        return errors

    if group_path:
        try:
            group_df = pd.read_csv(group_path)
            control_group_vals = group_df[group_df["group"] == control_label][
                "median_months"
            ]
            if control_group_vals.empty:
                errors.append(f"{group_file}: no '{control_label}' group row found")
                print(f"[ERROR] '{control_label}' group missing in {group_file}")
            else:
                control_group = control_group_vals.values[0]
                control_survival = survival_df["Median DFS (control)"].values[0]
                if medians_consistent(control_group, control_survival):
                    print(
                        f"[OK] Median DFS matches for {control_label} "
                        f"({control_group} vs {control_survival})"
                    )
                else:
                    errors.append(
                        f"{survival_file}: median mismatch ({control_group} vs {control_survival})"
                    )
                    print(
                        f"[ERROR] Median mismatch for {control_label}: "
                        f"{control_group} vs {control_survival}"
                    )
        except Exception as exc:
            errors.append(f"{group_file}: error during consistency check ({exc})")
            print(f"[ERROR] Failed group consistency check for {group_file}: {exc}")
    else:
        print(
            f"[WARN] No group summary found for {survival_file}; skipping median cross-check."
        )

    treatment_median = survival_df.get("Median DFS (treatment)")
    control_median = survival_df.get("Median DFS (control)")
    if treatment_median is not None and control_median is not None:
        print(
            f"[INFO] Medians -> {treatment_label}: {treatment_median.iloc[0]}, "
            f"{control_label}: {control_median.iloc[0]}"
        )
    return errors


def run() -> int:
    errors: List[str] = []

    required = [
        "README.md",
        "DATA_SOURCES.md",
        "demo_quickstart.py",
        "requirements.txt",
        "calibration_summary.csv",
        "calibration_curve_example.png",
        "synthetic_survival_summary.csv",
        "kaplan_meier_example.png",
        "test_demo_quickstart.py",
    ]

    _, missing = check_required_files(required)
    if missing:
        errors.extend([f"Missing required file: {name}" for name in missing])

    print("\nLoading calibration summary...")
    try:
        calib_path = find_file("calibration_summary.csv")
        if not calib_path:
            raise FileNotFoundError(
                "calibration_summary.csv not found in expected locations."
            )
        calib = pd.read_csv(calib_path)
        print(
            f"[OK] Loaded calibration summary with {len(calib)} models from {calib_path}"
        )

        ece_values = calib["ECE (10-bin)"].to_numpy()
        if all(ece < 0.05 for ece in ece_values):
            print("[OK] All ECE values below 0.05")
        else:
            max_ece = float(np.max(ece_values))
            print(f"[WARN] Some ECE values exceed 0.05 (max observed: {max_ece:.3f})")

        platt_rows = calib[calib["Model"].str.contains("Platt", case=False)]
        if not platt_rows.empty:
            platt_slope = float(platt_rows["Calibration Slope"].values[0])
            if abs(platt_slope - 1.0) > 2.0:
                print(
                    f"[WARN] Platt slope={platt_slope:.2f} deviates strongly from 1.0"
                )
            else:
                print("[OK] Platt slope within expected bounds.")
    except Exception as exc:
        errors.append(f"Calibration summary check failed: {exc}")
        print(f"[ERROR] Calibration check failed: {exc}")

    print("\nChecking survival data consistency...")
    survival_errors = []
    any_survival = False

    tcga_errs = check_survival_dataset(
        group_file="tcga_group_summary.csv",
        survival_file="tcga_survival_summary.csv",
        treatment_label="Node-positive",
        control_label="Node-negative",
    )
    if find_file("tcga_survival_summary.csv"):
        any_survival = True
    survival_errors.extend(tcga_errs)

    synthetic_errs = check_survival_dataset(
        group_file="synthetic_gastric_cancer_group_summary.csv",
        survival_file="synthetic_gastric_cancer_survival_summary.csv",
        treatment_label="Node-positive",
        control_label="Node-negative",
    )
    if find_file("synthetic_gastric_cancer_survival_summary.csv"):
        any_survival = True
    survival_errors.extend(synthetic_errs)

    if not any_survival:
        survival_errors.append("No gastric cancer survival summaries found.")
        print(
            "[ERROR] No survival summaries found for TCGA or synthetic gastric cancer."
        )

    errors.extend(survival_errors)

    print("\n" + "=" * 80)
    if errors:
        print("Validation failed with the following issues:")
        for err in errors:
            print(f" - {err}")
        print("=" * 80)
        return 1
    else:
        print("All checks passed - repository appears consistent.")
        print("=" * 80)
        return 0


if __name__ == "__main__":
    sys.exit(run())
