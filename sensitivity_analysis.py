"""
Sensitivity analysis: assess result stability across random seeds.

Demonstrates that findings are not artifacts of seed selection by running
the survival analysis pipeline multiple times with different random seeds
and summarizing the variability in key metrics.

Usage:
    python sensitivity_analysis.py
    python sensitivity_analysis.py --seeds 42 123 456 789 1011
    python sensitivity_analysis.py --n-patients 500 --output outputs/sensitivity_large.csv
"""

from __future__ import annotations

import argparse
import logging
import sys
import tempfile
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

import demo_quickstart as dq

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def run_sensitivity_analysis(
    seeds: List[int] | None = None,
    n_patients: int = 240,
    output_path: Path | None = None,
) -> pd.DataFrame:
    """Run survival analysis across multiple seeds and summarize variability.

    Parameters
    ----------
    seeds : List[int] | None
        List of random seeds to test. Defaults to [42, 123, 456, 789, 1011].
    n_patients : int
        Number of patients per cohort.
    output_path : Path | None
        Path to save results CSV. Defaults to outputs/sensitivity_analysis.csv.

    Returns
    -------
    pd.DataFrame
        DataFrame with results for each seed and summary statistics.
    """
    if seeds is None:
        seeds = [42, 123, 456, 789, 1011]
    
    if output_path is None:
        output_path = Path("outputs/sensitivity_analysis.csv")
    
    logger.info(f"Running sensitivity analysis with {len(seeds)} seeds")
    logger.info(f"Seeds: {seeds}")
    logger.info(f"Patients per cohort: {n_patients}")
    
    results = []
    
    for i, seed in enumerate(seeds, 1):
        logger.info(f"Processing seed {seed} ({i}/{len(seeds)})")
        
        try:
            # Generate synthetic cohort with this seed
            df = dq.build_survival_dataset(n_patients=n_patients, seed=seed)
            
            # Use a temporary file for the figure we won't keep
            with tempfile.TemporaryDirectory() as tmpdir:
                temp_fig_path = Path(tmpdir) / f"km_seed{seed}.png"
                survival_results = dq.analyse_survival(df, temp_fig_path)
            
            results.append({
                "seed": seed,
                "hazard_ratio": survival_results.hazard_ratio,
                "hr_ci_lower": survival_results.hazard_ratio_ci[0],
                "hr_ci_upper": survival_results.hazard_ratio_ci[1],
                "logrank_p": survival_results.logrank_p_value,
                "c_index": survival_results.c_index,
                "median_treatment": survival_results.median_treatment,
                "median_control": survival_results.median_control,
                "ph_global_p": survival_results.ph_global_p_value,
            })
            
        except Exception as e:
            logger.error(f"Failed for seed {seed}: {e}")
            results.append({
                "seed": seed,
                "hazard_ratio": np.nan,
                "hr_ci_lower": np.nan,
                "hr_ci_upper": np.nan,
                "logrank_p": np.nan,
                "c_index": np.nan,
                "median_treatment": np.nan,
                "median_control": np.nan,
                "ph_global_p": np.nan,
            })
    
    df_results = pd.DataFrame(results)
    
    # Calculate summary statistics
    print("\n" + "=" * 60)
    print("SENSITIVITY ANALYSIS SUMMARY")
    print("=" * 60)
    print(f"\nSeeds tested: {seeds}")
    print(f"Patients per cohort: {n_patients}")
    print("-" * 60)
    
    # Hazard ratio summary
    hr_mean = df_results["hazard_ratio"].mean()
    hr_std = df_results["hazard_ratio"].std()
    hr_min = df_results["hazard_ratio"].min()
    hr_max = df_results["hazard_ratio"].max()
    print(f"\nHazard Ratio:")
    print(f"  Mean: {hr_mean:.3f} (SD: {hr_std:.3f})")
    print(f"  Range: [{hr_min:.3f}, {hr_max:.3f}]")
    print(f"  All HR < 1.0 (treatment benefit): {(df_results['hazard_ratio'] < 1.0).all()}")
    
    # C-index summary
    ci_mean = df_results["c_index"].mean()
    ci_std = df_results["c_index"].std()
    ci_min = df_results["c_index"].min()
    ci_max = df_results["c_index"].max()
    print(f"\nConcordance Index:")
    print(f"  Mean: {ci_mean:.3f} (SD: {ci_std:.3f})")
    print(f"  Range: [{ci_min:.3f}, {ci_max:.3f}]")
    
    # Statistical significance consistency
    sig_count = (df_results["logrank_p"] < 0.05).sum()
    total_count = len(df_results)
    print(f"\nStatistical Significance (p < 0.05):")
    print(f"  {sig_count}/{total_count} seeds achieved significance")
    print(f"  All p-values < 0.05: {(df_results['logrank_p'] < 0.05).all()}")
    
    # Median survival times
    print(f"\nMedian Survival Times:")
    print(f"  Treatment: {df_results['median_treatment'].mean():.1f} months "
          f"(SD: {df_results['median_treatment'].std():.1f})")
    print(f"  Control: {df_results['median_control'].mean():.1f} months "
          f"(SD: {df_results['median_control'].std():.1f})")
    
    # Proportional hazards assumption
    ph_valid = (df_results["ph_global_p"] > 0.05).sum()
    print(f"\nProportional Hazards Assumption:")
    print(f"  Valid (p > 0.05): {ph_valid}/{total_count} seeds")
    
    print("-" * 60)
    
    # Save results
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_results.to_csv(output_path, index=False)
    logger.info(f"Detailed results saved to: {output_path}")
    
    # Add summary row
    summary_row = {
        "seed": "SUMMARY",
        "hazard_ratio": f"{hr_mean:.3f} (SD={hr_std:.3f})",
        "hr_ci_lower": f"{df_results['hr_ci_lower'].mean():.3f}",
        "hr_ci_upper": f"{df_results['hr_ci_upper'].mean():.3f}",
        "logrank_p": f"{sig_count}/{total_count} sig.",
        "c_index": f"{ci_mean:.3f} (SD={ci_std:.3f})",
        "median_treatment": f"{df_results['median_treatment'].mean():.1f}",
        "median_control": f"{df_results['median_control'].mean():.1f}",
        "ph_global_p": f"{ph_valid}/{total_count} valid",
    }
    
    print("\n" + "=" * 60)
    print("CONCLUSION")
    print("=" * 60)
    
    # Assess overall stability
    hr_cv = (hr_std / hr_mean) * 100 if hr_mean != 0 else float("inf")
    ci_cv = (ci_std / ci_mean) * 100 if ci_mean != 0 else float("inf")
    
    stability_issues = []
    if hr_cv > 20:
        stability_issues.append(f"HR coefficient of variation ({hr_cv:.1f}%) exceeds 20%")
    if ci_cv > 10:
        stability_issues.append(f"C-index CV ({ci_cv:.1f}%) exceeds 10%")
    if sig_count < total_count:
        stability_issues.append(f"Only {sig_count}/{total_count} seeds achieved significance")
    if ph_valid < total_count * 0.8:
        stability_issues.append(f"PH assumption violated in {total_count - ph_valid}/{total_count} seeds")
    
    if not stability_issues:
        print("\n[OK] Results are STABLE across random seeds.")
        print("  - Hazard ratio consistent (CV < 20%)")
        print("  - C-index consistent (CV < 10%)")
        print("  - All seeds achieved statistical significance")
        print("  - PH assumption valid in majority of runs")
    else:
        print("\n[WARN] Results show VARIABILITY across random seeds:")
        for issue in stability_issues:
            print(f"  - {issue}")
        print("\nRecommendation: Consider increasing sample size or investigating sources of instability.")
    
    print("=" * 60 + "\n")
    
    return df_results


def parse_args(args: List[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Assess result stability across random seeds",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python sensitivity_analysis.py
    python sensitivity_analysis.py --seeds 42 123 456 789 1011
    python sensitivity_analysis.py --n-patients 500 --output outputs/sensitivity_large.csv
        """,
    )
    parser.add_argument(
        "--seeds",
        type=int,
        nargs="+",
        default=[42, 123, 456, 789, 1011],
        help="Random seeds to test (default: 42 123 456 789 1011)",
    )
    parser.add_argument(
        "--n-patients",
        type=int,
        default=240,
        help="Number of patients per cohort (default: 240)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/sensitivity_analysis.csv"),
        help="Output CSV path (default: outputs/sensitivity_analysis.csv)",
    )
    return parser.parse_args(args)


def main() -> int:
    """Main entry point."""
    args = parse_args()
    
    try:
        run_sensitivity_analysis(
            seeds=args.seeds,
            n_patients=args.n_patients,
            output_path=args.output,
        )
        return 0
    except Exception as e:
        logger.error(f"Sensitivity analysis failed: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
