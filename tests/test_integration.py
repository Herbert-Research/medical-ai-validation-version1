"""
Integration tests for the complete validation pipeline.

These tests verify end-to-end execution without mocking.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import demo_quickstart as dq


class TestEndToEndPipeline:
    """Integration tests for demo_quickstart.py full execution."""

    @pytest.fixture
    def temp_output_dir(self, tmp_path: Path) -> Path:
        """Create temporary directory for test outputs."""
        return tmp_path

    def test_full_pipeline_executes_without_error(self, temp_output_dir: Path) -> None:
        """Verify complete pipeline runs successfully with default config."""
        config = dq.RunConfig(
            n_patients=50,  # Reduced for speed
            metastasis_samples=200,
            seed=42,
            output_dir=temp_output_dir,
            skip_tcga=True,  # Skip external data dependency
        )

        dq.run(config)

        # Verify expected outputs exist
        assert (temp_output_dir / "kaplan_meier_example.png").exists()
        assert (temp_output_dir / "calibration_curve_example.png").exists()
        assert (temp_output_dir / "synthetic_survival_summary.csv").exists()
        assert (temp_output_dir / "calibration_summary.csv").exists()

    def test_pipeline_with_tcga_data(self, temp_output_dir: Path) -> None:
        """Verify pipeline with TCGA data if available."""
        tcga_path = Path("data/tcga_2018_clinical_data.tsv")
        if not tcga_path.exists():
            pytest.skip("TCGA data not available")

        config = dq.RunConfig(
            n_patients=50,
            metastasis_samples=200,
            seed=42,
            output_dir=temp_output_dir,
            tcga_path=tcga_path,
            skip_tcga=False,
        )

        dq.run(config)

        assert (temp_output_dir / "tcga_survival_summary.csv").exists()
        assert (temp_output_dir / "tcga_group_summary.csv").exists()

    def test_output_csv_schema_validation(self, temp_output_dir: Path) -> None:
        """Verify output CSVs have expected columns."""
        import pandas as pd

        config = dq.RunConfig(
            n_patients=50,
            metastasis_samples=200,
            output_dir=temp_output_dir,
            skip_tcga=True,
        )
        dq.run(config)

        # Validate survival summary schema
        survival_df = pd.read_csv(temp_output_dir / "synthetic_survival_summary.csv")
        expected_cols = {
            "Cohort",
            "Treatment label",
            "Control label",
            "Median DFS (treatment)",
            "Median DFS (control)",
            "Log-rank p-value",
            "Hazard ratio",
            "C-index",
        }
        assert expected_cols.issubset(set(survival_df.columns))

        # Validate calibration summary schema
        calib_df = pd.read_csv(temp_output_dir / "calibration_summary.csv")
        assert "Model" in calib_df.columns
        assert "Brier Score" in calib_df.columns
        assert len(calib_df) == 3  # 3 calibration methods


class TestCLIInterface:
    """Tests for command-line argument parsing."""

    def test_all_cli_options_accepted(self) -> None:
        """Verify all documented CLI options are parsed."""
        config = dq.parse_cli_args(
            [
                "--n-patients",
                "100",
                "--metastasis-samples",
                "500",
                "--seed",
                "123",
                "--output-dir",
                "/tmp/test",
                "--calibration-strategy",
                "quantile",
                "--skip-tcga",
            ]
        )

        assert config.n_patients == 100
        assert config.metastasis_samples == 500
        assert config.seed == 123
        assert config.calibration_bin_strategy == "quantile"
        assert config.skip_tcga is True

    def test_invalid_calibration_strategy_rejected(self) -> None:
        """Verify invalid calibration strategy raises error."""
        with pytest.raises(SystemExit):
            dq.parse_cli_args(["--calibration-strategy", "invalid"])
