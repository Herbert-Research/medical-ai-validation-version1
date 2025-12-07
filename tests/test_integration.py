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


class TestNotebookExecution:
    """Tests for Jupyter notebook validity and structure."""

    def test_notebook_syntax_valid(self) -> None:
        """Verify notebook is valid JSON and has expected structure."""
        import json

        notebook_path = Path("survival_and_calibration_enhanced.ipynb")
        if not notebook_path.exists():
            pytest.skip("Notebook not present")

        with open(notebook_path, "r", encoding="utf-8") as f:
            nb = json.load(f)

        # Verify required top-level keys exist
        assert "cells" in nb, "Notebook must have 'cells' key"
        assert "metadata" in nb, "Notebook must have 'metadata' key"
        assert "nbformat" in nb, "Notebook must have 'nbformat' key"
        assert len(nb["cells"]) > 0, "Notebook must have at least one cell"

    def test_notebook_has_code_cells(self) -> None:
        """Verify notebook contains executable code cells."""
        import json

        notebook_path = Path("survival_and_calibration_enhanced.ipynb")
        if not notebook_path.exists():
            pytest.skip("Notebook not present")

        with open(notebook_path, "r", encoding="utf-8") as f:
            nb = json.load(f)

        code_cells = [c for c in nb["cells"] if c.get("cell_type") == "code"]
        assert len(code_cells) > 0, "Notebook should have code cells"

    def test_notebook_cells_have_valid_structure(self) -> None:
        """Verify each cell has required structure fields."""
        import json

        notebook_path = Path("survival_and_calibration_enhanced.ipynb")
        if not notebook_path.exists():
            pytest.skip("Notebook not present")

        with open(notebook_path, "r", encoding="utf-8") as f:
            nb = json.load(f)

        for i, cell in enumerate(nb["cells"]):
            assert "cell_type" in cell, f"Cell {i} missing 'cell_type'"
            assert "source" in cell, f"Cell {i} missing 'source'"
            assert cell["cell_type"] in {
                "code",
                "markdown",
                "raw",
            }, f"Cell {i} has invalid cell_type: {cell.get('cell_type')}"

    def test_notebook_metadata_has_kernel_info(self) -> None:
        """Verify notebook metadata includes kernel specification."""
        import json

        notebook_path = Path("survival_and_calibration_enhanced.ipynb")
        if not notebook_path.exists():
            pytest.skip("Notebook not present")

        with open(notebook_path, "r", encoding="utf-8") as f:
            nb = json.load(f)

        metadata = nb.get("metadata", {})
        # Check for either kernelspec or language_info
        has_kernel_info = "kernelspec" in metadata or "language_info" in metadata
        assert has_kernel_info, "Notebook metadata should include kernel information"


class TestSensitivityAnalysis:
    """Tests for sensitivity analysis script."""

    def test_sensitivity_analysis_runs(self, tmp_path: Path) -> None:
        """Verify sensitivity analysis executes without error."""
        import sensitivity_analysis as sa

        # Run with minimal seeds and patients for speed
        output_path = tmp_path / "sensitivity_test.csv"
        results = sa.run_sensitivity_analysis(
            seeds=[42, 123],
            n_patients=50,
            output_path=output_path,
        )

        assert len(results) == 2
        assert "hazard_ratio" in results.columns
        assert "c_index" in results.columns
        assert output_path.exists()

    def test_sensitivity_results_schema(self, tmp_path: Path) -> None:
        """Verify sensitivity analysis output has expected columns."""
        import sensitivity_analysis as sa

        output_path = tmp_path / "sensitivity_schema.csv"
        results = sa.run_sensitivity_analysis(
            seeds=[42],
            n_patients=50,
            output_path=output_path,
        )

        expected_columns = {
            "seed",
            "hazard_ratio",
            "hr_ci_lower",
            "hr_ci_upper",
            "logrank_p",
            "c_index",
            "median_treatment",
            "median_control",
            "ph_global_p",
        }
        assert expected_columns.issubset(set(results.columns))
