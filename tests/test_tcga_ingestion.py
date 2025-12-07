"""
Tests for TCGA clinical data ingestion and preprocessing.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import demo_quickstart as dq


class TestStatusFlagParsing:
    """Tests for parse_status_flag function."""

    @pytest.mark.parametrize(
        "input_val,expected",
        [
            ("0:LIVING", 0.0),
            ("1:DECEASED", 1.0),
            ("0:CENSORED", 0.0),
            ("1:RECURRED", 1.0),
            ("invalid", np.nan),
            ("", np.nan),
            (None, np.nan),
            (123, np.nan),
        ],
    )
    def test_status_flag_parsing(self, input_val, expected) -> None:
        result = dq.parse_status_flag(input_val)
        if np.isnan(expected):
            assert np.isnan(result)
        else:
            assert result == expected


class TestStageMapping:
    """Tests for AJCC stage code parsing."""

    @pytest.mark.parametrize(
        "input_val,expected",
        [
            ("T1", 2.0),
            ("T2A", 4.0),
            ("T3", 8.0),
            ("t4b", 11.0),  # Case insensitive
            ("  T2  ", 5.0),  # Whitespace handling
            ("TX", np.nan),  # Unknown stage
            ("", np.nan),
            (None, np.nan),
        ],
    )
    def test_stage_code_to_numeric(self, input_val, expected) -> None:
        result = dq.stage_code_to_numeric(input_val)
        if np.isnan(expected):
            assert np.isnan(result)
        else:
            assert result == expected


class TestGradeMapping:
    """Tests for histologic grade parsing."""

    @pytest.mark.parametrize(
        "input_val,expected",
        [
            ("G1", 1.0),
            ("G2", 2.0),
            ("G3", 3.0),
            ("G4", 4.0),
            ("g2", 2.0),  # Case insensitive
            ("GX", np.nan),
            (None, np.nan),
        ],
    )
    def test_grade_to_numeric(self, input_val, expected) -> None:
        result = dq.grade_to_numeric(input_val)
        if np.isnan(expected):
            assert np.isnan(result)
        else:
            assert result == expected


class TestTCGACohortPreparation:
    """Tests for prepare_tcga_survival_cohort function."""

    def test_missing_required_column_raises_keyerror(self) -> None:
        """Verify missing columns are detected."""
        incomplete_df = pd.DataFrame(
            {
                "Patient ID": ["P1", "P2"],
                "Progress Free Survival (Months)": [10.0, 20.0],
                # Missing other required columns
            }
        )

        with pytest.raises(KeyError, match="missing required columns"):
            dq.prepare_tcga_survival_cohort(incomplete_df)

    def test_invalid_n_stage_filtered(self) -> None:
        """Verify invalid N-stage values are excluded."""
        df = pd.DataFrame(
            {
                "Patient ID": ["P1", "P2", "P3"],
                "Progress Free Survival (Months)": [10.0, 20.0, 30.0],
                "Disease Free (Months)": [None, None, None],
                "Progression Free Status": ["1:RECURRED", "0:CENSORED", "1:RECURRED"],
                "Disease Free Status": [None, None, None],
                "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code": [
                    "N0",
                    "NX",
                    "N1",
                ],  # NX should be filtered
                "American Joint Committee on Cancer Tumor Stage Code": [
                    "T1",
                    "T2",
                    "T3",
                ],
                "Neoplasm Histologic Grade": ["G1", "G2", "G3"],
                "Diagnosis Age": [60, 65, 70],
            }
        )

        cohort = dq.prepare_tcga_survival_cohort(df)

        # Only N0 and N1 should remain
        assert len(cohort) == 2
        assert "NX" not in cohort["group"].values
