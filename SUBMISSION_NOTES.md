# Submission Notes - PhD Application Portfolio

## Repository Purpose

This repository demonstrates my technical competency in:
- Survival analysis methodology (Kaplan-Meier, Cox proportional hazards)
- Probabilistic calibration assessment and remediation
- Reproducible research practices
- Clinical AI validation workflows

## Key Points for Reviewers

1. All data is synthetic - Generated using clinically informed parameters to demonstrate methodology without requiring protected health information. See DATA_SOURCES.md for complete provenance.
2. Calibration challenges are intentional - The current results show realistic calibration difficulties that motivate my proposed PhD research on improving calibration for surgical AI systems.
3. Reproducibility - All results can be regenerated from scratch using provided scripts with fixed random seeds.
4. Time to execute - Complete workflow runs in under 2 minutes on standard hardware.

## Files Overview

**Core Analysis:**
- `demo_quickstart.py` - Main analysis script
- `outputs/calibration_summary.csv` - Calibration metrics
- `outputs/synthetic_survival_summary.csv` - Survival analysis results
- `outputs/kaplan_meier_example.png` - Survival curves (Standard vs AI-guided)
- `outputs/calibration_curve_example.png` - Calibration visualization

**Documentation:**
- `README.md` - Primary documentation
- `DATA_SOURCES.md` - Data provenance and methodology
- `LITERATURE_COMPARISON.md` - Contextualization within field
- `requirements.txt` - Python dependencies

**Testing and Validation:**
- `tests/test_demo_quickstart.py` - Unit tests
- `validate_repository.py` - Data consistency checks

**Supplementary:**
- `statistical_power_analysis.py` - Sample size calculations
- `power_analysis.png` - Power curve visualizations
- `outputs/tcga_survival_summary.csv` and `outputs/tcga_group_summary.csv` - TCGA pilot summaries (when data available)

## Execution Instructions

```bash
# 1. Set up environment
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 2. Run main analysis
python3 demo_quickstart.py

# 3. Run tests
python3 -m pytest tests/test_demo_quickstart.py

# 4. Validate consistency
python3 validate_repository.py
```
