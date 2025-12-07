# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-11-11

### Added
- Synthetic survival cohort generation with KLASS-02-informed clinical parameters
- Kaplan-Meier and Cox proportional hazards analysis pipeline with hazard ratio reporting
- Calibration assessment framework computing ECE, Brier score, and calibration slope
- Three calibration methods comparison: uncalibrated, Platt scaling, isotonic regression
- TCGA-STAD external validation pilot workflow for institutional readiness
- Statistical power analysis script for prospective study design (`statistical_power_analysis.py`)
- Docker containerization for reproducible execution across environments
- Comprehensive test suite with 43 tests covering unit and integration scenarios
- Repository self-validation script (`validate_repository.py`) for artifact consistency
- CLI interface with configurable parameters (`--n-patients`, `--seed`, `--output-dir`)
- Literature comparison documentation benchmarking results against published trials
- Data provenance documentation for synthetic and TCGA clinical streams
- API reference and validation playbook linking PhD aims to code artifacts

### Technical Details
- Python 3.9+ with type annotations throughout
- Fixed random seed (42) for deterministic reproducibility
- Pinned dependencies via `requirements-lock.txt` (192 transitive packages)
- PEP8-compliant code structure with modular design

### Known Limitations
- Expected Calibration Error exceeds clinical threshold (best achieved: 0.054 vs target ≤0.05)
- Platt scaling inappropriate for this dataset (calibration slope = 3.52, should be ≈1.0)
- Synthetic hazard ratio (0.52) represents optimistic upper bound vs published trials (0.85-0.92)
- Single synthetic cohort; cross-validation sensitivity analysis pending
- Sample size (n=240) insufficient for deployment; power analysis indicates n≥2000 needed

## [Unreleased]

### Planned
- Cross-validation framework for calibration stability assessment
- Sensitivity analysis across multiple random seeds
- Additional calibration methods (temperature scaling, beta calibration)
- Multi-institutional external validation beyond TCGA pilot
- Prospective data ingestion pipeline for IRB-approved institutional cohorts
- Interactive dashboard for real-time calibration monitoring
