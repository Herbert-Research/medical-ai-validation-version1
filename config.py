"""
Centralized configuration constants for the medical AI validation pipeline.

These values are derived from clinical literature and KLASS-02 trial parameters.
Modifications should be documented with literature references.
"""

from __future__ import annotations

# =============================================================================
# Clinical Thresholds (Literature-Derived)
# =============================================================================

# Expected Calibration Error threshold for clinical deployment
# Reference: Van Calster et al., BMC Medicine, 2019
ECE_CLINICAL_THRESHOLD: float = 0.05

# Acceptable ECE with appropriate risk communication
ECE_ACCEPTABLE_THRESHOLD: float = 0.08

# =============================================================================
# Synthetic Data Generation Parameters
# =============================================================================

# Base monthly hazard rate (corresponds to ~38-month median DFS for standard care)
# Reference: Kim et al., KLASS-02 trial, 2022
BASE_HAZARD_MONTHLY: float = 0.018

# Censoring rate for 5-year surgical oncology studies
CENSORING_RATE: float = 0.012

# Default random seed for reproducibility
DEFAULT_SEED: int = 42

# Default cohort sizes
DEFAULT_N_PATIENTS: int = 240
DEFAULT_METASTASIS_SAMPLES: int = 1000

# =============================================================================
# Clinical Feature Distributions
# =============================================================================

# Lymph node positivity rate in gastric cancer (epidemiological reference)
LYMPH_POSITIVE_RATE: float = 0.35

# Nerve invasion rate for advanced gastric cancer
NERVE_INVASION_RATE: float = 0.28

# Vascular invasion rate
VASCULAR_INVASION_RATE: float = 0.22

# Margin positivity rate
MARGIN_POSITIVE_RATE: float = 0.12

# =============================================================================
# Cox Model Coefficients (Clinically Informed)
# =============================================================================

# Hazard ratio coefficients for synthetic data generation
COEF_TUMOR_STAGE: float = 0.42
COEF_LYMPH_POSITIVE: float = 0.55
COEF_NERVE_INVASION: float = 0.35
COEF_VASCULAR_INVASION: float = 0.28
COEF_MARGIN_POSITIVE: float = 0.42
COEF_AGE_PER_YEAR: float = 0.014
COEF_ALBUMIN_PER_UNIT: float = -0.25
COEF_AI_BENEFIT_BASE: float = -0.6

# =============================================================================
# Calibration Parameters
# =============================================================================

# Number of bins for ECE calculation
ECE_N_BINS: int = 10

# Gradient Boosting hyperparameters
GB_N_ESTIMATORS: int = 250
GB_LEARNING_RATE: float = 0.05
