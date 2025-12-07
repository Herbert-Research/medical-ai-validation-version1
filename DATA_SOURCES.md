# Data Sources and Provenance

## Overview

This repository uses synthetic data to demonstrate survival analysis and calibration methodology for medical AI validation in gastric cancer surgical planning.

## Synthetic Data Generation

### Rationale for Synthetic Data

All data in this repository is **synthetic (simulated)** and generated using clinically informed parameters. Synthetic data was chosen to:

1. **Demonstrate methodology** without requiring access to protected health information
2. **Enable reproducibility** - anyone can regenerate the exact same results
3. **Avoid patient privacy concerns** - no real patient data is used or required
4. **Provide educational value** - parameter choices reflect real-world clinical relationships

### Data Generation Parameters

#### 1. Synthetic Survival Analysis (Standard vs AI-guided)

File: `demo_quickstart.py` -> Function: `build_survival_dataset()`

Clinical Basis:
- Cohort size: 240 patients (120 per arm)
- Tumor stage: I-III (realistic distribution)
- Lymph node positivity: 35% (matches gastric cancer epidemiology)
- Nerve invasion: 25% (typical for advanced gastric cancer)

Survival Model:
- Base monthly hazard: 0.018 (corresponds to median DFS ~38 months for standard care)
- AI-guided hazard ratio: ~0.52 (48% risk reduction - optimistic but plausible)
- Censoring rate: ~30% (typical for 5-year surgical oncology studies)

Results:
- Median DFS (AI-guided): ~28 months
- Median DFS (Standard): ~10 months
- Log-rank p-value: < 0.001
- C-index: ~0.70

Clinical Plausibility:
A 48% risk reduction is optimistic but within range of practice-changing surgical innovations (for example, sentinel node mapping or neo-adjuvant therapy). This represents an upper-bound effect size for validation methodology demonstration.

#### 2. Synthetic Gastric Cancer Node Status Analysis

File: `generate_synthetic_gastric_cancer_data.py`

Clinical Basis:
- Cohort size: 394 patients (matching TCGA-STAD approximate size)
- Node-negative: ~31% (realistic for surgical cohorts)
- Node-positive: ~69%
- Age: mean 63 years (SD 11) - matches gastric cancer demographics
- Additional features: tumor size, Lauren classification

Survival Model:
- Base monthly hazard: 0.015 (calibrated for stage II-III gastric cancer)
- Node-positive hazard ratio: ~1.7-1.9 (matches published literature)
- Median survival: 35-45 months (node-negative), 25-35 months (node-positive)
- Censoring rate: ~35%

Results:
- Hazard ratio: ~1.8 (95% CI: 1.2-2.9)
- Log-rank p-value: < 0.01
- C-index: ~0.63

Clinical Plausibility:
These parameters closely match published TCGA-STAD analyses and recent gastric cancer surgical trials. The hazard ratio is consistent with meta-analyses of node-positive versus node-negative disease.

Reference Ranges from Literature:
- Node-positive HR: 1.5-2.3 (various studies)
- 5-year DFS: 40-60% (node-negative), 20-40% (node-positive)
- Median DFS: 30-50 months overall for stage II-III

#### 3. Metastasis Risk Calibration Dataset

File: `demo_quickstart.py` -> Function: `simulate_metastasis_dataset()`

Clinical Basis:
- Sample size: 1000 stations
- Features: tumor stage, lymph node count, station distance, histology score, BMI
- Base metastasis risk: sigmoid function calibrated to ~25% overall metastasis rate

Results:
- Class balance: ~25% positive (metastasis present)
- Brier scores: 0.14 (reasonable for binary classification)
- Feature relationships: match clinical intuition (higher stage -> higher risk, and similar relationships for other features)

## Relationship to Real-World Data

While this data is entirely synthetic, the generative process was informed by:

1. **TCGA-STAD Database:** Public gastric adenocarcinoma dataset with survival and clinical annotations
   - Used for parameter calibration only, not as source data
   - Available at: https://portal.gdc.cancer.gov/

2. **Published Literature:**
   - Kim et al. (2022) - KLASS-02 trial results
   - Japanese Gastric Cancer Association treatment guidelines
   - AJCC 8th edition staging system

3. **Clinical Consultation:**
   - Parameters reviewed for plausibility with surgical oncology domain experts
   - Event rates, hazard ratios, and survival times within published ranges

## Intended Use

Appropriate uses:
- Demonstrating statistical methodology
- Teaching survival analysis and calibration concepts
- Validating software workflows
- PhD application portfolio demonstration

Not for:
- Clinical decision-making
- Representing actual patient outcomes
- Publication as novel clinical findings

## Reproducibility

All data can be regenerated exactly using the provided scripts with fixed random seeds:

```bash
python demo_quickstart.py  # Generates calibration and AI-guided survival data
python generate_synthetic_gastric_cancer_data.py  # Generates node status analysis
```

## For PhD Research Application

As part of this PhD research, the validated methodologies will be applied to:

1. Real clinical datasets with appropriate IRB approval and data use agreements
2. Multi-institutional cohorts for external validation
3. Prospective surgical trial data from ongoing AI-guided gastrectomy studies

The current synthetic data demonstrates the ability to:
- Implement appropriate statistical methods
- Generate clinically plausible test cases
- Validate computational workflows
- Understand survival analysis principles

## Questions or Clarifications

For questions about data generation parameters or clinical plausibility, please contact maximilian.dressler@stud.uni-heidelberg.de.

Last updated: 2025-11-22
