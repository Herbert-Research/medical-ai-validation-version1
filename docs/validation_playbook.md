# Validation Playbook

This repository is the code counterpart to my KLASS-focused PhD research statement. Each aim below references the components reviewers see in `demo_quickstart.py`, the notebooks, and exported artefacts so they can trace intent to implementation.

## Aim 1 – Station-Specific Risk Prediction
- **Synthetic cohorts (`build_survival_dataset`, `simulate_metastasis_dataset`)** encode the BIS/NOL/TOF-linked covariates I managed in the ENiMoN RCT. Coefficients mirror intraoperative heuristics described in the CV to show I can translate clinical priors into model-ready features.
- **Calibration stack (`evaluate_calibration_models`, `calibration_slope_intercept`)** demonstrates the decision-focused calibration workflow (raw vs Platt vs isotonic) promised in the research statement. Outputs are persisted to `calibration_summary.csv` for auditability.
- **`synthetic_survival_summary.csv` + `kaplan_meier_example.png`** provide the “validation-first” evidence package that would be shared with KLASS surgeons before activating a real console overlay.

## Aim 2 – Outcome-Linked Quality Metrics
- **TCGA ingestion (`prepare_tcga_survival_cohort`)** shows rigorous data governance: explicit parsing, AJCC mapping, and censoring logic aligned with KLASS-02-QC documentation.
- **Survival diagnostics (`analyse_survival`, Schoenfeld tests, Cox PH)** quantify the disease-free survival linkage required for Aim 2. Results stream to `tcga_survival_summary.csv`, `tcga_group_summary.csv`, and `tcga_kaplan_meier.png`.
- These artefacts are the template for the stepped-wedge trial reporting pipeline described in the statement.

## Aim 3 – Multimodal Fusion Pilot
- While full multimodal fusion will live in future modules, `plot_calibration_curves` and the saved PNG outputs document my latency-aware plotting and reporting workflow.
- The outlined structure leaves hooks for RGB/ICG inputs (see `simulate_metastasis_dataset` placeholders), matching the Aim 3 feasibility plan.

## Governance & CV Alignment
- **Clinical insight**: Parameters and acceptance thresholds copy directly from procedures observed at SNUH (97 cases) and ENiMoN trial SOPs listed in the CV.
- **Regulatory mindset**: Reproducible CSV artefacts, path-logged figures, and pre-registered style summaries mirror the GDPR/MFDS planning section of the statement.
- **Next steps**: Future commits will extend this playbook with the multimodal inference code and DAQ scripts referenced in the research statement’s Year 2–3 roadmap.
