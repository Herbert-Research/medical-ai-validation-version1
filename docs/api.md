# API Reference

## Core Modules

### `demo_quickstart`

Main entry point for the validation pipeline.

#### Data Generation

| Function | Description |
|----------|-------------|
| `build_survival_dataset(n_patients, seed)` | Generate synthetic survival cohort |
| `build_synthetic_gastric_node_dataset(n_patients, seed)` | Generate gastric cancer node-status cohort |
| `simulate_metastasis_dataset(n_samples, seed)` | Generate metastasis prediction dataset |

#### Analysis Functions

| Function | Description |
|----------|-------------|
| `analyse_survival(data, output_path, ...)` | Kaplan-Meier + Cox PH analysis |
| `evaluate_calibration_models(X_train, X_test, y_train, y_test, ...)` | Compare calibration methods |
| `expected_calibration_error(y_true, y_prob, n_bins, strategy)` | Compute ECE metric |
| `calibration_slope_intercept(y_true, y_prob)` | Compute calibration slope/intercept |

#### Data Classes

| Class | Description |
|-------|-------------|
| `RunConfig` | Pipeline configuration parameters |
| `SurvivalResults` | Container for survival analysis outputs |

### `validate_repository`

Repository consistency validation.

| Function | Description |
|----------|-------------|
| `run()` | Execute all validation checks |
| `check_required_files(required)` | Verify file existence |
| `check_survival_dataset(...)` | Validate survival summary consistency |

### `config`

Centralized configuration constants.

See [Configuration Constants](#configuration-constants) for full reference.

## Configuration Constants

| Constant | Value | Reference |
|----------|-------|-----------|
| `ECE_CLINICAL_THRESHOLD` | 0.05 | Van Calster et al., 2019 |
| `BASE_HAZARD_MONTHLY` | 0.018 | KLASS-02 trial |
| `DEFAULT_SEED` | 42 | - |

See `config.py` for complete listing.
