# Contributing to Medical AI Validation

Thank you for your interest in contributing to this research project.

## For Collaborators

This repository supports PhD research on AI-guided surgical validation. Contributions are welcome in the following areas:

1. **Methodology improvements** - Enhanced calibration techniques, alternative survival models
2. **External validation** - Application to additional clinical datasets
3. **Documentation** - Clarifications, additional clinical context

## Development Setup

```bash
# Clone and setup
git clone https://github.com/Herbert-Research/medical-ai-validation.git
cd medical-ai-validation
python -m venv .venv
source .venv/bin/activate  # or .venv\Scripts\Activate.ps1 on Windows
pip install -r requirements.txt
pip install -r requirements-dev.txt

# Verify setup
pytest tests/ -v
black --check *.py tests/*.py
```

## Pre-commit Hooks

This repository uses pre-commit hooks to enforce code quality:

```bash
pip install pre-commit
pre-commit install
pre-commit install --hook-type pre-push
```

Hooks run automatically on `git commit`. To run manually:

```bash
pre-commit run --all-files
```
