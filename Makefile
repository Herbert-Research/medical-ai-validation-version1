# Makefile for Medical AI Validation
# ===================================
# Standardized development workflow commands
#
# Usage: make <target>
# Run 'make help' to see all available targets

.PHONY: help install install-dev test test-cov lint format run clean \
        docker-build docker-run docker-test validate pre-commit \
        check-all notebook

# Default Python and pip commands
PYTHON := python
PIP := pip
PYTEST := pytest

# Project directories
OUTPUT_DIR := outputs
TEST_DIR := tests

# Colors for terminal output (works on Unix/macOS, ignored on Windows)
BLUE := \033[36m
GREEN := \033[32m
YELLOW := \033[33m
RESET := \033[0m

# ============================================================================
# HELP
# ============================================================================

help:  ## Show this help message
	@echo "Medical AI Validation - Development Commands"
	@echo "============================================="
	@echo ""
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "$(BLUE)%-20s$(RESET) %s\n", $$1, $$2}'

# ============================================================================
# INSTALLATION
# ============================================================================

install:  ## Install production dependencies
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt

install-dev:  ## Install all dependencies (production + development)
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt
	$(PIP) install -r requirements-dev.txt
	pre-commit install
	@echo "$(GREEN)Development environment ready. Pre-commit hooks installed.$(RESET)"

install-lock:  ## Install exact locked versions for reproducibility
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements-lock.txt

# ============================================================================
# TESTING
# ============================================================================

test:  ## Run test suite
	$(PYTEST) $(TEST_DIR)/ -v --tb=short

test-cov:  ## Run tests with coverage report
	$(PYTEST) $(TEST_DIR)/ -v --tb=short --cov=. --cov-report=term-missing --cov-report=html
	@echo "$(GREEN)Coverage report generated in htmlcov/index.html$(RESET)"

test-fast:  ## Run tests without slow integration tests
	$(PYTEST) $(TEST_DIR)/ -v --tb=short -m "not slow"

# ============================================================================
# CODE QUALITY
# ============================================================================

lint:  ## Run all linters (black, isort, mypy)
	@echo "Checking code formatting with black..."
	black --check *.py $(TEST_DIR)/*.py
	@echo "Checking import order with isort..."
	isort --check-only --profile black *.py $(TEST_DIR)/*.py
	@echo "Running type checker (mypy)..."
	mypy demo_quickstart.py config.py validate_repository.py --ignore-missing-imports
	@echo "$(GREEN)All lint checks passed!$(RESET)"

format:  ## Auto-format code (black + isort)
	black *.py $(TEST_DIR)/*.py
	isort --profile black *.py $(TEST_DIR)/*.py
	@echo "$(GREEN)Code formatted successfully.$(RESET)"

pre-commit:  ## Run pre-commit hooks on all files
	pre-commit run --all-files

check-all:  ## Run all quality checks (lint + test + validate)
	@echo "$(YELLOW)Running comprehensive quality checks...$(RESET)"
	$(MAKE) lint
	$(MAKE) test
	$(MAKE) validate
	@echo "$(GREEN)All checks passed!$(RESET)"

# ============================================================================
# EXECUTION
# ============================================================================

run:  ## Run main analysis pipeline
	$(PYTHON) demo_quickstart.py --output-dir $(OUTPUT_DIR)/

run-quick:  ## Run pipeline with reduced sample size (faster)
	$(PYTHON) demo_quickstart.py --output-dir $(OUTPUT_DIR)/ --n-patients 100

validate:  ## Validate repository consistency
	$(PYTHON) validate_repository.py

power-analysis:  ## Run statistical power analysis
	$(PYTHON) statistical_power_analysis.py

# ============================================================================
# JUPYTER NOTEBOOK
# ============================================================================

notebook:  ## Launch Jupyter notebook server
	jupyter notebook survival_and_calibration_enhanced.ipynb

# ============================================================================
# DOCKER
# ============================================================================

docker-build:  ## Build Docker image
	docker build -t medical-ai-validation .

docker-run:  ## Run analysis pipeline in Docker container
	docker run -v $(PWD)/$(OUTPUT_DIR):/app/outputs medical-ai-validation

docker-test:  ## Test Docker image (import check)
	docker run --rm medical-ai-validation python -c "import demo_quickstart; print('Docker OK')"

docker-shell:  ## Open interactive shell in Docker container
	docker run -it --rm -v $(PWD):/app medical-ai-validation /bin/bash

# ============================================================================
# CLEANUP
# ============================================================================

clean:  ## Remove generated artifacts and caches
	rm -rf __pycache__ $(TEST_DIR)/__pycache__
	rm -rf .pytest_cache .mypy_cache
	rm -rf .coverage htmlcov/
	rm -rf $(OUTPUT_DIR)/.mplconfig
	rm -rf *.egg-info dist/ build/
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -delete
	@echo "$(GREEN)Cleanup complete.$(RESET)"

clean-outputs:  ## Remove all generated output files
	rm -rf $(OUTPUT_DIR)/*.png $(OUTPUT_DIR)/*.csv
	@echo "$(GREEN)Output files removed.$(RESET)"

clean-all: clean clean-outputs  ## Remove all generated files and caches
	@echo "$(GREEN)Full cleanup complete.$(RESET)"

# ============================================================================
# DOCUMENTATION
# ============================================================================

docs-serve:  ## Serve documentation locally (if using mkdocs/sphinx)
	@echo "Documentation server not configured. See docs/ directory."

# ============================================================================
# RELEASE
# ============================================================================

checksums:  ## Generate checksums for release artifacts
	$(PYTHON) -c "from validate_repository import generate_checksums; generate_checksums()"
