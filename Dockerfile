# Medical AI Validation - Reproducible Environment
# Python 3.9 base image for compatibility with lifelines and scikit-learn

FROM python:3.9-slim

LABEL maintainer="Maximilian Herbert Dressler"
LABEL description="Reproducible environment for medical AI survival and calibration analysis"

# Set working directory
WORKDIR /app

# Install system dependencies for matplotlib
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    libffi-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for layer caching
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy project files
COPY . .

# Create outputs directory
RUN mkdir -p outputs

# Default command: run the demo
CMD ["python", "demo_quickstart.py", "--output-dir", "outputs/"]
