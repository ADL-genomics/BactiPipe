#!/bin/bash

# Name of the conda environment (from environment.yml)
ENV_NAME="bactipipe"

echo "📦 Creating conda environment: $ENV_NAME"
conda env create -f environment.yml

echo "✅ Environment created. Activating..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

echo "🔍 Verifying bactipipe installation..."
if command -v bactipipe >/dev/null 2>&1; then
    echo "✅ bactipipe is available as a CLI command."
else
    echo "❌ ERROR: bactipipe command not found."
    echo "Did pip install -e . fail? Check your environment.yml"
    exit 1
fi

echo "🚀 Optionally checking for updates..."
bactipipe check-updates

echo "🎉 Setup complete. You're ready to use bactipipe!"
