#!/usr/bin/env bash
set -e

# Check if conda is installed
check_conda_installed() {
    command -v conda >/dev/null 2>&1
}

# Get the latest Anaconda Linux x86_64 installer URL
get_latest_anaconda_url() {
    archive_url="https://repo.anaconda.com/archive/"

    # Debug message to stderr (won't be captured by command substitution)
    echo "Retrieving latest Anaconda installer URL..." >&2

    latest_filename=$(curl -s "$archive_url" | \
        grep -o 'Anaconda3-[0-9.]\+-[0-9]\+-Linux-x86_64.sh' | \
        sort -V | \
        tail -n 1)

    if [[ -z "$latest_filename" ]]; then
        echo "❌ Could not find the latest Anaconda installer." >&2
        return 1
    fi

    # Only echo the final URL (captured by caller)
    echo "${archive_url}${latest_filename}"
}

# Download Anaconda installer
download_anaconda() {
    installer_url=$(get_latest_anaconda_url)
    
    if [[ -z "$installer_url" ]]; then
        echo "❌ Installer URL is empty. Aborting download."
        exit 1
    fi

    echo "✅ Downloading installer from: $installer_url"
    wget -O ~/anaconda.sh "$installer_url"
}

# Install Anaconda
install_anaconda() {
    echo "Installing Anaconda..."
    echo "Follow the prompts to complete the installation."
    echo "Press Ctrl+C to cancel the installation if you don't want to proceed."
    echo "You can also run the installer manually: bash ~/anaconda.sh"
    bash ~/anaconda.sh -b -p "$HOME/anaconda3"
    echo "Initializing conda..."
    "$HOME/anaconda3/bin/conda" init
    echo "✅ Anaconda installed. Restart your terminal or run: source ~/.bashrc"
}

# Main
if check_conda_installed; then
    echo "✅ Anaconda is already installed."
else
    echo "⚠️  Anaconda not found. Installing latest version..."
    download_anaconda
    install_anaconda
fi

# Clean up
echo "Cleaning up..."
# rm ~/anaconda.sh
echo "✅ Installation script completed."

echo ""
echo "🔧 Setting up conda environment for bactipipe..."
echo ""
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
