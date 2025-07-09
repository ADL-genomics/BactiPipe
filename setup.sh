#!/usr/bin/env bash
set -e

# 1) Detect host OS and set the "platform" string used on repo.anaconda.com
detect_platform() {
  local os
  os=$(uname -s)
  case "$os" in
    Linux)   PLATFORM="Linux"   ;;
    Darwin)  PLATFORM="MacOSX"  ;;  # macOS installers use “MacOSX”
    *) 
      echo "❌ Unsupported OS: $os" >&2
      exit 1
      ;;
  esac
}

# 2) Get the latest Anaconda x86_64 installer URL for that platform
get_latest_anaconda_url() {
  archive_url="https://repo.anaconda.com/archive/"
  echo "🔍 Looking for latest Anaconda3-x86_64 for ${PLATFORM}..." >&2

  latest_filename=$(
    curl -s "$archive_url" |
    grep -oE "Anaconda3-[0-9]+\.[0-9]+-[0-9]+-${PLATFORM}-x86_64\.sh" |
    sort -V |
    tail -n1
  )

  if [[ -z "$latest_filename" ]]; then
    echo "❌ Could not find an Anaconda installer for ${PLATFORM}-x86_64." >&2
    return 1
  fi

  echo "${archive_url}${latest_filename}"
}

# 3) Download + install
download_anaconda() {
  installer_url=$(get_latest_anaconda_url)
  echo "✅ Downloading: $installer_url"
  wget -O ~/anaconda.sh "$installer_url"
}

install_anaconda() {
  echo "🚀 Installing Anaconda into \$HOME/anaconda3…"
  bash ~/anaconda.sh -b -p "$HOME/anaconda3"
  echo "🔧 Initializing conda…"
  "$HOME/anaconda3/bin/conda" init
  echo "✅ Installed and initialized. Restart your shell or run: source ~/.bashrc"
}

# ——— Main ———
detect_platform

if command -v conda >/dev/null 2>&1; then
  echo "✅ conda already installed."
else
  echo "⚠️  conda not found; installing latest Anaconda3-x86_64 for ${PLATFORM}…"
  download_anaconda
  install_anaconda
fi

# Clean up
echo "Cleaning up..."
# rm ~/anaconda.sh
echo "✅ Conda Installation script completed."

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

# Install basespace-cli
echo "🔧 Installing basespace-cli..."
if [ "$PLATFORM" == "Linux" ]; then
    wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O $HOME/bin/bs
elif [ "$PLATFORM" == "MacOSX" ]; then
    wget "https://launch.basespace.illumina.com/CLI/latest/amd64-osx/bs" -O $HOME/bin/bs
fi

chmod +x $HOME/bin/bs
echo "🎉 Setup complete. You're ready to use bactipipe!"
