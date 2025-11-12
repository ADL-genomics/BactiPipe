#!/usr/bin/env bash
set -eo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${REPO_ROOT}/envs"

# ---------- helpers ----------
msg() { echo -e "\033[1;36m$*\033[0m"; }
ok()  { echo -e "\033[1;32m$*\033[0m"; }
err() { echo -e "\033[1;31m$*\033[0m" 1>&2; }

detect_platform() {
  case "$(uname -s)" in
    Linux)   PLATFORM="Linux" ;;
    Darwin)  PLATFORM="MacOSX" ;;
    *) err "Unsupported OS: $(uname -s)"; exit 1 ;;
  esac
}

ensure_conda() {
  if command -v conda >/dev/null 2>&1; then
    ok "conda already installed."
    return
  fi
  msg "conda not found; installing latest Anaconda3 for ${PLATFORM}…"
  ARCHIVE="https://repo.anaconda.com/archive/"
  FILENAME=$(curl -s "${ARCHIVE}" | grep -oE "Anaconda3-[0-9]+\.[0-9]+-[0-9]+-${PLATFORM}-x86_64\.sh" | sort -V | tail -n1)
  [[ -n "${FILENAME}" ]] || { err "Could not resolve Anaconda installer."; exit 1; }
  wget -O "${HOME}/anaconda.sh" "${ARCHIVE}${FILENAME}"
  bash "${HOME}/anaconda.sh" -b -p "${HOME}/anaconda3"
  "${HOME}/anaconda3/bin/conda" init
  ok "Anaconda installed. Open a new shell or: source ~/.bashrc"
  # Re-exec shell so 'conda' is on PATH
  source "${HOME}/anaconda3/etc/profile.d/conda.sh"
}

ensure_mamba() {
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate base
  if ! command -v mamba >/dev/null 2>&1; then
    msg "Installing mamba in base…"
    conda install -y -n base -c conda-forge mamba
  fi
  ok "mamba ready."
}

create_env() {
  local yml="$1"
  [[ -f "${yml}" ]] || { err "Missing environment file: ${yml}"; exit 1; }
  msg "Creating env from ${yml}…"
  mamba env create -f "${yml}" || {
    msg "Env exists; updating instead…"
    mamba env update -f "${yml}" --prune
  }
}

pip_editable_install() {
  # Install BactiPipe in editable mode into the 'bactipipe' env
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate bactipipe
  msg "Installing BactiPipe in editable mode (-e .)…"
  pip install -U pip
  pip install -e "${REPO_ROOT}"
  ok "Editable install complete."
}

install_basespace_cli() {
  msg "Installing BaseSpace CLI (bs)…"
  mkdir -p "${HOME}/bin"
  if [[ "${PLATFORM}" == "Linux" ]]; then
    wget -q "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O "${HOME}/bin/bs"
  else
    wget -q "https://launch.basespace.illumina.com/CLI/latest/amd64-osx/bs" -O "${HOME}/bin/bs"
  fi
  chmod +x "${HOME}/bin/bs"
  ok "BaseSpace CLI installed to ${HOME}/bin/bs (ensure \$HOME/bin is on PATH)."
}

post_install_checks() {
  source "$(conda info --base)/etc/profile.d/conda.sh"

  msg "Verifying core CLI…"
  conda activate bactipipe
  if command -v bactipipe >/dev/null 2>&1; then ok "bactipipe CLI found."; else err "bactipipe CLI missing"; exit 1; fi

  msg "Checking representative tools…"
  set +e
  conda run -n bactipipe fastp -h >/dev/null 2>&1 && ok "fastp ok" || err "fastp missing"
  conda run -n bactipipe flye --version 2>/dev/null && ok "flye ok" || err "flye missing"
  conda run -n genepid abricate --version 2>/dev/null && ok "abricate ok" || err "abricate missing"
  conda run -n viramr amrfinder -h >/dev/null 2>&1 && ok "amrfinderplus ok" || err "amrfinderplus missing"
  set -e

  ok "Sanity checks complete."
}

# ---------- main ----------
detect_platform
ensure_conda
ensure_mamba

msg "Creating environments…"
create_env "${ENV_DIR}/bactipipe.yml"
create_env "${ENV_DIR}/genepid.yml"
create_env "${ENV_DIR}/viramr.yml"

pip_editable_install
install_basespace_cli
post_install_checks

echo
ok "All set! Common usage:"
echo "  - Activate core:    conda activate bactipipe"
echo "  - Run external tools without switching envs:"
echo "        conda run -n genepid abricate --help"
echo "        conda run -n viramr amrfinder -h"
