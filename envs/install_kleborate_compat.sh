#!/usr/bin/env bash
set -euo pipefail

environment_name="${1:-kleborate}"

# Kleborate 3.2.4 imports ``kaptive.database``. Kaptive 3.3 removed that
# module, so replace the solver-selected release with the last compatible
# Kaptive 3.2 release without changing the rest of this isolated environment.
conda install -y -n "$environment_name" --no-deps \
  https://conda.anaconda.org/bioconda/noarch/kaptive-3.2.2-pyhdfd78af_0.conda
