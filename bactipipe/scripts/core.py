import os

def env_cmd(tool: str) -> list:
    """
    Return 'conda run -n <env> <tool>' for each external tool,
    with optional overrides via env vars:
      VIRAMR_ENV  (default: viramr)      -> amrfinder, virulencefinder
      ABRICATE_ENV (default: abricate)   -> abricate
    """
    viramr_env = os.environ.get("VIRAMR_ENV", "viramr")
    abricate_env = os.environ.get("ABRICATE_ENV", "abricate")

    t = tool.lower()
    if t.startswith("amrfinder"):
        return ["conda", "run", "-n", viramr_env, "amrfinder"]
    if t.startswith("virulencefinder"):
        # prefer script name used in bioconda; adjust if your install uses 'virulencefinder'
        return ["conda", "run", "-n", viramr_env, "virulencefinder.py"]
    if t.startswith("abricate"):
        return ["conda", "run", "-n", abricate_env, "abricate"]
    return [tool]

def set_tmp_env(tmp_dir: str) -> dict:
    """Return a copy of os.environ with TMPDIR/TEMP/TMP set to tmp_dir."""
    env = os.environ.copy()
    env["TMPDIR"] = tmp_dir
    env["TEMP"]  = tmp_dir
    env["TMP"]   = tmp_dir
    return env
