# bactipipe/scripts/traits_db.py
"""
traits_db.py
------------
DB/tool discovery & version capture for the traits-detect pipeline.

✅ Behavior (as requested):
- No guessing of "common paths".
- Accept a single root path (--db-dir or $VIRULENCEFINDER_DB) that contains
  subdirectories for various databases.
- The VirulenceFinder database of interest is always:
      <db_root>/virulencefinder_db
"""

from __future__ import annotations

import os
import subprocess
from typing import Dict, Tuple, Optional, List
import json, shlex

from bactipipe.scripts.core import env_cmd


# -----------------------------------------------------------
# Generic subprocess helper
# -----------------------------------------------------------
def _run(cmd: List[str]) -> Tuple[str, str]:
    """Run command and return (stdout, stderr)."""
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return p.stdout.strip(), p.stderr.strip()


# -----------------------------------------------------------
# ABRicate DB listing
# -----------------------------------------------------------
def _abricate_list() -> Dict[str, dict]:
    """Call 'abricate --list' (through its conda env) and parse the table."""
    out, _ = _run(env_cmd("abricate") + ["--list"])
    dbs: Dict[str, dict] = {}
    if not out:
        return dbs
    for line in out.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 3:
            name, count, path = parts[0], parts[1], parts[-1]
            dbs[name] = {"entries": count, "path": path}
    return dbs


# -----------------------------------------------------------
# VirulenceFinder DB resolution
# -----------------------------------------------------------
def find_virulencefinder_db(db_root: Optional[str] = None) -> Optional[str]:
    """
    Resolve the path to the actual VirulenceFinder database directory.

    Parameters
    ----------
    db_root : str or None
        Root directory passed from CLI as --db-dir.
        If None, uses environment variable VIRULENCEFINDER_DB.

    Returns
    -------
    str or None
        Full path to <db_root>/virulencefinder_db if it exists, else None.

    Example
    -------
    If --db-dir ~/data/databases  →  returns ~/data/databases/virulencefinder_db
    """
    root = db_root or os.getenv("VIRULENCEFINDER_DB")
    if not root:
        return None

    vf_path = os.path.join(os.path.expanduser(root), "virulencefinder_db")
    if os.path.isdir(vf_path):
        return vf_path
    return None


# -----------------------------------------------------------
# Version / metadata collection
# -----------------------------------------------------------
def _run_bash(cmd: str, env: Optional[dict] = None) -> str:
    cp = subprocess.run(["bash", "-lc", cmd], env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if cp.returncode != 0:
        raise subprocess.CalledProcessError(cp.returncode, cmd, cp.stdout, cp.stderr)
    return cp.stdout.strip()

def _format_v(s: str) -> str:
    """Ensure version strings look like vX.Y.Z."""
    s = s.strip()
    if not s:
        return "v?"
    if s.startswith("v"):
        return s
    return f"v{s}"

def _detect_amrfinder_versions(viramr_env: str) -> Dict[str, str]:
    # amrfinder -V prints multiple lines; we parse software & db versions
    out = _run_bash(f"conda run -n {shlex.quote(viramr_env)} amrfinder -V || true")
    tool_v, db_v = "", ""
    for line in out.splitlines():
        line = line.strip()
        if line.startswith("Software version:"):
            tool_v = line.split(":", 1)[1].strip()
        elif line.startswith("Database version:"):
            db_v = line.split(":", 1)[1].strip()
    return {
        "tool": _format_v(tool_v) if tool_v else "v?",
        "database": db_v or "",
    }

def _detect_abricate_tool(abricate_env: str) -> str:
    out = _run_bash(f"conda run -n {shlex.quote(abricate_env)} abricate --version || true")
    # Examples: "abricate 1.0.1"  -> take 2nd token
    toks = out.strip().split()
    ver = toks[-1] if toks else ""
    return _format_v(ver) if ver else "v?"

def _detect_abricate_dbs(abricate_env: str) -> Dict[str, str]:
    """
    Parse 'abricate --list' to extract DB dates (last column).
    Return canonical keys to avoid colliding with CGE tool DBs:
      resfinder   -> 'abricate_resfinder'
      card        -> 'abricate_card'
      vfdb        -> 'abricate_vfdb'
    """
    out = _run_bash(f"conda run -n {shlex.quote(abricate_env)} abricate --list || true")
    db_dates: Dict[str, str] = {}
    for line in out.splitlines():
        line = line.strip()
        if not line or line.lower().startswith("database"):
            continue
        parts = line.split()
        db = parts[0].lower()
        date = parts[-1] if len(parts) >= 2 else ""
        if db == "resfinder":
            db_dates["abricate_resfinder"] = date
        elif db == "card":
            db_dates["abricate_card"] = date
        elif db == "vfdb":
            db_dates["abricate_vfdb"] = date
    return db_dates

def _detect_virulencefinder_version_via_python(viramr_env: str) -> str:
    """
    Try to read VirulenceFinder's version from the installed package.
    Fallback to v3.2.0 if not importable (per your install).
    """
    py = (
        "python - <<'PY'\n"
        "try:\n"
        "  import importlib\n"
        "  try:\n"
        "    import importlib.metadata as md\n"
        "    v = md.version('virulencefinder')\n"
        "  except Exception:\n"
        "    vf = importlib.import_module('virulencefinder')\n"
        "    v = getattr(vf, '__version__', '')\n"
        "  print(v)\n"
        "except Exception:\n"
        "  print('')\n"
        "PY"
    )
    out = _run_bash(f"conda run -n {shlex.quote(viramr_env)} {py} || true").strip()
    return _format_v(out) if out else "v3.2.0" # default fallback (manually set here)

def _detect_vfinder_db_version_from_file(db_root: Optional[str]) -> str:
    """
    VirulenceFinder DB version comes from a one-line file:
      <db_root>/virulencefinder_db/VERSION
    Returns the stripped line (e.g., '2025-07-16.1') or '' if missing.
    """
    if not db_root:
        return ""
    vf_dir = os.path.join(db_root, "virulencefinder_db")
    version_fp = os.path.join(vf_dir, "VERSION")
    try:
        with open(version_fp, "r", encoding="utf-8") as fh:
            line = fh.readline().strip()
            return line
    except Exception:
        return ""

def _detect_resfinder_version_via_python(viramr_env: str) -> str:
    """
    Try to read ResFinder's version from the installed package in the viramr env.
    Fallback to 'v?' if not found.
    """
    py = (
        "python - <<'PY'\n"
        "try:\n"
        "  import importlib\n"
        "  try:\n"
        "    import importlib.metadata as md\n"
        "    v = md.version('resfinder')\n"
        "  except Exception:\n"
        "    rf = importlib.import_module('resfinder')\n"
        "    v = getattr(rf, '__version__', '')\n"
        "  print(v)\n"
        "except Exception:\n"
        "  print('')\n"
        "PY"
    )
    out = _run_bash(f"conda run -n {shlex.quote(viramr_env)} {py} || true").strip()
    return _format_v(out) if out else "v?"
    
def _detect_resfinder_db_version_from_file(db_root: Optional[str]) -> str:
    """
    ResFinder DB version from a one-line VERSION file:
      <db_root>/resfinder_db/VERSION
    """
    if not db_root:
        return ""
    rf_dir = os.path.join(db_root, "resfinder_db")
    version_fp = os.path.join(rf_dir, "VERSION")
    try:
        with open(version_fp, "r", encoding="utf-8") as fh:
            return fh.readline().strip()
    except Exception:
        return ""

def ensure_and_collect_versions(
    *,
    want_amr: bool,
    want_vf: bool,
    want_resfinder: bool,
    want_card: bool,
    want_vfinder: bool,
    want_vfdb: bool,
    db_root: Optional[str] = None,
) -> Dict:
    """
    Collect tool & database 'versions' for header. Tools always 'vX.Y.Z'.
    DBs use date strings where available.
    Honors conda env separation via env vars (default names):
      VIRAMR_ENV (amrfinder & virulencefinder), ABRICATE_ENV (abricate)
    """
    viramr_env = os.getenv("VIRAMR_ENV", "viramr")
    abricate_env = os.getenv("ABRICATE_ENV", "abricate")
    tools, dbs = {}, {}

    # --- AMR track ---
    if want_amr:
        # AMRFinderPlus tool + DB
        try:
            amrinfo = _detect_amrfinder_versions(viramr_env)
            tools["AMRFinderPlus"] = amrinfo.get("tool", "v?")
            if amrinfo.get("database"):
                dbs["amrfinder_db"] = amrinfo["database"]
        except Exception:
            tools["AMRFinderPlus"] = "v?"

        # ResFinder (CGE Python tool) — detect if installed
        try:
            rf_ver = _detect_resfinder_version_via_python(viramr_env)
            if rf_ver != "v?":
                tools["ResFinder"] = rf_ver
        except Exception:
            pass

        # ResFinder DB from VERSION file (if your --db-dir has it)
        rf_db_ver = _detect_resfinder_db_version_from_file(db_root)
        if rf_db_ver:
            dbs["resfinder_db"] = rf_db_ver

        # ABRicate tool + DB dates (separately named to avoid collision)
        if want_resfinder or want_card or want_vfdb:
            try:
                tools["ABRicate"] = _detect_abricate_tool(abricate_env)
            except Exception:
                tools["ABRicate"] = tools.get("ABRicate", "v?")
            try:
                abdbs = _detect_abricate_dbs(abricate_env)
                dbs.update(abdbs)  # adds 'abricate_resfinder', 'abricate_card', 'abricate_vfdb'
            except Exception:
                pass

    # --- Virulence track ---
    if want_vf:
        if want_vfinder:
            try:
                tools["VirulenceFinder"] = _detect_virulencefinder_version_via_python(viramr_env)
            except Exception:
                tools["VirulenceFinder"] = "v3.2.0"
            vf_db_ver = _detect_vfinder_db_version_from_file(db_root)
            if vf_db_ver:
                dbs["virulencefinder_db"] = vf_db_ver

        if want_vfdb:
            # ABRicate DB date already captured above via _detect_abricate_dbs (key 'abricate_vfdb')
            pass

    return {"tools": tools, "databases": dbs, "db_root": db_root or ""}
