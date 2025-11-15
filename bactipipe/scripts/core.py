import os
import sys
import shlex
from shutil import which
import subprocess

def env_cmd(tool: str) -> list:
    """
    Return 'conda run -n <env> <tool>' for each external tool,
    with optional overrides via env vars:
      VIRAMR_ENV  (default: viramr)      -> amrfinder, virulencefinder
      ABRICATE_ENV (default: abricate)   -> abricate
    """
    viramr_env = os.environ.get("VIRAMR_ENV", "viramr")
    abricate_env = os.environ.get("ABRICATE_ENV", "genepid")

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

############CLI TOOL MANAGEMENT HELPERS ############

# DB root: 
DB_ROOT = os.environ.get("BACTIPIPE_DB_DIR", os.path.expanduser("~/src/database/BactiPipe"))

DB_SPECS = {
    "cgmlstfinder": {
        "label": "cgmlstfinder_db",
        "path": lambda: os.path.join(DB_ROOT, "cgmlstfinder_db"),
        "type": "git+INSTALL.py",
    },
    "kmerfinder": {
        "label": "kmerfinder_db",
        "path": lambda: os.path.join(DB_ROOT, "kmerfinder_db"),
        "type": "git+INSTALL.sh",
    },
    "mlst": {
        "label": "mlst_db",
        "path": lambda: os.path.join(DB_ROOT, "mlst_db"),
        "type": "git-only",
    },
    "serotypefinder": {
        "label": "serotypefinder_db",
        "path": lambda: os.path.join(DB_ROOT, "serotypefinder_db"),
        "type": "git-only",
    },
    "virulencefinder": {
        "label": "virulencefinder_db",
        "path": lambda: os.path.join(DB_ROOT, "virulencefinder_db"),
        "type": "git-only",
    },
    "amrfinder": {
        "label": "amrfinder",
        "path": lambda: os.path.expanduser("~/.amrfinderplus"),  # default DB location
        "type": "tool-managed",
    },
    "abricate": {
        "label": "abricate",
        "path": lambda: None,  # abricate has multiple DB dirs; query via CLI
        "type": "tool-managed",
    },
}

def _run(cmd, cwd=None):
    """Run a shell command with printing and error propagation."""
    if isinstance(cmd, str):
        printable = cmd
        cmd = shlex.split(cmd)
    else:
        printable = " ".join(cmd)
    print(f"[bactipipe] $ {printable}", flush=True)
    subprocess.run(cmd, check=True, cwd=cwd)
def list_databases():
    """List known databases and whether they appear to be installed."""
    print(f"Database root: {DB_ROOT}\n")

    lines = []
    for key, spec in DB_SPECS.items():
        label = spec["label"]
        db_type = spec["type"]
        path = spec["path"]()
        status = "unknown"

        if key in {"cgmlstfinder", "kmerfinder", "mlst", "serotypefinder", "virulencefinder"}:
            if path and os.path.isdir(path) and os.listdir(path):
                status = "installed"
            else:
                status = "missing"
        elif key == "amrfinder":
            # consider DB installed if default DB dir exists
            if path and os.path.isdir(path) and os.listdir(path):
                status = "installed"
            else:
                status = "missing"
        elif key == "abricate":
            try:
                result = subprocess.run(
                    ["conda", "run", "-n", "genepid", "abricate", "--list"],
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                # more than header line → DBs configured
                num_lines = len([ln for ln in result.stdout.splitlines() if ln.strip()])
                status = "installed" if num_lines > 1 else "missing"
            except Exception:
                status = "missing"

        location = path or "(managed by tool)"
        lines.append((label, key, db_type, status, location))

    # pretty-ish table
    colnames = ("Label", "Key", "Type", "Status", "Location")
    print("{:<18} {:<14} {:<15} {:<10} {}".format(*colnames))
    print("-" * 80)
    for label, key, db_type, status, location in lines:
        print("{:<18} {:<14} {:<15} {:<10} {}".format(label, key, db_type, status, location))
    
def _parse_db_names(args):
    if not args:
        raise SystemExit("ERROR: 'update-db' requires at least one database name or 'all'.")

    raw = " ".join(args)
    tokens = []
    for part in raw.replace(",", " ").split():
        tokens.append(part.strip())

    if not tokens:
        raise SystemExit("ERROR: No valid database names provided.")

    if any(t.lower() == "all" for t in tokens):
        return list(DB_SPECS.keys())

    # Allow either the "key" or "<label>"/"<label>_db" forms
    normalized = []
    for t in tokens:
        t_low = t.lower().rstrip("/")
        t_low = t_low.replace("_db", "")  # virulencefinder_db → virulencefinder
        if t_low in DB_SPECS:
            normalized.append(t_low)
        else:
            # try matching by label
            for key, spec in DB_SPECS.items():
                if t_low == spec["label"].lower().replace("_db", ""):
                    normalized.append(key)
                    break
            else:
                raise SystemExit(f"ERROR: Unknown database name: {t}")

    # preserve order but drop duplicates
    seen = set()
    result = []
    for n in normalized:
        if n not in seen:
            seen.add(n)
            result.append(n)
    return result

def _ensure_git():
    if which("git") is None:
        raise SystemExit("ERROR: 'git' not found on PATH. Please install git and retry.")


def _update_cgmlstfinder():
    _ensure_git()
    db_dir = DB_SPECS["cgmlstfinder"]["path"]()
    os.makedirs(DB_ROOT, exist_ok=True)
    if os.path.isdir(os.path.join(db_dir, ".git")):
        print(f"[bactipipe] Updating cgmlstfinder_db in {db_dir}")
        _run(["git", "pull", "--ff-only"], cwd=db_dir)
    else:
        print(f"[bactipipe] Cloning cgmlstfinder_db into {db_dir}")
        _run(["git", "clone", "https://bitbucket.org/genomicepidemiology/cgmlstfinder_db.git", db_dir])

    # Install databases via INSTALL.py in genepid env
    _run(["conda", "run", "-n", "genepid", "python", "INSTALL.py"], cwd=db_dir)


def _update_kmerfinder():
    _ensure_git()
    db_dir = DB_SPECS["kmerfinder"]["path"]()
    os.makedirs(DB_ROOT, exist_ok=True)
    if os.path.isdir(os.path.join(db_dir, ".git")):
        print(f"[bactipipe] Updating kmerfinder_db in {db_dir}")
        _run(["git", "pull", "--ff-only"], cwd=db_dir)
    else:
        print(f"[bactipipe] Cloning kmerfinder_db into {db_dir}")
        _run(["git", "clone", "https://bitbucket.org/genomicepidemiology/kmerfinder_db.git", db_dir])

    # Install only bacteria DB from latest update
    _run(["bash", "INSTALL.sh", db_dir, "bacteria", "latest"], cwd=db_dir)


def _update_simple_git_db(key, url):
    _ensure_git()
    db_dir = DB_SPECS[key]["path"]()
    os.makedirs(DB_ROOT, exist_ok=True)
    if os.path.isdir(os.path.join(db_dir, ".git")):
        print(f"[bactipipe] Updating {DB_SPECS[key]['label']} in {db_dir}")
        _run(["git", "pull", "--ff-only"], cwd=db_dir)
    else:
        print(f"[bactipipe] Cloning {DB_SPECS[key]['label']} into {db_dir}")
        _run(["git", "clone", url, db_dir])

def update_databases(args):
    dbs = _parse_db_names(args)
    print(f"Database root: {DB_ROOT}")
    os.makedirs(DB_ROOT, exist_ok=True)

    for db in dbs:
        print(f"\n=== Updating database: {db} ===")
        if db == "cgmlstfinder":
            _update_cgmlstfinder()
        elif db == "kmerfinder":
            _update_kmerfinder()
        elif db == "mlst":
            _update_simple_git_db("mlst", "https://bitbucket.org/genomicepidemiology/mlst_db.git")
        elif db == "serotypefinder":
            _update_simple_git_db("serotypefinder", "https://bitbucket.org/genomicepidemiology/serotypefinder_db.git")
        elif db == "virulencefinder":
            _update_simple_git_db("virulencefinder", "https://bitbucket.org/genomicepidemiology/virulencefinder_db.git")
        elif db == "amrfinder":
            print("[bactipipe] Updating AMRFinderPlus database (viramr env)…")
            _run(["conda", "run", "-n", "viramr", "amrfinder", "-u"])
        elif db == "abricate":
            print("[bactipipe] Updating ABRicate databases (genepid env)…")
            _run(["conda", "run", "-n", "genepid", "abricate", "--setupdb"])
        else:
            print(f"[bactipipe] WARNING: unknown db key {db}, skipping.")

    print("\n✅ Database update complete.")

def check_updates():
    """Check for environment and database updates without installing them."""

    print("=== Checking Conda environment updates ===\n")

    for env in ENV_NAMES:
        try:
            result = subprocess.run(
                ["conda", "search", "--outdated", "--name", env],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            lines = [ln for ln in result.stdout.splitlines() if ln.strip()]
            if len(lines) > 1:  # header + entries
                print(f"[{env}] {len(lines)-1} updates available")
            else:
                print(f"[{env}] up-to-date")
        except Exception:
            print(f"[{env}] unable to check (error)")

    print("\n=== Checking databases ===\n")

    for key, spec in DB_SPECS.items():
        label = spec["label"]
        db_type = spec["type"]
        path = spec["path"]()

        if db_type == "git-only" or db_type.startswith("git"):
            if path and os.path.isdir(os.path.join(path, ".git")):
                try:
                    result = subprocess.run(
                        ["git", "fetch"], cwd=path, stdout=subprocess.PIPE
                    )
                    local = subprocess.check_output(
                        ["git", "rev-parse", "HEAD"], cwd=path, text=True
                    ).strip()
                    remote = subprocess.check_output(
                        ["git", "rev-parse", "@{u}"], cwd=path, text=True
                    ).strip()
                    if local != remote:
                        print(f"{label:18} update available")
                    else:
                        print(f"{label:18} up-to-date")
                except Exception:
                    print(f"{label:18} cannot determine (error)")
            else:
                print(f"{label:18} not installed")

        elif key == "amrfinder":
            print(f"{label:18} run: 'bactipipe update-db amrfinder' to update")

        elif key == "abricate":
            print(f"{label:18} run: 'bactipipe update-db abricate' to update")

    print("\nDone.\n")






