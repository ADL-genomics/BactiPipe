import os
import sys
import shlex
from shutil import which
import subprocess
import datetime as dt

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
#Package management
ENV_NAMES = ("bactipipe", "genepid", "viramr")

def _run(cmd, **kwargs):
    """Run a shell command with pretty printing and error bubbling."""
    if isinstance(cmd, str):
        printable = cmd
        cmd = shlex.split(cmd)
    else:
        printable = " ".join(cmd)
    print(f"[bactipipe:update] $ {printable}", flush=True)
    subprocess.run(cmd, check=True, **kwargs)

def _which(cmd):
    from shutil import which
    return which(cmd)

def update_packages():
    """Update Python + conda packages in all envs, then pip-upgrade outdated pip packages."""
    mgr = "mamba" if _which("mamba") else "conda"

    # 1) Update all conda-managed packages (includes python) per env
    for env in ENV_NAMES:
        print(f"\n=== Updating conda packages in env: {env} ===", flush=True)
        _run([mgr, "update", "-n", env, "--all", "-y"])

    # 2) Upgrade pip-managed packages that are outdated in each env
    for env in ENV_NAMES:
        print(f"\n=== Upgrading pip packages in env: {env} ===", flush=True)
        # list outdated pip packages (name==version pins are not upgraded unless compatible)
        try:
            # capture the list of outdated in 'pkg==ver' lines
            result = subprocess.run(
                ["conda", "run", "-n", env, "python", "-m", "pip", "list", "--outdated", "--format=freeze"],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            # pip may not exist in some envs (e.g., viramr if no pip installs) -> skip
            if "No module named pip" in e.stderr or "No module named pip" in e.stdout:
                print(f"[bactipipe:update] pip not present in env {env}; skipping pip upgrades.")
                continue
            # other errors bubble up
            raise

        lines = [ln.strip() for ln in result.stdout.splitlines() if ln.strip()]
        if not lines:
            print("[bactipipe:update] No outdated pip packages.")
            continue

        # Each line looks like: 'package==current_version'
        pkgs = [ln.split("==", 1)[0] for ln in lines if "==" in ln]
        for pkg in pkgs:
            _run(["conda", "run", "-n", env, "python", "-m", "pip", "install", "-U", pkg])

    print("\n✅ Update complete. Tip: if you maintain lockfiles, regenerate them now.", flush=True)


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
    # NOTE: tool-managed inside envs, no DB_ROOT path
    "amrfinder": {
        "label": "amrfinder",
        "path": lambda: None,
        "type": "tool-managed",
        "env": "viramr",
    },
    "abricate": {
        "label": "abricate",
        "path": lambda: None,
        "type": "tool-managed",
        "env": "genepid",
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

    rows = []
    abricate_dbs = []   # we'll fill this from `abricate --list`

    for key, spec in DB_SPECS.items():
        label = spec["label"]
        db_type = spec["type"]
        path = spec["path"]()
        status = "unknown"
        location = ""

        # --- CGE git-based DBs ---
        if db_type.startswith("git"):
            if path and os.path.isdir(path) and os.listdir(path):
                status = "installed"
            else:
                status = "missing"
            location = path

        # --- AMRFinder (tool-managed in viramr env) ---
        elif key == "amrfinder":
            env_name = spec.get("env", "viramr")
            try:
                # Check that amrfinder runs; version output may include DB info.
                result = subprocess.run(
                    ["conda", "run", "-n", env_name, "amrfinder", "-V"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
                status = "installed"
                db_dir_line = next(
                    (ln for ln in result.stdout.splitlines()
                     if "Database directory:" in ln),
                    None,
                )
                if db_dir_line:
                    location = db_dir_line.split(":", 1)[1].strip()
                else:
                    location = f"(inside {env_name} environment)"
            except Exception:
                status = "missing"
                location = f"(inside {env_name} environment)"

        # --- ABRicate (tool-managed in genepid env) ---
        elif key == "abricate":
            env_name = spec.get("env", "genepid")
            try:
                result = subprocess.run(
                    ["conda", "run", "-n", env_name, "abricate", "--list"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
                lines = [ln for ln in result.stdout.splitlines() if ln.strip()]
                status = "installed" if len(lines) > 1 else "missing"
                location = f"(inside {env_name} environment)"

                # Parse the table: header + rows like:
                # DATABASE    SEQUENCES   DBTYPE   DATE
                if len(lines) > 1:
                    header = lines[0].split()
                    for ln in lines[1:]:
                        parts = ln.split()
                        if len(parts) < 4:
                            continue
                        dbname = parts[0]
                        sequences = parts[1]
                        dbtype = parts[2]
                        date = " ".join(parts[3:])
                        abricate_dbs.append(
                            {
                                "name": dbname,
                                "sequences": sequences,
                                "dbtype": dbtype,
                                "date": date,
                            }
                        )
            except Exception:
                status = "missing"
                location = f"(inside {env_name} environment)"

        else:
            location = path or ""

        rows.append((label, key, db_type, status, location or "-"))

    # Main summary table
    headers = ("Label", "Key", "Type", "Status", "Location")
    print("{:<18} {:<14} {:<15} {:<10} {}".format(*headers))
    print("-" * 80)
    for label, key, db_type, status, location in rows:
        print("{:<18} {:<14} {:<15} {:<10} {}".format(label, key, db_type, status, location))

    # Extra: ABRicate databases
    if abricate_dbs:
        print("\nABRicate databases (env: genepid):")
        print("  {:<15} {:>10}   {:<6}   {}".format("DATABASE", "SEQUENCES", "TYPE", "DATE"))
        print("  " + "-" * 55)
        for db in abricate_dbs:
            print(
                "  {:<15} {:>10}   {:<6}   {}".format(
                    db["name"],
                    db["sequences"],
                    db["dbtype"],
                    db["date"],
                )
            )
    
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
    print(f"\nDatabase root: {DB_ROOT}")
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

    import datetime as dt

    print("\n=== Checking Conda environment updates ===\n")

    # -----------------------------------------
    # Check outdated conda packages per env
    # -----------------------------------------
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

    today = dt.date.today()
    stale_days = 183  # ~6 months

    for key, spec in DB_SPECS.items():
        label = spec["label"]
        db_type = spec["type"]
        path = spec["path"]()

        # -------------------------
        # CGE git-based DBs
        # -------------------------
        if db_type == "git-only" or db_type.startswith("git"):
            if path and os.path.isdir(os.path.join(path, ".git")):
                try:
                    subprocess.run(
                        ["git", "fetch"],
                        cwd=path,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        check=False,
                        text=True,
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

        # -------------------------
        # AMRFinderPlus DB (tool-managed in viramr env)
        # SAFE: show installed DB version, do NOT update
        # -------------------------
        elif key == "amrfinder":
            env_name = spec.get("env", "viramr")
            try:
                result = subprocess.run(
                    ["conda", "run", "-n", env_name, "amrfinder", "--database_version"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
                db_version = None
                db_dir = None
                for ln in result.stdout.splitlines():
                    if "Database version:" in ln:
                        db_version = ln.split(":", 1)[1].strip()
                    if "Database directory:" in ln:
                        db_dir = ln.split(":", 1)[1].strip().strip("'\"")

                if db_version:
                    msg = f"installed DB version {db_version}"
                    if db_dir:
                        msg += f" at {db_dir}"
                    msg += " (not checked against NCBI; run 'bactipipe update-db amrfinder' to refresh)"
                else:
                    msg = "installed (database version unknown)"

                print(f"{label:18} {msg}")
            except Exception as e:
                print(f"{label:18} unable to check (error: {e})")

        # -------------------------
        # ABRicate DBs (tool-managed in genepid env)
        # classify by age based on DATE column
        # -------------------------
        elif key == "abricate":
            env_name = spec.get("env", "genepid")
            try:
                result = subprocess.run(
                    ["conda", "run", "-n", env_name, "abricate", "--list"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
                lines = [ln for ln in result.stdout.splitlines() if ln.strip()]
                if len(lines) <= 1:
                    print(f"{label:18} no databases installed")
                    continue

                stale = 0
                total = 0
                oldest = None
                newest = None

                # Expected format:
                # DATABASE    SEQUENCES   DBTYPE   DATE
                for ln in lines[1:]:
                    parts = ln.split()
                    if len(parts) < 4:
                        continue
                    date_str = parts[3]
                    total += 1
                    try:
                        db_date = dt.datetime.strptime(date_str, "%Y-%b-%d").date()
                    except ValueError:
                        # skip unparsable dates
                        continue

                    if oldest is None or db_date < oldest:
                        oldest = db_date
                    if newest is None or db_date > newest:
                        newest = db_date

                    if (today - db_date).days > stale_days:
                        stale += 1

                if total == 0 or oldest is None:
                    print(f"{label:18} installed (dates unavailable)")
                else:
                    if stale == 0:
                        print(
                            f"{label:18} all {total} DBs downloaded within 6 months "
                            f"(newest: {newest.isoformat()})"
                        )
                    else:
                        print(
                            f"{label:18} {stale}/{total} DBs older than 6 months "
                            f"(oldest: {oldest.isoformat()})"
                        )
            except Exception as e:
                print(f"{label:18} unable to check (error: {e})")

        else:
            # For any future DB types not explicitly handled
            print(f"{label:18} check not implemented")

    print("\nDone.\n")







