import subprocess
import os
import sys
import json
from datetime import datetime
from bactipipe.scripts.downloader import get_kmerfinder_latest_version_ftp, read_local_kmerfinder_version, read_local_version
def get_timestamp():
    return datetime.now().strftime("%Y%m%d_%H%M%S")

def prompt_user(msg):
    response = input(f"{msg} [y/N]: ").lower()
    return response == "y"

def check_outdated_pip():
    print("\n🔍 Checking for outdated Python packages...")

    result = subprocess.run(
        ["pip", "list", "--outdated", "--format=json"],
        capture_output=True,
        text=True
    )

    try:
        packages = json.loads(result.stdout)
    except json.JSONDecodeError:
        print("⚠️ Could not parse pip output.")
        print("⚠️ Proceeding with package updates despite pip upgrade failure.")

    if not packages:
        print("\n✅ All required Python packages are up to date.\n")
        return

    print("🔧 Outdated Python packages:")
    for pkg in packages:
        name = pkg["name"]
        current = pkg["version"]
        latest = pkg["latest_version"]
        print(f" - {name} ({current} → {latest})")

def check_outdated_conda():
    print("\n🔍 Checking for outdated conda packages...\n")

    result = subprocess.run(
        ["conda", "update", "--all", "--dry-run", "--json"],
        capture_output=True,
        text=True
    )

    try:
        info = json.loads(result.stdout)
    except json.JSONDecodeError:
        print("⚠️ Could not parse conda output.")
        return

    updates = []
    for pkg in info.get("actions", {}).get("LINK", []):
        name = pkg["name"]
        version = pkg["version"]
        channel = pkg["channel"]
        updates.append(f"{name} → {version} ({channel})")

    if not updates:
        print("✅ All required conda packages are up to date.\n")
        return

    print("🔧 Outdated conda packages:")
    for update in updates:
        print(f" - {update}")

def update_pip_packages():
    print("\n🔧 Checking for outdated Python packages...")
    max_retries = 3
    for attempt in range(max_retries):
        result = subprocess.run(["pip", "install", "--upgrade", "pip"], capture_output=True, text=True)
        if result.returncode == 0:
            break
        print(f"⚠️ Attempt {attempt + 1} to upgrade pip failed: {result.stderr}")
    else:
        print("❌ All attempts to upgrade pip have failed. Please try upgrading pip manually.")

        print("⚠️ Proceeding with package updates anyway.")
    
    output = subprocess.check_output(["pip", "list", "--outdated", "--format=json"]).decode()
    packages = json.loads(output)

    if not packages:
        print("\n✅ All required Python packages are up to date.")
        return

    print("\n🔧 Outdated Python packages:")
    for pkg in packages:
        name = pkg["name"]
        current = pkg["version"]
        latest = pkg["latest_version"]
        print(f" - {name} ({current} → {latest})")

    if prompt_user("\nProceed to update all listed Python packages?"):
        for pkg in packages:
            name = pkg["name"]
            subprocess.run(["pip", "install", "--upgrade", name])
        print("✅ All selected Python packages have been updated.")
    else:
        print("⏭️ Skipped pip package updates.")


def update_conda_packages():
    print("\n🔧 Checking for outdated conda packages...")

    # Run dry-run with JSON output
    result = subprocess.run(
        ["conda", "update", "--all", "--dry-run", "--json"],
        capture_output=True,
        text=True
    )

    try:
        info = json.loads(result.stdout)
    except json.JSONDecodeError:
        print("⚠️ Could not parse conda output.")
        return

    # Check what packages will be updated
    updates = []
    for pkg in info.get("actions", {}).get("LINK", []):
        name = pkg["name"]
        version = pkg["version"]
        channel = pkg["channel"]
        updates.append(f"{name} → {version} ({channel})")

    if not updates:
        print("\n✅ All required conda packages are up to date.\n")
        return

    print("\n🔧 Outdated conda packages:")
    for line in updates:
        print(f" - {line}")

    if prompt_user("\nProceed to update all listed conda packages?"):
        subprocess.run(["conda", "update", "--all", "--yes"])
        print("✅ All selected conda packages have been updated.")
    else:
        print("⏭️ Skipped conda package updates.")

def check_database_versions():
    print("\n🧬 Checking KmerFinder database version...")

    local_kf_version = read_local_kmerfinder_version()
    remote_kf_version = get_kmerfinder_latest_version_ftp()

    if not remote_kf_version:
        print("⚠️ Could not determine the latest remote version.")
    else:
        print(f"📂 Local version:  {local_kf_version or 'unknown'}")
        print(f"🌐 Remote version: {remote_kf_version}")

        if local_kf_version != remote_kf_version:
            print("🔔 KmerFinder database is outdated.")
        else:
            print("✅ KmerFinder database is up to date.")

    print("\n🧬 Checking CheckM database version...")

    CHECKM_DB_LATEST_VERSION = "checkm_data_2015_01_16"
    checkm_db = CHECKM_DB_LATEST_VERSION
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
    checkm_db_dir = os.path.join(base_dir, checkm_db)
    local_checkm_version = read_local_version(checkm_db_dir)

    print(f"📂 Local version:  {local_checkm_version or 'unknown'}")
    print(f"🌐 Remote version: {CHECKM_DB_LATEST_VERSION}")

    if local_checkm_version != CHECKM_DB_LATEST_VERSION:
        print("🔔 CheckM database is outdated.")
    else:
        print("✅ CheckM database is up to date.\n")
    