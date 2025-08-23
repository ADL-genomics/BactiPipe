# scripts/downloader.py
import os
import subprocess
import shutil
import urllib.request
import tarfile
import sys
from datetime import datetime
from ftplib import FTP
from .config import DATA_DIR


CHECKM_DB_LATEST_VERSION = "checkm_data_2015_01_16"
checkm_db = f"{CHECKM_DB_LATEST_VERSION}"
def get_timestamp():
    return datetime.now().strftime("%Y%m%d_%H%M%S")

def get_version_file(db_dir):
    return os.path.join(db_dir, "VERSION.txt")

def read_local_version(db_dir):
    version_file = get_version_file(db_dir)
    if os.path.exists(version_file):
        with open(version_file, "r") as f:
            return f.read().strip()
    return None

def write_version_file(db_dir, version_str):
    with open(get_version_file(db_dir), "w") as f:
        f.write(version_str)

def prompt_user(msg):
    response = input(f"{msg} [y/N]: ").lower()
    return response == "y"

def backup_directory(path):
    backup_path = f"{path}_backup_{get_timestamp()}"
    print(f"📦 Backing up current database to: {backup_path}")
    shutil.move(path, backup_path)
    return backup_path

def run_command(command, cwd=None):
    print(f"Running: {' '.join(command)}")
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        print("Error:", result.stderr)
        raise RuntimeError(f"Command failed: {' '.join(command)}")
    else:
        print(result.stdout)

def download_with_progress(url, dest):
    def show_progress(block_num, block_size, total_size):
        downloaded = block_num * block_size
        percent = int(downloaded * 100 / total_size) if total_size > 0 else 0
        bar = f"[{'=' * (percent // 2)}{' ' * (50 - percent // 2)}] {percent}%"
        sys.stdout.write(f"\rDownloading {os.path.basename(dest)} {bar}")
        sys.stdout.flush()

    print(f"Downloading {url} ...")
    urllib.request.urlretrieve(url, dest, reporthook=show_progress)
    print("\n✅ Download complete.")

def get_kmerfinder_latest_version_ftp():
    try:
        ftp = FTP("ftp.cbs.dtu.dk")
        ftp.login()
        ftp.cwd("/public/CGE/databases/KmerFinder/version")

        entries = []
        ftp.retrlines('LIST', entries.append)
        ftp.quit()

        for line in entries:
            if 'latest ->' in line:
                parts = line.split('->')
                if len(parts) == 2:
                    return parts[1].strip()
    except Exception as e:
        print(f"⚠️ Could not connect to KmerFinder FTP: {e}")

    return None

def read_local_kmerfinder_version():
    version_file = os.path.join(DATA_DIR,"kmerfinder_db", "bacteria", "VERSION.txt")
    if os.path.exists(version_file):
        with open(version_file) as f:
            return f.read().strip()
    return None

def setup_kmerfinder_database(dest_dir, remote_version=None, notifications = "on"):
    repo_url = "https://bitbucket.org/genomicepidemiology/kmerfinder_db.git"
    repo_path = os.path.join(dest_dir, "tmp_kmerfinder_repo")
    final_db_path = os.path.join(dest_dir, "bacteria")
    if notifications == "on":
        print("\nSetting up KmerFinder database...\n")

    if os.path.exists(os.path.join(final_db_path, "bacteria.ATG.seq.b")):
        if notifications == "on":
            print("✅ KmerFinder database already installed. Skipping.")
        return

    os.makedirs(dest_dir, exist_ok=True)

    print("Cloning KmerFinder database repository...")
    if not os.path.exists(repo_path):
        run_command(["git", "clone", repo_url, repo_path])

    print("Installing KmerFinder database...")
    print("This may take a while.The database is > 30 GB.\n")
    run_command(["bash", "INSTALL.sh", dest_dir, "bacteria", "latest"], cwd=repo_path)

    shutil.rmtree(repo_path)
    print("✅ KmerFinder database installed successfully.")
    write_version_file(final_db_path, get_kmerfinder_latest_version_ftp())

def setup_checkm_database(dest_dir, checkm_db="checkm_data_2015_01_16", notifications = "on"):
    url = f"https://data.ace.uq.edu.au/public/CheckM_databases/{checkm_db}.tar.gz"
    filename = os.path.join(dest_dir, os.path.basename(url))
    extracted_dir = os.path.join(dest_dir, checkm_db)
    if notifications == "on":
        print("\nSetting up CheckM database...\n")

    if os.path.exists(os.path.join(extracted_dir, "hmms/checkm.hmm")):
        if notifications == "on":
            print("✅ CheckM database already exists. Skipping.")
        return

    os.makedirs(dest_dir, exist_ok=True)

    if not os.path.exists(filename):
        download_with_progress(url, filename)

    print("Extracting CheckM database...")
    with tarfile.open(filename, "r:gz") as tar:
        tar.extractall(path=extracted_dir)

    print("✅ CheckM database installed successfully.\n")
    write_version_file(extracted_dir, CHECKM_DB_LATEST_VERSION)
    os.remove(filename)
    
def setup_databases(notifications="on"):
    base_dir = DATA_DIR
    remote_version = get_kmerfinder_latest_version_ftp()
    setup_checkm_database(base_dir, notifications=notifications)
    setup_kmerfinder_database(os.path.join(base_dir, "kmerfinder_db"), remote_version, notifications=notifications)
    
def update_databases():
    print("\n🔄 Checking for new versions of databases...\n")
    base_dir = DATA_DIR

    # --- CheckM ---
    checkm_db_dir = os.path.join(base_dir, checkm_db)
    checkm_hmm_path = os.path.join(checkm_db_dir, "hmms", "checkm.hmm")

    if os.path.exists(checkm_hmm_path):
        local_version = read_local_version(checkm_db_dir)
        print(f"📂 CheckM db local version: {local_version or 'unknown'}")
        print(f"🌐 CheckM db remote version: {CHECKM_DB_LATEST_VERSION}")

        if local_version == CHECKM_DB_LATEST_VERSION:
            print("✅ CheckM database is up to date.")
        else:
            if prompt_user("Do you want to update the CheckM database?"):
                if os.path.exists(checkm_db_dir):
                    backup_directory(checkm_db_dir)
                setup_checkm_database(base_dir)
    else:
        print("📁 CheckM DB not found — installing fresh.")
        setup_checkm_database(base_dir)

    # --- KmerFinder ---
    kf_dir = os.path.join(base_dir, "kmerfinder_db")
    kf_bacteria_dir = os.path.join(kf_dir, "bacteria")
    remote_version = get_kmerfinder_latest_version_ftp()
    local_version = read_local_version(kf_bacteria_dir)

    if remote_version:
        print(f"\n📂 Local KmerFinder db version: {local_version or 'unknown'}")
        print(f"🌐 Remote KmerFinder db version: {remote_version}")
        
        if local_version is None:
            print("⚠️ Local KmerFinder database not found.")
            setup_kmerfinder_database(kf_dir, remote_version)
        elif remote_version != local_version:
            if prompt_user("Update KmerFinder database to the latest version?"):
                if os.path.exists(kf_dir):
                    backup_directory(kf_dir)
                setup_kmerfinder_database(kf_dir, remote_version)
        else:
            print("✅ KmerFinder database is up to date.")
    else:
        print("⚠️ Could not retrieve remote KmerFinder version.")
        # Optional: install from scratch if not already present
        if not os.path.exists(kf_bacteria_dir):
            print("📁 KmerFinder DB not found — installing fresh.")
            setup_kmerfinder_database(kf_dir)

    print("\n✅ All selected databases updated.\n")

    