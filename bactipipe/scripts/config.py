import os
import json
from datetime import datetime, timedelta

CONFIG_DIR = os.path.expanduser("~/.bactipipe")
CONFIG_PATH = os.path.join(CONFIG_DIR, "config.json")

# Fallback to the “data” folder inside the package
PACKAGE_DIR      = os.path.dirname(__file__)            # .../bactipipe/scripts
DEFAULT_DATA_DIR = os.path.normpath(os.path.join(PACKAGE_DIR, "..", "data"))

# In Docker it will be /opt/bactipipe/data)
# otherwise fall back to the shipped data folder.
DATA_DIR = os.environ.get("BACTIPIPE_DB_DIR", DEFAULT_DATA_DIR)

# CHECKM_DATA_DIR      = os.path.join(DATA_DIR, "checkm_data")
# KMERFINDER_DB_PARENT = os.path.join(DATA_DIR, "kmerfinder_db")


DEFAULT_CONFIG = {
    "last_update_check": None,
    "update_check_interval_days": 60
}

def load_config():
    if not os.path.exists(CONFIG_PATH):
        return DEFAULT_CONFIG.copy()
    with open(CONFIG_PATH, "r") as f:
        return json.load(f)

def save_config(config):
    os.makedirs(CONFIG_DIR, exist_ok=True)
    with open(CONFIG_PATH, "w") as f:
        json.dump(config, f, indent=2)

def should_check_for_updates():
    config = load_config()
    last_check = config.get("last_update_check")
    interval_days = config.get("update_check_interval_days", 60)

    if last_check:
        last_check_time = datetime.strptime(last_check, "%Y-%m-%d %H:%M:%S")
        days_since_last_check = (datetime.now() - last_check_time).days
        if days_since_last_check < interval_days:
            return False
    if last_check:
        print(f"\nThe program is designed to automatically check for updates to the pipeline packages and databases every 60 days. Last check: {last_check} ({days_since_last_check} days)\n")
    else: 
        print("\nThe program must first configure automatic updates, which will occur every 60 days.\n")
    print("\n🕒 Auto-checking for updates...")
    return True

def update_last_checked():
    config = load_config()
    config["last_update_check"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    save_config(config)
