#!/usr/bin/env python3

import sys
import subprocess
import os
from bactipipe.scripts.downloader import setup_databases, update_databases
from bactipipe.scripts import updater, config

def main():
    # if config.should_check_for_updates():
    #     updater.check_outdated_pip()
    #     updater.check_outdated_conda()
    #     updater.check_database_versions()
    #     config.update_last_checked() 
    
    usage = '''\nUsage: bactipipe.py <command> [options]
    \nPackage/Database management commands:
      Check-updates    :  Check for updates to pipeline packages.
      Update-packages  :  Update all Python/Conda packages.
      Update-databases :  Re-download the latest databases.
    \nPipeline commands:
      qc-illumina      :  Run Illumina QC pipeline.
      qc-nanopore      :  Run Nanopore QC pipeline.
      relate           :  Type assembled genomes and compute relatedness.
      detect           :  Detect virulence and antimicrobial resistance genes in assemblies.

    \nUse 'bactipipe <command> -h' to see help for the selected pipeline. [Only for pipeline commands.]\n'''

    if len(sys.argv) < 2 or sys.argv[1] in ["-h", "--help"]:
        print(usage)
        sys.exit(1)
    
    command = sys.argv[1].lower()
    
    tool_manage_commands = ["check-updates", "update-packages", "update-databases"]

    # Map run commands to the corresponding script
    script_dir = os.path.join(os.path.dirname(__file__), "scripts")
    script_map = {
        "qc-illumina": os.path.join(script_dir, "run_illumina_qc.py"),
        "qc-nanopore": os.path.join(script_dir, "run_nanopore_qc.py"),
        "relate": os.path.join(script_dir, "type_genomes.py"),
        "detect": os.path.join(script_dir, "run_traits.py"),
    }

    if command not in script_map and command not in tool_manage_commands:
        print(f"\nError: Invalid command '{command}'.")
        print(usage)
        sys.exit(1)

    elif command in tool_manage_commands:
        if command == "check-updates":
            updater.check_outdated_pip()
            updater.check_outdated_conda()
            updater.check_database_versions()
            config.update_last_checked()
        elif command == "update-packages":
            updater.update_pip_packages()
            updater.update_conda_packages()
        elif command == "update-databases":
            update_databases()
        sys.exit(0)

    elif command in script_map:
        # Ensure databases are present before running any pipeline
        setup_databases(notifications="off")
        script_to_run = script_map[command]

    # Forward all remaining arguments to the appropriate script
    run_command = ["python3", script_to_run] + sys.argv[2:]

    try:
        subprocess.run(run_command, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(e.returncode)

if __name__ == "__main__":
    main()
