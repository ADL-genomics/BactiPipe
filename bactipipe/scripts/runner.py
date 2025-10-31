# bactipipe/scripts/runner.py
from __future__ import annotations
from pathlib import Path
from typing import List, Optional, Tuple
import logging
import shutil
import subprocess

def compose_tool_cmd(cmd: List[str], env_name: Optional[str]) -> List[str]:
    """Wrap command to run inside a conda/mamba/micromamba env if requested."""
    if env_name:
        if shutil.which("conda"):
            return ["conda", "run", "--no-capture-output", "-n", env_name] + cmd
        if shutil.which("micromamba"):
            return ["micromamba", "run", "-n", env_name] + cmd
        if shutil.which("mamba"):
            return ["mamba", "run", "-n", env_name] + cmd
    return cmd

def run(cmd, cwd, log, env_name=None):
    cmd_exec = compose_tool_cmd(cmd, env_name)
    # normalize everything to str for safe logging/subprocess
    if any(isinstance(x, bool) for x in cmd_exec):
        log.error(f"BUG: command contains boolean(s): {cmd_exec!r}")
    cmd_exec = [str(x) for x in cmd_exec]
    log.debug(f"RUN: {' '.join(cmd_exec)}")
    try:
        p = subprocess.run(
            cmd_exec,
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        return p.returncode, p.stdout
    except Exception as e:
        return 127, f"[runner] Failed to execute: {e}"

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)
