#!/usr/bin/env python3
"""
Capture reproducibility metadata for pipeline runs.

This script generates files containing:
- Command-line invocation
- Snakemake version
- Python version
- OS information
- Dependency snapshots (conda env export and pip freeze)
"""

import json
import os
import platform
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def get_snakemake_version():
    """Get Snakemake version."""
    try:
        result = subprocess.run(
            ["snakemake", "--version"],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def get_python_version():
    """Get Python version."""
    return sys.version


def get_os_info():
    """Get OS information."""
    return {
        "system": platform.system(),
        "release": platform.release(),
        "version": platform.version(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "platform": platform.platform(),
    }


def get_conda_env_export():
    """Get conda environment export."""
    try:
        # Try to get the active conda environment
        conda_env = os.environ.get("CONDA_DEFAULT_ENV", None)
        if conda_env:
            result = subprocess.run(
                ["conda", "env", "export", "--name", conda_env, "--no-builds"],
                capture_output=True,
                text=True,
                check=True
            )
            return result.stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    # Fallback: try to export current environment
    try:
        result = subprocess.run(
            ["conda", "env", "export", "--no-builds"],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def get_pip_freeze():
    """Get pip freeze output."""
    try:
        result = subprocess.run(
            [sys.executable, "-m", "pip", "freeze"],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def get_command_line_invocation():
    """Get command-line invocation from environment or arguments."""
    # Snakemake stores the command line in SNAKEMAKE_COMMANDLINE if available
    cmd = os.environ.get("SNAKEMAKE_COMMANDLINE", None)
    if cmd:
        return cmd

    # Fallback: try to reconstruct from sys.argv
    # This won't be perfect but better than nothing
    if len(sys.argv) > 1:
        return " ".join(sys.argv[1:])

    return "unknown"


def main():
    """Generate reproducibility files."""
    if len(sys.argv) < 2:
        print("Usage: capture_reproducibility.py <output_dir>", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(sys.argv[1])
    output_dir.mkdir(parents=True, exist_ok=True)

    # Collect all metadata
    metadata = {
        "timestamp": datetime.now().isoformat(),
        "command_line": get_command_line_invocation(),
        "snakemake_version": get_snakemake_version(),
        "python_version": get_python_version(),
        "os_info": get_os_info(),
        "conda_env": os.environ.get("CONDA_DEFAULT_ENV", None),
        "python_executable": sys.executable,
    }

    # Write metadata JSON
    metadata_file = output_dir / "metadata.json"
    with open(metadata_file, "w") as f:
        json.dump(metadata, f, indent=2)

    # Write human-readable metadata
    metadata_txt = output_dir / "metadata.txt"
    with open(metadata_txt, "w") as f:
        f.write("FUSILLI Pipeline Reproducibility Metadata\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Timestamp: {metadata['timestamp']}\n")
        f.write(f"\nCommand Line Invocation:\n")
        f.write(f"  {metadata['command_line']}\n")
        f.write(f"\nSnakemake Version:\n")
        f.write(f"  {metadata['snakemake_version']}\n")
        f.write(f"\nPython Version:\n")
        f.write(f"  {metadata['python_version']}\n")
        f.write(f"\nPython Executable:\n")
        f.write(f"  {metadata['python_executable']}\n")
        f.write(f"\nConda Environment:\n")
        f.write(f"  {metadata['conda_env'] or 'Not in conda environment'}\n")
        f.write(f"\nOS Information:\n")
        for key, value in metadata['os_info'].items():
            f.write(f"  {key}: {value}\n")

    # Write conda environment export (always create file, even if empty)
    conda_env_export = get_conda_env_export()
    conda_file = output_dir / "conda-env.yaml"
    with open(conda_file, "w") as f:
        if conda_env_export:
            f.write(conda_env_export)
        else:
            f.write("# Conda environment export not available\n")
            f.write("# This may occur if conda is not installed or not in use\n")

    # Write pip freeze (always create file, even if empty)
    pip_freeze = get_pip_freeze()
    pip_file = output_dir / "pip-freeze.txt"
    with open(pip_file, "w") as f:
        if pip_freeze:
            f.write(pip_freeze)
        else:
            f.write("# Pip freeze not available\n")
            f.write("# This may occur if pip is not installed or not in use\n")

    print(f"Reproducibility metadata written to {output_dir}")


if __name__ == "__main__":
    main()

