#!/usr/bin/env python3
"""
SLURM job status checker for Snakemake.

This script queries SLURM to determine the status of a job.
Returns one of: 'running', 'success', 'failed', 'canceled', 'pending'

Usage:
    python cluster-status.py <jobid>
"""

import subprocess
import sys
import time
from enum import Enum


class JobStatus(Enum):
    RUNNING = "running"
    SUCCESS = "success"
    FAILED = "failed"
    CANCELED = "canceled"
    PENDING = "pending"
    UNKNOWN = "unknown"


def get_job_status(job_id: str) -> JobStatus:
    """
    Query SLURM for job status using sacct.

    Falls back to squeue if sacct fails or returns incomplete data.
    """
    try:
        result = subprocess.run(
            [
                "sacct",
                "-j",
                job_id,
                "--format=JobID,State,ExitCode",
                "--noheader",
                "--parsable2",
            ],
            capture_output=True,
            text=True,
            timeout=30,
        )

        if result.returncode != 0:
            return get_status_via_squeue(job_id)

        lines = result.stdout.strip().split("\n")
        if not lines:
            return get_status_via_squeue(job_id)

        primary_job = None
        for line in lines:
            if not line:
                continue
            parts = line.split("|")
            if len(parts) >= 2:
                jobid_part = parts[0]
                state = parts[1].strip()

                if jobid_part == job_id or jobid_part.startswith(job_id + "."):
                    if primary_job is None:
                        primary_job = state
                    elif state in ["RUNNING", "PD", "PENDING", "COMPLETING"]:
                        primary_job = state

        if primary_job:
            return parse_state(primary_job)

        return get_status_via_squeue(job_id)

    except (subprocess.TimeoutExpired, subprocess.SubprocessError, FileNotFoundError):
        return get_status_via_squeue(job_id)


def get_status_via_squeue(job_id: str) -> JobStatus:
    """Fallback to squeue when sacct is unavailable or times out."""
    try:
        result = subprocess.run(
            ["squeue", "-j", job_id, "--format=%T", "--noheader"],
            capture_output=True,
            text=True,
            timeout=10,
        )

        if result.returncode == 0 and result.stdout.strip():
            state = result.stdout.strip()
            return parse_state(state)

        return JobStatus.UNKNOWN

    except (subprocess.TimeoutExpired, subprocess.SubprocessError, FileNotFoundError):
        return JobStatus.UNKNOWN


def parse_state(state: str) -> JobStatus:
    """Map SLURM state strings to JobStatus enum."""
    state_upper = state.upper()

    if state_upper in ("CD", "COMPLETED", "NF", "NO_REASON"):
        return JobStatus.SUCCESS
    elif state_upper in (
        "F",
        "FAILED",
        "DL",
        "DEADLINE",
        "OM",
        "OOM",
        "OUT_OF_MEMORY",
        "PREEMPTED",
        "TO",
        "TIMEOUT",
    ):
        return JobStatus.FAILED
    elif state_upper in ("CG", "COMPLETING", "RV", "REVOKED"):
        return JobStatus.RUNNING
    elif state_upper in (
        "PD",
        "PENDING",
        "PR",
        "PREEMPTED_SUSPENDED",
        "S",
        "SUSPENDED",
        "SI",
        "SIGNALING",
    ):
        return JobStatus.PENDING
    elif state_upper in ("R", "RUNNING", "BF", "BOOT_FAIL"):
        return JobStatus.RUNNING
    elif state_upper in ("CA", "CANCELLED", "NC", "NODE_FAIL"):
        return JobStatus.CANCELED
    else:
        return JobStatus.UNKNOWN


def main():
    if len(sys.argv) < 2:
        print("UNKNOWN")
        sys.exit(1)

    job_id = sys.argv[1]

    max_retries = 3
    retry_delay = 1

    for attempt in range(max_retries):
        status = get_job_status(job_id)

        if status != JobStatus.UNKNOWN:
            print(status.value)
            sys.exit(0)

        if attempt < max_retries - 1:
            time.sleep(retry_delay)

    print(JobStatus.UNKNOWN.value)
    sys.exit(1)


if __name__ == "__main__":
    main()
