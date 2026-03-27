#!/usr/bin/env python3
"""
Tests for SLURM cluster profile configuration.

These tests validate the profile structure and configuration without
requiring access to an actual SLURM cluster.
"""

import os
import pytest
import tempfile
import yaml
from pathlib import Path


class TestSLURMProfileConfig:
    """Tests for SLURM profile configuration file."""

    @pytest.fixture
    def profile_dir(self):
        """Get path to SLURM profile directory."""
        return Path(__file__).parent.parent.parent / "profiles" / "slurm"

    def test_profile_config_exists(self, profile_dir):
        """Profile config.yaml should exist."""
        config_path = profile_dir / "config.yaml"
        assert config_path.exists(), f"Profile config not found at {config_path}"

    def test_cluster_status_script_exists(self, profile_dir):
        """cluster-status.py should exist."""
        script_path = profile_dir / "cluster-status.py"
        assert script_path.exists(), f"Status script not found at {script_path}"

    def test_profile_config_valid_yaml(self, profile_dir):
        """Profile config should be valid YAML."""
        config_path = profile_dir / "config.yaml"
        with open(config_path) as f:
            config = yaml.safe_load(f)

        assert config is not None
        assert isinstance(config, dict)

    def test_profile_has_executor(self, profile_dir):
        """Profile should specify executor."""
        config_path = profile_dir / "config.yaml"
        with open(config_path) as f:
            config = yaml.safe_load(f)

        assert "executor" in config
        assert config["executor"] == "cluster-generic"

    def test_profile_has_submit_cmd(self, profile_dir):
        """Profile should have cluster-generic-submit-cmd."""
        config_path = profile_dir / "config.yaml"
        with open(config_path) as f:
            config = yaml.safe_load(f)

        assert "cluster-generic-submit-cmd" in config
        submit_cmd = config["cluster-generic-submit-cmd"]
        assert "sbatch" in submit_cmd

    def test_profile_has_status_cmd(self, profile_dir):
        """Profile should have cluster-generic-status-cmd."""
        config_path = profile_dir / "config.yaml"
        with open(config_path) as f:
            config = yaml.safe_load(f)

        assert "cluster-generic-status-cmd" in config
        assert "cluster-status.py" in config["cluster-generic-status-cmd"]

    def test_profile_has_cancel_cmd(self, profile_dir):
        """Profile should have cluster-generic-cancel-cmd."""
        config_path = profile_dir / "config.yaml"
        with open(config_path) as f:
            config = yaml.safe_load(f)

        assert "cluster-generic-cancel-cmd" in config
        assert "scancel" in config["cluster-generic-cancel-cmd"]

    def test_profile_has_job_limits(self, profile_dir):
        """Profile should specify job limits."""
        config_path = profile_dir / "config.yaml"
        with open(config_path) as f:
            config = yaml.safe_load(f)

        assert "jobs" in config
        assert "max-jobs-per-second" in config
        assert "restart-times" in config


class TestClusterStatusScript:
    """Tests for SLURM status monitoring script."""

    @pytest.fixture
    def script_path(self):
        """Get path to cluster-status.py."""
        return (
            Path(__file__).parent.parent.parent
            / "profiles"
            / "slurm"
            / "cluster-status.py"
        )

    def test_script_exists(self, script_path):
        """Status script should exist."""
        assert script_path.exists()

    def test_script_is_executable(self, script_path):
        """Script should have execute permission or be runnable."""
        with open(script_path) as f:
            first_line = f.readline()
        assert first_line.startswith("#!")

    def test_script_has_main(self, script_path):
        """Script should have main function."""
        content = script_path.read_text()
        assert "def main(" in content

    def test_script_handles_no_args(self, script_path):
        """Script should handle missing arguments gracefully."""
        import subprocess

        result = subprocess.run(
            ["python3", str(script_path)],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0

    def test_script_imports(self, script_path):
        """Script should import required modules."""
        content = script_path.read_text()
        assert "import subprocess" in content
        assert "import sys" in content


class TestGenericProfile:
    """Tests for generic local execution profile."""

    @pytest.fixture
    def profile_dir(self):
        """Get path to generic profile directory."""
        return Path(__file__).parent.parent.parent / "profiles" / "generic"

    def test_profile_config_exists(self, profile_dir):
        """Profile config.yaml should exist."""
        config_path = profile_dir / "config.yaml"
        assert config_path.exists()

    def test_profile_config_valid_yaml(self, profile_dir):
        """Profile config should be valid YAML."""
        config_path = profile_dir / "config.yaml"
        with open(config_path) as f:
            config = yaml.safe_load(f)

        assert config is not None
        assert isinstance(config, dict)

    def test_generic_has_local_cores(self, profile_dir):
        """Generic profile should specify local-cores."""
        config_path = profile_dir / "config.yaml"
        with open(config_path) as f:
            config = yaml.safe_load(f)

        assert "local-cores" in config


class TestProfileReadme:
    """Tests for profile documentation."""

    @pytest.fixture
    def readme_path(self):
        """Get path to profiles README."""
        return Path(__file__).parent.parent.parent / "profiles" / "README.md"

    def test_readme_exists(self, readme_path):
        """README.md should exist."""
        assert readme_path.exists()

    def test_readme_has_slurm_info(self, readme_path):
        """README should document SLURM profile."""
        content = readme_path.read_text()
        assert "slurm" in content.lower()
        assert "sbatch" in content
