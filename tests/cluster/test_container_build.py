#!/usr/bin/env python3
"""
Tests for container definitions.

These tests validate container configuration files without
requiring Docker or Singularity to be installed.
"""

import pytest
from pathlib import Path


class TestDockerfile:
    """Tests for Docker container definition."""

    @pytest.fixture
    def dockerfile_path(self):
        """Get path to Dockerfile."""
        return Path(__file__).parent.parent.parent / "containers" / "Dockerfile"

    def test_dockerfile_exists(self, dockerfile_path):
        """Dockerfile should exist."""
        assert dockerfile_path.exists()

    def test_dockerfile_has_from(self, dockerfile_path):
        """Dockerfile should have FROM instruction."""
        content = dockerfile_path.read_text()
        assert content.startswith("FROM")

    def test_dockerfile_uses_micromamba(self, dockerfile_path):
        """Dockerfile should use micromamba base image."""
        content = dockerfile_path.read_text()
        assert "micromamba" in content.lower()

    def test_dockerfile_has_workdir(self, dockerfile_path):
        """Dockerfile should set working directory."""
        content = dockerfile_path.read_text()
        assert "WORKDIR" in content

    def test_dockerfile_copies_env_file(self, dockerfile_path):
        """Dockerfile should copy environment file."""
        content = dockerfile_path.read_text()
        assert "COPY" in content
        assert "fusilli_env.yaml" in content

    def test_dockerfile_installs_packages(self, dockerfile_path):
        """Dockerfile should install packages."""
        content = dockerfile_path.read_text()
        assert "micromamba install" in content

    def test_dockerfile_has_cmd(self, dockerfile_path):
        """Dockerfile should have CMD instruction."""
        content = dockerfile_path.read_text()
        assert "CMD" in content


class TestSingularityDefinition:
    """Tests for Singularity/Apptainer definition file."""

    @pytest.fixture
    def def_file(self):
        """Get path to Singularity definition file."""
        return Path(__file__).parent.parent.parent / "containers" / "fusilli.def"

    def test_def_file_exists(self, def_file):
        """Singularity definition file should exist."""
        assert def_file.exists()

    def test_def_uses_docker_bootstrap(self, def_file):
        """Definition should use Docker bootstrap."""
        content = def_file.read_text()
        assert "Bootstrap: docker" in content
        assert "From:" in content

    def test_def_has_micromamba(self, def_file):
        """Definition should use micromamba."""
        content = def_file.read_text()
        assert "micromamba" in content.lower()

    def test_def_has_files_section(self, def_file):
        """Definition should have %files section."""
        content = def_file.read_text()
        assert "%files" in content

    def test_def_has_post_section(self, def_file):
        """Definition should have %post section."""
        content = def_file.read_text()
        assert "%post" in content

    def test_def_has_environment_section(self, def_file):
        """Definition should have %environment section."""
        content = def_file.read_text()
        assert "%environment" in content

    def test_def_has_runscript(self, def_file):
        """Definition should have %runscript section."""
        content = def_file.read_text()
        assert "%runscript" in content


class TestContainerBuildScripts:
    """Tests for container build scripts."""

    def test_docker_build_script_exists(self):
        """Docker build script should exist."""
        script_path = (
            Path(__file__).parent.parent.parent / "containers" / "build_docker.sh"
        )
        assert script_path.exists()

    def test_docker_build_script_is_executable(self):
        """Docker build script should be executable."""
        script_path = (
            Path(__file__).parent.parent.parent / "containers" / "build_docker.sh"
        )
        content = script_path.read_text()
        assert content.startswith("#!")
        assert "docker build" in content

    def test_singularity_build_script_exists(self):
        """Singularity build script should exist."""
        script_path = (
            Path(__file__).parent.parent.parent / "containers" / "build_singularity.sh"
        )
        assert script_path.exists()

    def test_singularity_build_script_is_executable(self):
        """Singularity build script should be executable."""
        script_path = (
            Path(__file__).parent.parent.parent / "containers" / "build_singularity.sh"
        )
        content = script_path.read_text()
        assert content.startswith("#!")
        assert "singularity build" in content


class TestContainerReadme:
    """Tests for container documentation."""

    @pytest.fixture
    def readme_path(self):
        """Get path to containers README."""
        return Path(__file__).parent.parent.parent / "containers" / "README.md"

    def test_readme_exists(self, readme_path):
        """README.md should exist."""
        assert readme_path.exists()

    def test_readme_documents_docker(self, readme_path):
        """README should document Docker usage."""
        content = readme_path.read_text()
        assert "docker" in content.lower()
        assert "Dockerfile" in content

    def test_readme_documents_singularity(self, readme_path):
        """README should document Singularity usage."""
        content = readme_path.read_text()
        assert "singularity" in content.lower() or "apptainer" in content.lower()

    def test_readme_has_quick_start(self, readme_path):
        """README should have quick start section."""
        content = readme_path.read_text()
        assert "Quick Start" in content or "quick start" in content.lower()


class TestContainerEnvironmentFile:
    """Tests for container environment file reference."""

    def test_env_file_referenced_in_dockerfile(self):
        """Dockerfile should reference fusilli_env.yaml."""
        dockerfile_path = (
            Path(__file__).parent.parent.parent / "containers" / "Dockerfile"
        )
        content = dockerfile_path.read_text()
        assert "fusilli_env.yaml" in content

    def test_env_file_referenced_in_singularity(self):
        """Singularity definition should reference fusilli_env.yaml."""
        def_path = Path(__file__).parent.parent.parent / "containers" / "fusilli.def"
        content = def_path.read_text()
        assert "fusilli_env.yaml" in content
