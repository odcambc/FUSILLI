#!/usr/bin/env python3
"""
Tests for resource allocation configuration.

These tests validate that resource specifications are correctly
defined for all rules and can be overridden per-rule.
"""

import pytest
import yaml
from pathlib import Path


class TestResourceSpecifications:
    """Tests for resource specification in rules."""

    @pytest.fixture
    def rules_dir(self):
        """Get path to rules directory."""
        return Path(__file__).parent.parent.parent / "workflow" / "rules"

    def test_common_smk_has_resource_helpers(self, rules_dir):
        """common.smk should define resource helper functions."""
        common_path = rules_dir / "common.smk"
        content = common_path.read_text()

        assert "get_partition" in content
        assert "get_runtime" in content
        assert "get_slurm_extra" in content

    def test_common_smk_has_default_runtimes(self, rules_dir):
        """common.smk should define default runtimes for rules."""
        common_path = rules_dir / "common.smk"
        content = common_path.read_text()

        assert "get_default_runtime" in content
        assert "runtimes = {" in content or "runtimes=" in content

    def test_common_smk_has_default_partition(self, rules_dir):
        """common.smk should define default SLURM partition."""
        common_path = rules_dir / "common.smk"
        content = common_path.read_text()

        assert "DEFAULT_PARTITION" in content

    def test_qc_rules_have_resources(self, rules_dir):
        """QC rules should define resources."""
        qc_path = rules_dir / "qc.smk"
        content = qc_path.read_text()

        assert "resources:" in content
        assert "mem_mb=" in content
        assert "runtime=" in content
        assert "partition=" in content

    def test_filter_rules_have_resources(self, rules_dir):
        """Filter rules should define resources."""
        filter_path = rules_dir / "filter_paired.smk"
        content = filter_path.read_text()

        assert "resources:" in content
        assert "mem_mb=" in content
        assert "runtime=" in content
        assert "partition=" in content

    def test_process_rules_have_resources(self, rules_dir):
        """String processing rules should define resources."""
        process_path = rules_dir / "process_strings.smk"
        content = process_path.read_text()

        assert "resources:" in content
        assert "mem_mb=" in content
        assert "runtime=" in content
        assert "partition=" in content

    def test_subsample_rules_have_resources(self, rules_dir):
        """Subsample rules should define resources."""
        subsample_path = rules_dir / "subsample.smk"
        content = subsample_path.read_text()

        assert "resources:" in content
        assert "mem_mb=" in content
        assert "runtime=" in content
        assert "partition=" in content

    def test_reproducibility_rules_have_resources(self, rules_dir):
        """Reproducibility rules should define resources."""
        repro_path = rules_dir / "reproducibility.smk"
        content = repro_path.read_text()

        assert "resources:" in content
        assert "mem_mb=" in content
        assert "runtime=" in content
        assert "partition=" in content


class TestConfigSchemaResources:
    """Tests for resource configuration in schema."""

    @pytest.fixture
    def schema_path(self):
        """Get path to config schema."""
        return (
            Path(__file__).parent.parent.parent
            / "workflow"
            / "schemas"
            / "config.schema.yaml"
        )

    def test_schema_has_resources_section(self, schema_path):
        """Schema should define resources section."""
        with open(schema_path) as f:
            schema = yaml.safe_load(f)

        assert "properties" in schema
        assert "resources" in schema["properties"]

    def test_schema_has_memory_mb(self, schema_path):
        """Schema should define memory_mb property."""
        with open(schema_path) as f:
            schema = yaml.safe_load(f)

        resources = schema["properties"]["resources"]
        assert "memory_mb" in resources["properties"]

    def test_schema_has_threads(self, schema_path):
        """Schema should define threads property."""
        with open(schema_path) as f:
            schema = yaml.safe_load(f)

        resources = schema["properties"]["resources"]
        assert "threads" in resources["properties"]

    def test_schema_has_rule_resources(self, schema_path):
        """Schema should define rule_resources for per-rule overrides."""
        with open(schema_path) as f:
            schema = yaml.safe_load(f)

        resources = schema["properties"]["resources"]
        assert "rule_resources" in resources["properties"]

    def test_schema_has_cluster_section(self, schema_path):
        """Schema should define cluster section."""
        with open(schema_path) as f:
            schema = yaml.safe_load(f)

        assert "cluster" in schema["properties"]
        cluster = schema["properties"]["cluster"]
        assert "slurm" in cluster["properties"]

    def test_schema_has_container_section(self, schema_path):
        """Schema should define container section."""
        with open(schema_path) as f:
            schema = yaml.safe_load(f)

        assert "container" in schema["properties"]
        container = schema["properties"]["container"]
        assert "runtime" in container["properties"]


class TestConfigYamlResources:
    """Tests for resource configuration in example config."""

    @pytest.fixture
    def config_path(self):
        """Get path to config.yaml."""
        return Path(__file__).parent.parent.parent / "config" / "config.yaml"

    def test_config_has_resources_section(self, config_path):
        """Config should define resources section."""
        with open(config_path) as f:
            config = yaml.safe_load(f)

        assert "resources" in config
        assert "memory_mb" in config["resources"]
        assert "threads" in config["resources"]

    def test_config_has_cluster_section(self, config_path):
        """Config should have cluster section (commented or active)."""
        content = config_path.read_text()
        assert "cluster:" in content

    def test_config_has_container_section(self, config_path):
        """Config should have container section (commented or active)."""
        content = config_path.read_text()
        assert "container:" in content

    def test_config_has_rule_resources_example(self, config_path):
        """Config should have rule_resources examples."""
        content = config_path.read_text()
        assert "rule_resources:" in content


class TestDefaultRuntimes:
    """Tests for default runtime specifications."""

    @pytest.fixture
    def common_path(self):
        """Get path to common.smk."""
        return Path(__file__).parent.parent.parent / "workflow" / "rules" / "common.smk"

    def test_fastqc_has_runtime(self, common_path):
        """FastQC should have runtime defined."""
        content = common_path.read_text()
        assert '"fastqc":' in content or "'fastqc':" in content

    def test_trim_adapters_has_runtime(self, common_path):
        """Trim adapters should have runtime defined."""
        content = common_path.read_text()
        assert '"trim_adapters":' in content or "'trim_adapters':" in content

    def test_merge_reads_has_runtime(self, common_path):
        """Merge reads should have runtime defined."""
        content = common_path.read_text()
        assert '"merge_reads":' in content or "'merge_reads':" in content

    def test_detect_fusions_has_runtime(self, common_path):
        """Detect fusions should have runtime defined."""
        content = common_path.read_text()
        assert (
            '"detect_fusions_string":' in content
            or "'detect_fusions_string':" in content
        )
