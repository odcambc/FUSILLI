#!/usr/bin/env python3
"""
Pytest configuration and shared fixtures for FUSILLI tests.
"""

import pytest
import sys
from pathlib import Path

# Add workflow/scripts to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))


# =============================================================================
# SHARED TEST DATA
# =============================================================================

# Standard test sequences used across multiple tests
TEST_ANCHOR_SEQ = (
    "ATGAAAAAGAGAAAGCAAATTAAAGATCTGGGCAGTGAATTAGTTCGCTACGATGCAAGA"
    "GTACACACTCCTCATTTGGATAGGCTTGTAAGTGCCCGAAGTGTAAGCCCAACTACAGAA"
    "ATGGTTTCAAATGAATCTGTAGACTACCGAGCTACTTTTCCAGAAGATCAGTTTCCTAAT"
    "TCATCTCAGAACGGTTCATGCCGACAAGTGCAGTATCCTCTGACAGACATGTCCCCCATC"
)  # 240 nt - simulates kinase domain (anchor)

TEST_PARTNER_SEQ = (
    "ATGGCGGCGGTGTTGCAGCAAGTCCTGGAGCGCACGGAGCTGAACAAGCTGCCCAAGTCT"
    "GTCCAGAACAAACTTGAAAAGTTCCTTGCTGATCAGCAATCCGAGATCGATGGCCTGAAG"
    "GGGCGGCATGAGAAATTTAAGGTGGAGAGCGAACAACAGTATTTTGAAATAGAAAAGAGG"
)  # 180 nt - simulates fusion partner (truncated upstream)

TEST_LINKER = "GGGAGC"  # Standard GS linker


# =============================================================================
# FIXTURES
# =============================================================================

@pytest.fixture
def test_sequences():
    """Standard test sequences dict."""
    return {
        'Met_WT': TEST_ANCHOR_SEQ,
        'TPR': TEST_PARTNER_SEQ,
    }


@pytest.fixture
def test_partners():
    """Standard test partners configuration."""
    return {
        'TPR': {
            'include': True,
            'sequence_length': len(TEST_PARTNER_SEQ),
            'description': 'Test TPR partner'
        },
    }


@pytest.fixture
def test_config():
    """Standard test fusion library configuration."""
    from fusion_sequences import FusionLibraryConfig

    return FusionLibraryConfig(
        anchor_name='Met_WT',
        anchor_position='downstream',
        linker_sequence=TEST_LINKER,
        breakpoint_window=12,
        maintain_frame=True,
        kmer_size=15
    )


@pytest.fixture
def simple_sequences():
    """Minimal test sequences for unit tests."""
    return {
        'Anchor': 'TTTCCCAAAGGG',  # 12 nt
        'Partner': 'ATGCCCAAAGGG',  # 12 nt
    }


@pytest.fixture
def simple_partners():
    """Minimal partners config for unit tests."""
    return {
        'Partner': {
            'include': True,
            'sequence_length': 12,
            'description': 'Simple test partner'
        },
    }


# =============================================================================
# TEST FILE FIXTURES
# =============================================================================

@pytest.fixture
def fasta_file(tmp_path, test_sequences):
    """Create a temporary FASTA file with test sequences."""
    fasta_path = tmp_path / "test.fasta"
    with open(fasta_path, 'w') as f:
        for name, seq in test_sequences.items():
            f.write(f">{name}\n{seq}\n")
    return fasta_path


@pytest.fixture
def partners_csv(tmp_path, test_partners):
    """Create a temporary partners CSV file."""
    csv_path = tmp_path / "partners.csv"
    with open(csv_path, 'w') as f:
        f.write("partner_name,sequence_length,include,description\n")
        for name, config in test_partners.items():
            include_str = 'true' if config['include'] else 'false'
            f.write(f"{name},{config['sequence_length']},{include_str},{config['description']}\n")
    return csv_path


@pytest.fixture
def fastq_file(tmp_path):
    """Create an empty temporary FASTQ file."""
    fastq_path = tmp_path / "test.fastq"
    fastq_path.touch()
    return fastq_path


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def create_fastq_with_reads(filepath, reads):
    """
    Create a FASTQ file with specified reads.

    Args:
        filepath: Path to write FASTQ file
        reads: List of sequence strings
    """
    with open(filepath, 'w') as f:
        for i, seq in enumerate(reads):
            f.write(f"@read{i}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write("I" * len(seq) + "\n")


def create_gzipped_fastq(filepath, reads):
    """Create a gzipped FASTQ file with specified reads."""
    import gzip

    with gzip.open(filepath, 'wt') as f:
        for i, seq in enumerate(reads):
            f.write(f"@read{i}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write("I" * len(seq) + "\n")


# =============================================================================
# PYTEST CONFIGURATION
# =============================================================================

def pytest_configure(config):
    """Configure custom pytest markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "integration: marks integration tests"
    )


# =============================================================================
# TEST COLLECTION HOOKS
# =============================================================================

def pytest_collection_modifyitems(config, items):
    """
    Modify test collection to add markers based on test names/locations.
    """
    for item in items:
        # Mark tests in test_integration.py as integration tests
        if "integration" in item.fspath.basename:
            item.add_marker(pytest.mark.integration)
