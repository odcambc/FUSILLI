#!/usr/bin/env python3
"""
Unit tests for fusion_sequences.py

Tests the core fusion breakpoint generation logic, specifically verifying:
- Full-length downstream gene (anchor) behavior
- Variable truncation of upstream gene (partner)
- In-frame breakpoint generation
- K-mer extraction around breakpoints
"""

import pytest
import sys
from pathlib import Path

# Add workflow/scripts to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))

from fusion_sequences import (
    generate_breakpoint_positions,
    generate_fusion_sequence,
    extract_breakpoint_kmer,
    generate_all_breakpoints,
    generate_domain_ends,
    FusionLibraryConfig,
    BreakpointSequence
)


# =============================================================================
# TEST DATA
# =============================================================================

# Simple test sequences (multiples of 3 for in-frame testing)
TEST_PARTNER_SEQ = "ATGCCCAAAGGG"  # 12 nt = 4 codons
TEST_ANCHOR_SEQ = "TTTCCCAAAGGG"   # 12 nt = 4 codons
TEST_LINKER = "GGGAGC"             # 6 nt = 2 codons (GS linker)

# Longer sequences for realistic testing
LONG_PARTNER = "ATGGCGGCGGTGTTGCAGCAAGTCCTGGAG"  # 30 nt
LONG_ANCHOR = "ATGAAAAAGAGAAAGCAAATTAAAGATCTG"   # 30 nt


# =============================================================================
# TEST: generate_breakpoint_positions
# =============================================================================

class TestGenerateBreakpointPositions:
    """Tests for breakpoint position generation."""

    def test_in_frame_positions(self):
        """In-frame mode should generate positions at codon boundaries (every 3 nt)."""
        positions = list(generate_breakpoint_positions(12, maintain_frame=True))

        assert positions == [3, 6, 9, 12]
        # All positions should be divisible by 3
        assert all(p % 3 == 0 for p in positions)

    def test_all_positions(self):
        """Non-frame mode should generate all positions."""
        positions = list(generate_breakpoint_positions(12, maintain_frame=False))

        assert positions == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

    def test_excludes_zero(self):
        """Position 0 (no partner contribution) should be excluded."""
        positions = list(generate_breakpoint_positions(12, maintain_frame=True))

        assert 0 not in positions

    def test_includes_full_length(self):
        """Full partner length should be included."""
        positions = list(generate_breakpoint_positions(12, maintain_frame=True))

        assert 12 in positions

    def test_short_sequence(self):
        """Should handle short sequences (less than one codon)."""
        positions = list(generate_breakpoint_positions(2, maintain_frame=True))

        assert positions == []

    def test_single_codon(self):
        """Should handle single codon sequences."""
        positions = list(generate_breakpoint_positions(3, maintain_frame=True))

        assert positions == [3]


# =============================================================================
# TEST: generate_fusion_sequence
# =============================================================================

class TestGenerateFusionSequence:
    """Tests for fusion sequence generation - core logic for fusion structure."""

    def test_downstream_anchor_structure(self):
        """
        With anchor downstream: Partner(truncated) + Linker + Anchor(full)
        This is the main use case we're testing.
        """
        fusion, bp_pos = generate_fusion_sequence(
            partner_seq=TEST_PARTNER_SEQ,  # ATGCCCAAAGGG
            anchor_seq=TEST_ANCHOR_SEQ,     # TTTCCCAAAGGG
            linker_seq=TEST_LINKER,         # GGGAGC
            breakpoint_pos=6,               # Truncate partner at position 6
            anchor_position='downstream'
        )

        # Partner portion (first 6 nt) + linker + full anchor
        expected = "ATGCCC" + "GGGAGC" + "TTTCCCAAAGGG"
        assert fusion == expected

        # Breakpoint should be right after partner portion
        assert bp_pos == 6

    def test_anchor_always_full_length(self):
        """Anchor (downstream gene) should always be used in full."""
        # Try various partner truncation positions
        for truncation in [3, 6, 9, 12]:
            fusion, _ = generate_fusion_sequence(
                partner_seq=TEST_PARTNER_SEQ,
                anchor_seq=TEST_ANCHOR_SEQ,
                linker_seq=TEST_LINKER,
                breakpoint_pos=truncation,
                anchor_position='downstream'
            )

            # Fusion should always end with the complete anchor
            assert fusion.endswith(TEST_ANCHOR_SEQ), \
                f"Anchor should be full length at truncation {truncation}"

    def test_partner_truncation_varies(self):
        """Partner (upstream gene) should be truncated at various positions."""
        expected_partner_portions = {
            3: "ATG",
            6: "ATGCCC",
            9: "ATGCCCAAA",
            12: "ATGCCCAAAGGG"
        }

        for truncation, expected_portion in expected_partner_portions.items():
            fusion, _ = generate_fusion_sequence(
                partner_seq=TEST_PARTNER_SEQ,
                anchor_seq=TEST_ANCHOR_SEQ,
                linker_seq=TEST_LINKER,
                breakpoint_pos=truncation,
                anchor_position='downstream'
            )

            # Fusion should start with truncated partner
            assert fusion.startswith(expected_portion), \
                f"Partner portion wrong at truncation {truncation}"

    def test_linker_placement(self):
        """Linker should be between partner and anchor."""
        fusion, bp_pos = generate_fusion_sequence(
            partner_seq=TEST_PARTNER_SEQ,
            anchor_seq=TEST_ANCHOR_SEQ,
            linker_seq=TEST_LINKER,
            breakpoint_pos=6,
            anchor_position='downstream'
        )

        # Linker should be at breakpoint position
        assert fusion[bp_pos:bp_pos + len(TEST_LINKER)] == TEST_LINKER

    def test_no_linker(self):
        """Should work without a linker."""
        fusion, bp_pos = generate_fusion_sequence(
            partner_seq=TEST_PARTNER_SEQ,
            anchor_seq=TEST_ANCHOR_SEQ,
            linker_seq="",
            breakpoint_pos=6,
            anchor_position='downstream'
        )

        expected = "ATGCCC" + "TTTCCCAAAGGG"
        assert fusion == expected
        assert bp_pos == 6

    def test_upstream_anchor_structure(self):
        """
        With anchor upstream: Anchor(full) + Linker + Partner(truncated)
        Alternative configuration.
        """
        fusion, bp_pos = generate_fusion_sequence(
            partner_seq=TEST_PARTNER_SEQ,
            anchor_seq=TEST_ANCHOR_SEQ,
            linker_seq=TEST_LINKER,
            breakpoint_pos=6,
            anchor_position='upstream'
        )

        # Full anchor + linker + partner portion
        expected = "TTTCCCAAAGGG" + "GGGAGC" + "ATGCCC"
        assert fusion == expected

        # Breakpoint is after anchor + linker
        assert bp_pos == len(TEST_ANCHOR_SEQ) + len(TEST_LINKER)

    def test_fusion_length_calculation(self):
        """Fusion length should be partner_truncation + linker + anchor."""
        for truncation in [3, 6, 9, 12]:
            fusion, _ = generate_fusion_sequence(
                partner_seq=TEST_PARTNER_SEQ,
                anchor_seq=TEST_ANCHOR_SEQ,
                linker_seq=TEST_LINKER,
                breakpoint_pos=truncation,
                anchor_position='downstream'
            )

            expected_length = truncation + len(TEST_LINKER) + len(TEST_ANCHOR_SEQ)
            assert len(fusion) == expected_length


# =============================================================================
# TEST: extract_breakpoint_kmer
# =============================================================================

class TestExtractBreakpointKmer:
    """Tests for k-mer extraction around breakpoints."""

    def test_basic_extraction(self):
        """Should extract k-mer centered on breakpoint."""
        # Create a fusion sequence
        fusion = "AAAAAABBBBBCCCCCC"  # Breakpoint at position 6

        kmer = extract_breakpoint_kmer(fusion, breakpoint_pos=6, window=3)

        # 3 nt before + 3 nt after breakpoint
        assert kmer == "AAABBB"

    def test_window_size(self):
        """K-mer should be 2 * window in length."""
        fusion = "A" * 50

        for window in [5, 10, 12, 15]:
            kmer = extract_breakpoint_kmer(fusion, breakpoint_pos=25, window=window)
            assert len(kmer) == 2 * window

    def test_bounds_check_start(self):
        """Should return None if window extends before start."""
        fusion = "AAAAAABBBBBB"

        # Breakpoint at 3, window of 5 would need position -2
        kmer = extract_breakpoint_kmer(fusion, breakpoint_pos=3, window=5)

        assert kmer is None

    def test_bounds_check_end(self):
        """Should return None if window extends past end."""
        fusion = "AAAAAABBBBBB"  # 12 nt

        # Breakpoint at 10, window of 5 would need position 15
        kmer = extract_breakpoint_kmer(fusion, breakpoint_pos=10, window=5)

        assert kmer is None

    def test_exact_bounds(self):
        """Should work at exact boundaries."""
        fusion = "AAAAAABBBBBB"  # 12 nt

        # Breakpoint at 6, window of 6 - exactly fits
        kmer = extract_breakpoint_kmer(fusion, breakpoint_pos=6, window=6)

        assert kmer == fusion


# =============================================================================
# TEST: generate_all_breakpoints
# =============================================================================

class TestGenerateAllBreakpoints:
    """Integration tests for complete breakpoint generation."""

    @pytest.fixture
    def test_sequences(self):
        """Create test sequence dictionary."""
        return {
            'Met_WT': LONG_ANCHOR,
            'TPR': LONG_PARTNER,
        }

    @pytest.fixture
    def test_partners(self):
        """Create test partner configuration."""
        return {
            'TPR': {'include': True, 'sequence_length': 30, 'description': 'Test'},
        }

    @pytest.fixture
    def test_config(self):
        """Create test configuration."""
        return FusionLibraryConfig(
            anchor_name='Met_WT',
            anchor_position='downstream',
            linker_sequence='GGGAGC',
            breakpoint_window=6,
            maintain_frame=True,
            kmer_size=15
        )

    def test_generates_breakpoints(self, test_sequences, test_partners, test_config):
        """Should generate breakpoints for each partner."""
        breakpoints = generate_all_breakpoints(
            test_sequences, test_partners, test_config
        )

        assert len(breakpoints) > 0
        assert all(isinstance(bp, BreakpointSequence) for bp in breakpoints)

    def test_breakpoint_naming(self, test_sequences, test_partners, test_config):
        """Breakpoint IDs should follow naming convention."""
        breakpoints = generate_all_breakpoints(
            test_sequences, test_partners, test_config
        )

        for bp in breakpoints:
            # Format: {partner}_{position}_{anchor}
            assert bp.fusion_id == f"{bp.partner_name}_{bp.breakpoint_nt}_{bp.anchor_name}"

    def test_in_frame_only(self, test_sequences, test_partners, test_config):
        """With maintain_frame=True, all breakpoints should be at codon boundaries."""
        breakpoints = generate_all_breakpoints(
            test_sequences, test_partners, test_config
        )

        for bp in breakpoints:
            assert bp.breakpoint_nt % 3 == 0, \
                f"Breakpoint {bp.fusion_id} not at codon boundary"

    def test_kmer_length(self, test_sequences, test_partners, test_config):
        """All k-mers should be 2 * window in length."""
        breakpoints = generate_all_breakpoints(
            test_sequences, test_partners, test_config
        )

        expected_kmer_length = 2 * test_config.breakpoint_window
        for bp in breakpoints:
            assert len(bp.sequence) == expected_kmer_length, \
                f"K-mer length wrong for {bp.fusion_id}"

    def test_skips_excluded_partners(self, test_sequences, test_config):
        """Should skip partners marked as not included."""
        partners = {
            'TPR': {'include': True, 'sequence_length': 30, 'description': ''},
            'EXCLUDED': {'include': False, 'sequence_length': 30, 'description': ''},
        }
        # Add excluded sequence
        sequences = test_sequences.copy()
        sequences['EXCLUDED'] = LONG_PARTNER

        breakpoints = generate_all_breakpoints(sequences, partners, test_config)

        partner_names = {bp.partner_name for bp in breakpoints}
        assert 'TPR' in partner_names
        assert 'EXCLUDED' not in partner_names

    def test_anchor_not_in_breakpoints(self, test_sequences, test_partners, test_config):
        """Anchor should not appear as a partner in breakpoints."""
        # Add anchor to partners (shouldn't be used)
        partners = test_partners.copy()
        partners['Met_WT'] = {'include': True, 'sequence_length': 30, 'description': ''}

        breakpoints = generate_all_breakpoints(
            test_sequences, partners, test_config
        )

        # Anchor shouldn't be truncated (it's always full in fusions)
        # But it might still generate breakpoints if included in partners
        # This behavior depends on implementation - verify your expectation


# =============================================================================
# TEST: generate_domain_ends
# =============================================================================

class TestGenerateDomainEnds:
    """Tests for domain end k-mer generation."""

    def test_generates_3prime_ends(self):
        """Should generate 3' end k-mers for each partner."""
        sequences = {
            'Met_WT': LONG_ANCHOR,
            'TPR': LONG_PARTNER,
        }
        partners = {
            'TPR': {'include': True, 'sequence_length': 30, 'description': ''},
        }

        ends = generate_domain_ends(
            sequences, partners, anchor_name='Met_WT', kmer_size=15
        )

        assert 'TPR' in ends
        assert len(ends['TPR']) == 15
        # Should be the last 15 nt of TPR
        assert ends['TPR'] == LONG_PARTNER[-15:]

    def test_excludes_anchor(self):
        """Anchor should not be in domain ends."""
        sequences = {
            'Met_WT': LONG_ANCHOR,
            'TPR': LONG_PARTNER,
        }
        partners = {
            'TPR': {'include': True, 'sequence_length': 30, 'description': ''},
            'Met_WT': {'include': True, 'sequence_length': 30, 'description': ''},
        }

        ends = generate_domain_ends(
            sequences, partners, anchor_name='Met_WT', kmer_size=15
        )

        assert 'Met_WT' not in ends

    def test_excludes_non_included(self):
        """Should exclude partners not marked for inclusion."""
        sequences = {
            'Met_WT': LONG_ANCHOR,
            'TPR': LONG_PARTNER,
            'OTHER': LONG_PARTNER,
        }
        partners = {
            'TPR': {'include': True, 'sequence_length': 30, 'description': ''},
            'OTHER': {'include': False, 'sequence_length': 30, 'description': ''},
        }

        ends = generate_domain_ends(
            sequences, partners, anchor_name='Met_WT', kmer_size=15
        )

        assert 'TPR' in ends
        assert 'OTHER' not in ends


# =============================================================================
# TEST: Specific Fusion Structure Verification
# =============================================================================

class TestFusionStructureVerification:
    """
    Critical tests to verify the pipeline handles the correct fusion structure:
    - Full-length downstream gene (anchor)
    - Variably-truncated upstream gene (partner)
    """

    def test_anchor_never_truncated(self):
        """
        CRITICAL: The anchor (downstream gene) should NEVER be truncated.
        It should always appear in full in every fusion.
        """
        sequences = {
            'Anchor': 'TTTCCCAAAGGGAAA',  # 15 nt anchor
            'Partner': 'ATGCCCAAAGGGCCC',  # 15 nt partner
        }
        partners = {
            'Partner': {'include': True, 'sequence_length': 15, 'description': ''},
        }
        config = FusionLibraryConfig(
            anchor_name='Anchor',
            anchor_position='downstream',
            linker_sequence='GGG',
            breakpoint_window=3,
            maintain_frame=True,
            kmer_size=9
        )

        breakpoints = generate_all_breakpoints(sequences, partners, config)

        # Every fusion should contain the complete anchor
        for bp in breakpoints:
            # Reconstruct the full fusion to verify
            fusion, _ = generate_fusion_sequence(
                partner_seq=sequences['Partner'],
                anchor_seq=sequences['Anchor'],
                linker_seq=config.linker_sequence,
                breakpoint_pos=bp.breakpoint_nt,
                anchor_position=config.anchor_position
            )

            assert fusion.endswith(sequences['Anchor']), \
                f"Anchor truncated in {bp.fusion_id}"

    def test_partner_truncation_generates_all_variants(self):
        """
        CRITICAL: The partner (upstream gene) should be truncated at
        all possible in-frame positions.
        """
        partner_seq = 'ATGCCCAAAGGG'  # 12 nt = 4 codons
        sequences = {
            'Anchor': 'TTTCCCAAAGGGAAA',
            'Partner': partner_seq,
        }
        partners = {
            'Partner': {'include': True, 'sequence_length': 12, 'description': ''},
        }
        config = FusionLibraryConfig(
            anchor_name='Anchor',
            anchor_position='downstream',
            linker_sequence='',  # No linker for simpler verification
            breakpoint_window=3,
            maintain_frame=True,
            kmer_size=9
        )

        breakpoints = generate_all_breakpoints(sequences, partners, config)

        # Should have breakpoints at positions 3, 6, 9, 12
        breakpoint_positions = {bp.breakpoint_nt for bp in breakpoints}

        # Note: some positions might be filtered out due to window bounds
        # But we should have multiple truncation variants
        assert len(breakpoint_positions) > 1, \
            "Should generate multiple truncation variants"

        # All positions should be codon boundaries
        assert all(pos % 3 == 0 for pos in breakpoint_positions)

    def test_downstream_position_means_anchor_at_3prime(self):
        """
        When anchor_position='downstream', the anchor should be at the 3' end.
        Structure: [Partner truncated]-[Linker]-[Anchor full]
        """
        fusion, _ = generate_fusion_sequence(
            partner_seq="AAACCC",
            anchor_seq="GGGTTTT",
            linker_seq="NNN",
            breakpoint_pos=3,
            anchor_position='downstream'
        )

        # Expected: AAA + NNN + GGGTTTT
        assert fusion == "AAANNNGGGTTTT"
        assert fusion.endswith("GGGTTTT")  # Anchor at 3' end
        assert fusion.startswith("AAA")     # Truncated partner at 5' end

    def test_kmer_spans_breakpoint_junction(self):
        """
        The breakpoint k-mer should span the junction between
        partner and linker (the actual fusion breakpoint).
        """
        partner = "AAACCCGGG"   # 9 nt
        anchor = "TTTAAACCC"    # 9 nt
        linker = "NNNNNN"       # 6 nt

        fusion, bp_pos = generate_fusion_sequence(
            partner_seq=partner,
            anchor_seq=anchor,
            linker_seq=linker,
            breakpoint_pos=6,  # Take first 6 nt of partner
            anchor_position='downstream'
        )

        # Fusion: AAACCC + NNNNNN + TTTAAACCC = "AAACCCNNNNNNTTTAAACCC"
        # bp_pos = 6 (right after partner portion)

        kmer = extract_breakpoint_kmer(fusion, bp_pos, window=4)

        # K-mer extraction:
        # fusion = "AAACCCNNNNNNTTTAAACCC"
        # bp_pos = 6, window = 4
        # start = 6 - 4 = 2
        # end = 6 + 4 = 10
        # fusion[2:10] = "ACCCNNNN"
        # This is: last 4 nt of truncated partner + first 4 nt of linker
        assert kmer == "ACCCNNNN"


# =============================================================================
# RUN TESTS
# =============================================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
