#!/usr/bin/env python3
"""
Unit tests for string_matcher.py

Tests the fusion detection logic via string matching.
"""

import pytest
import sys
import tempfile
import gzip
from pathlib import Path

# Add workflow/scripts to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))

from string_matcher import (
    find_matches_in_read,
    find_unfused_matches_in_read,
    load_breakpoint_sequences,
    load_domain_ends,
    load_unfused_kmers,
    count_fusion_matches,
    parse_fastq_python,
    write_counts_csv,
)


# =============================================================================
# TEST DATA
# =============================================================================

# Simulated breakpoint sequences
TEST_BREAKPOINTS = {
    'TPR': {
        'TPR_126_Met_WT': 'CTGAAGGGGCGGGGGAGCATGAAA',
        'TPR_129_Met_WT': 'AAGGGGCGGCATGGGAGCATGAAA',
    },
    'CCDC6': {
        'CCDC6_300_Met_WT': 'GCCAGCGTGACCGGGAGCATGAAA',
        'CCDC6_303_Met_WT': 'AGCGTGACCATCGGGAGCATGAAA',
    }
}

# Domain end k-mers (3' ends of partners)
TEST_DOMAIN_ENDS = {
    'TPR': 'GAATACTTAACA',  # Last 12 nt of TPR
    'CCDC6': 'GTGACCATC',   # Last 9 nt of CCDC6
}


# =============================================================================
# TEST: find_matches_in_read
# =============================================================================

class TestFindMatchesInRead:
    """Tests for the core matching function."""

    def test_exact_match(self):
        """Should find exact breakpoint match."""
        # Read contains TPR_126 breakpoint sequence
        read = "AAAA" + TEST_BREAKPOINTS['TPR']['TPR_126_Met_WT'] + "TTTT"

        # Need domain end for pre-filter
        domain_ends = {'TPR': TEST_BREAKPOINTS['TPR']['TPR_126_Met_WT'][-12:]}

        matches = find_matches_in_read(read, domain_ends, TEST_BREAKPOINTS)

        assert 'TPR_126_Met_WT' in matches

    def test_no_match(self):
        """Should return empty list when no match."""
        read = "AAAAAAAAAAAATTTTTTTTTTTT"

        matches = find_matches_in_read(read, TEST_DOMAIN_ENDS, TEST_BREAKPOINTS)

        assert matches == []

    def test_prefilter_gates_matching(self):
        """Pre-filter should prevent unnecessary breakpoint searches."""
        # Read contains breakpoint but NOT the domain end (pre-filter fails)
        breakpoint_seq = TEST_BREAKPOINTS['TPR']['TPR_126_Met_WT']
        # Use a domain end that doesn't appear in the read
        domain_ends = {'TPR': 'XXXXXXXXXXXX'}

        read = "AAAA" + breakpoint_seq + "TTTT"

        matches = find_matches_in_read(read, domain_ends, TEST_BREAKPOINTS)

        # Should not find match because pre-filter fails
        assert matches == []

    def test_multiple_matches(self):
        """Should find multiple breakpoints in same read if present."""
        # Unlikely in real data, but the function should handle it
        bp1 = TEST_BREAKPOINTS['TPR']['TPR_126_Met_WT']
        bp2 = TEST_BREAKPOINTS['TPR']['TPR_129_Met_WT']

        # Construct read with both (with proper domain end for pre-filter)
        read = bp1 + "NNNN" + bp2
        domain_ends = {'TPR': 'GGGAGCATGAAA'}  # Common suffix

        matches = find_matches_in_read(read, domain_ends, TEST_BREAKPOINTS)

        assert 'TPR_126_Met_WT' in matches
        assert 'TPR_129_Met_WT' in matches

    def test_match_from_different_partners(self):
        """Should correctly identify matches from different partners."""
        # Read with TPR breakpoint
        tpr_read = "AAAA" + TEST_BREAKPOINTS['TPR']['TPR_126_Met_WT'] + "TTTT"
        domain_ends = {'TPR': 'GGGAGCATGAAA'}

        matches = find_matches_in_read(tpr_read, domain_ends, TEST_BREAKPOINTS)

        assert any('TPR' in m for m in matches)
        assert not any('CCDC6' in m for m in matches)


# =============================================================================
# TEST: load_breakpoint_sequences
# =============================================================================

class TestLoadBreakpointSequences:
    """Tests for loading breakpoint CSV files."""

    def test_loads_csv_correctly(self, tmp_path):
        """Should parse breakpoint CSV into nested dict."""
        csv_content = """fusion_id,partner_name,anchor_name,breakpoint_nt,breakpoint_aa,breakpoint_sequence,full_fusion_length
TPR_126_Met_WT,TPR,Met_WT,126,42,CTGAAGGGGCGGGGGAGCATGAAA,372
TPR_129_Met_WT,TPR,Met_WT,129,43,AAGGGGCGGCATGGGAGCATGAAA,375
CCDC6_300_Met_WT,CCDC6,Met_WT,300,100,GCCAGCGTGACCGGGAGCATGAAA,546
"""
        csv_file = tmp_path / "breakpoints.csv"
        csv_file.write_text(csv_content)

        breakpoints = load_breakpoint_sequences(str(csv_file))

        assert 'TPR' in breakpoints
        assert 'CCDC6' in breakpoints
        assert 'TPR_126_Met_WT' in breakpoints['TPR']
        assert breakpoints['TPR']['TPR_126_Met_WT'] == 'CTGAAGGGGCGGGGGAGCATGAAA'

    def test_organizes_by_partner(self, tmp_path):
        """Should organize breakpoints by partner name for fast lookup."""
        csv_content = """fusion_id,partner_name,anchor_name,breakpoint_nt,breakpoint_aa,breakpoint_sequence,full_fusion_length
A_1_X,A,X,1,1,AAA,10
A_2_X,A,X,2,1,BBB,10
B_1_X,B,X,1,1,CCC,10
"""
        csv_file = tmp_path / "breakpoints.csv"
        csv_file.write_text(csv_content)

        breakpoints = load_breakpoint_sequences(str(csv_file))

        assert len(breakpoints['A']) == 2
        assert len(breakpoints['B']) == 1
        assert set(breakpoints.keys()) == {'A', 'B'}


# =============================================================================
# TEST: load_domain_ends
# =============================================================================

class TestLoadDomainEnds:
    """Tests for loading domain ends CSV files."""

    def test_loads_csv_correctly(self, tmp_path):
        """Should parse domain ends CSV into dict."""
        csv_content = """domain_name,end_kmer
TPR,GAATACTTAACA
CCDC6,GTGACCATC
"""
        csv_file = tmp_path / "ends.csv"
        csv_file.write_text(csv_content)

        ends = load_domain_ends(str(csv_file))

        assert ends['TPR'] == 'GAATACTTAACA'
        assert ends['CCDC6'] == 'GTGACCATC'


# =============================================================================
# TEST: load_unfused_kmers
# =============================================================================

class TestLoadUnfusedKmers:
    """Tests for loading unfused k-mer CSV files."""

    def test_loads_kmers_grouped_by_length(self, tmp_path):
        """Should group k-mers by length and map to sequence names."""
        csv_content = """sequence_name,kmer,position_in_sequence,sequence_length
KRAS,AAAAA,10,100
KRAS,TTTT,20,100
NRAS,AAAAA,15,90
"""
        csv_file = tmp_path / "unfused.csv"
        csv_file.write_text(csv_content)

        kmers_by_len = load_unfused_kmers(str(csv_file))

        assert 5 in kmers_by_len
        assert 4 in kmers_by_len
        assert kmers_by_len[5]["AAAAA"] == ["KRAS", "NRAS"]
        assert kmers_by_len[4]["TTTT"] == ["KRAS"]


# =============================================================================
# TEST: find_unfused_matches_in_read
# =============================================================================

class TestFindUnfusedMatchesInRead:
    """Tests for unfused k-mer matching in reads."""

    def test_detects_unfused_kmer(self):
        """Should detect unfused sequence when k-mer is present."""
        kmers_by_len = {3: {"AAA": ["KRAS"]}}
        matches = find_unfused_matches_in_read("GGGAAATTT", kmers_by_len)
        assert matches == {"KRAS"}

    def test_detects_reverse_complement(self):
        """Should detect reverse-complement matches when enabled."""
        kmers_by_len = {3: {"AAA": ["KRAS"]}}
        matches = find_unfused_matches_in_read("TTT", kmers_by_len, orientation_check=True)
        assert matches == {"KRAS"}


# =============================================================================
# TEST: write_counts_csv
# =============================================================================

class TestWriteCountsCsv:
    """Tests for writing fusion and unfused counts."""

    def test_writes_type_column_with_unfused(self, tmp_path):
        """Should include type column when unfused counts are provided."""
        out_file = tmp_path / "counts.csv"
        write_counts_csv(
            {"TPR_1_Met_WT": 3},
            str(out_file),
            unfused_counts={"KRAS": 2},
            include_type=True
        )

        rows = out_file.read_text().strip().splitlines()
        assert rows[0] == "fusion_id,type,count"
        assert "TPR_1_Met_WT,fusion,3" in rows
        assert "KRAS,unfused,2" in rows

# =============================================================================
# TEST: FASTQ Parsing
# =============================================================================

class TestFastqParsing:
    """Tests for FASTQ file parsing."""

    def test_parse_plain_fastq(self, tmp_path):
        """Should parse plain text FASTQ."""
        fastq_content = """@read1
ATGCATGCATGC
+
IIIIIIIIIIII
@read2
GCTAGCTAGCTA
+
HHHHHHHHHHHH
"""
        fastq_file = tmp_path / "test.fastq"
        fastq_file.write_text(fastq_content)

        reads = list(parse_fastq_python(str(fastq_file)))

        assert len(reads) == 2
        assert reads[0][0] == 'read1'
        assert reads[0][1] == 'ATGCATGCATGC'
        assert reads[1][0] == 'read2'
        assert reads[1][1] == 'GCTAGCTAGCTA'

    def test_parse_gzipped_fastq(self, tmp_path):
        """Should parse gzipped FASTQ."""
        fastq_content = """@read1
ATGCATGCATGC
+
IIIIIIIIIIII
"""
        fastq_file = tmp_path / "test.fastq.gz"
        with gzip.open(fastq_file, 'wt') as f:
            f.write(fastq_content)

        reads = list(parse_fastq_python(str(fastq_file)))

        assert len(reads) == 1
        assert reads[0][1] == 'ATGCATGCATGC'

    def test_uppercase_conversion(self, tmp_path):
        """Should convert sequences to uppercase."""
        fastq_content = """@read1
atgcatgcatgc
+
IIIIIIIIIIII
"""
        fastq_file = tmp_path / "test.fastq"
        fastq_file.write_text(fastq_content)

        reads = list(parse_fastq_python(str(fastq_file)))

        assert reads[0][1] == 'ATGCATGCATGC'


# =============================================================================
# TEST: count_fusion_matches (Integration)
# =============================================================================

class TestCountFusionMatches:
    """Integration tests for full counting workflow."""

    def test_counts_matches_in_file(self, tmp_path):
        """Should count fusion matches across a FASTQ file."""
        # Create test breakpoints
        breakpoints = {
            'TPR': {
                'TPR_test': 'ATGCATGCATGC',
            }
        }
        domain_ends = {
            'TPR': 'ATGCATGC',
        }

        # Create FASTQ with known matches
        fastq_content = """@read1
AAAATGCATGCATGCAAA
+
IIIIIIIIIIIIIIIIII
@read2
AAAATGCATGCATGCAAA
+
IIIIIIIIIIIIIIIIII
@read3
GGGGGGGGGGGGGGGGGG
+
IIIIIIIIIIIIIIIIII
"""
        fastq_file = tmp_path / "test.fastq"
        fastq_file.write_text(fastq_content)

        counts = count_fusion_matches(
            str(fastq_file),
            breakpoints,
            domain_ends,
            show_progress=False
        )

        assert counts.get('TPR_test', 0) == 2

    def test_empty_file_returns_empty_counts(self, tmp_path):
        """Should handle empty/minimal files gracefully."""
        # pyfastx throws RuntimeError on truly empty files, so use a minimal valid FASTQ
        # with a read that won't match anything
        fastq_file = tmp_path / "empty.fastq"
        fastq_file.write_text("@empty\nN\n+\nI\n")

        counts = count_fusion_matches(
            str(fastq_file),
            TEST_BREAKPOINTS,
            TEST_DOMAIN_ENDS,
            show_progress=False
        )

        assert counts == {}

    def test_prefilter_fallback_detects_breakpoint(self, tmp_path):
        """Prefilter fallback should detect breakpoint when domain end is absent."""
        breakpoints = {
            'TPR': {
                'TPR_no_end': "AATTT",
            }
        }
        # Domain end does not appear in the read
        domain_ends = {'TPR': 'CCCC'}

        fastq_file = tmp_path / "nolinker.fastq"
        fastq_file.write_text("@r1\nCCAATTTGG\n+\nIIIIIIIII\n")

        counts = count_fusion_matches(
            str(fastq_file),
            breakpoints,
            domain_ends,
            show_progress=False,
            prefilter_fallback=True
        )

        assert counts.get('TPR_no_end', 0) == 1

    def test_prefilter_fallback_disabled(self, tmp_path):
        """Without fallback, missing domain end should prevent detection."""
        breakpoints = {
            'TPR': {
                'TPR_no_end': "AATTT",
            }
        }
        domain_ends = {'TPR': 'CCCC'}

        fastq_file = tmp_path / "nofallback.fastq"
        fastq_file.write_text("@r1\nCCAATTTGG\n+\nIIIIIIIII\n")

        counts = count_fusion_matches(
            str(fastq_file),
            breakpoints,
            domain_ends,
            show_progress=False,
            prefilter_fallback=False
        )

        assert counts.get('TPR_no_end', 0) == 0


# =============================================================================
# RUN TESTS
# =============================================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
