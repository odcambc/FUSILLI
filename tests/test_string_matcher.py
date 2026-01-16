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
    find_partner_hits,
    load_breakpoint_sequences,
    load_domain_ends,
    load_unfused_kmers,
    count_fusion_matches,
    parse_fastq_python,
    write_counts_csv,
    build_domain_ends_automaton,
    build_breakpoints_automaton,
    build_partner_breakpoints_automata,
    build_unfused_kmers_automata,
    find_matches_aho,
    reverse_complement,
    get_reverse_complement_cached,
    HAS_AHOCORASICK,
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
            include_type=True,
            expected_fusions=["TPR_1_Met_WT", "EGFR_2_Met_WT"],
            expected_unfused=["KRAS", "NRAS"]
        )

        rows = out_file.read_text().strip().splitlines()
        assert rows[0] == "fusion_id,type,count"
        assert "TPR_1_Met_WT,fusion,3" in rows
        assert "EGFR_2_Met_WT,fusion,0" in rows
        assert "KRAS,unfused,2" in rows
        assert "NRAS,unfused,0" in rows

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
# TEST: Unmerged Read Processing
# =============================================================================

class TestUnmergedReadProcessing:
    """Tests for processing unmerged paired-end reads (R1 and R2 separately)."""

    def test_process_unmerged_r1_separately(self, tmp_path):
        """Should process R1 unmerged reads separately and detect breakpoints."""
        breakpoints = {
            'TPR': {
                'TPR_126_Met_WT': 'CTGAAGGGGCGGGGGAGCATGAAA',
            }
        }
        # Use a domain end that appears in the breakpoint sequence
        domain_ends = {'TPR': 'GGGAGCATGAAA'}

        # Create R1 FASTQ with known breakpoint
        r1_fastq = tmp_path / "sample_R1.unmerged.fastq.gz"
        with gzip.open(r1_fastq, 'wt') as f:
            f.write("@read1\n")
            f.write("AAAA" + breakpoints['TPR']['TPR_126_Met_WT'] + "TTTT\n")
            f.write("+\n")
            f.write("I" * (4 + len(breakpoints['TPR']['TPR_126_Met_WT']) + 4) + "\n")

        counts = count_fusion_matches(
            str(r1_fastq),
            breakpoints,
            domain_ends,
            show_progress=False
        )

        assert counts.get('TPR_126_Met_WT', 0) == 1

    def test_process_unmerged_r2_separately(self, tmp_path):
        """Should process R2 unmerged reads separately and detect breakpoints."""
        breakpoints = {
            'TPR': {
                'TPR_126_Met_WT': 'CTGAAGGGGCGGGGGAGCATGAAA',
            }
        }
        domain_ends = {'TPR': 'GGGAGCATGAAA'}

        # Create R2 FASTQ with known breakpoint
        r2_fastq = tmp_path / "sample_R2.unmerged.fastq.gz"
        with gzip.open(r2_fastq, 'wt') as f:
            f.write("@read1\n")
            f.write("CCCC" + breakpoints['TPR']['TPR_126_Met_WT'] + "GGGG\n")
            f.write("+\n")
            f.write("I" * (4 + len(breakpoints['TPR']['TPR_126_Met_WT']) + 4) + "\n")

        counts = count_fusion_matches(
            str(r2_fastq),
            breakpoints,
            domain_ends,
            show_progress=False
        )

        assert counts.get('TPR_126_Met_WT', 0) == 1

    def test_unmerged_counts_distinct_per_mate(self, tmp_path):
        """Verify R1 and R2 unmerged counts are separate and distinct."""
        breakpoints = {
            'TPR': {
                'TPR_126_Met_WT': 'CTGAAGGGGCGGGGGAGCATGAAA',
                'TPR_129_Met_WT': 'AAGGGGCGGCATGGGAGCATGAAA',
            }
        }
        domain_ends = {'TPR': 'GGGAGCATGAAA'}

        # R1 has one breakpoint
        r1_fastq = tmp_path / "sample_R1.unmerged.fastq.gz"
        with gzip.open(r1_fastq, 'wt') as f:
            f.write("@read1\n")
            f.write("AAAA" + breakpoints['TPR']['TPR_126_Met_WT'] + "TTTT\n")
            f.write("+\n")
            f.write("I" * (4 + len(breakpoints['TPR']['TPR_126_Met_WT']) + 4) + "\n")

        # R2 has a different breakpoint
        r2_fastq = tmp_path / "sample_R2.unmerged.fastq.gz"
        with gzip.open(r2_fastq, 'wt') as f:
            f.write("@read1\n")
            f.write("CCCC" + breakpoints['TPR']['TPR_129_Met_WT'] + "GGGG\n")
            f.write("+\n")
            f.write("I" * (4 + len(breakpoints['TPR']['TPR_129_Met_WT']) + 4) + "\n")

        r1_counts = count_fusion_matches(
            str(r1_fastq),
            breakpoints,
            domain_ends,
            show_progress=False
        )

        r2_counts = count_fusion_matches(
            str(r2_fastq),
            breakpoints,
            domain_ends,
            show_progress=False
        )

        # R1 should only have TPR_126_Met_WT
        assert r1_counts.get('TPR_126_Met_WT', 0) == 1
        assert r1_counts.get('TPR_129_Met_WT', 0) == 0

        # R2 should only have TPR_129_Met_WT
        assert r2_counts.get('TPR_126_Met_WT', 0) == 0
        assert r2_counts.get('TPR_129_Met_WT', 0) == 1

    def test_empty_unmerged_file_handled_gracefully(self, tmp_path):
        """Should handle empty unmerged FASTQ files without errors."""
        breakpoints = {
            'TPR': {
                'TPR_126_Met_WT': 'CTGAAGGGGCGGGGGAGCATGAAA',
            }
        }
        domain_ends = {'TPR': 'GGGAGCATGAAA'}

        # Create empty gzipped FASTQ (just header, no reads)
        empty_fastq = tmp_path / "empty.unmerged.fastq.gz"
        with gzip.open(empty_fastq, 'wt') as f:
            pass  # Empty file

        counts = count_fusion_matches(
            str(empty_fastq),
            breakpoints,
            domain_ends,
            show_progress=False
        )

        # Should return empty counts, not raise error
        assert counts == {}

    def test_unmerged_file_with_no_matches(self, tmp_path):
        """Should return empty counts when no breakpoints are found."""
        breakpoints = {
            'TPR': {
                'TPR_126_Met_WT': 'CTGAAGGGGCGGGGGAGCATGAAA',
            }
        }
        domain_ends = {'TPR': 'GGGAGCATGAAA'}

        # Create FASTQ with random sequence (no breakpoint)
        no_match_fastq = tmp_path / "no_match.unmerged.fastq.gz"
        with gzip.open(no_match_fastq, 'wt') as f:
            f.write("@read1\n")
            f.write("ATGCATGCATGCATGCATGCATGCATGC\n")
            f.write("+\n")
            f.write("I" * 28 + "\n")

        counts = count_fusion_matches(
            str(no_match_fastq),
            breakpoints,
            domain_ends,
            show_progress=False
        )

        assert counts == {}

    def test_unmerged_r1_and_r2_same_fusion_detected_separately(self, tmp_path):
        """Same fusion in R1 and R2 should be counted separately."""
        breakpoints = {
            'TPR': {
                'TPR_126_Met_WT': 'CTGAAGGGGCGGGGGAGCATGAAA',
            }
        }
        domain_ends = {'TPR': 'GGGAGCATGAAA'}

        # Both R1 and R2 have the same breakpoint
        r1_fastq = tmp_path / "sample_R1.unmerged.fastq.gz"
        with gzip.open(r1_fastq, 'wt') as f:
            f.write("@read1\n")
            f.write("AAAA" + breakpoints['TPR']['TPR_126_Met_WT'] + "TTTT\n")
            f.write("+\n")
            f.write("I" * (4 + len(breakpoints['TPR']['TPR_126_Met_WT']) + 4) + "\n")

        r2_fastq = tmp_path / "sample_R2.unmerged.fastq.gz"
        with gzip.open(r2_fastq, 'wt') as f:
            f.write("@read1\n")
            f.write("CCCC" + breakpoints['TPR']['TPR_126_Met_WT'] + "GGGG\n")
            f.write("+\n")
            f.write("I" * (4 + len(breakpoints['TPR']['TPR_126_Met_WT']) + 4) + "\n")

        r1_counts = count_fusion_matches(
            str(r1_fastq),
            breakpoints,
            domain_ends,
            show_progress=False
        )

        r2_counts = count_fusion_matches(
            str(r2_fastq),
            breakpoints,
            domain_ends,
            show_progress=False
        )

        # Both should detect the same fusion, but counts are separate
        assert r1_counts.get('TPR_126_Met_WT', 0) == 1
        assert r2_counts.get('TPR_126_Met_WT', 0) == 1

        # Total would be 2 if combined, but they're kept separate
        assert r1_counts.get('TPR_126_Met_WT', 0) + r2_counts.get('TPR_126_Met_WT', 0) == 2


# =============================================================================
# TEST: Aho-Corasick Optimizations
# =============================================================================

class TestAhoCorasickOptimizations:
    """Tests for Aho-Corasick automaton-based optimizations."""

    @pytest.mark.skipif(not HAS_AHOCORASICK, reason="pyahocorasick not available")
    def test_build_domain_ends_automaton(self):
        """Should build automaton for domain ends."""
        automaton = build_domain_ends_automaton(TEST_DOMAIN_ENDS)
        assert automaton is not None

        # Test matching
        matches = find_matches_aho("AAA" + TEST_DOMAIN_ENDS['TPR'] + "TTT", automaton)
        assert 'TPR' in matches
        assert 'CCDC6' not in matches

    @pytest.mark.skipif(not HAS_AHOCORASICK, reason="pyahocorasick not available")
    def test_build_breakpoints_automaton(self):
        """Should build automaton for breakpoint sequences."""
        automaton = build_breakpoints_automaton(TEST_BREAKPOINTS)
        assert automaton is not None

        # Test matching
        bp_seq = TEST_BREAKPOINTS['TPR']['TPR_126_Met_WT']
        matches = find_matches_aho("AAA" + bp_seq + "TTT", automaton)
        # Matches should be tuples of (partner, fusion_id)
        assert any(m[0] == 'TPR' and m[1] == 'TPR_126_Met_WT' for m in matches if isinstance(m, tuple))

    @pytest.mark.skipif(not HAS_AHOCORASICK, reason="pyahocorasick not available")
    def test_build_partner_breakpoints_automata(self):
        """Should build separate automata per partner."""
        automata = build_partner_breakpoints_automata(TEST_BREAKPOINTS)
        assert 'TPR' in automata
        assert 'CCDC6' in automata

        # Test TPR automaton
        bp_seq = TEST_BREAKPOINTS['TPR']['TPR_126_Met_WT']
        matches = find_matches_aho("AAA" + bp_seq + "TTT", automata['TPR'])
        assert 'TPR_126_Met_WT' in matches

    @pytest.mark.skipif(not HAS_AHOCORASICK, reason="pyahocorasick not available")
    def test_build_unfused_kmers_automata(self):
        """Should build automata for unfused k-mers grouped by length."""
        unfused_kmers = {
            3: {"AAA": ["KRAS"], "TTT": ["NRAS"]},
            4: {"AAAA": ["KRAS", "NRAS"]}
        }
        automata = build_unfused_kmers_automata(unfused_kmers)
        assert 3 in automata
        assert 4 in automata

        # Test matching
        matches = find_matches_aho("GGGAAATTT", automata[3])
        # Matches should be tuples of sequence names
        assert any("KRAS" in m for m in matches if isinstance(m, tuple))
        assert any("NRAS" in m for m in matches if isinstance(m, tuple))

    @pytest.mark.skipif(not HAS_AHOCORASICK, reason="pyahocorasick not available")
    def test_find_matches_aho_returns_correct_types(self):
        """Should return correct value types from automata."""
        # Domain ends return strings
        domain_automaton = build_domain_ends_automaton(TEST_DOMAIN_ENDS)
        matches = find_matches_aho(TEST_DOMAIN_ENDS['TPR'], domain_automaton)
        assert all(isinstance(m, str) for m in matches)

        # Breakpoints return tuples
        bp_automaton = build_breakpoints_automaton(TEST_BREAKPOINTS)
        bp_seq = TEST_BREAKPOINTS['TPR']['TPR_126_Met_WT']
        matches = find_matches_aho(bp_seq, bp_automaton)
        assert any(isinstance(m, tuple) and len(m) == 2 for m in matches)

    @pytest.mark.skipif(not HAS_AHOCORASICK, reason="pyahocorasick not available")
    def test_aho_matches_same_as_original(self):
        """Aho-Corasick implementation should produce same results as original."""
        read = "AAAA" + TEST_BREAKPOINTS['TPR']['TPR_126_Met_WT'] + "TTTT"
        domain_ends = {'TPR': TEST_BREAKPOINTS['TPR']['TPR_126_Met_WT'][-12:]}

        # Original implementation
        original_matches = find_matches_in_read(
            read, domain_ends, TEST_BREAKPOINTS,
            domain_ends_automaton=None  # Force original
        )

        # Aho-Corasick implementation
        domain_automaton = build_domain_ends_automaton(domain_ends)
        partner_automata = build_partner_breakpoints_automata(TEST_BREAKPOINTS)
        aho_matches = find_matches_in_read(
            read, domain_ends, TEST_BREAKPOINTS,
            domain_ends_automaton=domain_automaton,
            partner_breakpoints_automata=partner_automata
        )

        assert set(original_matches) == set(aho_matches)

    @pytest.mark.skipif(not HAS_AHOCORASICK, reason="pyahocorasick not available")
    def test_aho_unfused_matches_same_as_original(self):
        """Aho-Corasick unfused matching should match original."""
        kmers_by_len = {3: {"AAA": ["KRAS"], "TTT": ["NRAS"]}}
        read = "GGGAAATTTCCC"

        # Original implementation
        original_matches = find_unfused_matches_in_read(
            read, kmers_by_len, unfused_automata=None
        )

        # Aho-Corasick implementation
        automata = build_unfused_kmers_automata(kmers_by_len)
        aho_matches = find_unfused_matches_in_read(
            read, kmers_by_len, unfused_automata=automata
        )

        assert original_matches == aho_matches

    @pytest.mark.skipif(not HAS_AHOCORASICK, reason="pyahocorasick not available")
    def test_aho_partner_hits_same_as_original(self):
        """Aho-Corasick partner hits should match original."""
        read = "AAA" + TEST_DOMAIN_ENDS['TPR'] + "TTT"
        linker = "GGG"

        # Original implementation
        original_hits, original_linker = find_partner_hits(
            read, TEST_DOMAIN_ENDS, linker_sequence=linker,
            domain_ends_automaton=None
        )

        # Aho-Corasick implementation
        automaton = build_domain_ends_automaton(TEST_DOMAIN_ENDS)
        aho_hits, aho_linker = find_partner_hits(
            read, TEST_DOMAIN_ENDS, linker_sequence=linker,
            domain_ends_automaton=automaton
        )

        assert original_hits == aho_hits
        assert original_linker == aho_linker

    @pytest.mark.skipif(not HAS_AHOCORASICK, reason="pyahocorasick not available")
    def test_aho_with_orientation_check(self):
        """Aho-Corasick should work with orientation checking."""
        read = "AAAA" + TEST_BREAKPOINTS['TPR']['TPR_126_Met_WT'] + "TTTT"
        domain_ends = {'TPR': TEST_BREAKPOINTS['TPR']['TPR_126_Met_WT'][-12:]}

        # Build reverse complement data
        rc_domain_ends = {k: reverse_complement(v) for k, v in domain_ends.items()}
        rc_breakpoints = {
            partner: {fid: reverse_complement(seq) for fid, seq in bp_dict.items()}
            for partner, bp_dict in TEST_BREAKPOINTS.items()
        }

        # Original
        original_matches, orig_f, orig_rc = find_matches_in_read(
            read, domain_ends, TEST_BREAKPOINTS,
            orientation_check=True,
            rc_domain_ends=rc_domain_ends,
            rc_breakpoints=rc_breakpoints,
            return_orientation=True,
            domain_ends_automaton=None
        )

        # Aho-Corasick
        domain_automaton = build_domain_ends_automaton(domain_ends)
        rc_domain_automaton = build_domain_ends_automaton(rc_domain_ends)
        partner_automata = build_partner_breakpoints_automata(TEST_BREAKPOINTS)
        rc_partner_automata = build_partner_breakpoints_automata(rc_breakpoints)

        aho_matches, aho_f, aho_rc = find_matches_in_read(
            read, domain_ends, TEST_BREAKPOINTS,
            orientation_check=True,
            rc_domain_ends=rc_domain_ends,
            rc_breakpoints=rc_breakpoints,
            return_orientation=True,
            domain_ends_automaton=domain_automaton,
            rc_domain_ends_automaton=rc_domain_automaton,
            partner_breakpoints_automata=partner_automata,
            rc_partner_breakpoints_automata=rc_partner_automata
        )

        assert set(original_matches) == set(aho_matches)
        assert orig_f == aho_f
        assert orig_rc == aho_rc

    @pytest.mark.skipif(not HAS_AHOCORASICK, reason="pyahocorasick not available")
    def test_aho_with_prefilter_fallback(self):
        """Aho-Corasick should work with prefilter fallback."""
        breakpoints = {'TPR': {'TPR_test': 'AATTT'}}
        domain_ends = {'TPR': 'CCCC'}  # Won't match

        read = "CCAATTTGG"

        # Original
        original_matches = find_matches_in_read(
            read, domain_ends, breakpoints,
            prefilter_fallback=True,
            domain_ends_automaton=None
        )

        # Aho-Corasick
        domain_automaton = build_domain_ends_automaton(domain_ends)
        breakpoints_automaton = build_breakpoints_automaton(breakpoints)
        partner_automata = build_partner_breakpoints_automata(breakpoints)

        aho_matches = find_matches_in_read(
            read, domain_ends, breakpoints,
            prefilter_fallback=True,
            domain_ends_automaton=domain_automaton,
            breakpoints_automaton=breakpoints_automaton,
            partner_breakpoints_automata=partner_automata
        )

        assert set(original_matches) == set(aho_matches)

    def test_fallback_when_ahocorasick_unavailable(self, monkeypatch):
        """Should fall back to original implementation when pyahocorasick unavailable."""
        # Mock HAS_AHOCORASICK to be False
        import string_matcher
        original_has_aho = string_matcher.HAS_AHOCORASICK
        string_matcher.HAS_AHOCORASICK = False

        try:
            read = "AAAA" + TEST_BREAKPOINTS['TPR']['TPR_126_Met_WT'] + "TTTT"
            domain_ends = {'TPR': TEST_BREAKPOINTS['TPR']['TPR_126_Met_WT'][-12:]}

            # Should still work, using original implementation
            matches = find_matches_in_read(
                read, domain_ends, TEST_BREAKPOINTS,
                domain_ends_automaton=None  # None because HAS_AHOCORASICK is False
            )

            assert 'TPR_126_Met_WT' in matches
        finally:
            # Restore original value
            string_matcher.HAS_AHOCORASICK = original_has_aho

    @pytest.mark.skipif(not HAS_AHOCORASICK, reason="pyahocorasick not available")
    def test_count_fusion_matches_uses_automata(self, tmp_path):
        """count_fusion_matches should use automata when available."""
        breakpoints = {
            'TPR': {
                'TPR_test': 'ATGCATGCATGC',
            }
        }
        domain_ends = {'TPR': 'ATGCATGC'}

        fastq_content = """@read1
AAAATGCATGCATGCAAA
+
IIIIIIIIIIIIIIIIII
@read2
AAAATGCATGCATGCAAA
+
IIIIIIIIIIIIIIIIII
"""
        fastq_file = tmp_path / "test.fastq"
        fastq_file.write_text(fastq_content)

        # Should use automata internally (no error means it worked)
        counts = count_fusion_matches(
            str(fastq_file),
            breakpoints,
            domain_ends,
            show_progress=False
        )

        assert counts.get('TPR_test', 0) == 2

    @pytest.mark.skipif(not HAS_AHOCORASICK, reason="pyahocorasick not available")
    def test_empty_automata_handled_gracefully(self):
        """Should handle empty automata gracefully."""
        # Empty domain ends
        empty_automaton = build_domain_ends_automaton({})
        assert empty_automaton is None or find_matches_aho("ATGC", empty_automaton) == set()

        # Empty breakpoints
        empty_bp_automaton = build_breakpoints_automaton({})
        assert empty_bp_automaton is None or find_matches_aho("ATGC", empty_bp_automaton) == set()

        # Empty unfused
        empty_unfused = build_unfused_kmers_automata({})
        assert empty_unfused == {}


# =============================================================================
# RUN TESTS
# =============================================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
