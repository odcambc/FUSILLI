#!/usr/bin/env python3
"""
Integration tests for FUSILLI pipeline.

These tests verify the end-to-end behavior of the pipeline, specifically:
    1. Full-length upstream gene (partner) in all fusions
    2. Variable truncation of downstream gene (anchor)
3. Correct detection of fusion breakpoints in reads
"""

import pytest
import sys
import gzip
import csv
import shutil
import subprocess
import os
from pathlib import Path

# Add workflow/scripts to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))

from fusion_sequences import (
    generate_fusion_sequence,
    generate_all_breakpoints,
    generate_domain_ends,
    FusionLibraryConfig,
    run_generation,
)
from string_matcher import (
    find_matches_in_read,
    count_fusion_matches,
    run_matching,
)


# =============================================================================
# TEST DATA - Realistic Sequences
# =============================================================================

# Simulated sequences based on real kinase fusion structure
MET_KINASE = (
    "ATGAAAAAGAGAAAGCAAATTAAAGATCTGGGCAGTGAATTAGTTCGCTACGATGCAAGA"
    "GTACACACTCCTCATTTGGATAGGCTTGTAAGTGCCCGAAGTGTAAGCCCAACTACAGAA"
    "ATGGTTTCAAATGAATCTGTAGACTACCGAGCTACTTTTCCAGAAGATCAGTTTCCTAAT"
    "TCATCTCAGAACGGTTCATGCCGACAAGTGCAGTATCCTCTGACAGACATGTCCCCCATC"
)  # 240 nt

TPR_PARTNER = (
    "ATGGCGGCGGTGTTGCAGCAAGTCCTGGAGCGCACGGAGCTGAACAAGCTGCCCAAGTCT"
    "GTCCAGAACAAACTTGAAAAGTTCCTTGCTGATCAGCAATCCGAGATCGATGGCCTGAAG"
    "GGGCGGCATGAGAAATTTAAGGTGGAGAGCGAACAACAGTATTTTGAAATAGAAAAGAGG"
)  # 180 nt

LINKER = "GGGAGC"  # GS linker


def _write_fastq(path: Path, reads: list[str]) -> None:
    """Write a minimal FASTQ file for testing."""
    with path.open("w") as fh:
        for i, seq in enumerate(reads, 1):
            fh.write(f"@read{i}\n{seq}\n+\n{'I' * len(seq)}\n")


def _write_fastq_gz(path: Path, reads: list[str]) -> None:
    """Write a minimal gzipped FASTQ file for testing."""
    with gzip.open(path, "wt") as fh:
        for i, seq in enumerate(reads, 1):
            fh.write(f"@read{i}\n{seq}\n+\n{'I' * len(seq)}\n")


# =============================================================================
# TEST: End-to-End Fusion Generation and Detection
# =============================================================================

class TestEndToEndFusionDetection:
    """
    Integration tests that verify the full pipeline correctly:
    1. Generates breakpoint sequences for full partner + truncated anchor
    2. Detects those breakpoints in simulated reads
    """

    @pytest.fixture
    def sequences(self):
        """Test sequences dict."""
        return {
            'Met_WT': MET_KINASE,
            'TPR': TPR_PARTNER,
        }

    @pytest.fixture
    def partners(self):
        """Test partners config."""
        return {
            'TPR': {'include': True, 'sequence_length': 180, 'description': ''},
        }

    @pytest.fixture
    def config(self):
        """Test fusion library config."""
        return FusionLibraryConfig(
            anchor_name='Met_WT',
            anchor_position='downstream',
            truncated_component='anchor',
            linker_sequence=LINKER,
            breakpoint_window=13,
            maintain_frame=True,
            kmer_size=15
        )

    def test_generated_fusions_have_full_partner(self, sequences, partners, config):
        """
        CRITICAL: Every generated fusion should contain the complete partner.
        """
        breakpoints = generate_all_breakpoints(sequences, partners, config)

        for bp in breakpoints:
            # Reconstruct the full fusion
            fusion, _ = generate_fusion_sequence(
                partner_seq=sequences['TPR'],
                anchor_seq=sequences['Met_WT'],
                linker_seq=config.linker_sequence,
                breakpoint_pos=bp.breakpoint_nt,
                anchor_position=config.anchor_position,
                truncated_component=config.truncated_component
            )

            # Verify partner is complete at 5' end
            assert fusion.startswith(sequences['TPR']), \
                f"Fusion {bp.fusion_id} does not start with complete partner"
            # Verify anchor truncation matches breakpoint position
            expected_anchor = sequences['Met_WT'][bp.breakpoint_nt:]
            assert fusion.endswith(expected_anchor), \
                f"Fusion {bp.fusion_id} does not end with truncated anchor"

    def test_different_truncations_produce_different_breakpoints(self, sequences, partners, config):
        """Each partner truncation should produce a unique breakpoint sequence."""
        breakpoints = generate_all_breakpoints(sequences, partners, config)

        # Get unique breakpoint sequences
        bp_sequences = [bp.sequence for bp in breakpoints]
        unique_sequences = set(bp_sequences)

        # Most should be unique (some edge cases might not be)
        assert len(unique_sequences) == len(bp_sequences), \
            "Breakpoint sequences should be unique for each truncation"

    def test_breakpoint_kmer_contains_junction(self, sequences, partners, config):
        """
        The breakpoint k-mer should span the junction between
        linker and anchor.
        """
        breakpoints = generate_all_breakpoints(sequences, partners, config)

        for bp in breakpoints:
            # The k-mer should contain part of the partner AND part of linker/anchor
            kmer = bp.sequence

            # With window=13, the kmer is 26 nt
            # First 13 nt from end of partner+linker
            # Last 13 nt from anchor

            # Check that the linker appears in the junction region
            # The linker "GGGAGC" should be in the k-mer for most breakpoints
            upstream_portion = sequences['TPR'] + LINKER

            if len(upstream_portion) >= config.breakpoint_window:
                # K-mer should contain the last part of partner+linker
                upstream_end = upstream_portion[-config.breakpoint_window:]
                assert kmer.startswith(upstream_end), \
                    f"K-mer for {bp.fusion_id} doesn't start with upstream end"
                # K-mer should contain the first part of truncated anchor
                anchor_start = sequences['Met_WT'][bp.breakpoint_nt:bp.breakpoint_nt + config.breakpoint_window]
                assert kmer.endswith(anchor_start), \
                    f"K-mer for {bp.fusion_id} doesn't end with anchor start"

    def test_detection_finds_correct_fusion(self, sequences, partners, config):
        """
        A simulated read containing a specific breakpoint should
        be detected correctly.
        """
        breakpoints = generate_all_breakpoints(sequences, partners, config)

        # Convert breakpoints to detection format
        bp_dict = {}
        for bp in breakpoints:
            if bp.partner_name not in bp_dict:
                bp_dict[bp.partner_name] = {}
            bp_dict[bp.partner_name][bp.fusion_id] = bp.sequence

        # Pick a specific breakpoint to test
        test_bp = breakpoints[len(breakpoints) // 2]  # Middle one

        # For the pre-filter to work, the domain end must be present in the read.
        # The domain end is the 3' end of the partner, which may not be in a
        # mid-truncation breakpoint k-mer. Instead, we'll use a domain end
        # that IS present in the breakpoint sequence.
        # Use the breakpoint sequence suffix as the "domain end" for this test.
        domain_ends = {test_bp.partner_name: test_bp.sequence[-config.kmer_size:]}

        # Create a "read" that contains this breakpoint
        simulated_read = "AAAA" + test_bp.sequence + "TTTT"

        # Find matches
        matches = find_matches_in_read(simulated_read, domain_ends, bp_dict)

        assert test_bp.fusion_id in matches, \
            f"Failed to detect {test_bp.fusion_id} in simulated read"


# =============================================================================
# TEST: Fusion Structure Verification
# =============================================================================

class TestFusionStructure:
    """
    Tests to verify the expected fusion structure:
    [Truncated Partner] -- [Linker] -- [Full Anchor]
    """

    def test_structure_downstream_anchor(self):
        """
        With anchor downstream, structure should be:
        Partner(truncated) + Linker + Anchor(full)
        """
        partner = "AAACCCGGG"    # 9 nt partner
        anchor = "TTTAAACCC"     # 9 nt anchor
        linker = "NNN"           # 3 nt linker

        # Truncate partner at position 6
        fusion, bp_pos = generate_fusion_sequence(
            partner_seq=partner,
            anchor_seq=anchor,
            linker_seq=linker,
            breakpoint_pos=6,
            anchor_position='downstream'
        )

        # Expected: AAACCC + NNN + TTTAAACCC
        assert fusion == "AAACCCNNNTTTAAACCC"

        # Breakpoint is after partner + linker
        assert bp_pos == 6 + len(linker)

        # Verify structure components
        assert fusion[:6] == "AAACCC"     # Truncated partner
        assert fusion[6:9] == "NNN"       # Linker
        assert fusion[9:] == anchor       # Full anchor

    def test_anchor_length_constant_across_truncations(self):
        """
        The anchor portion should be the same length regardless
        of partner truncation position.
        """
        partner = "AAACCCGGGTTTTTT"  # 15 nt
        anchor = "XXXYYYZZZAAA"       # 12 nt
        linker = "NNN"

        anchor_lengths = []
        for truncation in [3, 6, 9, 12, 15]:
            fusion, bp_pos = generate_fusion_sequence(
                partner_seq=partner,
                anchor_seq=anchor,
                linker_seq=linker,
                breakpoint_pos=truncation,
                anchor_position='downstream'
            )

            # Extract anchor portion (everything after linker)
            anchor_portion = fusion[bp_pos:]

            anchor_lengths.append(len(anchor_portion))

            # Verify it's the complete anchor
            assert anchor_portion == anchor

        # All anchor lengths should be identical
        assert len(set(anchor_lengths)) == 1

    def test_partner_length_varies_with_truncation(self):
        """Partner length should vary based on truncation position."""
        partner = "AAACCCGGGTTTTTT"  # 15 nt
        anchor = "XXX"
        linker = ""

        for truncation in [3, 6, 9, 12, 15]:
            fusion, _ = generate_fusion_sequence(
                partner_seq=partner,
                anchor_seq=anchor,
                linker_seq=linker,
                breakpoint_pos=truncation,
                anchor_position='downstream'
            )

            # Fusion length = truncation + anchor
            expected_length = truncation + len(anchor)
            assert len(fusion) == expected_length


# =============================================================================
# TEST: Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_minimal_truncation(self):
        """Should handle minimal (3 nt / 1 codon) truncation."""
        fusion, bp_pos = generate_fusion_sequence(
            partner_seq="AAACCCGGG",
            anchor_seq="TTTAAACCC",
            linker_seq="NNN",
            breakpoint_pos=3,  # Just one codon
            anchor_position='downstream'
        )

        assert fusion == "AAANNNTTTAAACCC"
        assert bp_pos == 3 + len("NNN")

    def test_full_partner_inclusion(self):
        """Should handle full partner (no truncation)."""
        partner = "AAACCCGGG"
        anchor = "TTTAAACCC"
        linker = "NNN"

        fusion, bp_pos = generate_fusion_sequence(
            partner_seq=partner,
            anchor_seq=anchor,
            linker_seq=linker,
            breakpoint_pos=len(partner),  # Full partner
            anchor_position='downstream'
        )

        expected = partner + linker + anchor
        assert fusion == expected

    def test_no_linker(self):
        """Should work without a linker."""
        fusion, bp_pos = generate_fusion_sequence(
            partner_seq="AAACCC",
            anchor_seq="TTTAAA",
            linker_seq="",
            breakpoint_pos=3,
            anchor_position='downstream'
        )

        assert fusion == "AAATTTAAA"
        assert bp_pos == 3

    def test_long_linker(self):
        """Should handle longer linkers."""
        long_linker = "GGGAGCGGGAGCGGGAGC"  # 18 nt (6x GS)

        fusion, bp_pos = generate_fusion_sequence(
            partner_seq="AAACCC",
            anchor_seq="TTTAAA",
            linker_seq=long_linker,
            breakpoint_pos=3,
            anchor_position='downstream'
        )

        expected = "AAA" + long_linker + "TTTAAA"
        assert fusion == expected
        assert bp_pos == 3 + len(long_linker)


# =============================================================================
# TEST: Detection Accuracy
# =============================================================================

class TestDetectionAccuracy:
    """Tests for detection accuracy with various read scenarios."""

    @pytest.fixture
    def setup_detection(self):
        """Set up detection components."""
        sequences = {
            'Anchor': 'TTTAAACCCGGGTTT',
            'Partner': 'AAACCCGGGAAACCC',
        }
        partners = {
            'Partner': {'include': True, 'sequence_length': 15, 'description': ''},
        }
        config = FusionLibraryConfig(
            anchor_name='Anchor',
            anchor_position='downstream',
            linker_sequence='NNN',
            breakpoint_window=6,
            maintain_frame=True,
            kmer_size=9
        )

        breakpoints = generate_all_breakpoints(sequences, partners, config)

        # For testing, create domain ends from breakpoint sequences themselves
        # This ensures the pre-filter will pass when the breakpoint is present
        # In a real scenario, the domain ends come from the full partner sequence
        bp_dict = {}
        domain_ends = {}
        for bp in breakpoints:
            if bp.partner_name not in bp_dict:
                bp_dict[bp.partner_name] = {}
            bp_dict[bp.partner_name][bp.fusion_id] = bp.sequence
            # Use a portion of the breakpoint as the domain end for testing
            if bp.partner_name not in domain_ends:
                domain_ends[bp.partner_name] = bp.sequence[-config.kmer_size:]

        return breakpoints, domain_ends, bp_dict

    def test_detects_exact_match(self, setup_detection):
        """Should detect exact breakpoint sequence."""
        breakpoints, domain_ends, bp_dict = setup_detection

        # Use exact breakpoint sequence as read
        test_bp = breakpoints[0]
        read = test_bp.sequence

        matches = find_matches_in_read(read, domain_ends, bp_dict)

        assert test_bp.fusion_id in matches

    def test_detects_with_flanking_sequence(self, setup_detection):
        """Should detect breakpoint within longer read."""
        breakpoints, domain_ends, bp_dict = setup_detection

        test_bp = breakpoints[0]
        read = "GGGGGG" + test_bp.sequence + "CCCCCC"

        matches = find_matches_in_read(read, domain_ends, bp_dict)

        assert test_bp.fusion_id in matches

    def test_no_false_positives_on_partial(self, setup_detection):
        """Should not match partial breakpoint sequences."""
        breakpoints, domain_ends, bp_dict = setup_detection

        test_bp = breakpoints[0]
        # Only half of the breakpoint sequence
        partial_read = test_bp.sequence[:len(test_bp.sequence)//2]

        matches = find_matches_in_read(partial_read, domain_ends, bp_dict)

        assert test_bp.fusion_id not in matches

    def test_no_match_on_random_sequence(self, setup_detection):
        """Should not match random sequences."""
        _, domain_ends, bp_dict = setup_detection

        random_read = "ATATATATATATATATATAT"

        matches = find_matches_in_read(random_read, domain_ends, bp_dict)

        assert len(matches) == 0


# =============================================================================
# TEST: File-Based Integration
# =============================================================================

class TestFileBasedIntegration:
    """Integration tests using actual file I/O."""

    def test_full_workflow(self, tmp_path):
        """Test complete workflow from sequences to counts."""
        # Set up test data
        sequences = {
            'Anchor': 'TTTAAACCCGGGTTTAAA',
            'Partner': 'AAACCCGGGAAACCCGGG',
        }
        partners = {
            'Partner': {'include': True, 'sequence_length': 18, 'description': ''},
        }
        config = FusionLibraryConfig(
            anchor_name='Anchor',
            anchor_position='downstream',
            linker_sequence='NNN',
            breakpoint_window=6,
            maintain_frame=True,
            kmer_size=9
        )

        # Generate breakpoints
        breakpoints = generate_all_breakpoints(sequences, partners, config)

        # Pick a breakpoint to test
        test_bp = breakpoints[len(breakpoints) // 2]

        # For detection to work, domain ends must be present in reads
        # Use the suffix of the test breakpoint sequence as domain end
        domain_ends = {'Partner': test_bp.sequence[-config.kmer_size:]}

        # Create test FASTQ with reads containing this breakpoint
        fastq_content = f"""@read1
{test_bp.sequence}AAAAAA
+
{'I' * (len(test_bp.sequence) + 6)}
@read2
{test_bp.sequence}TTTTTT
+
{'I' * (len(test_bp.sequence) + 6)}
@read3
GGGGGGGGGGGGGGGGGGGG
+
IIIIIIIIIIIIIIIIIIII
"""
        fastq_file = tmp_path / "test.fastq"
        fastq_file.write_text(fastq_content)

        # Convert breakpoints to detection format
        bp_dict = {}
        for bp in breakpoints:
            if bp.partner_name not in bp_dict:
                bp_dict[bp.partner_name] = {}
            bp_dict[bp.partner_name][bp.fusion_id] = bp.sequence

        # Count matches
        counts = count_fusion_matches(
            str(fastq_file),
            bp_dict,
            domain_ends,
            show_progress=False
        )

        # Should have 2 counts for our test breakpoint
        assert counts.get(test_bp.fusion_id, 0) == 2


# =============================================================================
# TEST: Full Run with Known Output
# =============================================================================

def test_full_run_small_dataset(tmp_path):
    """Full run from references to counts with a known expected outcome."""
    sequences = {
        "Anchor": MET_KINASE[:36],
        "Partner": TPR_PARTNER[:30],
        "Unfused": MET_KINASE[60:90],
    }

    fasta_path = tmp_path / "sequences.fasta"
    with fasta_path.open("w") as fh:
        for name, seq in sequences.items():
            fh.write(f">{name}\n{seq}\n")

    partners_path = tmp_path / "partners.csv"
    with partners_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["partner_name", "sequence_length", "include", "description"])
        writer.writerow(["Partner", len(sequences["Partner"]), "true", ""])

    unfused_path = tmp_path / "unfused.csv"
    with unfused_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sequence_name", "sequence_length", "include", "exclude_overlap", "description"])
        writer.writerow(["Unfused", len(sequences["Unfused"]), "true", "false", "control"])

    refs_dir = tmp_path / "refs"
    breakpoints_path = refs_dir / "breakpoints.csv"
    ends_path = refs_dir / "domain_ends.csv"
    unfused_kmers_path = refs_dir / "unfused_sequences.csv"
    variants_path = refs_dir / "variant_catalog.csv"
    expected_counts_path = refs_dir / "expected_counts_template.csv"

    run_generation(
        sequences_file=str(fasta_path),
        partners_file=str(partners_path),
        output_breakpoints=str(breakpoints_path),
        output_ends=str(ends_path),
        anchor_name="Anchor",
        anchor_position="downstream",
        truncated_component="anchor",
        linker_sequence="GGGAGC",
        breakpoint_window=15,
        maintain_frame=True,
        kmer_size=9,
        unfused_sequences_file=str(unfused_path),
        output_unfused=str(unfused_kmers_path),
        output_variants=str(variants_path),
        output_expected_counts=str(expected_counts_path),
        logger=None,
    )

    with breakpoints_path.open() as fh:
        reader = csv.DictReader(fh)
        bp_rows = list(reader)

    seq_counts: dict[str, int] = {}
    for row in bp_rows:
        seq = row["breakpoint_sequence"]
        seq_counts[seq] = seq_counts.get(seq, 0) + 1

    target_row = next(
        (row for row in bp_rows if seq_counts[row["breakpoint_sequence"]] == 1),
        bp_rows[0],
    )
    target_id = target_row["fusion_id"]
    target_seq = target_row["breakpoint_sequence"]

    reads = [f"AAAA{target_seq}TTTT"] * 3 + [
        "ACGT" * 10,
        "TTTTGGGGCCCCAAAATTTT",
    ]
    fastq_path = tmp_path / "reads.fastq"
    _write_fastq(fastq_path, reads)

    counts_path = tmp_path / "counts.csv"
    run_matching(
        fastq_file=str(fastq_path),
        breakpoints_file=str(breakpoints_path),
        ends_file=str(ends_path),
        unfused_file=str(unfused_kmers_path),
        output_file=str(counts_path),
        show_progress=False,
        progress_interval=50,
        logger=None,
        prefilter_fallback=False,
    )

    with counts_path.open() as fh:
        counts_rows = list(csv.DictReader(fh))

    counts = {row["fusion_id"]: int(row["count"]) for row in counts_rows}
    count_types = {row["fusion_id"]: row.get("type", "") for row in counts_rows}

    assert counts.get(target_id, 0) == 3
    nonzero = {fid for fid, count in counts.items() if count > 0}
    assert nonzero == {target_id}

    with variants_path.open() as fh:
        variants_rows = list(csv.DictReader(fh))
    with expected_counts_path.open() as fh:
        expected_rows = list(csv.DictReader(fh))

    variants_set = {(row["fusion_id"], row["type"]) for row in variants_rows}
    expected_set = {(row["fusion_id"], row["type"]) for row in expected_rows}
    counts_set = {(fid, count_types.get(fid, "")) for fid in counts.keys()}

    assert variants_set == expected_set == counts_set
    assert all(int(row["count"]) == 0 for row in expected_rows)


def test_snakemake_counts_only(tmp_path):
    """Run Snakemake on a tiny dataset and confirm expected counts."""
    if shutil.which("snakemake") is None:
        pytest.skip("snakemake not available in PATH")

    repo_root = Path(__file__).resolve().parents[1]
    snakefile = repo_root / "workflow" / "Snakefile"

    experiment = "smoke_test"
    sample = "sample1"
    data_dir = tmp_path / "data"
    ref_dir = tmp_path / "refs"
    data_dir.mkdir()
    ref_dir.mkdir()

    sequences = {
        "Anchor": MET_KINASE[:60],
        "Partner": TPR_PARTNER[:45],
        "Unfused": MET_KINASE[60:90],
    }

    fasta_path = ref_dir / "sequences.fasta"
    with fasta_path.open("w") as fh:
        for name, seq in sequences.items():
            fh.write(f">{name}\n{seq}\n")

    partners_path = tmp_path / "partners.csv"
    with partners_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["partner_name", "sequence_length", "include", "description"])
        writer.writerow(["Partner", len(sequences["Partner"]), "true", ""])

    unfused_path = tmp_path / "unfused.csv"
    with unfused_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sequence_name", "sequence_length", "include", "exclude_overlap", "description"])
        writer.writerow(["Unfused", len(sequences["Unfused"]), "true", "false", "control"])

    samples_path = tmp_path / "samples.csv"
    with samples_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample", "condition", "file"])
        writer.writerow([sample, "baseline", "S1"])

    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "\n".join([
            f"experiment: '{experiment}'",
            f"data_dir: '{data_dir}'",
            f"ref_dir: '{ref_dir.name}'",
            f"samples_file: '{samples_path.name}'",
            "fusion_library:",
            "  anchor:",
            "    name: 'Anchor'",
            "    position: 'downstream'",
            "    truncated_component: 'anchor'",
            "  linker_sequence: 'GGGAGC'",
            f"  partners_file: '{partners_path.name}'",
            f"  sequences_file: '{fasta_path.name}'",
            f"  unfused_sequences_file: '{unfused_path.name}'",
            "detection:",
            "  method: 'string'",
            "  breakpoint_window: 14",
            "  maintain_frame: true",
            "  kmer_size: 8",
            "  orientation_check: false",
            "  prefilter_fallback: false",
            "sequencing:",
            "  paired: true",
            "qc:",
            "  run_qc: false",
            "mem_fastqc: 1000",
            "pipeline:",
            "  show_progress: false",
            "resources:",
            "  threads: 1",
        ])
    )

    config = FusionLibraryConfig(
        anchor_name="Anchor",
        anchor_position="downstream",
        truncated_component="anchor",
        linker_sequence="GGGAGC",
        breakpoint_window=14,
        maintain_frame=True,
        kmer_size=8,
    )
    partners = {"Partner": {"include": True, "sequence_length": len(sequences["Partner"]), "description": ""}}
    breakpoints = generate_all_breakpoints(sequences, partners, config)

    seq_counts: dict[str, int] = {}
    for bp in breakpoints:
        seq_counts[bp.sequence] = seq_counts.get(bp.sequence, 0) + 1
    target_bp = next((bp for bp in breakpoints if seq_counts[bp.sequence] == 1), breakpoints[0])

    merged_path = tmp_path / "results" / experiment / "merged" / f"{sample}_merged.fastq.gz"
    merged_path.parent.mkdir(parents=True, exist_ok=True)
    _write_fastq_gz(merged_path, [f"AAAA{target_bp.sequence}TTTT"] * 2)

    env = os.environ.copy()
    env["XDG_CACHE_HOME"] = str(tmp_path / "cache")

    result = subprocess.run(
        [
            "snakemake",
            "-s",
            str(snakefile),
            "--configfile",
            str(config_path),
            "--directory",
            str(tmp_path),
            "--cores",
            "1",
            "--shared-fs-usage",
            "none",
            "--",
            "counts_only",
        ],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr

    counts_path = tmp_path / "results" / experiment / "counts" / f"{sample}.fusion_counts.csv"
    assert counts_path.exists()

    with counts_path.open() as fh:
        rows = list(csv.DictReader(fh))
    counts = {row["fusion_id"]: int(row["count"]) for row in rows}
    nonzero = {fid for fid, count in counts.items() if count > 0}

    assert counts.get(target_bp.fusion_id, 0) == 2
    assert nonzero == {target_bp.fusion_id}

    variants_path = tmp_path / "results" / experiment / "references" / "variant_catalog.csv"
    expected_counts_path = tmp_path / "results" / experiment / "references" / "expected_counts_template.csv"
    assert variants_path.exists()
    assert expected_counts_path.exists()


# =============================================================================
# RUN TESTS
# =============================================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
