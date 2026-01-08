#!/usr/bin/env python3
"""
Synthetic-read detection test.

This test generates breakpoint references and a tiny synthetic FASTQ on the fly,
then runs the string-matching detector to verify counts match expectations. It
serves as a true regression test (not just a smoke run) for the detection
pipeline on controlled data.
"""

import csv
import gzip
from pathlib import Path

from workflow.scripts.fusion_sequences import (
    FusionLibraryConfig,
    generate_all_breakpoints,
    generate_domain_ends,
    write_breakpoints_csv,
    write_domain_ends_csv,
)
from workflow.scripts.string_matcher import run_matching


# Reuse realistic sequences from integration tests (shortened for speed)
MET_KINASE = (
    "ATGAAAAAGAGAAAGCAAATTAAAGATCTGGGCAGTGAATTAGTTCGCTACGATGCAAGA"
    "GTACACACTCCTCATTTGGATAGGCTTGTAAGTGCCCGAAGTGTAAGCCCAACTACAGAA"
    "ATGGTTTCAAATGAATCTGTAGACTACCGAGCTACTTTTCCAGAAGATCAGTTTCCTAAT"
    "TCATCTCAGAACGGTTCATGCCGACAAGTGCAGTATCCTCTGACAGACATGTCCCCCATC"
)

TPR_PARTNER = (
    "ATGGCGGCGGTGTTGCAGCAAGTCCTGGAGCGCACGGAGCTGAACAAGCTGCCCAAGTCT"
    "GTCCAGAACAAACTTGAAAAGTTCCTTGCTGATCAGCAATCCGAGATCGATGGCCTGAAG"
    "GGGCGGCATGAGAAATTTAAGGTGGAGAGCGAACAACAGTATTTTGAAATAGAAAAGAGG"
)


def _write_fastq(path: Path, reads: list[str]) -> Path:
    """Write a gzipped FASTQ with constant qualities."""
    with gzip.open(path, "wt") as fh:
        for i, seq in enumerate(reads):
            fh.write(f"@read{i}\n")
            fh.write(f"{seq}\n")
            fh.write("+\n")
            fh.write(f"{'I' * len(seq)}\n")
    return path


def test_detects_synthetic_breakpoint(tmp_path: Path):
    # Config: keep kmer_size == breakpoint_window so domain end is present in k-mer
    cfg = FusionLibraryConfig(
        anchor_name="Met_WT",
        anchor_position="downstream",
        linker_sequence="GGGAGC",
        breakpoint_window=12,
        maintain_frame=True,
        kmer_size=12,
    )

    sequences = {"Met_WT": MET_KINASE, "TPR": TPR_PARTNER}
    partners = {"TPR": {"include": True, "sequence_length": len(TPR_PARTNER), "description": ""}}

    # Generate references
    breakpoints = generate_all_breakpoints(sequences, partners, cfg)
    domain_ends = generate_domain_ends(sequences, partners, cfg.anchor_name, kmer_size=cfg.kmer_size)

    # Choose a representative breakpoint
    target_bp = breakpoints[len(breakpoints) // 3]

    # Paths
    bp_csv = tmp_path / "breakpoints.csv"
    ends_csv = tmp_path / "domain_ends.csv"
    fastq_path = tmp_path / "synthetic.fastq.gz"
    counts_path = tmp_path / "counts.csv"

    # Write reference CSVs
    write_breakpoints_csv(breakpoints, bp_csv)
    write_domain_ends_csv(domain_ends, ends_csv)

    # Build synthetic reads: multiple positives + negatives
    positive_read = f"AAAA{target_bp.sequence}TTTT"
    reads = [positive_read for _ in range(5)] + [
        "ACGT" * 10,  # negative read
        "TTTTGGGGCCCCAAAATTTT",  # negative read
    ]
    _write_fastq(fastq_path, reads)

    # Run detection
    run_matching(
        fastq_file=str(fastq_path),
        breakpoints_file=str(bp_csv),
        ends_file=str(ends_csv),
        unfused_file=None,
        output_file=str(counts_path),
        show_progress=False,
        progress_interval=50,
        logger=None,
        prefilter_fallback=True,
    )

    # Validate counts
    with counts_path.open() as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)

    fusion_ids = [row["fusion_id"] for row in rows]
    counts = {row["fusion_id"]: int(row["count"]) for row in rows}

    assert target_bp.fusion_id in fusion_ids
    assert counts.get(target_bp.fusion_id, 0) == 5, "Expected exact count for synthetic positives"

    # Ensure no false positives beyond the target breakpoint
    nonzero = {fid for fid, count in counts.items() if count > 0}
    assert nonzero == {target_bp.fusion_id}, "Unexpected extra fusion detections"
