#!/usr/bin/env python3
"""
Detect fusion breakpoints via string matching in FASTQ files.

This script searches for fusion breakpoint sequences in sequencing reads
using a two-stage approach:

1. PRE-FILTER: Check if read contains any partner domain 3' end k-mer
2. MATCH: If pre-filter passes, search for specific breakpoint sequences

USAGE:
    # Via Snakemake (automatic parameter passing)
    script: "scripts/string_matcher.py"

    # Standalone
    python string_matcher.py \
        -i reads.fastq.gz \
        -b breakpoint_sequences.csv \
        -e domain_ends.csv \
        -o fusion_counts.csv \
        --progress
"""

import argparse
import csv
import gzip
import os
from collections import defaultdict
from pathlib import Path
from typing import Iterator

# Handle both Snakemake and standalone execution
try:
    from utils import setup_logging, ProgressReporter
except ImportError:
    try:
        from workflow.scripts.utils import setup_logging, ProgressReporter
    except ImportError:
        # Fallback for direct execution
        import logging
        def setup_logging(log_file=None, level='INFO', name='fusilli'):
            logger = logging.getLogger(name)
            logger.setLevel(getattr(logging, level))
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter('%(asctime)s | %(levelname)s | %(message)s'))
            logger.addHandler(handler)
            return logger

        class ProgressReporter:
            def __init__(self, total, desc='', interval_pct=1, enabled=True, logger=None):
                self.total = total
                self.count = 0
                self.enabled = enabled
            def update(self, n=1):
                self.count += n
            def finish(self):
                pass


# Try to import pyfastx for fast FASTQ parsing
try:
    import pyfastx
    HAS_PYFASTX = True
except ImportError:
    HAS_PYFASTX = False


# =============================================================================
# FASTQ PARSING
# =============================================================================

def parse_fastq_pyfastx(filepath: str) -> Iterator[tuple[str, str, str]]:
    """
    Parse FASTQ using pyfastx (fast C library).

    Yields:
        Tuples of (name, sequence, quality)
    """
    fq = pyfastx.Fastx(filepath, uppercase=True)
    for name, seq, qual in fq:
        yield name, seq, qual


def parse_fastq_python(filepath: str) -> Iterator[tuple[str, str, str]]:
    """
    Parse FASTQ using pure Python (fallback).

    Yields:
        Tuples of (name, sequence, quality)
    """
    opener = gzip.open if filepath.endswith('.gz') else open

    with opener(filepath, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break

            seq = f.readline().strip().upper()
            _ = f.readline()
            qual = f.readline().strip()

            name = header[1:].split()[0]
            yield name, seq, qual


def parse_fastq(filepath: str) -> Iterator[tuple[str, str, str]]:
    """
    Parse FASTQ file, using fastest available method.
    """
    if HAS_PYFASTX:
        yield from parse_fastq_pyfastx(filepath)
    else:
        yield from parse_fastq_python(filepath)


def estimate_read_count(filepath: str) -> int:
    """Estimate number of reads in FASTQ file based on file size."""
    file_size = os.path.getsize(filepath)

    if filepath.endswith('.gz'):
        return max(1, file_size // 60)
    else:
        return max(1, file_size // 400)


# =============================================================================
# DATA LOADING
# =============================================================================

def load_breakpoint_sequences(filepath: str) -> dict[str, dict]:
    """
    Load breakpoint sequences from CSV, organized by partner for fast lookup.
    """
    breakpoints = defaultdict(dict)

    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            partner = row['partner_name']
            fusion_id = row['fusion_id']
            sequence = row['breakpoint_sequence']
            breakpoints[partner][fusion_id] = sequence

    return dict(breakpoints)


def load_domain_ends(filepath: str) -> dict[str, str]:
    """
    Load domain end k-mers from CSV.
    """
    ends = {}

    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            ends[row['domain_name']] = row['end_kmer']

    return ends


# =============================================================================
# UNFUSED LOADING
# =============================================================================

def load_unfused_kmers(filepath: str | None) -> dict[int, dict[str, list[str]]]:
    """
    Load unfused sequence k-mers from CSV, grouped by k-mer length.

    Returns:
        Dict of {kmer_length: {kmer: [sequence_name, ...]}}
    """
    if not filepath:
        return {}

    path = Path(filepath)
    if not path.exists():
        return {}

    kmers_by_len: dict[int, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))

    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            return {}
        for row in reader:
            seq_name = row.get('sequence_name', '').strip()
            kmer = row.get('kmer', '').strip()
            if not seq_name or not kmer:
                continue
            kmers_by_len[len(kmer)][kmer].append(seq_name)

    return {k: dict(v) for k, v in kmers_by_len.items()}


def load_unfused_sequence_names(filepath: str | None) -> list[str]:
    """Load unique unfused sequence names from a k-mer CSV file."""
    if not filepath:
        return []

    path = Path(filepath)
    if not path.exists():
        return []

    sequence_names: set[str] = set()
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            return []
        for row in reader:
            seq_name = row.get('sequence_name', '').strip()
            if seq_name:
                sequence_names.add(seq_name)

    return sorted(sequence_names)


def collect_expected_fusion_ids(breakpoints: dict[str, dict]) -> list[str]:
    """Collect all fusion IDs from breakpoint definitions."""
    fusion_ids: set[str] = set()
    for bp_map in breakpoints.values():
        fusion_ids.update(bp_map.keys())
    return sorted(fusion_ids)


# =============================================================================
# MATCHING
# =============================================================================

def reverse_complement(seq: str) -> str:
    """Return reverse complement of a nucleotide sequence."""
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


def find_matches_in_read(
    sequence: str,
    domain_ends: dict[str, str],
    breakpoints: dict[str, dict],
    orientation_check: bool = False,
    rc_domain_ends: dict[str, str] | None = None,
    rc_breakpoints: dict[str, dict] | None = None,
    return_orientation: bool = False,
    prefilter_fallback: bool = False,
) -> list[str] | tuple[list[str], bool, bool]:
    """
    Find all fusion breakpoint matches in a single read.

    Two-stage approach:
    1. Check if any domain end k-mer is present (fast pre-filter)
    2. If found, search for specific breakpoint sequences
    """
    matches: list[str] = []
    forward_hit = False
    rc_hit = False

    matched_domains = []
    for domain_name, end_kmer in domain_ends.items():
        if end_kmer in sequence:
            matched_domains.append(domain_name)

    if orientation_check and rc_domain_ends:
        rc_seq = reverse_complement(sequence)
        for domain_name, end_kmer in rc_domain_ends.items():
            if end_kmer in rc_seq:
                matched_domains.append(domain_name)
                rc_hit = True

    if matched_domains:
        forward_hit = not rc_hit
    else:
        if prefilter_fallback:
            for domain, bp_map in breakpoints.items():
                for fusion_id, bp_sequence in bp_map.items():
                    if bp_sequence in sequence:
                        matches.append(fusion_id)
                        forward_hit = True

                if orientation_check and rc_breakpoints:
                    rc_seq = rc_seq if "rc_seq" in locals() else reverse_complement(sequence)
                    for fusion_id, bp_sequence in rc_breakpoints.get(domain, {}).items():
                        if bp_sequence in rc_seq:
                            matches.append(fusion_id)
                            rc_hit = True

        if return_orientation:
            return matches, forward_hit, rc_hit
        return matches

    for domain in matched_domains:
        if domain not in breakpoints:
            continue

        for fusion_id, bp_sequence in breakpoints[domain].items():
            if bp_sequence in sequence:
                matches.append(fusion_id)
                forward_hit = True

        if orientation_check and rc_breakpoints:
            rc_seq = rc_seq if "rc_seq" in locals() else reverse_complement(sequence)
            for fusion_id, bp_sequence in rc_breakpoints.get(domain, {}).items():
                if bp_sequence in rc_seq:
                    matches.append(fusion_id)
                    rc_hit = True

    if return_orientation:
        return matches, forward_hit, rc_hit
    return matches


def find_unfused_matches_in_read(
    sequence: str,
    unfused_kmers_by_len: dict[int, dict[str, list[str]]],
    orientation_check: bool = False,
    return_orientation: bool = False
) -> set[str] | tuple[set[str], bool, bool]:
    """
    Find unfused sequence matches in a single read.

    Returns a set of unfused sequence names with at least one k-mer hit.
    """
    matches: set[str] = set()
    forward_hit = False
    rc_hit = False

    if not unfused_kmers_by_len:
        if return_orientation:
            return matches, forward_hit, rc_hit
        return matches

    sequences_to_check = [(sequence, False)]
    if orientation_check:
        sequences_to_check.append((reverse_complement(sequence), True))

    for seq, is_rc in sequences_to_check:
        for kmer_len, kmer_map in unfused_kmers_by_len.items():
            if len(seq) < kmer_len:
                continue
            for i in range(0, len(seq) - kmer_len + 1):
                window = seq[i:i + kmer_len]
                if window in kmer_map:
                    for seq_name in kmer_map[window]:
                        matches.add(seq_name)
                    if is_rc:
                        rc_hit = True
                    else:
                        forward_hit = True

        # If we've matched everything possible in this orientation, continue
        # to allow rc matches to contribute orientation metrics.

    if return_orientation:
        return matches, forward_hit, rc_hit
    return matches


def count_fusion_matches(
    fastq_file: str,
    breakpoints: dict[str, dict],
    domain_ends: dict[str, str],
    show_progress: bool = True,
    progress_interval: float = 1.0,
    logger=None,
    orientation_check: bool = False,
    return_metrics: bool = False,
    prefilter_fallback: bool = False,
) -> dict[str, int] | tuple[dict[str, int], dict[str, int]]:
    """
    Count fusion breakpoint matches in a FASTQ file.
    """
    counts = defaultdict(int)
    metrics = {
        "reads_processed": 0,
        "prefilter_pass_reads": 0,
        "matched_reads": 0,
        "total_matches": 0,
        "unique_fusions_detected": 0,
        "forward_match_events": 0,
        "rc_match_events": 0,
    }

    estimated_total = estimate_read_count(fastq_file)

    if logger:
        logger.info(f"Processing: {fastq_file}")
        logger.info(f"Estimated reads: {estimated_total:,}")
        logger.info(f"Partners to search: {len(domain_ends)}")
        total_breakpoints = sum(len(bps) for bps in breakpoints.values())
        logger.info(f"Breakpoint sequences: {total_breakpoints:,}")

    rc_domain_ends = None
    rc_breakpoints = None
    if orientation_check:
        rc_domain_ends = {k: reverse_complement(v) for k, v in domain_ends.items()}
        rc_breakpoints = {
            partner: {fid: reverse_complement(seq) for fid, seq in bp_dict.items()}
            for partner, bp_dict in breakpoints.items()
        }

    progress = ProgressReporter(
        total=estimated_total,
        desc="Matching reads",
        interval_pct=progress_interval,
        enabled=show_progress,
        logger=logger
    )

    read_count = 0
    match_count = 0

    for name, seq, qual in parse_fastq(fastq_file):
        read_count += 1
        metrics["reads_processed"] = read_count

        matches, f_hit, rc_hit = find_matches_in_read(
            seq,
            domain_ends,
            breakpoints,
            orientation_check=orientation_check,
            rc_domain_ends=rc_domain_ends,
            rc_breakpoints=rc_breakpoints,
            return_orientation=True,
            prefilter_fallback=prefilter_fallback,
        )

        if matches:
            match_count += len(matches)
            metrics["prefilter_pass_reads"] += 1
            metrics["matched_reads"] += 1
            for fusion_id in matches:
                counts[fusion_id] += 1
            metrics["forward_match_events"] += int(f_hit) * len(matches)
            metrics["rc_match_events"] += int(rc_hit) * len(matches)

        progress.update()

    progress.finish()

    metrics["total_matches"] = match_count
    metrics["unique_fusions_detected"] = len(counts)

    if logger:
        logger.info(f"Processed {read_count:,} reads")
        logger.info(f"Found {match_count:,} total matches")
        logger.info(f"Unique fusions detected: {len(counts)}")

    if return_metrics:
        return dict(counts), metrics
    return dict(counts)


def count_all_matches(
    fastq_file: str,
    breakpoints: dict[str, dict],
    domain_ends: dict[str, str],
    unfused_kmers_by_len: dict[int, dict[str, list[str]]] | None,
    show_progress: bool = True,
    progress_interval: float = 1.0,
    logger=None,
    orientation_check: bool = False,
    prefilter_fallback: bool = False,
) -> tuple[dict[str, int], dict[str, int], dict[str, int]]:
    """
    Count fusion breakpoint matches and unfused sequence matches in a FASTQ file.
    """
    fusion_counts = defaultdict(int)
    unfused_counts = defaultdict(int)
    metrics = {
        "reads_processed": 0,
        "prefilter_pass_reads": 0,
        "matched_reads": 0,
        "total_matches": 0,
        "unique_fusions_detected": 0,
        "forward_match_events": 0,
        "rc_match_events": 0,
        "unfused_reads": 0,
        "unfused_total_matches": 0,
        "unique_unfused_detected": 0,
    }

    estimated_total = estimate_read_count(fastq_file)

    if logger:
        logger.info(f"Processing: {fastq_file}")
        logger.info(f"Estimated reads: {estimated_total:,}")
        logger.info(f"Partners to search: {len(domain_ends)}")
        total_breakpoints = sum(len(bps) for bps in breakpoints.values())
        logger.info(f"Breakpoint sequences: {total_breakpoints:,}")
        unfused_total = sum(len(kmers) for kmers in (unfused_kmers_by_len or {}).values())
        logger.info(f"Unfused k-mers: {unfused_total:,}")

    rc_domain_ends = None
    rc_breakpoints = None
    if orientation_check:
        rc_domain_ends = {k: reverse_complement(v) for k, v in domain_ends.items()}
        rc_breakpoints = {
            partner: {fid: reverse_complement(seq) for fid, seq in bp_dict.items()}
            for partner, bp_dict in breakpoints.items()
        }

    progress = ProgressReporter(
        total=estimated_total,
        desc="Matching reads",
        interval_pct=progress_interval,
        enabled=show_progress,
        logger=logger
    )

    read_count = 0
    match_count = 0
    unfused_match_count = 0

    for _, seq, _ in parse_fastq(fastq_file):
        read_count += 1
        metrics["reads_processed"] = read_count

        matches, f_hit, rc_hit = find_matches_in_read(
            seq,
            domain_ends,
            breakpoints,
            orientation_check=orientation_check,
            rc_domain_ends=rc_domain_ends,
            rc_breakpoints=rc_breakpoints,
            return_orientation=True,
            prefilter_fallback=prefilter_fallback,
        )

        if matches:
            match_count += len(matches)
            metrics["prefilter_pass_reads"] += 1
            metrics["matched_reads"] += 1
            for fusion_id in matches:
                fusion_counts[fusion_id] += 1
            metrics["forward_match_events"] += int(f_hit) * len(matches)
            metrics["rc_match_events"] += int(rc_hit) * len(matches)

        if unfused_kmers_by_len:
            unfused_matches, uf_hit, uf_rc_hit = find_unfused_matches_in_read(
                seq,
                unfused_kmers_by_len,
                orientation_check=orientation_check,
                return_orientation=True
            )
            if unfused_matches:
                unfused_match_count += len(unfused_matches)
                metrics["unfused_reads"] += 1
                for seq_name in unfused_matches:
                    unfused_counts[seq_name] += 1
                metrics["forward_match_events"] += int(uf_hit) * len(unfused_matches)
                metrics["rc_match_events"] += int(uf_rc_hit) * len(unfused_matches)

        progress.update()

    progress.finish()

    metrics["total_matches"] = match_count
    metrics["unique_fusions_detected"] = len(fusion_counts)
    metrics["unfused_total_matches"] = unfused_match_count
    metrics["unique_unfused_detected"] = len(unfused_counts)

    if logger:
        logger.info(f"Processed {read_count:,} reads")
        logger.info(f"Found {match_count:,} total fusion matches")
        logger.info(f"Unique fusions detected: {len(fusion_counts)}")
        if unfused_kmers_by_len:
            logger.info(f"Found {unfused_match_count:,} unfused matches")
            logger.info(f"Unique unfused detected: {len(unfused_counts)}")

    return dict(fusion_counts), dict(unfused_counts), metrics


# =============================================================================
# OUTPUT
# =============================================================================

def write_counts_csv(
    counts: dict[str, int],
    filepath: str,
    unfused_counts: dict[str, int] | None = None,
    include_type: bool = False,
    expected_fusions: list[str] | None = None,
    expected_unfused: list[str] | None = None,
) -> None:
    """
    Write fusion counts to CSV file.
    """
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)

    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        if include_type:
            writer.writerow(['fusion_id', 'type', 'count'])
        else:
            writer.writerow(['fusion_id', 'count'])

        def write_rows(
            row_counts: dict[str, int],
            row_type: str | None = None,
            expected_ids: list[str] | None = None
        ) -> None:
            merged_counts = dict(row_counts)
            if expected_ids:
                for expected_id in expected_ids:
                    merged_counts.setdefault(expected_id, 0)
            for fusion_id, count in sorted(
                merged_counts.items(),
                key=lambda x: (-x[1], x[0])
            ):
                if include_type and row_type:
                    writer.writerow([fusion_id, row_type, count])
                else:
                    writer.writerow([fusion_id, count])

        write_rows(
            counts,
            'fusion' if include_type else None,
            expected_ids=expected_fusions
        )
        if unfused_counts is not None:
            write_rows(
                unfused_counts,
                'unfused',
                expected_ids=expected_unfused
            )


def write_metrics_json(metrics: dict, filepath: str) -> None:
    """Write metrics dictionary to JSON."""
    import json

    Path(filepath).parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, 'w') as f:
        json.dump(metrics, f, indent=2)


# =============================================================================
# MAIN
# =============================================================================

def main_snakemake(snakemake) -> None:
    """Entry point when called from Snakemake."""
    fastq_file = snakemake.input.fastq
    breakpoints_file = snakemake.input.breakpoints
    ends_file = snakemake.input.ends
    unfused_file = getattr(snakemake.input, "unfused", None)

    output_file = snakemake.output.counts
    metrics_file = getattr(snakemake.output, "metrics", None)

    show_progress = snakemake.params.get('show_progress', True)
    progress_interval = snakemake.params.get('progress_interval', 1)
    orientation_check = snakemake.params.get('orientation_check', False)
    prefilter_fallback = snakemake.params.get('prefilter_fallback', False)

    log_file = snakemake.log[0] if snakemake.log else None
    logger = setup_logging(log_file, name='string_matcher')

    run_matching(
        fastq_file=fastq_file,
        breakpoints_file=breakpoints_file,
        ends_file=ends_file,
        output_file=output_file,
        show_progress=show_progress,
        progress_interval=progress_interval,
        logger=logger,
        metrics_file=metrics_file,
        orientation_check=orientation_check,
        prefilter_fallback=prefilter_fallback,
        unfused_file=unfused_file,
    )


def main_cli() -> None:
    """Entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description='Detect fusion breakpoints via string matching',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('-i', '--input', required=True, help='Input FASTQ file (gzipped or plain)')
    parser.add_argument('-b', '--breakpoints', required=True, help='Breakpoint sequences CSV file')
    parser.add_argument('-e', '--ends', required=True, help='Domain ends CSV file')
    parser.add_argument('-o', '--output', required=True, help='Output counts CSV file')
    parser.add_argument('--metrics', help='Optional metrics JSON output')
    parser.add_argument('-u', '--unfused', help='Optional unfused sequences CSV file')
    parser.add_argument('--progress', action='store_true', default=True, help='Show progress (default: on)')
    parser.add_argument('--no-progress', action='store_true', help='Disable progress display')
    parser.add_argument('--progress-interval', type=float, default=1.0, help='Progress update interval in % (default: 1)')
    parser.add_argument('--orientation-check', action='store_true', help='Also search reverse complement to gauge orientation issues')
    parser.add_argument('--log', help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')

    args = parser.parse_args()

    level = 'DEBUG' if args.verbose else 'INFO'
    logger = setup_logging(args.log, level=level, name='string_matcher')

    show_progress = args.progress and not args.no_progress

    run_matching(
        fastq_file=args.input,
        breakpoints_file=args.breakpoints,
        ends_file=args.ends,
        unfused_file=args.unfused,
        output_file=args.output,
        show_progress=show_progress,
        progress_interval=args.progress_interval,
        logger=logger,
        metrics_file=args.metrics,
        orientation_check=args.orientation_check,
    )


def run_matching(
    fastq_file: str,
    breakpoints_file: str,
    ends_file: str,
    unfused_file: str | None,
    output_file: str,
    show_progress: bool,
    progress_interval: float,
    logger=None,
    metrics_file: str | None = None,
    orientation_check: bool = False,
    prefilter_fallback: bool = False,
) -> None:
    """Core matching logic."""
    estimated_reads = estimate_read_count(fastq_file)
    if prefilter_fallback and logger and estimated_reads >= 1_000_000:
        logger.warning(
            "prefilter_fallback is enabled for a large input; this can be extremely slow. "
            "Consider disabling it for production runs."
        )
    if logger:
        logger.info(f"Loading breakpoint sequences from: {breakpoints_file}")
    breakpoints = load_breakpoint_sequences(breakpoints_file)

    if logger:
        logger.info(f"Loading domain ends from: {ends_file}")
    domain_ends = load_domain_ends(ends_file)

    unfused_kmers_by_len = load_unfused_kmers(unfused_file)
    include_type = bool(unfused_file)
    expected_fusions = collect_expected_fusion_ids(breakpoints)
    expected_unfused = (
        load_unfused_sequence_names(unfused_file) if include_type else None
    )

    if not HAS_PYFASTX and logger:
        logger.warning(
            "pyfastx not available - using slower pure Python FASTQ parser. "
            "Install pyfastx for better performance: pip install pyfastx"
        )

    fusion_counts, unfused_counts, metrics = count_all_matches(
        fastq_file=fastq_file,
        breakpoints=breakpoints,
        domain_ends=domain_ends,
        unfused_kmers_by_len=unfused_kmers_by_len,
        show_progress=show_progress,
        progress_interval=progress_interval,
        logger=logger,
        orientation_check=orientation_check,
        prefilter_fallback=prefilter_fallback,
    )

    if logger:
        logger.info(f"Writing counts to: {output_file}")
    write_counts_csv(
        fusion_counts,
        output_file,
        unfused_counts=unfused_counts if include_type else None,
        include_type=include_type,
        expected_fusions=expected_fusions,
        expected_unfused=expected_unfused
    )

    if metrics_file:
        if logger:
            logger.info(f"Writing metrics to: {metrics_file}")
        metrics["orientation_check_enabled"] = orientation_check
        if fusion_counts:
            top_fusion_id, top_count = max(fusion_counts.items(), key=lambda x: x[1])
            metrics["top_fusion_id"] = top_fusion_id
            metrics["top_fusion_count"] = int(top_count)
            metrics["top_fusion_frac_reads"] = (
                top_count / metrics["matched_reads"] if metrics["matched_reads"] else 0.0
            )
        write_metrics_json(metrics, metrics_file)

    if logger:
        logger.info("Done!")


if __name__ == '__main__':
    try:
        snakemake
        main_snakemake(snakemake)
    except NameError:
        main_cli()
