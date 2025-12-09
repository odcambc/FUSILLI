#!/usr/bin/env python3
"""
Detect fusion breakpoints via string matching in FASTQ files.

This script searches for fusion breakpoint sequences in sequencing reads
using an optimized two-stage approach:

1. PRE-FILTER: Check if read contains any partner domain 3' end k-mer
2. MATCH: If pre-filter passes, search for specific breakpoint sequences

This approach is much faster than searching all breakpoints against all reads.

USAGE:
    # Via Snakemake (automatic parameter passing)
    script: "scripts/string_matcher.py"

    # Standalone
    python string_matcher.py \\
        -i reads.fastq.gz \\
        -b breakpoint_sequences.csv \\
        -e domain_ends.csv \\
        -o fusion_counts.csv \\
        --progress
"""

import argparse
import csv
import gzip
import os
import sys
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
            # Read 4 lines at a time
            header = f.readline().strip()
            if not header:
                break

            seq = f.readline().strip().upper()
            _ = f.readline()  # + line
            qual = f.readline().strip()

            # Extract name from header (remove @)
            name = header[1:].split()[0]

            yield name, seq, qual


def parse_fastq(filepath: str) -> Iterator[tuple[str, str, str]]:
    """
    Parse FASTQ file, using fastest available method.

    Yields:
        Tuples of (name, sequence, quality)
    """
    if HAS_PYFASTX:
        yield from parse_fastq_pyfastx(filepath)
    else:
        yield from parse_fastq_python(filepath)


def estimate_read_count(filepath: str) -> int:
    """
    Estimate number of reads in FASTQ file based on file size.

    Args:
        filepath: Path to FASTQ file

    Returns:
        Estimated read count
    """
    file_size = os.path.getsize(filepath)

    if filepath.endswith('.gz'):
        # Compressed: assume ~60 bytes per read after compression
        return file_size // 60
    else:
        # Uncompressed: assume ~400 bytes per read (4 lines × ~100 bytes)
        return file_size // 400


# =============================================================================
# DATA LOADING
# =============================================================================

def load_breakpoint_sequences(filepath: str) -> dict[str, dict]:
    """
    Load breakpoint sequences from CSV, organized by partner for fast lookup.

    Args:
        filepath: Path to breakpoint sequences CSV

    Returns:
        Nested dict: {partner_name: {fusion_id: sequence}}
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

    Args:
        filepath: Path to domain ends CSV

    Returns:
        Dict mapping domain names to end k-mers
    """
    ends = {}

    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            ends[row['domain_name']] = row['end_kmer']

    return ends


# =============================================================================
# MATCHING
# =============================================================================

def find_matches_in_read(
    sequence: str,
    domain_ends: dict[str, str],
    breakpoints: dict[str, dict]
) -> list[str]:
    """
    Find all fusion breakpoint matches in a single read.

    Uses two-stage approach:
    1. Check if any domain end k-mer is present (fast pre-filter)
    2. If found, search for specific breakpoint sequences

    Args:
        sequence: Read sequence
        domain_ends: Dict of domain -> end k-mer
        breakpoints: Nested dict of partner -> fusion_id -> breakpoint sequence

    Returns:
        List of matching fusion IDs
    """
    matches = []

    # Stage 1: Pre-filter by domain ends
    # Using 'in' operator is faster than .find() for presence check
    matched_domains = []
    for domain_name, end_kmer in domain_ends.items():
        if end_kmer in sequence:
            matched_domains.append(domain_name)

    if not matched_domains:
        return matches

    # Stage 2: Search for specific breakpoints
    for domain in matched_domains:
        if domain not in breakpoints:
            continue

        for fusion_id, bp_sequence in breakpoints[domain].items():
            if bp_sequence in sequence:
                matches.append(fusion_id)
                # Don't break - a read might match multiple breakpoints
                # (though this is rare and might indicate ambiguity)

    return matches


def count_fusion_matches(
    fastq_file: str,
    breakpoints: dict[str, dict],
    domain_ends: dict[str, str],
    show_progress: bool = True,
    progress_interval: float = 1.0,
    logger=None
) -> dict[str, int]:
    """
    Count fusion breakpoint matches in a FASTQ file.

    Args:
        fastq_file: Path to FASTQ file
        breakpoints: Breakpoint sequences organized by partner
        domain_ends: Domain end k-mers for pre-filtering
        show_progress: Display progress updates
        progress_interval: Progress update interval (percentage)
        logger: Logger instance

    Returns:
        Dict mapping fusion IDs to counts
    """
    counts = defaultdict(int)

    # Estimate total reads for progress
    estimated_total = estimate_read_count(fastq_file)

    if logger:
        logger.info(f"Processing: {fastq_file}")
        logger.info(f"Estimated reads: {estimated_total:,}")
        logger.info(f"Partners to search: {len(domain_ends)}")
        total_breakpoints = sum(len(bps) for bps in breakpoints.values())
        logger.info(f"Breakpoint sequences: {total_breakpoints:,}")

    # Pre-cache for faster lookup
    domain_ends_items = list(domain_ends.items())

    # Progress reporter
    progress = ProgressReporter(
        total=estimated_total,
        desc="Matching reads",
        interval_pct=progress_interval,
        enabled=show_progress,
        logger=logger
    )

    # Process reads
    read_count = 0
    match_count = 0

    for name, seq, qual in parse_fastq(fastq_file):
        read_count += 1

        # Find matches
        matches = find_matches_in_read(seq, domain_ends, breakpoints)

        if matches:
            match_count += len(matches)
            for fusion_id in matches:
                counts[fusion_id] += 1

        progress.update()

    progress.finish()

    if logger:
        logger.info(f"Processed {read_count:,} reads")
        logger.info(f"Found {match_count:,} total matches")
        logger.info(f"Unique fusions detected: {len(counts)}")

    return dict(counts)


# =============================================================================
# OUTPUT
# =============================================================================

def write_counts_csv(counts: dict[str, int], filepath: str) -> None:
    """
    Write fusion counts to CSV file.

    Args:
        counts: Dict mapping fusion IDs to counts
        filepath: Output file path
    """
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)

    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['fusion_id', 'count'])

        # Sort by count (descending) then by name
        for fusion_id, count in sorted(
            counts.items(),
            key=lambda x: (-x[1], x[0])
        ):
            writer.writerow([fusion_id, count])


# =============================================================================
# MAIN
# =============================================================================

def main_snakemake(snakemake) -> None:
    """Entry point when called from Snakemake."""
    # Input files
    fastq_file = snakemake.input.fastq
    breakpoints_file = snakemake.input.breakpoints
    ends_file = snakemake.input.ends

    # Output
    output_file = snakemake.output[0]

    # Parameters
    show_progress = snakemake.params.get('show_progress', True)
    progress_interval = snakemake.params.get('progress_interval', 1)

    # Logging
    log_file = snakemake.log[0] if snakemake.log else None
    logger = setup_logging(log_file, name='string_matcher')

    # Run matching
    run_matching(
        fastq_file=fastq_file,
        breakpoints_file=breakpoints_file,
        ends_file=ends_file,
        output_file=output_file,
        show_progress=show_progress,
        progress_interval=progress_interval,
        logger=logger
    )


def main_cli() -> None:
    """Entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description='Detect fusion breakpoints via string matching',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Input files
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input FASTQ file (gzipped or plain)'
    )
    parser.add_argument(
        '-b', '--breakpoints',
        required=True,
        help='Breakpoint sequences CSV file'
    )
    parser.add_argument(
        '-e', '--ends',
        required=True,
        help='Domain ends CSV file'
    )

    # Output
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output counts CSV file'
    )

    # Options
    parser.add_argument(
        '--progress',
        action='store_true',
        default=True,
        help='Show progress (default: on)'
    )
    parser.add_argument(
        '--no-progress',
        action='store_true',
        help='Disable progress display'
    )
    parser.add_argument(
        '--progress-interval',
        type=float,
        default=1.0,
        help='Progress update interval in %% (default: 1)'
    )

    # Logging
    parser.add_argument(
        '--log',
        help='Log file path'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Setup logging
    level = 'DEBUG' if args.verbose else 'INFO'
    logger = setup_logging(args.log, level=level, name='string_matcher')

    # Determine progress setting
    show_progress = args.progress and not args.no_progress

    # Run matching
    run_matching(
        fastq_file=args.input,
        breakpoints_file=args.breakpoints,
        ends_file=args.ends,
        output_file=args.output,
        show_progress=show_progress,
        progress_interval=args.progress_interval,
        logger=logger
    )


def run_matching(
    fastq_file: str,
    breakpoints_file: str,
    ends_file: str,
    output_file: str,
    show_progress: bool,
    progress_interval: float,
    logger=None
) -> None:
    """Core matching logic."""
    # Load reference data
    if logger:
        logger.info(f"Loading breakpoint sequences from: {breakpoints_file}")
    breakpoints = load_breakpoint_sequences(breakpoints_file)

    if logger:
        logger.info(f"Loading domain ends from: {ends_file}")
    domain_ends = load_domain_ends(ends_file)

    # Check for pyfastx
    if not HAS_PYFASTX:
        if logger:
            logger.warning(
                "pyfastx not available - using slower pure Python FASTQ parser. "
                "Install pyfastx for better performance: pip install pyfastx"
            )

    # Run matching
    counts = count_fusion_matches(
        fastq_file=fastq_file,
        breakpoints=breakpoints,
        domain_ends=domain_ends,
        show_progress=show_progress,
        progress_interval=progress_interval,
        logger=logger
    )

    # Write output
    if logger:
        logger.info(f"Writing counts to: {output_file}")
    write_counts_csv(counts, output_file)

    if logger:
        logger.info("Done!")


# Entry point
if __name__ == '__main__':
    # Check if running under Snakemake
    try:
        snakemake
        main_snakemake(snakemake)
    except NameError:
        main_cli()

