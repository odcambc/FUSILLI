#!/usr/bin/env python3
"""
Generate fusion breakpoint sequences for detection.

This script generates all possible breakpoint sequences for a fusion library,
producing the reference sequences needed for string-matching detection.

FUSION STRUCTURE:
    A typical kinase fusion has the structure:

    [Fusion Partner N-terminus]--[Linker]--[Kinase Domain]

    For example, TPR-MET fusion:
    [TPR (variable truncation)]--[GS linker]--[MET kinase domain (full)]

    The "breakpoint" is where the partner truncates. We generate all possible
    in-frame truncation points and create k-mer sequences spanning each breakpoint.

OUTPUT FILES:
    1. breakpoint_sequences.csv - All breakpoint k-mers for detection
    2. domain_ends.csv - K-mer ends of each domain for pre-filtering

USAGE:
    # Via Snakemake (automatic parameter passing)
    script: "scripts/fusion_sequences.py"

    # Standalone
    python fusion_sequences.py \\
        --sequences references/kinase_sequences.fasta \\
        --partners config/fusion_partners.csv \\
        --anchor Met_WT \\
        --anchor-position downstream \\
        --linker GGGAGC \\
        --window 12 \\
        --output-breakpoints results/breakpoint_sequences.csv \\
        --output-ends results/domain_ends.csv
"""

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

# Handle both Snakemake and standalone execution
try:
    from utils import (
        parse_fasta,
        parse_partners_csv,
        validate_nucleotide_sequence,
        validate_sequences_match_config,
        setup_logging,
        ProgressReporter
    )
except ImportError:
    from workflow.scripts.utils import (
        parse_fasta,
        parse_partners_csv,
        validate_nucleotide_sequence,
        validate_sequences_match_config,
        setup_logging,
        ProgressReporter
    )


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class BreakpointSequence:
    """Represents a single fusion breakpoint sequence."""
    fusion_id: str          # Unique identifier: {partner}_{bp_pos}_{anchor}
    partner_name: str       # Name of fusion partner
    anchor_name: str        # Name of anchor domain
    breakpoint_nt: int      # Breakpoint position in partner (nucleotides from start)
    breakpoint_aa: int      # Breakpoint position (amino acids/codons)
    sequence: str           # The actual breakpoint k-mer for detection
    full_fusion_length: int # Length of full fusion construct

    def to_dict(self) -> dict:
        """Convert to dictionary for CSV output."""
        return {
            'fusion_id': self.fusion_id,
            'partner_name': self.partner_name,
            'anchor_name': self.anchor_name,
            'breakpoint_nt': self.breakpoint_nt,
            'breakpoint_aa': self.breakpoint_aa,
            'breakpoint_sequence': self.sequence,
            'full_fusion_length': self.full_fusion_length
        }


@dataclass
class FusionLibraryConfig:
    """Configuration for a fusion library."""
    anchor_name: str
    anchor_position: str  # 'upstream' or 'downstream'
    linker_sequence: str
    breakpoint_window: int
    maintain_frame: bool
    kmer_size: int


# =============================================================================
# FUSION GENERATION
# =============================================================================

def generate_breakpoint_positions(
    partner_length: int,
    maintain_frame: bool = True
) -> Iterator[int]:
    """
    Generate valid breakpoint positions for a partner sequence.

    Args:
        partner_length: Length of partner sequence in nucleotides
        maintain_frame: Only generate in-frame (codon boundary) positions

    Yields:
        Valid breakpoint positions (nucleotide index)
    """
    step = 3 if maintain_frame else 1

    # Start from 0 (complete partner removal) to end (full partner)
    # Typically we exclude 0 as that would be no partner contribution
    for pos in range(step, partner_length + 1, step):
        yield pos


def generate_fusion_sequence(
    partner_seq: str,
    anchor_seq: str,
    linker_seq: str,
    breakpoint_pos: int,
    anchor_position: str
) -> tuple[str, int]:
    """
    Generate a complete fusion sequence at a specific breakpoint.

    Args:
        partner_seq: Full partner sequence
        anchor_seq: Full anchor sequence
        linker_seq: Linker sequence
        breakpoint_pos: Position to truncate partner
        anchor_position: Where anchor appears ('upstream' or 'downstream')

    Returns:
        Tuple of (fusion_sequence, breakpoint_absolute_position)
        where breakpoint_absolute_position is the index in the fusion
        where the junction occurs
    """
    partner_portion = partner_seq[:breakpoint_pos]

    if anchor_position == 'downstream':
        # Partner--Linker--Anchor  (most common for kinase fusions)
        fusion = partner_portion + linker_seq + anchor_seq
        # Breakpoint is right after partner portion
        breakpoint_absolute = len(partner_portion)
    else:
        # Anchor--Linker--Partner
        fusion = anchor_seq + linker_seq + partner_portion
        # Breakpoint is after anchor + linker
        breakpoint_absolute = len(anchor_seq) + len(linker_seq)

    return fusion, breakpoint_absolute


def extract_breakpoint_kmer(
    fusion_seq: str,
    breakpoint_pos: int,
    window: int
) -> str | None:
    """
    Extract a k-mer sequence centered on the breakpoint.

    Args:
        fusion_seq: Complete fusion sequence
        breakpoint_pos: Position of breakpoint in fusion
        window: Nucleotides on each side of breakpoint

    Returns:
        K-mer sequence, or None if window extends beyond sequence
    """
    start = breakpoint_pos - window
    end = breakpoint_pos + window

    # Check bounds
    if start < 0 or end > len(fusion_seq):
        return None

    return fusion_seq[start:end]


def generate_all_breakpoints(
    sequences: dict[str, str],
    partners: dict[str, dict],
    config: FusionLibraryConfig,
    logger=None
) -> list[BreakpointSequence]:
    """
    Generate all breakpoint sequences for a fusion library.

    Args:
        sequences: Dictionary of domain sequences
        partners: Partner configuration
        config: Fusion library configuration
        logger: Optional logger

    Returns:
        List of BreakpointSequence objects
    """
    anchor_seq = sequences[config.anchor_name]
    breakpoints = []

    # Count total breakpoints for progress reporting
    active_partners = [p for p, cfg in partners.items() if cfg['include']]

    if logger:
        logger.info(f"Generating breakpoints for {len(active_partners)} partners")
        logger.info(f"Anchor: {config.anchor_name} ({len(anchor_seq)} nt)")
        logger.info(f"Linker: {config.linker_sequence or '(none)'}")
        logger.info(f"Window: ±{config.breakpoint_window} nt")

    for partner_name in active_partners:
        if partner_name not in sequences:
            if logger:
                logger.warning(f"Skipping {partner_name}: not found in sequences")
            continue

        partner_seq = sequences[partner_name]

        for bp_pos in generate_breakpoint_positions(
            len(partner_seq),
            config.maintain_frame
        ):
            # Generate the fusion sequence
            fusion_seq, bp_absolute = generate_fusion_sequence(
                partner_seq=partner_seq,
                anchor_seq=anchor_seq,
                linker_seq=config.linker_sequence,
                breakpoint_pos=bp_pos,
                anchor_position=config.anchor_position
            )

            # Extract the breakpoint k-mer centered on the junction
            kmer = extract_breakpoint_kmer(
                fusion_seq,
                bp_absolute,  # Junction between partner and linker
                config.breakpoint_window
            )

            if kmer is None:
                continue

            # Create the breakpoint record
            fusion_id = f"{partner_name}_{bp_pos}_{config.anchor_name}"

            breakpoints.append(BreakpointSequence(
                fusion_id=fusion_id,
                partner_name=partner_name,
                anchor_name=config.anchor_name,
                breakpoint_nt=bp_pos,
                breakpoint_aa=bp_pos // 3,
                sequence=kmer,
                full_fusion_length=len(fusion_seq)
            ))

    if logger:
        logger.info(f"Generated {len(breakpoints)} breakpoint sequences")

    return breakpoints


def generate_domain_ends(
    sequences: dict[str, str],
    partners: dict[str, dict],
    anchor_name: str,
    kmer_size: int = 15
) -> dict[str, str]:
    """
    Generate k-mer ends for each domain (used for pre-filtering).

    The 3' end of each partner domain is used to quickly identify
    which fusion partner might be present in a read before doing
    full breakpoint matching.

    Args:
        sequences: Domain sequences
        partners: Partner configuration
        anchor_name: Name of anchor (excluded from ends)
        kmer_size: Size of k-mer to extract

    Returns:
        Dictionary mapping domain names to their 3' end k-mers
    """
    ends = {}

    for partner_name, config in partners.items():
        if not config['include']:
            continue
        if partner_name == anchor_name:
            continue
        if partner_name not in sequences:
            continue

        seq = sequences[partner_name]
        if len(seq) >= kmer_size:
            ends[partner_name] = seq[-kmer_size:]

    return ends


# =============================================================================
# OUTPUT
# =============================================================================

def write_breakpoints_csv(
    breakpoints: list[BreakpointSequence],
    filepath: Path
) -> None:
    """Write breakpoint sequences to CSV."""
    filepath.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        'fusion_id', 'partner_name', 'anchor_name',
        'breakpoint_nt', 'breakpoint_aa', 'breakpoint_sequence',
        'full_fusion_length'
    ]

    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for bp in breakpoints:
            writer.writerow(bp.to_dict())


def write_domain_ends_csv(
    ends: dict[str, str],
    filepath: Path
) -> None:
    """Write domain end k-mers to CSV."""
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['domain_name', 'end_kmer'])
        for name, kmer in sorted(ends.items()):
            writer.writerow([name, kmer])


# =============================================================================
# MAIN
# =============================================================================

def main_snakemake(snakemake) -> None:
    """Entry point when called from Snakemake."""
    # Extract parameters from snakemake object
    sequences_file = snakemake.input.sequences
    partners_file = snakemake.params.partners_file

    output_breakpoints = snakemake.output.breakpoints
    output_ends = snakemake.output.ends

    # Get config values
    anchor_name = snakemake.params.anchor_name
    anchor_position = snakemake.params.get('anchor_position', 'downstream')
    linker_sequence = snakemake.params.get('linker_sequence', '')
    breakpoint_window = snakemake.params.get('breakpoint_window', 12)
    maintain_frame = snakemake.params.get('maintain_frame', True)
    kmer_size = snakemake.params.get('kmer_size', 15)

    # Setup logging
    log_file = snakemake.log[0] if snakemake.log else None
    logger = setup_logging(log_file, name='fusion_sequences')

    # Run generation
    run_generation(
        sequences_file=sequences_file,
        partners_file=partners_file,
        output_breakpoints=output_breakpoints,
        output_ends=output_ends,
        anchor_name=anchor_name,
        anchor_position=anchor_position,
        linker_sequence=linker_sequence,
        breakpoint_window=breakpoint_window,
        maintain_frame=maintain_frame,
        kmer_size=kmer_size,
        logger=logger
    )


def main_cli() -> None:
    """Entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description='Generate fusion breakpoint sequences for detection',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Input files
    parser.add_argument(
        '-s', '--sequences',
        required=True,
        help='FASTA file containing domain sequences'
    )
    parser.add_argument(
        '-p', '--partners',
        required=True,
        help='CSV file defining fusion partners'
    )

    # Fusion library configuration
    parser.add_argument(
        '-a', '--anchor',
        required=True,
        help='Name of anchor domain (e.g., Met_WT)'
    )
    parser.add_argument(
        '--anchor-position',
        choices=['upstream', 'downstream'],
        default='downstream',
        help='Position of anchor in fusion (default: downstream/3\')'
    )
    parser.add_argument(
        '-l', '--linker',
        default='',
        help='Linker sequence (default: none)'
    )

    # Detection parameters
    parser.add_argument(
        '-w', '--window',
        type=int,
        default=12,
        help='Breakpoint window size in nt (default: 12)'
    )
    parser.add_argument(
        '--no-frame',
        action='store_true',
        help='Generate all breakpoints, not just in-frame'
    )
    parser.add_argument(
        '-k', '--kmer-size',
        type=int,
        default=15,
        help='K-mer size for domain ends (default: 15)'
    )

    # Output files
    parser.add_argument(
        '-o', '--output-breakpoints',
        required=True,
        help='Output CSV for breakpoint sequences'
    )
    parser.add_argument(
        '-e', '--output-ends',
        required=True,
        help='Output CSV for domain end k-mers'
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
    logger = setup_logging(args.log, level=level, name='fusion_sequences')

    # Run generation
    run_generation(
        sequences_file=args.sequences,
        partners_file=args.partners,
        output_breakpoints=args.output_breakpoints,
        output_ends=args.output_ends,
        anchor_name=args.anchor,
        anchor_position=args.anchor_position,
        linker_sequence=args.linker,
        breakpoint_window=args.window,
        maintain_frame=not args.no_frame,
        kmer_size=args.kmer_size,
        logger=logger
    )


def run_generation(
    sequences_file: str,
    partners_file: str,
    output_breakpoints: str,
    output_ends: str,
    anchor_name: str,
    anchor_position: str,
    linker_sequence: str,
    breakpoint_window: int,
    maintain_frame: bool,
    kmer_size: int,
    logger=None
) -> None:
    """Core generation logic."""
    # Parse inputs
    if logger:
        logger.info(f"Loading sequences from: {sequences_file}")
    sequences = parse_fasta(sequences_file)

    if logger:
        logger.info(f"Loading partners from: {partners_file}")
    partners = parse_partners_csv(partners_file)

    # Validate
    messages = validate_sequences_match_config(sequences, partners, anchor_name)
    for msg in messages:
        if logger:
            if msg.startswith('ERROR'):
                logger.error(msg)
            else:
                logger.warning(msg)
        else:
            print(msg, file=sys.stderr)

    if any(m.startswith('ERROR') for m in messages):
        sys.exit(1)

    # Create config
    config = FusionLibraryConfig(
        anchor_name=anchor_name,
        anchor_position=anchor_position,
        linker_sequence=linker_sequence,
        breakpoint_window=breakpoint_window,
        maintain_frame=maintain_frame,
        kmer_size=kmer_size
    )

    # Generate breakpoints
    breakpoints = generate_all_breakpoints(sequences, partners, config, logger)

    # Generate domain ends
    domain_ends = generate_domain_ends(
        sequences, partners, anchor_name, kmer_size
    )

    # Write outputs
    if logger:
        logger.info(f"Writing breakpoints to: {output_breakpoints}")
    write_breakpoints_csv(breakpoints, Path(output_breakpoints))

    if logger:
        logger.info(f"Writing domain ends to: {output_ends}")
    write_domain_ends_csv(domain_ends, Path(output_ends))

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

