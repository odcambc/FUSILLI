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
        parse_unfused_sequences_csv,
        parse_exon_partners_csv,
        validate_nucleotide_sequence,
        validate_sequences_match_config,
        setup_logging,
        ProgressReporter
    )
except ImportError:
    from workflow.scripts.utils import (
        parse_fasta,
        parse_partners_csv,
        parse_unfused_sequences_csv,
        parse_exon_partners_csv,
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
    """
    Represents a single fusion breakpoint sequence.

    The breakpoint position (breakpoint_nt) is relative to whichever component
    is truncated, as specified by the truncated_component configuration:
    - When truncated_component == 'partner': position is from partner start
    - When truncated_component == 'anchor': position is from anchor start
    """
    fusion_id: str          # Unique identifier: {partner}_{bp_pos}_{anchor}
    partner_name: str       # Name of fusion partner
    anchor_name: str        # Name of anchor domain
    breakpoint_nt: int      # Breakpoint position (nucleotides from start of truncated component)
    breakpoint_aa: int      # Breakpoint position (amino acids/codons)
    sequence: str           # The actual breakpoint k-mer for detection
    full_fusion_length: int  # Length of full fusion construct
    variant_anchor_name: str | None = None  # Name of variant anchor (None for regular fusions)

    def to_dict(self) -> dict:
        """Convert to dictionary for CSV output."""
        return {
            'fusion_id': self.fusion_id,
            'partner_name': self.partner_name,
            'anchor_name': self.anchor_name,
            'breakpoint_nt': self.breakpoint_nt,
            'breakpoint_aa': self.breakpoint_aa,
            'breakpoint_sequence': self.sequence,
            'full_fusion_length': self.full_fusion_length,
            'variant_anchor_name': self.variant_anchor_name or ''
        }


@dataclass
class UnfusedSequence:
    """Represents a k-mer from an unfused sequence."""
    sequence_name: str
    kmer: str
    position_in_sequence: int
    sequence_length: int

    def to_dict(self) -> dict:
        """Convert to dictionary for CSV output."""
        return {
            'sequence_name': self.sequence_name,
            'kmer': self.kmer,
            'position_in_sequence': self.position_in_sequence,
            'sequence_length': self.sequence_length
        }


@dataclass
class FusionLibraryConfig:
    """Configuration for a fusion library."""
    anchor_name: str
    anchor_position: str  # 'upstream' or 'downstream'
    truncated_component: str = 'partner'  # 'partner' (default) or 'anchor'
    linker_sequence: str = ''
    breakpoint_window: int = 12
    maintain_frame: bool = True
    kmer_size: int = 15


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
    anchor_position: str,
    truncated_component: str = 'partner'
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
        where the junction between upstream+linker and downstream occurs
    """
    # Determine which component is truncated
    if truncated_component == 'anchor':
        partner_portion = partner_seq  # full-length upstream component
        # Keep the 3' portion of the anchor (truncate N-terminal region).
        anchor_portion = anchor_seq[breakpoint_pos:]
    else:
        partner_portion = partner_seq[:breakpoint_pos]
        anchor_portion = anchor_seq

    if anchor_position == 'downstream':
        # Upstream--Linker--Downstream
        fusion = partner_portion + linker_seq + anchor_portion
        breakpoint_absolute = len(partner_portion) + len(linker_seq)
    else:
        # Downstream--Linker--Upstream
        fusion = anchor_portion + linker_seq + partner_portion
        # Breakpoint is after downstream portion + linker
        breakpoint_absolute = len(anchor_portion) + len(linker_seq)

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

        # Decide which component is truncated
        truncated_seq = (
            sequences[config.anchor_name]
            if config.truncated_component == 'anchor'
            else partner_seq
        )

        for bp_pos in generate_breakpoint_positions(
            len(truncated_seq),
            config.maintain_frame
        ):
            # Generate the fusion sequence
            fusion_seq, bp_absolute = generate_fusion_sequence(
                partner_seq=partner_seq,
                anchor_seq=anchor_seq,
                linker_seq=config.linker_sequence,
                breakpoint_pos=bp_pos,
                anchor_position=config.anchor_position,
                truncated_component=config.truncated_component
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
    full breakpoint matching. This remains true even when the anchor
    is truncated, since short reads covering the breakpoint are most
    likely to include the partner 3' end rather than the anchor 3' end.

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
# UNFUSED SEQUENCE GENERATION
# =============================================================================

def generate_unfused_kmers(
    sequences: dict[str, str],
    unfused_config: dict[str, dict],
    kmer_size: int = 15,
    spacing: int = 50,
    exclude_kmers: set[str] | None = None,
    exclude_names: set[str] | None = None,
    logger=None
) -> list[UnfusedSequence]:
    """
    Generate multiple k-mers distributed across unfused sequences.

    Args:
        sequences: Dictionary of domain sequences
        unfused_config: Unfused sequence configuration
        kmer_size: Size of k-mers to generate
        spacing: Spacing between k-mers in nucleotides (default: 50)
        exclude_kmers: Optional set of k-mers to exclude (e.g., overlap with fusion breakpoints)
        exclude_names: Optional set of sequence names to apply exclusion to
        logger: Optional logger

    Returns:
        List of UnfusedSequence objects
    """
    unfused_kmers = []
    active_sequences = [name for name, cfg in unfused_config.items() if cfg['include']]

    if logger:
        logger.info(f"Generating k-mers for {len(active_sequences)} unfused sequences")
        logger.info(f"K-mer size: {kmer_size} nt, spacing: {spacing} nt")

    for seq_name in active_sequences:
        if seq_name not in sequences:
            if logger:
                logger.warning(f"Skipping {seq_name}: not found in sequences")
            continue

        seq = sequences[seq_name]
        seq_length = len(seq)

        # Generate k-mers distributed across the sequence
        # Start from kmer_size/2 to ensure we can extract a full k-mer
        start_pos = kmer_size // 2
        # End before the last kmer_size/2 to avoid going past the end
        end_pos = seq_length - (kmer_size // 2)

        positions = []
        if end_pos > start_pos:
            # Generate positions evenly spaced
            for pos in range(start_pos, end_pos, spacing):
                positions.append(pos)
            # Always include the last position if we haven't already
            if positions and positions[-1] < end_pos - spacing:
                positions.append(end_pos - (kmer_size // 2))
        else:
            # Sequence is too short, just use the middle
            positions = [seq_length // 2] if seq_length >= kmer_size else []

        for pos in positions:
            # Extract k-mer centered at position
            kmer_start = max(0, pos - kmer_size // 2)
            kmer_end = min(seq_length, pos + (kmer_size + 1) // 2)
            kmer = seq[kmer_start:kmer_end]

            # Only add if we got a full k-mer
            if len(kmer) >= kmer_size:
                kmer = kmer[:kmer_size]
                if exclude_kmers and exclude_names and seq_name in exclude_names and kmer in exclude_kmers:
                    continue
                unfused_kmers.append(UnfusedSequence(
                    sequence_name=seq_name,
                    kmer=kmer,
                    position_in_sequence=pos,
                    sequence_length=seq_length
                ))

    if logger:
        logger.info(f"Generated {len(unfused_kmers)} unfused k-mers")

    return unfused_kmers


# =============================================================================
# VARIANT ANCHOR BREAKPOINT GENERATION
# =============================================================================

def generate_variant_breakpoints(
    sequences: dict[str, str],
    partners: dict[str, dict],
    variant_anchor_config: dict,
    config: FusionLibraryConfig,
    logger=None
) -> list[BreakpointSequence]:
    """
    Generate breakpoint sequences for variant anchor fusions.

    Only generates fusions for partners specified in variant_anchor_config['partners'].

    Args:
        sequences: Dictionary of domain sequences
        partners: Partner configuration
        variant_anchor_config: Variant anchor config with 'name' and 'partners' list
        config: Fusion library configuration
        logger: Optional logger

    Returns:
        List of BreakpointSequence objects
    """
    variant_anchor_name = variant_anchor_config['name']
    variant_partners = variant_anchor_config.get('partners', [])
    all_partners = variant_anchor_config.get('all_partners', False)

    if isinstance(variant_partners, str):
        variant_partners = [variant_partners]

    if all_partners or not variant_partners or any(p in ("*", "all") for p in variant_partners):
        variant_partners = [p for p, cfg in partners.items() if cfg.get('include')]

    if variant_anchor_name not in sequences:
        if logger:
            logger.warning(f"Variant anchor '{variant_anchor_name}' not found in sequences")
        return []

    variant_anchor_seq = sequences[variant_anchor_name]
    breakpoints = []

    if logger:
        logger.info(f"Generating variant breakpoints for anchor: {variant_anchor_name}")
        logger.info(f"Partners: {', '.join(variant_partners)}")

    for partner_name in variant_partners:
        if partner_name not in partners:
            if logger:
                logger.warning(f"Skipping variant partner '{partner_name}': not in partners config")
            continue

        if not partners[partner_name]['include']:
            if logger:
                logger.debug(f"Skipping variant partner '{partner_name}': not included")
            continue

        if partner_name not in sequences:
            if logger:
                logger.warning(f"Skipping variant partner '{partner_name}': not found in sequences")
            continue

        partner_seq = sequences[partner_name]

        # Decide which component is truncated (use same logic as regular fusions)
        truncated_seq = (
            variant_anchor_seq
            if config.truncated_component == 'anchor'
            else partner_seq
        )

        for bp_pos in generate_breakpoint_positions(
            len(truncated_seq),
            config.maintain_frame
        ):
            # Generate the fusion sequence using variant anchor
            fusion_seq, bp_absolute = generate_fusion_sequence(
                partner_seq=partner_seq,
                anchor_seq=variant_anchor_seq,
                linker_seq=config.linker_sequence,
                breakpoint_pos=bp_pos,
                anchor_position=config.anchor_position,
                truncated_component=config.truncated_component
            )

            # Extract the breakpoint k-mer centered on the junction
            kmer = extract_breakpoint_kmer(
                fusion_seq,
                bp_absolute,
                config.breakpoint_window
            )

            if kmer is None:
                continue

            # Create the breakpoint record with variant anchor name
            fusion_id = f"{partner_name}_{bp_pos}_{variant_anchor_name}"

            breakpoints.append(BreakpointSequence(
                fusion_id=fusion_id,
                partner_name=partner_name,
                anchor_name=variant_anchor_name,
                breakpoint_nt=bp_pos,
                breakpoint_aa=bp_pos // 3,
                sequence=kmer,
                full_fusion_length=len(fusion_seq),
                variant_anchor_name=variant_anchor_name
            ))

    if logger:
        logger.info(f"Generated {len(breakpoints)} variant breakpoint sequences")

    return breakpoints


def _reverse_complement(seq: str) -> str:
    """Return reverse complement of a nucleotide sequence."""
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


def build_exclusion_kmers(
    breakpoints: list[BreakpointSequence],
    kmer_size: int,
    include_reverse_complement: bool = True
) -> set[str]:
    """Build a set of k-mers to exclude based on breakpoint sequences."""
    exclude = set()
    for bp in breakpoints:
        seq = bp.sequence
        if len(seq) < kmer_size:
            continue
        for i in range(0, len(seq) - kmer_size + 1):
            exclude.add(seq[i:i + kmer_size])
    if include_reverse_complement:
        exclude |= {_reverse_complement(k) for k in exclude}
    return exclude


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
        'full_fusion_length', 'variant_anchor_name'
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


def write_unfused_sequences_csv(
    unfused_kmers: list[UnfusedSequence],
    filepath: Path
) -> None:
    """Write unfused sequence k-mers to CSV."""
    filepath.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = ['sequence_name', 'kmer', 'position_in_sequence', 'sequence_length']

    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for unfused in unfused_kmers:
            writer.writerow(unfused.to_dict())


def build_variant_rows(
    breakpoints: list[BreakpointSequence],
    partners: dict[str, dict],
    unfused_config: dict[str, dict]
) -> list[dict]:
    """Build variant metadata rows for dataframe-friendly CSV output."""
    rows: list[dict] = []

    for bp in breakpoints:
        partner_cfg = partners.get(bp.partner_name, {})
        rows.append({
            'fusion_id': bp.fusion_id,
            'type': 'fusion',
            'partner_name': bp.partner_name,
            'anchor_name': bp.anchor_name,
            'breakpoint_nt': bp.breakpoint_nt,
            'breakpoint_aa': bp.breakpoint_aa,
            'full_fusion_length': bp.full_fusion_length,
            'variant_anchor_name': bp.variant_anchor_name or '',
            'sequence_length': partner_cfg.get('sequence_length', ''),
            'description': partner_cfg.get('description', '')
        })

    for name, cfg in sorted(unfused_config.items()):
        if not cfg.get('include', False):
            continue
        rows.append({
            'fusion_id': name,
            'type': 'unfused',
            'partner_name': '',
            'anchor_name': '',
            'breakpoint_nt': '',
            'breakpoint_aa': '',
            'full_fusion_length': '',
            'variant_anchor_name': '',
            'sequence_length': cfg.get('sequence_length', ''),
            'description': cfg.get('description', '')
        })

    return rows


def write_variants_csv(rows: list[dict], filepath: Path) -> None:
    """Write variant metadata to a CSV file."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        'fusion_id',
        'type',
        'partner_name',
        'anchor_name',
        'breakpoint_nt',
        'breakpoint_aa',
        'full_fusion_length',
        'variant_anchor_name',
        'sequence_length',
        'description'
    ]

    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_expected_counts_csv(rows: list[dict], filepath: Path) -> None:
    """Write a zero-count template covering all expected variants."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ['fusion_id', 'type', 'count']

    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                'fusion_id': row.get('fusion_id', ''),
                'type': row.get('type', ''),
                'count': 0
            })


# =============================================================================
# MAIN
# =============================================================================

def main_snakemake(snakemake) -> None:
    """Entry point when called from Snakemake."""
    # Extract parameters from snakemake object
    sequences_file = snakemake.input.sequences
    partners_file = snakemake.params.partners_file
    exon_partners_file = snakemake.params.get('exon_partners_file', None)

    output_breakpoints = snakemake.output.breakpoints
    output_ends = snakemake.output.ends
    output_unfused = snakemake.output.get('unfused', None)
    output_variants = snakemake.output.get('variants', None)
    output_expected_counts = snakemake.output.get('expected_counts', None)

    # Get config values
    anchor_name = snakemake.params.anchor_name
    anchor_position = snakemake.params.get('anchor_position', 'downstream')
    truncated_component = snakemake.params.get('truncated_component', 'partner')
    linker_sequence = snakemake.params.get('linker_sequence', '')
    breakpoint_window = snakemake.params.get('breakpoint_window', 12)
    maintain_frame = snakemake.params.get('maintain_frame', True)
    kmer_size = snakemake.params.get('kmer_size', 15)

    # Optional: unfused sequences and variant anchors
    unfused_sequences_file = snakemake.params.get('unfused_sequences_file', None)
    variant_anchors = snakemake.params.get('variant_anchors', [])

    # Setup logging
    log_file = snakemake.log[0] if snakemake.log else None
    logger = setup_logging(log_file, name='fusion_sequences')

    # Run generation
    run_generation(
        sequences_file=sequences_file,
        partners_file=partners_file,
        output_breakpoints=output_breakpoints,
        output_ends=output_ends,
        output_unfused=output_unfused,
        output_variants=output_variants,
        output_expected_counts=output_expected_counts,
        anchor_name=anchor_name,
        anchor_position=anchor_position,
        truncated_component=truncated_component,
        linker_sequence=linker_sequence,
        breakpoint_window=breakpoint_window,
        maintain_frame=maintain_frame,
        kmer_size=kmer_size,
        unfused_sequences_file=unfused_sequences_file,
        exon_partners_file=exon_partners_file,
        variant_anchors=variant_anchors,
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
    parser.add_argument(
        '--exon-partners',
        help='Optional CSV file defining exon-based partners'
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
        '--truncated-component',
        choices=['partner', 'anchor'],
        default='partner',
        help='Which component is truncated: partner (default) or anchor'
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
    parser.add_argument(
        '--output-variants',
        help='Output CSV for variant metadata'
    )
    parser.add_argument(
        '--output-expected-counts',
        help='Output CSV for zero-count template'
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
        truncated_component=args.truncated_component,
        maintain_frame=not args.no_frame,
        kmer_size=args.kmer_size,
        exon_partners_file=args.exon_partners,
        output_variants=args.output_variants,
        output_expected_counts=args.output_expected_counts,
        logger=logger
    )


def run_generation(
    sequences_file: str,
    partners_file: str,
    output_breakpoints: str,
    output_ends: str,
    anchor_name: str,
    anchor_position: str,
    truncated_component: str,
    linker_sequence: str,
    breakpoint_window: int,
    maintain_frame: bool,
    kmer_size: int,
    unfused_sequences_file: str | None = None,
    exon_partners_file: str | None = None,
    variant_anchors: list[dict] | None = None,
    output_unfused: str | None = None,
    output_variants: str | None = None,
    output_expected_counts: str | None = None,
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

    if exon_partners_file:
        if logger:
            logger.info(f"Loading exon partners from: {exon_partners_file}")
        try:
            exon_partners = parse_exon_partners_csv(exon_partners_file, sequences)
            for name, cfg in exon_partners.items():
                partners.setdefault(name, cfg)
        except FileNotFoundError:
            if logger:
                logger.warning(f"Exon partners file not found: {exon_partners_file}")

    # Parse unfused sequences if provided
    unfused_config = {}
    if unfused_sequences_file:
        if logger:
            logger.info(f"Loading unfused sequences from: {unfused_sequences_file}")
        try:
            unfused_config = parse_unfused_sequences_csv(unfused_sequences_file)
        except FileNotFoundError:
            if logger:
                logger.warning(f"Unfused sequences file not found: {unfused_sequences_file}")
            unfused_config = {}

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
        truncated_component=truncated_component,
        linker_sequence=linker_sequence,
        breakpoint_window=breakpoint_window,
        maintain_frame=maintain_frame,
        kmer_size=kmer_size
    )

    # Generate regular breakpoints
    breakpoints = generate_all_breakpoints(sequences, partners, config, logger)

    # Generate variant anchor breakpoints if provided
    if variant_anchors:
        for variant_config in variant_anchors:
            variant_bps = generate_variant_breakpoints(
                sequences, partners, variant_config, config, logger
            )
            breakpoints.extend(variant_bps)

    # Generate unfused sequence k-mers if provided
    unfused_kmers = []
    if unfused_config:
        exclude_names = {
            name for name, cfg in unfused_config.items()
            if cfg.get('exclude_overlap')
        }
        exclude_kmers = (
            build_exclusion_kmers(breakpoints, kmer_size=kmer_size)
            if exclude_names else None
        )
        unfused_kmers = generate_unfused_kmers(
            sequences,
            unfused_config,
            kmer_size=kmer_size,
            exclude_kmers=exclude_kmers,
            exclude_names=exclude_names,
            logger=logger
        )

    # Generate domain ends
    domain_ends = generate_domain_ends(
        sequences,
        partners,
        anchor_name,
        kmer_size=kmer_size
    )

    # Write outputs
    if logger:
        logger.info(f"Writing breakpoints to: {output_breakpoints}")
    write_breakpoints_csv(breakpoints, Path(output_breakpoints))

    if logger:
        logger.info(f"Writing domain ends to: {output_ends}")
    write_domain_ends_csv(domain_ends, Path(output_ends))

    if output_variants:
        if logger:
            logger.info(f"Writing variant catalog to: {output_variants}")
        variant_rows = build_variant_rows(breakpoints, partners, unfused_config)
        write_variants_csv(variant_rows, Path(output_variants))
    else:
        variant_rows = build_variant_rows(breakpoints, partners, unfused_config)

    if output_expected_counts:
        if logger:
            logger.info(f"Writing expected counts template to: {output_expected_counts}")
        write_expected_counts_csv(variant_rows, Path(output_expected_counts))

    # Always write unfused sequences file (even if empty) to satisfy Snakemake
    if output_unfused:
        if unfused_kmers:
            if logger:
                logger.info(f"Writing unfused sequences to: {output_unfused}")
            write_unfused_sequences_csv(unfused_kmers, Path(output_unfused))
        else:
            # Write empty file with just header
            if logger:
                logger.info(f"Writing empty unfused sequences file: {output_unfused}")
            Path(output_unfused).parent.mkdir(parents=True, exist_ok=True)
            with open(output_unfused, 'w', newline='') as f:
                fieldnames = [
                    'sequence_name', 'kmer',
                    'position_in_sequence', 'sequence_length'
                ]
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()

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
