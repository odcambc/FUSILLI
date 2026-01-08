"""
Shared utility functions for FUSILLI pipeline.

This module provides common functionality used across the pipeline:
- Sequence I/O (FASTA, CSV)
- Validation functions
- Logging setup
- Progress reporting
"""

import csv
import logging
import sys
import time
from pathlib import Path
from typing import Iterator, TextIO


# =============================================================================
# CONSTANTS
# =============================================================================

VALID_NUCLEOTIDES = set("ACGT")
VALID_NUCLEOTIDES_AMBIGUOUS = set("ACGTNRYSWKMBDHV")


# =============================================================================
# LOGGING
# =============================================================================

def setup_logging(
    log_file: str | Path | None = None,
    level: str = "INFO",
    name: str = "fusilli"
) -> logging.Logger:
    """
    Configure logging for FUSILLI scripts.

    Args:
        log_file: Path to log file. If None, logs to stderr only.
        level: Logging level ('DEBUG', 'INFO', 'WARNING', 'ERROR')
        name: Logger name

    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper()))

    # Clear any existing handlers
    logger.handlers.clear()

    # Console handler (stderr)
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(logging.INFO)
    console_format = logging.Formatter(
        '%(asctime)s | %(levelname)-8s | %(message)s',
        datefmt='%H:%M:%S'
    )
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)

    # File handler (if specified)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_format = logging.Formatter(
            '%(asctime)s | %(levelname)-8s | %(name)s | %(message)s'
        )
        file_handler.setFormatter(file_format)
        logger.addHandler(file_handler)

    return logger


# =============================================================================
# SEQUENCE I/O
# =============================================================================

def parse_fasta(filepath: str | Path) -> dict[str, str]:
    """
    Parse a FASTA file into a dictionary.

    Args:
        filepath: Path to FASTA file

    Returns:
        Dictionary mapping sequence names to sequences (uppercase)

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file is malformed
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"FASTA file not found: {filepath}")

    sequences = {}
    current_name = None
    current_seq = []

    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()

            if not line:
                continue

            if line.startswith('>'):
                # Save previous sequence
                if current_name is not None:
                    sequences[current_name] = ''.join(current_seq).upper()

                # Start new sequence
                current_name = line[1:].split()[0]  # Take first word after >
                if not current_name:
                    raise ValueError(f"Empty sequence name at line {line_num}")
                if current_name in sequences:
                    raise ValueError(f"Duplicate sequence name: {current_name}")
                current_seq = []
            else:
                if current_name is None:
                    raise ValueError(f"Sequence data before header at line {line_num}")
                current_seq.append(line)

    # Don't forget the last sequence
    if current_name is not None:
        sequences[current_name] = ''.join(current_seq).upper()

    return sequences


def parse_partners_csv(filepath: str | Path) -> dict[str, dict]:
    """
    Parse the fusion partners CSV file.

    Args:
        filepath: Path to partners CSV file

    Returns:
        Dictionary mapping partner names to their properties:
        {
            'partner_name': {
                'sequence_length': int,
                'include': bool,
                'description': str
            }
        }
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Partners file not found: {filepath}")

    partners = {}

    with open(filepath, 'r') as f:
        # Skip comment lines
        lines = [l for l in f if not l.strip().startswith('#')]

    reader = csv.DictReader(lines)

    required_cols = {'partner_name', 'sequence_length', 'include'}
    if not required_cols.issubset(set(reader.fieldnames or [])):
        missing = required_cols - set(reader.fieldnames or [])
        raise ValueError(f"Missing required columns in partners file: {missing}")

    for row in reader:
        name = row['partner_name'].strip()
        if not name:
            continue

        partners[name] = {
            'sequence_length': int(row['sequence_length']),
            'include': row['include'].lower() in ('true', 'yes', '1'),
            'description': row.get('description', '').strip()
        }

    return partners


def parse_exon_partners_csv(
    filepath: str | Path,
    sequences: dict[str, str]
) -> dict[str, dict]:
    """
    Parse exon-based partner CSV file and normalize sequence names.

    Args:
        filepath: Path to exon partners CSV file
        sequences: Parsed FASTA sequences used to resolve lengths

    Returns:
        Dictionary mapping exon partner names to their properties.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Exon partners file not found: {filepath}")

    partners: dict[str, dict] = {}

    with open(filepath, 'r') as f:
        lines = [l for l in f if not l.strip().startswith('#')]

    reader = csv.DictReader(lines)
    required_cols = {'domain_exon_bp', 'domain', 'exon'}
    if not required_cols.issubset(set(reader.fieldnames or [])):
        missing = required_cols - set(reader.fieldnames or [])
        raise ValueError(f"Missing required columns in exon partners file: {missing}")

    for row in reader:
        raw_name = row['domain_exon_bp'].strip()
        if not raw_name:
            continue

        # Normalize exon partner names to match FASTA headers.
        normalized = raw_name.replace("_Rev_", "_")
        partner_name = normalized if normalized in sequences else raw_name

        seq_length = len(sequences.get(partner_name, "")) if partner_name in sequences else 0
        description = f"{row.get('domain', '').strip()} exon {row.get('exon', '').strip()}".strip()

        partners[partner_name] = {
            'sequence_length': seq_length,
            'include': True,
            'description': description
        }

    return partners


def parse_unfused_sequences_csv(filepath: str | Path) -> dict[str, dict]:
    """
    Parse the unfused sequences CSV file.

    Args:
        filepath: Path to unfused sequences CSV file

    Returns:
        Dictionary mapping sequence names to their properties:
        {
            'sequence_name': {
                'sequence_length': int,
                'include': bool,
                'description': str
            }
        }
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Unfused sequences file not found: {filepath}")

    sequences = {}

    with open(filepath, 'r') as f:
        # Skip comment lines
        lines = [l for l in f if not l.strip().startswith('#')]

    reader = csv.DictReader(lines)

    required_cols = {'sequence_name', 'sequence_length', 'include'}
    if not required_cols.issubset(set(reader.fieldnames or [])):
        missing = required_cols - set(reader.fieldnames or [])
        raise ValueError(f"Missing required columns in unfused sequences file: {missing}")

    for row in reader:
        name = row['sequence_name'].strip()
        if not name:
            continue

        sequences[name] = {
            'sequence_length': int(row['sequence_length']),
            'include': row['include'].lower() in ('true', 'yes', '1'),
            'exclude_overlap': row.get('exclude_overlap', '').lower() in ('true', 'yes', '1'),
            'description': row.get('description', '').strip()
        }

    return sequences


def parse_samples_csv(filepath: str | Path) -> dict[str, dict]:
    """
    Parse the samples CSV file.

    Args:
        filepath: Path to samples CSV file

    Returns:
        Dictionary mapping sample names to their properties
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Samples file not found: {filepath}")

    samples = {}

    with open(filepath, 'r') as f:
        # Skip comment lines
        lines = [l for l in f if not l.strip().startswith('#')]

    reader = csv.DictReader(lines)

    required_cols = {'sample', 'condition', 'file'}
    if not required_cols.issubset(set(reader.fieldnames or [])):
        missing = required_cols - set(reader.fieldnames or [])
        raise ValueError(f"Missing required columns in samples file: {missing}")

    for row in reader:
        name = row['sample'].strip()
        if not name:
            continue

        samples[name] = {
            'condition': row['condition'].strip(),
            'file': row['file'].strip(),
            'replicate': int(row.get('replicate', 1) or 1),
            'time': row.get('time', '0').strip(),
            'tile': row.get('tile', '1').strip()
        }

    return samples


# =============================================================================
# OUTPUT WRITERS
# =============================================================================

def write_csv(
    filepath: str | Path,
    data: list[dict],
    fieldnames: list[str],
    header_comment: str | None = None
) -> None:
    """
    Write data to a CSV file.

    Args:
        filepath: Output file path
        data: List of dictionaries to write
        fieldnames: Column names (determines order)
        header_comment: Optional comment to add at top of file
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with open(filepath, 'w', newline='') as f:
        if header_comment:
            for line in header_comment.strip().split('\n'):
                f.write(f"# {line}\n")

        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)


def write_fasta(
    filepath: str | Path,
    sequences: dict[str, str],
    line_width: int = 80
) -> None:
    """
    Write sequences to a FASTA file.

    Args:
        filepath: Output file path
        sequences: Dictionary mapping names to sequences
        line_width: Characters per line for sequence wrapping
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with open(filepath, 'w') as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n")
            # Wrap sequence
            for i in range(0, len(seq), line_width):
                f.write(f"{seq[i:i+line_width]}\n")


# =============================================================================
# VALIDATION
# =============================================================================

def validate_nucleotide_sequence(
    sequence: str,
    allow_ambiguous: bool = False
) -> tuple[bool, str | None]:
    """
    Validate that a sequence contains only valid nucleotides.

    Args:
        sequence: Nucleotide sequence to validate
        allow_ambiguous: Allow IUPAC ambiguity codes

    Returns:
        Tuple of (is_valid, error_message)
    """
    valid_chars = VALID_NUCLEOTIDES_AMBIGUOUS if allow_ambiguous else VALID_NUCLEOTIDES

    sequence = sequence.upper()
    invalid_chars = set(sequence) - valid_chars

    if invalid_chars:
        return False, f"Invalid characters: {', '.join(sorted(invalid_chars))}"

    return True, None


def validate_sequences_match_config(
    sequences: dict[str, str],
    partners: dict[str, dict],
    anchor_name: str
) -> list[str]:
    """
    Validate that sequence file matches configuration.

    Args:
        sequences: Parsed FASTA sequences
        partners: Partner configuration
        anchor_name: Name of anchor domain

    Returns:
        List of warning/error messages (empty if all valid)
    """
    messages = []

    # Check anchor exists
    if anchor_name not in sequences:
        messages.append(f"ERROR: Anchor domain '{anchor_name}' not found in sequences")

    # Check each included partner
    for partner_name, config in partners.items():
        if not config['include']:
            continue

        if partner_name not in sequences:
            messages.append(f"ERROR: Partner '{partner_name}' not found in sequences")
            continue

        actual_length = len(sequences[partner_name])
        expected_length = config['sequence_length']

        if actual_length != expected_length:
            messages.append(
                f"WARNING: Partner '{partner_name}' length mismatch: "
                f"config={expected_length}, actual={actual_length}"
            )

    return messages


# =============================================================================
# PROGRESS REPORTING
# =============================================================================

class ProgressReporter:
    """
    Report progress for long-running operations.

    Example:
        progress = ProgressReporter(total=1000000, desc="Processing reads")
        for item in items:
            process(item)
            progress.update()
        progress.finish()
    """

    def __init__(
        self,
        total: int,
        desc: str = "Processing",
        interval_pct: float = 1.0,
        enabled: bool = True,
        logger: logging.Logger | None = None,
        report_every_seconds: int = 60
    ):
        """
        Initialize progress reporter.

        Args:
            total: Total number of items
            desc: Description for progress messages
            interval_pct: Update interval as percentage
            enabled: Whether to show progress
            logger: Logger to use (defaults to print)
        """
        self.total = total
        self.desc = desc
        self.interval = max(1, int(total * interval_pct / 100))
        self.enabled = enabled
        self.logger = logger
        self.report_every_seconds = report_every_seconds

        self.count = 0
        self.start_time = time.time()
        self.last_report_pct = -1
        self.last_report_time = self.start_time

    def update(self, n: int = 1) -> None:
        """Update progress by n items."""
        self.count += n

        if not self.enabled:
            return

        now = time.time()
        if self.count % self.interval == 0 or self.count == self.total:
            self._report()
        elif now - self.last_report_time >= self.report_every_seconds:
            self._report()

    def _report(self) -> None:
        """Print progress report."""
        self.last_report_time = time.time()
        # Cap percentage at 100% if count exceeds total (estimate was too low)
        pct = min(100, int(100 * self.count / self.total)) if self.total > 0 else 0

        if pct == self.last_report_pct:
            return
        self.last_report_pct = pct

        elapsed = time.time() - self.start_time
        rate = self.count / elapsed if elapsed > 0 else 0

        # Calculate remaining time, handling cases where count exceeds total
        if rate > 0 and self.count < self.total:
            remaining = (self.total - self.count) / rate
            time_str = self._format_time(remaining)
            msg = f"{self.desc}: {pct}% ({self.count:,}/{self.total:,}) - {time_str} remaining"
        elif self.count >= self.total:
            # Estimate was too low - show actual count vs estimate
            msg = f"{self.desc}: {pct}% ({self.count:,}/{self.total:,} estimated) - processing..."
        else:
            msg = f"{self.desc}: {pct}% ({self.count:,}/{self.total:,})"

        if self.logger:
            self.logger.info(msg)
        else:
            print(msg, flush=True)

    def finish(self) -> None:
        """Report completion."""
        if not self.enabled:
            return

        elapsed = time.time() - self.start_time
        rate = self.count / elapsed if elapsed > 0 else 0

        msg = (
            f"{self.desc}: Complete! "
            f"Processed {self.count:,} items in {self._format_time(elapsed)} "
            f"({rate:,.0f}/sec)"
        )

        if self.logger:
            self.logger.info(msg)
        else:
            print(msg, flush=True)

    @staticmethod
    def _format_time(seconds: float) -> str:
        """Format seconds as human-readable string."""
        if seconds < 60:
            return f"{int(seconds)}s"
        elif seconds < 3600:
            mins = int(seconds / 60)
            secs = int(seconds % 60)
            return f"{mins}m {secs}s"
        else:
            hours = int(seconds / 3600)
            mins = int((seconds % 3600) / 60)
            return f"{hours}h {mins}m"


def progress_iterator(
    iterable,
    total: int | None = None,
    desc: str = "Processing",
    interval_pct: float = 1.0,
    enabled: bool = True
) -> Iterator:
    """
    Wrap an iterable with progress reporting.

    Args:
        iterable: Iterable to wrap
        total: Total count (if known)
        desc: Progress description
        interval_pct: Update interval percentage
        enabled: Whether to show progress

    Yields:
        Items from the iterable
    """
    if total is None:
        try:
            total = len(iterable)
        except TypeError:
            total = 0
            enabled = False

    progress = ProgressReporter(
        total=total,
        desc=desc,
        interval_pct=interval_pct,
        enabled=enabled
    )

    for item in iterable:
        yield item
        progress.update()

    progress.finish()
