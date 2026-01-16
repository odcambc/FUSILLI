# Technical Patterns & Conventions

## Code Style

### Python

**Formatting:**
- Use [Black](https://black.readthedocs.io/) for code formatting
- Maximum line length: 88 characters (Black default)
- Use `isort` for import sorting

**Naming Conventions:**
- Functions and variables: `snake_case`
- Classes: `PascalCase`
- Constants: `UPPER_CASE`
- Private functions: `_leading_underscore`

**Example:**
```python
# Constants
DEFAULT_WINDOW_SIZE = 12

# Functions
def parse_fasta(filepath: str) -> dict:
    """Parse FASTA file and return sequence dictionary."""
    pass

# Classes
class BreakpointGenerator:
    """Generate fusion breakpoint sequences."""
    pass
```

### Snakemake

**Rule Naming:**
- Use descriptive names: `generate_breakpoints`, not `gen_bp`
- Group related rules in modules
- Use consistent prefixes for related rules

**Rule Structure:**
```python
rule rule_name:
    """
    Docstring describing what this rule does.
    """
    input:
        # Input files
    output:
        # Output files
    params:
        # Parameters
    resources:
        # Resource limits
    script:
        # Script path
    log:
        # Log file
```

## Type Hints

**Required for:**
- All function parameters
- All function return values
- Class attributes (where applicable)

**Import from `typing`:**
```python
from typing import Dict, List, Optional, Tuple, Iterator

def process_reads(
    filepath: str,
    breakpoints: Dict[str, str],
    min_quality: Optional[int] = None
) -> List[Tuple[str, int]]:
    """Process reads and return matches."""
    pass
```

**Use `Optional[Type]` instead of `Type | None`:**
```python
# Good
def find_match(sequence: str, threshold: Optional[int] = None) -> Optional[str]:
    pass

# Avoid (for Python < 3.10 compatibility)
def find_match(sequence: str, threshold: int | None = None) -> str | None:
    pass
```

## Documentation

### Docstrings

**Use Google-style docstrings:**
```python
def generate_breakpoints(
    partner_seq: str,
    anchor_seq: str,
    linker: str,
    window: int = 12
) -> List[str]:
    """
    Generate all possible breakpoint k-mers for a fusion.

    Args:
        partner_seq: Partner sequence (nucleotides)
        anchor_seq: Anchor sequence (nucleotides)
        linker: Linker sequence between domains
        window: Number of nucleotides on each side of breakpoint

    Returns:
        List of breakpoint k-mer sequences

    Raises:
        ValueError: If sequences are invalid or window is too large
    """
    pass
```

### Inline Comments

**Use comments to explain "why", not "what":**
```python
# Good: Explains reasoning
# Use 12-nt window to balance specificity and sensitivity
# based on typical read length and error rates
window = 12

# Bad: States the obvious
# Set window to 12
window = 12
```

## Error Handling

### Exception Types

**Use specific exceptions:**
```python
# Good
if not os.path.exists(filepath):
    raise FileNotFoundError(f"Reference file not found: {filepath}")

if sequence_length < window * 2:
    raise ValueError(f"Sequence too short for window size {window}")

# Bad
if not os.path.exists(filepath):
    raise Exception("File not found")
```

### Error Messages

**Include context in error messages:**
```python
# Good
raise ValueError(
    f"Partner '{partner_name}' not found in sequences file. "
    f"Available partners: {', '.join(available_partners)}"
)

# Bad
raise ValueError("Partner not found")
```

### Logging

**Use structured logging:**
```python
import logging

logger = logging.getLogger(__name__)

def process_sample(sample_name: str):
    logger.info(f"Processing sample: {sample_name}")
    try:
        # Process
        logger.debug(f"Found {count} matches")
    except Exception as e:
        logger.error(f"Failed to process {sample_name}: {e}", exc_info=True)
        raise
```

## Testing

### Test Structure

**Follow pytest conventions:**
```python
# tests/test_fusion_sequences.py
import pytest
from workflow.scripts.fusion_sequences import generate_breakpoints

def test_generate_breakpoints_basic():
    """Test basic breakpoint generation."""
    partner = "ATGC" * 10
    anchor = "CGTA" * 10
    linker = "GG"
    breakpoints = generate_breakpoints(partner, anchor, linker, window=2)
    assert len(breakpoints) > 0
    assert all(len(bp) == 6 for bp in breakpoints)  # 2 + 2 + 2

def test_generate_breakpoints_invalid_window():
    """Test error handling for invalid window size."""
    with pytest.raises(ValueError, match="window size"):
        generate_breakpoints("ATGC", "CGTA", "", window=100)
```

### Test Fixtures

**Use `conftest.py` for shared fixtures:**
```python
# tests/conftest.py
import pytest
from pathlib import Path

@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return Path(__file__).parent / "data"

@pytest.fixture
def sample_fasta(tmp_path):
    """Create a temporary FASTA file for testing."""
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(">seq1\nATGC\n>seq2\nCGTA\n")
    return fasta_file
```

### Test Coverage

**Aim for:**
- 80%+ coverage for core functionality
- 100% coverage for error handling paths
- Integration tests for critical workflows

## Performance

### File I/O

**Use streaming for large files:**
```python
# Good: Stream processing
def process_fastq(filepath: str):
    with gzip.open(filepath, 'rt') as f:
        for record in parse_fastq(f):
            yield process_record(record)

# Bad: Load everything into memory
def process_fastq(filepath: str):
    records = list(parse_fastq(filepath))
    return [process_record(r) for r in records]
```

**Use `pyfastx` when available:**
```python
try:
    import pyfastx
    # Fast C-based parsing
    for name, seq, qual in pyfastx.Fastq(filepath):
        process_read(name, seq, qual)
except ImportError:
    # Fallback to pure Python
    for record in parse_fastq(filepath):
        process_record(record)
```

### Data Structures

**Choose appropriate data structures:**
```python
# Good: Set for membership testing
partner_ends = set(domain_ends)
if read_kmer in partner_ends:
    # Pre-filter passed

# Good: Dict for lookups
breakpoint_map = {bp_id: sequence for bp_id, sequence in breakpoints}

# Avoid: List for membership testing
if read_kmer in partner_ends_list:  # O(n) lookup
    pass
```

### Parallelization

**Let Snakemake handle parallelization:**
- Use wildcards for sample-level parallelism
- Configure resources appropriately
- Avoid manual threading in scripts unless necessary

## Common Patterns

### Configuration Loading

**In Snakemake rules:**
```python
# Access config via global variable (set in common.smk)
rule process_sample:
    input:
        fastq = get_input_fastq_r1(wildcards)
    output:
        counts = "results/{experiment}/counts/{sample}.csv"
    script:
        "scripts/process.py"
```

**In Python scripts:**
```python
# Handle both Snakemake and standalone execution
try:
    from snakemake import snakemake
    # Running via Snakemake
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    window = snakemake.params.window
except NameError:
    # Standalone execution
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--window', type=int, default=12)
    args = parser.parse_args()
    input_file = args.input
    output_file = args.output
    window = args.window
```

### Path Handling

**Use `pathlib.Path`:**
```python
from pathlib import Path

# Good
data_dir = Path(config["data_dir"])
fastq_file = data_dir / f"{sample}_R1_001.fastq.gz"
if not fastq_file.exists():
    raise FileNotFoundError(f"FASTQ not found: {fastq_file}")

# Avoid
import os
fastq_file = os.path.join(config["data_dir"], f"{sample}_R1_001.fastq.gz")
if not os.path.exists(fastq_file):
    raise FileNotFoundError(f"FASTQ not found: {fastq_file}")
```

### Progress Reporting

**Use `tqdm` or custom progress reporter:**
```python
from tqdm import tqdm

# For long operations
for read in tqdm(parse_fastq(filepath), desc="Processing reads"):
    process_read(read)
```

**Or use custom progress reporter (from `utils.py`):**
```python
from utils import ProgressReporter

reporter = ProgressReporter(
    total=total_reads,
    desc="Matching breakpoints",
    interval_pct=1,
    enabled=show_progress
)
for read in parse_fastq(filepath):
    process_read(read)
    reporter.update(1)
```

### CSV Handling

**Use `csv` module for proper handling:**
```python
import csv

# Reading
with open(csv_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        process_row(row)

# Writing
with open(output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['fusion_id', 'count'])
    writer.writeheader()
    for fusion_id, count in counts.items():
        writer.writerow({'fusion_id': fusion_id, 'count': count})
```

## Code Organization

### File Size

**Keep files under 300 lines:**
- Split large files into smaller modules
- Extract common functionality to utilities
- Use classes to group related functions

### Import Organization

**Follow isort conventions:**
```python
# Standard library
import csv
import gzip
from pathlib import Path
from typing import Dict, List

# Third-party
import pandas as pd
from Bio import SeqIO

# Local
from utils import parse_fasta, setup_logging
```

### Function Length

**Keep functions focused:**
- Single responsibility
- < 50 lines when possible
- Extract complex logic to helper functions

## Security

### Input Validation

**Validate all user inputs:**
```python
def process_config(config: dict):
    """Process and validate configuration."""
    # Validate required fields
    required = ['experiment', 'data_dir', 'samples_file']
    for field in required:
        if field not in config:
            raise ValueError(f"Missing required config field: {field}")

    # Validate paths
    data_dir = Path(config['data_dir'])
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")

    # Validate types
    if not isinstance(config.get('breakpoint_window', 12), int):
        raise TypeError("breakpoint_window must be an integer")
```

### File Path Safety

**Avoid path traversal:**
```python
# Good: Validate paths are within expected directory
def get_output_path(experiment: str, filename: str) -> Path:
    base_dir = Path("results") / experiment
    output_path = base_dir / filename
    # Ensure path is within base_dir
    if not str(output_path.resolve()).startswith(str(base_dir.resolve())):
        raise ValueError(f"Invalid output path: {filename}")
    return output_path
```

## Debugging

### Logging Levels

**Use appropriate log levels:**
- `DEBUG`: Detailed information for debugging
- `INFO`: General informational messages
- `WARNING`: Warning messages (non-fatal)
- `ERROR`: Error messages (recoverable)
- `CRITICAL`: Critical errors (may abort)

### Debugging Tools

**Use Python debugger when needed:**
```python
import pdb

def complex_function():
    # Set breakpoint
    pdb.set_trace()
    # Or use breakpoint() in Python 3.7+
    breakpoint()
    # Continue execution
```

**Use Snakemake debugging:**
```bash
# Dry run to see what would execute
snakemake -n

# Show detailed execution plan
snakemake --detailed-summary

# Run with debug output
snakemake --debug
```
