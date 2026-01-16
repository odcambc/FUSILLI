# fusilli-multiqc

Custom MultiQC modules for the FUSILLI fusion detection pipeline.

This package provides four custom MultiQC modules that visualize FUSILLI-specific quality control metrics:

- **Detection Metrics**: Detection efficiency, library coverage, and sensitivity analysis
- **Diversity Metrics**: Library diversity metrics (Shannon, Simpson, evenness)
- **Preprocessing Metrics**: Read retention and loss through preprocessing steps
- **Partner Detection**: Partner detection metrics and coverage across samples

## Installation

### Option 1: Local Development Installation

For development or when the package is in the same repository as the FUSILLI pipeline:

```bash
pip install -e fusilli-multiqc/
```

### Option 2: Installation from Git Repository

If the package is in a separate git repository:

```bash
pip install git+https://github.com/user/fusilli-multiqc.git@main
```

Or for a specific version/tag:

```bash
pip install git+https://github.com/user/fusilli-multiqc.git@v1.0.0
```

### Option 3: Installation from PyPI (Future)

Once published to PyPI:

```bash
pip install fusilli-multiqc
```

## Requirements

- Python >= 3.8
- MultiQC >= 1.0
- pandas >= 1.0
- numpy >= 1.0

## Usage

After installation, the modules will be automatically discovered by MultiQC when you run:

```bash
multiqc <results_directory>
```

The modules will automatically detect and parse the following FUSILLI output files:

- `fusion_qc_metrics.csv` - Detection efficiency and coverage metrics
- `sensitivity_metrics.csv` - Sensitivity analysis data
- `fusion_counts_summary.csv` - Fusion count data for diversity analysis
- `decay_metrics.csv` - Read decay through preprocessing steps
- `partner_counts_summary.csv` - Partner detection counts

## Module Descriptions

### Detection Metrics Module (`fusilli_detection`)

Visualizes:
- **Detection Efficiency**: Fraction of processed reads that matched fusion sequences
- **Prefilter Efficiency**: Fraction of reads containing partner domain 3' ends
- **Matching Efficiency**: Fraction of prefiltered reads that matched specific breakpoints
- **Library Coverage**: Variant, breakpoint, and partner coverage metrics
- **Sensitivity Analysis**: Sensitivity index vs expected detection fraction

### Diversity Metrics Module (`fusilli_diversity`)

Visualizes:
- **Shannon Diversity**: Measures both richness and evenness
- **Simpson Diversity**: Measures dominance
- **Evenness**: Pielou's evenness index
- **Top N Fractions**: Distribution of counts across top variants
- **Variant Count Distribution**: Histogram of variant abundances

### Preprocessing Metrics Module (`fusilli_preprocessing`)

Visualizes:
- **Read Decay**: Read counts at each preprocessing step (log scale)
- **Retention Rates**: Fraction of reads retained after each step
- **Step Loss Breakdown**: Cumulative loss of reads through preprocessing

### Partner Detection Module (`fusilli_partners`)

Visualizes:
- **Partner Detection Heatmap**: Binary detection matrix of partners across samples
- **Partner Coverage**: Number of samples with detection per partner
- **Partner End vs Linker**: Comparison of detection methods

## Development

### Package Structure

```
fusilli-multiqc/
├── fusilli_multiqc/           # Main package directory
│   ├── __init__.py
│   ├── utils.py               # Shared utilities (parse_csv_file, parse_json_file)
│   └── modules/               # MultiQC modules
│       ├── __init__.py
│       ├── detection.py       # Detection metrics module
│       ├── diversity.py       # Diversity metrics module
│       ├── preprocessing.py   # Preprocessing metrics module
│       └── partners.py        # Partner detection module
├── setup.py                   # Setuptools configuration
├── pyproject.toml             # Modern Python packaging
└── README.md                  # This file
```

### Building the Package

To build a distribution package:

```bash
python -m build
```

This will create `dist/` directory with source and wheel distributions.

### Testing Installation

After installation, verify that entry points are registered:

```python
import pkg_resources
entry_points = list(pkg_resources.iter_entry_points('multiqc.modules.v1'))
print([ep.name for ep in entry_points if 'fusilli' in ep.name])
```

You should see:
- `fusilli_detection`
- `fusilli_diversity`
- `fusilli_preprocessing`
- `fusilli_partners`

## Version Management

This package uses semantic versioning (MAJOR.MINOR.PATCH):
- **MAJOR**: Breaking changes to module interfaces or MultiQC compatibility
- **MINOR**: New features or modules (backward compatible)
- **PATCH**: Bug fixes and minor improvements

Current version: **1.0.0**

## License

MIT License - see LICENSE file for details.

## Author

Chris Macdonald

## Related Projects

- [FUSILLI Pipeline](https://github.com/user/fusilli) - Main fusion detection pipeline
- [MultiQC](https://multiqc.info/) - Aggregates bioinformatics analysis results
