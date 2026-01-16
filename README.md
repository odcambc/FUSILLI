# FUSILLI

## Fusion Utility for Scanning and Identification of Library Linked Interactions

A Snakemake pipeline for analyzing mutational scanning data from fusion protein libraries.

---

## Overview

FUSILLI is a tool for analyzing mutational scanning data from fusion protein libraries. It was designed for a project building off work studying MET kinase domain fusions and exon skipping in disrupting signalling pathways using deep mutational scanning (DMS) ([Estevam et al., 2024](https://elifesciences.org/articles/91619), [Estevam et al., 2025](https://elifesciences.org/articles/101882)). This tool was created to support the analysis of a library consisting of a variety of domains (TPR, CCDC6, etc.) fused to a variably truncated anchor domain (e.g, MET kinase).

FUSILLI processes short-read sequencing data from libraries of fusion constructs (e.g., kinase domain fusions) and produces counts of detected variants-specific fusion breakpoints. It's intended to be used in DMS-type experiments where counting the number and identity of specific variants is the primary goal.

### How It Works

```(markdown)
┌─────────────────────────────────────────────────────────────────────────┐
│                        FUSION LIBRARY STRUCTURE                         │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│   A typical kinase fusion:                                              │
│                                                                         │
│   [Partner N-term]────[Linker]────[Kinase Domain (truncated)]           │
│         ↑                              ↑                                │
│    Full length partner            Variable truncated anchor             │
│    (breakpoint)                                                         │
│                                                                         │
│   Example: TPR-MET fusion                                               │
│   [TPR]──[GS]──[MET kinase domain (truncated)]                          │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

The pipeline:

1. **Preprocesses** reads (adapter trimming, quality filtering, read merging)
2. **Identifies** all possible breakpoint sequences for your fusion library
3. **Detects** fusions by string matching breakpoint k-mers against reads
4. **Counts** occurrences of each fusion variant per sample
5. **Generates** a summary of the counts across all samples, statistics, and reproducibility metadata

---

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/odcambc/FUSILLI
cd FUSILLI

# Create conda environment
conda env create --file fusilli_env.yaml
conda activate FUSILLI
```

> **Note for Apple Silicon Macs:** It may be necessary to use Rosetta emulation to build the conda environment:
>
> ```bash
> CONDA_SUBDIR=osx-64 conda env create --file fusilli_env.yaml
> ```

### Configuration

1. **Edit the main config file** (`config/config.yaml`):
   - Set `experiment` name
   - Set `data_dir` to your FASTQ directory
   - Configure your fusion library (anchor, partners, linker)

2. **Create your samples file** (`config/samples.csv`):
   This file defines the sequencing samples to be processed. It should contain the following columns:

   ```csv
   sample,condition,replicate,file
   plasmid_rep1,baseline,1,pDNA_S1_L001
   treated_rep1,drug,1,treated_S2_L001
   ```

   - `sample`: Unique sample identifier
   - `condition`: Experimental condition
   - `replicate`: Replicate number (optional)
   - `time`: Timepoint (optional)
   - `tile`: Tile identifier for tiled amplicon sequencing (optional)
   - `file`: FASTQ filename prefix (without `_R1_001.fastq.gz`, `_R2_001.fastq.gz`)

3. **Create your partners file** (`config/fusion_partners.csv`):

   ```csv
   partner_name,sequence_length,include,description
   TPR,426,true,TPR-MET fusion partner
   CCDC6,303,true,RET fusion partner
   ```

   - `partner_name`: Must match FASTA header
   - `sequence_length`: Length in nucleotides
   - `include`: Whether to include the partner in the analysis (boolean)
   - `description`: Human-readable description (optional)

4. **Add reference sequences** to `references/`:
   - FASTA file with sequences for anchor and all partners
   - Sequence names must match `partner_name` in `fusion_partners.csv`

> **Note:** The pipeline expects paired-end data. Single-end data is currently not supported.

### Run

```bash
# Dry run (see what will execute)
snakemake -s workflow/Snakefile -n

# Full run
snakemake -s workflow/Snakefile --cores 16

# With conda environment management
snakemake -s workflow/Snakefile --software-deployment-method conda --cores 16
```

---

## Configuration Reference

### Main Config (`config/config.yaml`)

```yaml
# Unique experiment identifier
experiment: 'my_experiment'

# Data locations
data_dir: '/path/to/fastq/files'
ref_dir: 'references'
samples_file: 'config/samples.csv'

# Fusion library definition
fusion_library:
  anchor:
    name: 'Met_WT'           # Constant domain in all fusions
    position: 'downstream'    # 'upstream' or 'downstream'
  linker_sequence: 'GGGAGC'  # Linker between domains (or '')
  partners_file: 'config/fusion_partners.csv'
  sequences_file: 'kinase_sequences.fasta'

# Detection parameters
detection:
  method: 'string'           # Detection algorithm
  breakpoint_window: 12      # nt on each side of breakpoint
  maintain_frame: true       # Only in-frame breakpoints
  kmer_size: 15              # K-mer size for pre-filtering
  unmerged_detection: false # Also scan unmerged R1/R2 reads separately (optional)

# Sequencing settings
sequencing:
  paired: true
  min_quality: 30

# Progress reporting
pipeline:
  show_progress: true
  progress_interval: 1       # Update every 1%

# Quick-mode subsampling (optional)
quick:
  enabled: false        # Set true for fast sanity-check runs
  max_reads: 100000     # Cap per mate when enabled
  fraction: null        # Optional samplerate override (e.g., 0.01)
  seed: 1337            # Deterministic subsampling seed
```

### Samples File (`config/samples.csv`)

| Column      | Required | Description                                        |
| ----------- | -------- | -------------------------------------------------- |
| `sample`    | Yes      | Unique sample identifier                           |
| `condition` | Yes      | Experimental condition                             |
| `file`      | Yes      | FASTQ filename prefix (without `_R1_001.fastq.gz`) |
| `replicate` | No       | Replicate number                                   |
| `time`      | No       | Timepoint                                          |
| `tile`      | No       | Tile identifier (for tiled amplicon sequencing)    |

### Partners File (`config/fusion_partners.csv`)

| Column            | Required | Description                |
| ----------------- | -------- | -------------------------- |
| `partner_name`    | Yes      | Must match FASTA header    |
| `sequence_length` | Yes      | Length in nucleotides      |
| `include`         | Yes      | `true` or `false`          |
| `description`     | No       | Human-readable description |

---

## Output Files

```(markdown)
results/{experiment}/
├── references/
│   ├── breakpoint_sequences.csv   # All possible breakpoint k-mers
│   └── domain_ends.csv            # Partner 3' ends for pre-filtering
├── counts/
│   ├── {sample}.fusion_counts.csv # Per-sample fusion counts (merged reads)
│   ├── {sample}.R1.unmerged_fusion_counts.csv  # R1 unmerged counts (if enabled)
│   └── {sample}.R2.unmerged_fusion_counts.csv  # R2 unmerged counts (if enabled)
├── repro/                         # Reproducibility metadata
│   ├── metadata.json              # Machine-readable metadata
│   ├── metadata.txt               # Human-readable metadata
│   ├── conda-env.yaml             # Conda environment export (if available)
│   └── pip-freeze.txt             # Pip freeze output (if available)
├── fusion_counts_summary.csv      # Aggregated counts matrix (merged reads)
├── unmerged_counts_summary.csv    # Aggregated unmerged counts (if enabled)
└── unmerged_partner_counts_summary.csv  # Aggregated unmerged partner counts (if enabled)

stats/{experiment}/
├── {experiment}_multiqc.html       # MultiQC report
└── (tool name))/                   # Additional intermediate tool reports
```

### Fusion Counts Format

```(markdown)
fusion_id,count
TPR_426_Met_WT,15234
TPR_423_Met_WT,14892
CCDC6_303_Met_WT,8921
...
```

The `fusion_id` encodes:

- Partner name
- Breakpoint position (nucleotides from partner start)
- Anchor name

---

## Understanding Breakpoints

### Naming Convention

A fusion ID like `TPR_126_Met_WT` means:

- **Partner:** TPR
- **Breakpoint:** 126 nucleotides from TPR start (= amino acid 42)
- **Anchor:** Met_WT (full kinase domain)

### Breakpoint Sequences

For each breakpoint, the pipeline generates a k-mer spanning the junction:

```(markdown)
Partner sequence:  ...ATGCTAGCTAGC[BREAKPOINT]
Linker:                            GGGAGC
Anchor sequence:                          ATGAAAAAG...

Breakpoint k-mer (window=12):
              GCTAGCGGGAGCATGAAAAA
              ←─12nt─→←─12nt─→
```

---

## Reproducibility

The pipeline automatically captures reproducibility metadata for each run, stored in `results/{experiment}/repro/`. This includes:

- **metadata.json**: Machine-readable metadata (JSON format)
- **metadata.txt**: Human-readable metadata with command-line invocation, versions, and OS info
- **conda-env.yaml**: Conda environment export (if conda is available)
- **pip-freeze.txt**: Pip freeze output (if pip is available)

### Viewing Reproducibility Information

```bash
# View human-readable metadata
cat results/{experiment}/repro/metadata.txt

# View conda environment (if available)
cat results/{experiment}/repro/conda-env.yaml

# View pip packages (if available)
cat results/{experiment}/repro/pip-freeze.txt
```

### Regenerating Reproducibility Files

If you need to regenerate reproducibility files for an existing run:

```bash
# Run the reproducibility capture script directly
python workflow/scripts/capture_reproducibility.py results/{experiment}/repro
```

Or regenerate as part of the pipeline:

```bash
# Regenerate reproducibility files only
snakemake -s workflow/Snakefile results/{experiment}/repro/metadata.json --cores 1
```

---

## Advanced Usage

### Running Specific Targets

```bash
# Generate reference files only (validate config)
snakemake -s workflow/Snakefile references --cores 1

# Generate counts without QC
snakemake -s workflow/Snakefile counts_only --cores 16

# Generate summary only (if counts exist)
snakemake -s workflow/Snakefile summary --cores 1
```

### Standalone Script Usage

The core scripts can be used outside Snakemake:

```bash
# Generate breakpoint sequences
python workflow/scripts/fusion_sequences.py \
    --sequences references/kinase_sequences.fasta \
    --partners config/fusion_partners.csv \
    --anchor Met_WT \
    --linker GGGAGC \
    --window 12 \
    --output-breakpoints breakpoints.csv \
    --output-ends ends.csv

# Run string matching
python workflow/scripts/string_matcher.py \
    --input reads.fastq.gz \
    --breakpoints breakpoints.csv \
    --ends ends.csv \
    --output counts.csv \
    --progress
```

### Unmerged Read Processing

By default, FUSILLI processes merged (error-corrected) reads for fusion detection. However, you can also enable processing of unmerged reads separately. This is useful when:

- You want to capture fusions that may be missed in merged reads
- You need to analyze R1 and R2 reads independently
- You want to compare detection rates between merged and unmerged reads

To enable unmerged read processing, set in your `config/config.yaml`:

```yaml
detection:
  unmerged_detection: true  # Enable processing of unmerged R1/R2 reads
```

When enabled, the pipeline will:

1. **Process R1 and R2 separately**: Each unmerged mate is processed independently using the same string matching algorithm as merged reads
2. **Keep counts distinct**: Unmerged counts are stored in separate files with the naming pattern:
   - `{sample}.R1.unmerged_fusion_counts.csv` - Counts from R1 unmerged reads
   - `{sample}.R2.unmerged_fusion_counts.csv` - Counts from R2 unmerged reads
3. **Generate separate summaries**: Aggregated unmerged counts are written to:
   - `unmerged_counts_summary.csv` - Combined R1 and R2 counts per sample
   - `unmerged_partner_counts_summary.csv` - Partner-level unmerged counts

**Important notes:**

- Unmerged counts are **always kept separate** from merged counts - they are never combined
- Empty unmerged files (when all reads merge successfully) are handled gracefully
- The same detection parameters (`breakpoint_window`, `kmer_size`, etc.) are used for both merged and unmerged detection
- Processing unmerged reads increases computational time and storage requirements

### Performance Tuning

- **Memory:** Adjust `resources.memory_mb` in config
- **Threads:** Adjust `resources.threads` in config
- **Progress:** Disable with `pipeline.show_progress: false` for batch jobs
- **QC:** `qc.run_qc` defaults to `true`; adjust `mem_fastqc` if needed.

---

## Troubleshooting

### Common Issues

#### "Partner X not found in sequences"

- Ensure partner name in CSV exactly matches FASTA header
- Check for trailing whitespace

#### "Sequence length mismatch"

- Update `sequence_length` in partners CSV to match actual FASTA sequence

#### *Slow performance

- Install `pyfastx` for faster FASTQ parsing: `pip install pyfastx`
- Reduce `breakpoint_window` if appropriate for your read length

### Getting Help

For issues, please open a GitHub issue with:

1. Your config file (anonymized if needed)
2. The error message
3. Snakemake version (`snakemake --version`)

---

## License

MIT License - see [LICENSE](LICENSE) for details.

## Citation

If you use FUSILLI in your research, please cite:
> [Citation to be added upon publication]

## Documentation

Additional documentation is available in the `docs/` directory:

- **[Architecture](docs/ARCHITECTURE.md)** - System architecture and design decisions
- **[Technical Patterns](docs/TECHNICAL.md)** - Coding conventions and best practices
- **[Project Management](docs/PROJECT_MANAGEMENT.md)** - Branch management, task design, and workflow organization

## Contributing

Contributions welcome! Please submit issues or pull requests on GitHub.

Before contributing, please review:

- [Technical Patterns](docs/TECHNICAL.md) for coding conventions
- [Project Management Guide](docs/PROJECT_MANAGEMENT.md) for workflow and task organization
