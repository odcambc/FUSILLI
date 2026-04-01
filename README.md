# FUSILLI

## Fusion Utility for Scanning and Identification of Library Linked Interactions

A Snakemake pipeline for analyzing mutational scanning data from fusion protein libraries.

## Overview

FUSILLI is a tool for analyzing mutational scanning data from fusion protein libraries. It was designed for a project building off work studying MET kinase domain fusions and exon skipping in disrupting signalling pathways using deep mutational scanning (DMS) ([Estevam et al., 2024](https://elifesciences.org/articles/91619), [Estevam et al., 2025](https://elifesciences.org/articles/101882)). This tool was created to support the analysis of a library consisting of a variety of domains (TPR, CCDC6, etc.) fused to a variably truncated anchor domain (e.g, MET kinase).

FUSILLI processes short-read sequencing data from libraries of fusion constructs (e.g., kinase domain fusions) and produces counts of detected variants-specific fusion breakpoints. It's intended to be used in DMS-type experiments where counting the number and identity of specific variants is the primary goal.

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

## Quick Start

```bash
git clone https://github.com/odcambc/FUSILLI
cd FUSILLI
conda env create --file fusilli_env.yaml
conda activate FUSILLI
```

If the environment installed and activated properly, edit the configuration files in the config directory as needed. Then run the pipeline with:

```bash
snakemake -s workflow/Snakefile --software-deployment-method conda --cores 16
```

## Installation

### Install from GitHub

Download or fork this repository and edit the configuration files as needed.

### Install from Docker

A docker image is not currently available, but is planned for the future.

## Configuration

The details of an experiment are specified in a YAML configuration file plus a few supporting data files:

- A main config file describing paths, library structure, and analysis options
- A samples CSV describing experimental conditions and FASTQ file prefixes
- A fusion partners CSV listing the partner domains to include
- A reference FASTA containing the anchor and partner sequences

The active config schema is defined in `workflow/schemas/config.schema.yaml`.

Repository examples:

- `config/examples/config.yaml`: fully annotated template
- `config/examples/test.yaml`: compact minimal example
- `config/test.yaml`: repo-local smoke-test config using the fixtures under `tests/`

### Working directory structure

```(markdown)
├── workflow
│   ├── rules
│   ├── envs
│   ├── scripts
│   ├── schemas
│   │   └── config.schema.yaml
│   └── Snakefile
├── config
│   ├── config.yaml
│   ├── test.yaml
│   └── examples/
├── logs
│   └── ...
├── references
│   └── kinase_sequences_hov.fasta
├── results
│   └── ...
├── stats
│   └── ...
├── tests
│   └── ...
├── resources
│   ├── adapters.fa
│   ├── sequencing_artifacts.fa.gz
│   └── ...
```

### Configuration details

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
   partner_name,include,description
   TPR,true,TPR-MET fusion partner
   CCDC6,true,RET fusion partner
   ```

   - `partner_name`: Must match FASTA header (sequence length is derived from the reference)
   - `include`: Whether to include the partner in the analysis (boolean)
   - `description`: Human-readable description (optional)

4. **Add reference sequences** to `references/`:
   - FASTA file with sequences for anchor and all partners
   - Sequence names must match `partner_name` in `fusion_partners.csv`

### Run the pipeline

```bash
# Dry run (see what will execute)
snakemake -s workflow/Snakefile -n

# Full run with 16 cores
snakemake -s workflow/Snakefile --cores 16

# Full run using conda for environment management
snakemake -s workflow/Snakefile --software-deployment-method conda --cores 16
```

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
    truncated_component: 'partner'  # Breakpoint positions measured from 'partner' (default) or 'anchor'
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

# QC settings
qc:
  run_qc: false
  baseline_condition: 'baseline'
  mem_fastqc: 4000

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

| Column         | Required | Description                                           |
| --------------- | -------- | ----------------------------------------------------- |
| `partner_name` | Yes      | Must match FASTA header (length derived from reference) |
| `include`      | Yes      | `true` or `false`                                    |
| `description`  | No       | Human-readable description                            |

## Output Files

```(markdown)
benchmarks/{experiment}/           # Benchmarking: time and usage of each rule
└── ...

logs/{experiment}/                 # Logs from each rule
├── bbduk/
├── bbmerge/
└── ...

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

### QC Metrics Files

The pipeline generates QC metrics files in `results/{experiment}/`:

- **`fusion_qc_metrics.csv`**: Per-sample fusion detection metrics including:
  - Detection efficiency (reads matched / reads processed)
  - Library coverage (variant, breakpoint, partner coverage fractions)
  - Diversity indices (Shannon, Simpson, evenness)
  - Detection yield (detections per read, detections per million)
  - Top N fractions (top 1, top 10 variant fractions)

- **`sensitivity_metrics.csv`**: Sensitivity analysis metrics:
  - Read length statistics
  - Expected detection fraction (fraction of reads long enough for breakpoint detection)
  - Sensitivity index (actual detections / expected detections)
  - Overlap statistics from read merging

- **`decay_metrics.csv`**: Read retention through preprocessing steps:
  - Read and base counts at each step (raw, trimmed, cleaned, quality-filtered, merged)
  - Retention fractions at each step

- **`trim_metrics.csv`**, **`contam_metrics.csv`**, **`quality_metrics.csv`**: Step-specific preprocessing metrics

- **`partner_counts_summary.csv`**: Partner domain detection counts across samples

All metrics are automatically aggregated into a MultiQC report at `stats/{experiment}/{experiment}_multiqc.html`.

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
- Breakpoint position (nucleotides from the start of the truncated component)
- Anchor name

**Note:** The breakpoint position is relative to whichever component is truncated (as specified by `fusion_library.anchor.truncated_component` in the config). When `truncated_component: 'anchor'`, the breakpoint is measured from the anchor start. When `truncated_component: 'partner'`, it's measured from the partner start.

## Breakpoint Conventions

### Naming Convention

A fusion ID like `TPR_126_Met_WT` means:

- **Partner:** TPR
- **Breakpoint:** 126 nucleotides from anchor start (e.g., amino acids 1-42 of MET missing) when `truncated_component: 'anchor'`
- **Anchor:** Met_WT (kinase domain, variably truncated)

**Note:** The breakpoint position interpretation depends on the `truncated_component` configuration setting. In the default configuration (`truncated_component: 'anchor'`), the breakpoint represents nucleotides from the anchor start, indicating how much of the anchor's N-terminal region is truncated.

### Breakpoint Sequences

For each breakpoint, the pipeline generates a k-mer spanning the junction to look for in the reads:

```(markdown)
Partner sequence:  ...ATGCTAGCTAGC[BREAKPOINT]
Linker:                            GGGAGC
Anchor sequence:                          ATGAAAAAG...

Breakpoint k-mer (window=8):
              TGCTAGCGGGAGCATGAAAAA
              ←─8 nt─→     ←─8 nt─→
```

Changing the `breakpoint_window` configuration will change the size of the k-mer generated, with a larger window allowing for more specific breakpoint detection but potentially missing reads that don't span the full window. The length of the k-mer is twice the `breakpoint_window` value plus the length of the linker sequence. In the above example, with a `breakpoint_window` of 8 and a linker sequence of `GGGAGC`, the k-mer length is 16 + 6 = 22.

By default, the window is set to 12 nt. Depending on your particular sequencing data and library construction, you may need to adjust this value. 12 nt is a reasonable default for our libraries, but consider how your library is constructed.

In general, a larger window will increase specificity but decrease sensitivity. In our experience, within a reasonable range, the window size is not particularly critical to the accuracy of the results.

## Reproducibility

The pipeline automatically captures reproducibility metadata for each run, stored in `results/{experiment}/repro/`. This includes:

- **metadata.json**: Machine-readable metadata (JSON format)
- **metadata.txt**: Human-readable metadata with command-line invocation, versions, and OS info
- **conda-env.yaml**: Conda environment export (if conda is available)
- **pip-freeze.txt**: Pip freeze output (if pip is available)

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

## Optional Modes

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

### Performance Tuning

- **Memory:** Adjust `resources.memory_mb` in config
- **Threads:** Adjust `resources.threads` in config
- **Progress:** Disable with `pipeline.show_progress: false` for batch jobs
- **QC:** `qc.run_qc` defaults to `false`; adjust `qc.mem_fastqc` if needed.

## Limitations

There are currently some limitations that are worth noting:

- The pipeline is not designed for long-read data. We intend to support this in the future but it is not currently implemented.
- The pipeline expects paired-end data. Single-end short-read data is currently not supported.

## License

This is licensed under the MIT license. See the [LICENSE](LICENSE) file for details.

## Contributing

Contributions and feedback are welcome. Please submit an issue or pull request.

## Getting help

For any issues, please open an issue on the GitHub repository.
