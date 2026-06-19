# FUSILLI

## Fusion Utility for Scanning and Identification of Library Linked Interactions

A Snakemake pipeline for analyzing mutational scanning data from fusion protein libraries.

## Overview

FUSILLI is a tool for analyzing mutational scanning data from fusion protein libraries. It was designed for a project building off work studying MET kinase domain fusions and exon skipping in disrupting signalling pathways using deep mutational scanning (DMS) ([Estevam et al., 2024](https://elifesciences.org/articles/91619), [Estevam et al., 2025](https://elifesciences.org/articles/101882)). This tool was created to support the analysis of a library consisting of a variety of **fusion partner** domains (TPR, CCDC6, etc.) fused to a variably truncated **retained domain** (e.g., the MET kinase domain).

FUSILLI processes short-read sequencing data from libraries of fusion constructs (e.g., kinase domain fusions) and produces counts of detected variant-specific fusion junctions. It applies the deep mutational scanning (DMS) approach — a pooled library read out by NGS — where counting the number and identity of specific variants is the primary goal.

```(markdown)
┌─────────────────────────────────────────────────────────────────────────┐
│                        FUSION LIBRARY STRUCTURE                         │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│   A typical kinase fusion:                                              │
│                                                                         │
│   [ fusion partner ]──[ linker ]──[ retained domain ]                   │
│         ↑                              ↑                                │
│    variable (TPR, CCDC6)          constant (e.g. MET kinase)            │
│    (junction →)                                                         │
│                                                                         │
│   Example: TPR-MET fusion                                               │
│   [TPR]──[GS]──[MET kinase domain (truncated)]                          │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

## Quick Start

> **New to FUSILLI?** Follow the [hands-on tutorial](docs/tutorial.md) — it runs the bundled test dataset end to end (no real data needed) and shows you how to read the results.

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
- A reference FASTA containing the retained domain and the fusion partner sequences

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
   - Configure your fusion library (retained domain, fusion partners, linker)

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
   - FASTA file with sequences for the retained domain and all fusion partners
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

Only four top-level keys are required; everything else is optional with defaults:

```yaml
experiment: 'my_experiment'              # names output paths: results/{experiment}/
data_dir: '/path/to/fastq/files'
samples_file: 'config/samples.csv'

fusion_library:
  retained:                              # the retained / constant domain
    name: 'Met_WT'                       # must match a FASTA header
    truncated_component: 'retained'      # 'retained' (default; 3′ side) or 'partner' (5′ side) — see Breakpoint Conventions
  linker_sequence: 'GGGAGC'              # or '' for none
  partners_file: 'config/fusion_partners.csv'
  sequences_file: 'kinase_sequences.fasta'   # filename inside ref_dir
```

The optional sections — `detection` (window/frame/orientation), `preprocessing`,
`sequencing`, `qc`, `resources`, `quick` (subsampling), plus `variant_retained` and
`unfused_sequences_file` — all have sensible defaults. To avoid drift, this README does
**not** restate every option. Two sources are authoritative:

- **`config/examples/config.yaml`** — every option, annotated, ready to copy.
- **`workflow/schemas/config.schema.yaml`** — the schema (types, defaults, allowed values) the pipeline validates against at startup.

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
│   ├── junction_sequences.csv     # All possible junction k-mers
│   └── domain_ends.csv            # Partner 3' ends for pre-filtering
├── counts/
│   └── {sample}.fusion_counts.csv # Per-sample fusion counts (one per read pair)
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

- Fusion partner name
- Breakpoint coordinate (nucleotides from the start of the truncated component)
- Retained domain name

**Note:** The breakpoint coordinate is relative to whichever component is truncated (set by `fusion_library.retained.truncated_component`). With `truncated_component: 'retained'`, it is measured from the retained domain start; with `truncated_component: 'partner'`, from the fusion partner start.

## Breakpoint Conventions

### Naming Convention

A fusion ID like `TPR_126_Met_WT` means:

- **Fusion partner:** TPR
- **Breakpoint:** coordinate 126 — nucleotides from the start of the truncated component (see note)
- **Retained domain:** Met_WT (the MET kinase domain)

**Note:** The breakpoint coordinate is measured from the start of whichever component `truncated_component` truncates:

- `truncated_component: 'retained'` (the default; canonical for kinase fusions) — from the **retained domain** start, indicating how much of its N-terminus is removed (e.g. exon 14 skipping); the fusion partner is full-length.
- `truncated_component: 'partner'` — from the **fusion partner** start; the fusion partner is variably truncated and the retained domain is full-length.

Both modes detect correctly; choose the value that matches how your library is constructed.

### Junction k-mers

For each candidate junction, the pipeline generates a k-mer spanning the join to search for in the reads. The window is centered on the single junction between the 5′ side (fusion partner + linker) and the 3′ side (retained domain):

```(markdown)
fusion construct:  [ ...fusion partner... ][ linker ][ retained domain... ]
                                              ▲ junction (end of linker)

junction k-mer (window = w):   w nt ending at the junction | w nt after it
                               └────────── 2·w nt total ──────────┘
```

The linker is **not** added on top of the window — when `w` is at least the linker length, the linker sits entirely within the 5′ half. So each k-mer is exactly `2 × breakpoint_window` nt, independent of the linker length. For example, `breakpoint_window: 8` with the `GGGAGC` linker yields a 16 nt k-mer (its 3′ 6 nt being the linker); the default `breakpoint_window: 12` yields a 24 nt k-mer.

Changing `breakpoint_window` changes the k-mer size: a larger window is more specific but can miss reads that don't span the full window; a smaller window is more sensitive but less specific.

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

For any issues, please open an issue on the GitHub repository. For
questions or feedback, [email Chris](https://www.waymentsteelelab.org).