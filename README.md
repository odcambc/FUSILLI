# FUSILLI

**Fusion Utility for Scanning and Identification of Library Linked Interactions**

A Snakemake pipeline for analyzing mutational scanning data from fusion protein libraries.

---

## Overview

FUSILLI processes short-read sequencing data from libraries of fusion constructs (e.g., kinase domain fusions) and produces counts of detected fusion breakpoints. It's designed for experiments where you want to identify which fusion truncation variants are present in a sample.

### How It Works

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        FUSION LIBRARY STRUCTURE                         │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│   A typical kinase fusion:                                              │
│                                                                         │
│   [Partner N-term]────[Linker]────[Kinase Domain (full)]                │
│         ↑                              ↑                                │
│    Variable truncation            Constant anchor                       │
│    (breakpoint)                                                         │
│                                                                         │
│   Example: TPR-MET fusion                                               │
│   [TPR aa 1-142]──[GS]──[MET kinase domain]                            │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

The pipeline:
1. **Preprocesses** reads (adapter trimming, quality filtering, read merging)
2. **Generates** all possible breakpoint sequences for your fusion library
3. **Detects** fusions by string matching breakpoint k-mers against reads
4. **Counts** occurrences of each fusion variant per sample

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

> **Note for Apple Silicon Macs:** Use Rosetta emulation:
> ```bash
> CONDA_SUBDIR=osx-64 conda env create --file fusilli_env.yaml
> ```

### Configuration

1. **Edit the main config file** (`config/config.yaml`):
   - Set `experiment` name
   - Set `data_dir` to your FASTQ directory
   - Configure your fusion library (anchor, partners, linker)

2. **Create your samples file** (`config/samples.csv`):
   ```csv
   sample,condition,replicate,file
   plasmid_rep1,baseline,1,pDNA_S1_L001
   treated_rep1,drug,1,treated_S2_L001
   ```

3. **Create your partners file** (`config/fusion_partners.csv`):
   ```csv
   partner_name,sequence_length,include,description
   TPR,426,true,TPR-MET fusion partner
   CCDC6,303,true,RET fusion partner
   ```

4. **Add reference sequences** to `references/`:
   - FASTA file with sequences for anchor and all partners

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

# Sequencing settings
sequencing:
  paired: true
  min_quality: 30

# Progress reporting
pipeline:
  show_progress: true
  progress_interval: 1       # Update every 1%
```

### Samples File (`config/samples.csv`)

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Unique sample identifier |
| `condition` | Yes | Experimental condition |
| `file` | Yes | FASTQ filename prefix (without `_R1_001.fastq.gz`) |
| `replicate` | No | Replicate number |
| `time` | No | Timepoint |
| `tile` | No | Tile identifier |

### Partners File (`config/fusion_partners.csv`)

| Column | Required | Description |
|--------|----------|-------------|
| `partner_name` | Yes | Must match FASTA header |
| `sequence_length` | Yes | Length in nucleotides |
| `include` | Yes | `true` or `false` |
| `description` | No | Human-readable description |

---

## Output Files

```
results/{experiment}/
├── references/
│   ├── breakpoint_sequences.csv   # All possible breakpoint k-mers
│   └── domain_ends.csv            # Partner 3' ends for pre-filtering
├── counts/
│   └── {sample}.fusion_counts.csv # Per-sample fusion counts
└── fusion_counts_summary.csv      # Aggregated counts matrix
```

### Fusion Counts Format

```csv
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

```
Partner sequence:  ...ATGCTAGCTAGC[BREAKPOINT]
Linker:                            GGGAGC
Anchor sequence:                          ATGAAAAAG...

Breakpoint k-mer (window=12):
              GCTAGCGGGAGCATGAAAAA
              ←─12nt─→←─12nt─→
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

### Performance Tuning

- **Memory:** Adjust `resources.memory_mb` in config
- **Threads:** Adjust `resources.threads` in config
- **Progress:** Disable with `pipeline.show_progress: false` for batch jobs

---

## Troubleshooting

### Common Issues

**"Partner X not found in sequences"**
- Ensure partner name in CSV exactly matches FASTA header
- Check for trailing whitespace

**"Sequence length mismatch"**
- Update `sequence_length` in partners CSV to match actual FASTA sequence

**Slow performance**
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

## Contributing

Contributions welcome! Please submit issues or pull requests on GitHub.
