# Architecture Overview

## System Components

### Workflow Layer (Snakemake)

**Entry Point:** `workflow/Snakefile`
- Main pipeline orchestrator
- Loads configuration and includes rule modules
- Defines target rules (`all`, `counts_only`, `summary`, `references`)

**Rule Modules:** `workflow/rules/*.smk`
- `common.smk` - Configuration parsing, helper functions, sample loading
- `subsample.smk` - Quick-mode subsampling (optional)
- `filter_paired.smk` - Read preprocessing (bbduk, bbmerge)
- `process_strings.smk` - Fusion detection via string matching
- `qc.smk` - Quality control (FastQC, MultiQC)
- `reproducibility.smk` - Metadata capture

**Scripts:** `workflow/scripts/*.py`
- `fusion_sequences.py` - Generate breakpoint reference sequences
- `string_matcher.py` - Detect fusions in reads
- `capture_reproducibility.py` - Capture run metadata
- `utils.py` - Shared utilities

### Configuration Layer

**Main Config:** `config/config.yaml`
- Single source of truth for experiment parameters
- Validated against `workflow/schemas/config.schema.yaml`
- Defines experiment, samples, fusion library, detection parameters

**Supporting Config Files:**
- `config/samples.csv` - Sample metadata
- `config/fusion_partners.csv` - Fusion partner definitions
- `workflow/schemas/*.schema.yaml` - JSON schemas for validation

### Data Flow

```
┌─────────────────────────────────────────────────────────────┐
│                    INPUT DATA                               │
│  FASTQ files (paired-end) from data_dir/                    │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│              PREPROCESSING (filter_paired.smk)              │
│  • Adapter trimming (bbduk)                                 │
│  • Quality filtering                                        │
│  • Read merging (bbmerge)                                   │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│         REFERENCE GENERATION (process_strings.smk)          │
│  • Load partner and anchor sequences                        │
│  • Generate all breakpoint k-mers (fusion_sequences.py)     │
│  • Generate domain end k-mers for pre-filtering             │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│          DETECTION (process_strings.smk)                    │
│  • Pre-filter: Check for partner domain ends                │
│  • Match: Search for specific breakpoint sequences          │
│  • Count: Aggregate matches per fusion variant              │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│                    OUTPUT                                   │
│  • Per-sample fusion counts                                 │
│  • Aggregated summary matrix                                │
│  • QC metrics and reports                                   │
│  • Reproducibility metadata                                 │
└─────────────────────────────────────────────────────────────┘
```

## Key Design Decisions

### Why Snakemake?
- **Reproducibility:** Automatic dependency tracking and execution
- **Parallelization:** Built-in support for parallel execution
- **Resume capability:** Can resume from failures
- **Environment management:** Integration with conda/mamba
- **Scientific workflow standard:** Common in bioinformatics

### Why String Matching?
- **Speed:** Fast for known breakpoint sequences
- **Simplicity:** No alignment overhead
- **Sufficient:** Works well for library-based fusion detection
- **Deterministic:** Exact matches, no scoring ambiguity
- **Robust:** Low risk of false positives

### Why Modular Rules?
- **Maintainability:** Each module has a clear purpose
- **Testability:** Rules can be tested independently
- **Reusability:** Common patterns extracted to `common.smk`
- **Clarity:** Easy to understand pipeline structure

### Why Standalone Scripts?
- **Testability:** Scripts can be tested outside Snakemake
- **Reusability:** Can be used in other contexts
- **Debugging:** Easier to debug standalone than embedded code
- **Documentation:** Clear command-line interfaces

## Module Boundaries

### Configuration
- **Single source of truth:** `config/config.yaml`
- **Validation:** Automatic via Snakemake schema validation
- **Loading:** Centralized in `common.smk`
- **Access:** Via global variables (e.g., `EXPERIMENT`, `SAMPLES`)

### Data Processing
- **Input:** FASTQ files from `data_dir/`
- **Intermediate:** Stored in `results/{experiment}/`
- **Output:** Final counts and summaries in `results/{experiment}/`
- **Cleanup:** Utility rules for removing intermediates

### Statistics and logs
- **Statistics:** Stored in `stats/{experiment}/`
- **Logs:** Stored in `logs/{experiment}/`
- **Aggregation:** Intermediate tool and script outputs are aggregated with MultiQC and saved in `stats/{experiment}/`

### Scripts
- **Independence:** Scripts are standalone and importable
- **Parameters:** Accept command-line arguments
- **Snakemake integration:** Can access `snakemake` object when run via Snakemake
- **Utilities:** Shared code in `utils.py`

### Rules
- **No direct script imports:** Rules call scripts, don't import them
- **Helper functions:** Defined in `common.smk` for reuse
- **Wildcards:** Used for sample-level parallelization
- **Resources:** Configurable memory and thread limits

## Component Interactions

### Configuration → Rules
- Config loaded once in `Snakefile`
- Parsed in `common.smk`
- Available as global variables to all rules

### Rules → Scripts
- Rules call scripts via `script:` directive or `shell:` with Python
- Parameters passed via Snakemake wildcards and config
- Outputs defined in rule `output:`

### Scripts → Utils
- Scripts import from `utils.py`
- Utils provide common functionality (FASTA parsing, logging, etc.)
- Utils are pure functions where possible

### Rules → Rules
- Dependencies defined via `input:` and `output:`
- Snakemake automatically resolves dependencies
- Parallel execution where possible

## Extension Points

### Adding New Detection Methods
1. Create new script in `workflow/scripts/` (e.g., `bowtie_matcher.py`)
2. Add new rule module `workflow/rules/process_bowtie.smk`
3. Update config schema to include new method
4. Conditionally include rule module in `Snakefile`
5. Update `common.smk` to handle new method

### Adding New QC Tools
1. Add rule to `workflow/rules/qc.smk`
2. Update `get_all_targets()` in `common.smk` if needed
3. Add config options for enabling/disabling
4. Incorporate new QC tools into MultiQC report

### Adding New Preprocessing Steps
1. Add rules to `workflow/rules/filter_paired.smk` or new module
2. Update dependencies in downstream rules
3. Add config options if needed

## Data Structures

### Fusion ID Format
`{partner}_{breakpoint_nt}_{anchor}`
- Example: `TPR_126_Met_WT`
- Partner: TPR
- Breakpoint: 126 nucleotides from partner start
- Anchor: Met_WT

### Breakpoint Sequences
- K-mer spanning breakpoint junction
- Window size configurable (default: 12 nt each side)
- Includes linker sequence if present

### Count Files
- CSV format: `fusion_id,count`
- One file per sample
- Aggregated into summary matrix

## Performance Considerations

### Parallelization
- Sample-level parallelization via Snakemake wildcards
- Configurable thread and memory limits per rule
- Automatic resource management

### Memory Usage
- Pre-filtering reduces search space
- Streaming processing for large FASTQ files
- Configurable memory limits per rule

### I/O Optimization
- Compressed FASTQ files (gzip)
- Intermediate files kept for debugging
- Cleanup rules available

## Error Handling

### Configuration Errors
- Schema validation catches invalid config early
- Clear error messages for missing required fields
- Type checking via schema

### Runtime Errors
- Snakemake logs all rule execution
- Scripts log errors with context
- Failed rules don't block unrelated work

### Data Errors
- Validation of sample files
- Validation of partner files
- Sequence length checks

## Testing Strategy

### Unit Tests
- Test scripts independently (`tests/test_*.py`)
- Mock file I/O where appropriate
- Test edge cases and error conditions

### Integration Tests
- Test full pipeline with sample data
- Verify output format and content
- Test configuration variations

### Regression Tests
- Compare outputs across versions
- Track performance metrics
- Validate reproducibility
