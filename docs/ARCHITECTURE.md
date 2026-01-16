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

### Unmerged Read Processing Flow

When `detection.unmerged_detection` is enabled in the configuration, the pipeline processes unmerged reads separately from merged reads. This allows detection of fusions in reads that could not be merged by bbmerge (e.g., due to insufficient overlap or low quality).

**Configuration:**
- Controlled by `detection.unmerged_detection` (default: `false`)
- When enabled, processes both R1 and R2 unmerged reads independently

**Data Flow for Unmerged Processing:**

```
┌─────────────────────────────────────────────────────────────┐
│         PREPROCESSING (filter_paired.smk)                   │
│  • bbmerge outputs:                                          │
│    - {sample}_merged.fastq.gz (merged reads)                │
│    - {sample}_R1.unmerged.fastq.gz (unmerged R1)            │
│    - {sample}_R2.unmerged.fastq.gz (unmerged R2)           │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ├─────────────────────┬──────────────────┐
                     ▼                     ▼                  ▼
         ┌──────────────────┐  ┌──────────────────┐  ┌──────────────────┐
         │  MERGED DETECTION│  │ UNMERGED R1      │  │ UNMERGED R2      │
         │  (detect_fusions_│  │ (detect_fusions_ │  │ (detect_fusions_ │
         │   _string)       │  │  _unmerged_string│  │  _unmerged_string│
         └────────┬─────────┘  └────────┬─────────┘  └────────┬─────────┘
                  │                     │                     │
                  ▼                     ▼                     ▼
         {sample}.fusion_      {sample}.R1.unmerged_  {sample}.R2.unmerged_
         counts.csv            fusion_counts.csv       fusion_counts.csv
                  │                     │                     │
                  └─────────────────────┴─────────────────────┘
                                       │
                                       ▼
                           ┌───────────────────────────┐
                           │  AGGREGATION              │
                           │  • aggregate_counts        │
                           │    (merged only)           │
                           │  • aggregate_unmerged_     │
                           │    counts (R1 + R2)        │
                           └───────────────────────────┘
```

**Key Characteristics:**

1. **Separate Processing:**
   - R1 and R2 unmerged reads are processed independently
   - Each mate uses the same string matching algorithm as merged reads
   - Same breakpoint reference sequences and detection parameters

2. **Output File Naming:**
   - Merged counts: `{sample}.fusion_counts.csv`
   - Unmerged R1 counts: `{sample}.R1.unmerged_fusion_counts.csv`
   - Unmerged R2 counts: `{sample}.R2.unmerged_fusion_counts.csv`
   - Unmerged metrics: `{sample}.{mate}.unmerged_fusion_metrics.json`
   - Unmerged partner counts: `{sample}.{mate}.unmerged_partner_counts.csv`

3. **Count Separation:**
   - Unmerged counts are kept completely separate from merged counts
   - Aggregation creates distinct summary files:
     - `fusion_counts_summary.csv` (merged only)
     - `unmerged_counts_summary.csv` (R1 + R2 combined)
   - Same fusion detected in both merged and unmerged reads results in separate entries

4. **Aggregation:**
   - `aggregate_unmerged_counts` rule combines R1 and R2 counts per sample
   - Sample columns in unmerged summary use format: `{sample}.R1` and `{sample}.R2`
   - Produces separate QC metrics: `unmerged_qc_metrics.csv`

5. **Rule Execution:**
   - Unmerged detection rules are always defined but only execute when:
     - `detection.unmerged_detection: true` in config
     - Targets are requested (via `get_all_targets()` in `common.smk`)
   - Snakemake automatically handles conditional execution based on target dependencies

6. **Empty File Handling:**
   - When bbmerge successfully merges all reads, unmerged files may be empty
   - The string matcher script handles empty FASTQ files gracefully
   - Empty files produce empty count files (no errors)

**Use Cases:**
- When a significant fraction of reads fail to merge
- To compare detection sensitivity between merged and unmerged reads
- To capture fusions that may be present in only one mate
- For quality control and sensitivity analysis

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
- Breakpoint: 126 nucleotides from anchor start (when `truncated_component: 'anchor'`)
- Anchor: Met_WT

**Note:** The breakpoint position is relative to whichever component is truncated, as specified by `fusion_library.anchor.truncated_component` in the configuration. When `truncated_component: 'anchor'` (the default in the actual configuration), the breakpoint represents nucleotides from the anchor start, indicating how much of the anchor's N-terminal region is truncated.

### Breakpoint Sequences
- K-mer spanning breakpoint junction
- Window size configurable (default: 12 nt each side)
- Includes linker sequence if present

### Count Files
- CSV format: `fusion_id,count` (or `fusion_id,type,count` if type column present)
- Merged counts: `{sample}.fusion_counts.csv` (one per sample)
- Unmerged counts: `{sample}.{mate}.unmerged_fusion_counts.csv` (one per sample per mate when unmerged detection enabled)
- Aggregated into summary matrices:
  - `fusion_counts_summary.csv` (merged counts across all samples)
  - `unmerged_counts_summary.csv` (unmerged counts across all samples and mates)

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
