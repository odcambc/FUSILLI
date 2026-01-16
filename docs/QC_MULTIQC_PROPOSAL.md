# QC Metrics & MultiQC Integration Proposal

## Executive Summary

This proposal outlines an enhanced QC architecture for FUSILLI that leverages MultiQC to aggregate outputs from standard bioinformatics tools (FastQC, BBDuk, BBMerge) and integrates FUSILLI-specific fusion detection metrics. The architecture provides comprehensive quality assessment across read quality, preprocessing efficiency, library representation, and detection performance.

## Current State

### Existing QC Infrastructure

The pipeline currently generates:

1. **FastQC outputs**: Per-sample, per-read (R1/R2) HTML and ZIP reports
2. **BBDuk statistics**: Trim, contaminant removal, and quality filtering stats/logs
3. **BBMerge statistics**: Insert size histograms (ihist) and merge logs
4. **FUSILLI-specific metrics**: Fusion counts, detection metrics, sensitivity metrics

### Current MultiQC Integration

MultiQC is already configured to aggregate:
- FastQC ZIP files and HTML reports
- BBMerge ihist files (insert size histograms)
- BBDuk stats files (trim, contam, quality)
- BBDuk log files (trim, clean, quality)
- BBMerge log files
- FUSILLI-specific CSV and JSON metrics files

However, the current implementation has limitations:
- BBDuk histogram files (qhist, bhist, gchist, aqhist, lhist) are not included
- MultiQC configuration is minimal (basic wrapper usage)
- FUSILLI-specific metrics are passed as CSV/JSON but not visualized
- No custom MultiQC modules for fusion-specific visualizations

## Proposed QC Metrics Architecture

### 1. Standard Tool Metrics (MultiQC Native Support)

#### 1.1 FastQC Metrics
**Status**: Already aggregated, but can be enhanced

**Available Metrics**:
- **Basic Statistics**: Total sequences, sequence length, %GC
- **Per-base sequence quality**: Mean quality scores across read positions
- **Per-sequence quality scores**: Quality distribution per read
- **Per-base GC content**: GC% across read positions
- **Per-sequence GC content**: GC% distribution per read
- **Per-base N content**: N content across positions
- **Sequence length distribution**: Read length histogram
- **Sequence duplication levels**: Duplication rate and distribution
- **Overrepresented sequences**: Top overrepresented sequences
- **Adapter content**: Adapter contamination by position

**Enhancement Opportunities**:
- Extract and tabulate key statistics for cross-sample comparison
- Flag samples with quality issues (e.g., mean quality < Q20)
- Compare R1 vs R2 quality metrics

**Implementation**:
- **No additional processing required**: FastQC already generates all metrics in `fastqc_data.txt` within the ZIP files
- **MultiQC already parses**: MultiQC natively extracts all FastQC metrics from ZIP files
- **Configuration only**: Enhancements are achieved through MultiQC configuration:
  - Configure which metrics appear in summary tables (`table_columns_visible` in `multiqc_config.yaml`)
  - Set up custom table columns for specific statistics (mean quality, GC content, duplication rate)
  - Configure plot types and comparisons
  - Set quality thresholds for flagging samples (e.g., highlight samples with mean quality < Q20)
  - Enable R1 vs R2 comparison plots (MultiQC can automatically group by read direction from filename patterns)
- **No custom parsing needed**: All data is already extracted by MultiQC's FastQC module; we just configure what to display

#### 1.2 BBDuk Metrics
**Status**: Partially aggregated (stats and logs), histograms not included

**Available Outputs**:
- **Stats files** (`*.stats.txt`): Input/output read and base counts
- **Log files** (`*.log`): Detailed processing information
- **Histogram files** (currently not aggregated):
  - `*.qhist`: Quality score histogram
  - `*.bhist`: Base composition histogram
  - `*.gchist`: GC content histogram
  - `*.aqhist`: Average quality histogram
  - `*.lhist`: Read length histogram

**Metrics to Extract**:
- **Trim step**:
  - Read retention rate: `output_reads / input_reads`
  - Base retention rate: `output_bases / input_bases`
  - Mean read length before/after trimming
  - Adapter removal efficiency
- **Contaminant removal step**:
  - Contaminant read fraction: `(input_reads - output_reads) / input_reads`
  - PhiX detection rate (if tracked separately)
- **Quality filtering step**:
  - Quality filter retention: `output_reads / input_reads`
  - Quality threshold used (from config)
  - Mean quality before/after filtering

**Implementation**:
- MultiQC has native support for BBDuk stats files
- Histogram files can be parsed and visualized with custom parsing or existing modules
- Log files can be parsed for additional metrics

#### 1.3 BBMerge Metrics
**Status**: Partially aggregated (ihist and logs)

**Available Outputs**:
- **Log files** (`*.log`): Merge statistics (Pairs, Joined, Avg Insert)
- **Insert size histogram** (`*.ihist`): Distribution of insert sizes

**Metrics to Extract**:
- **Merge rate**: `merged_reads / total_pairs`
- **Unmerged fraction**: `unmerged_pairs / total_pairs`
- **Average insert size**: From log or calculated from ihist
- **Insert size distribution**: Median, mode, standard deviation from ihist
- **Overlap statistics**: Median/mean overlap length (if available)

**Implementation**:
- MultiQC has native support for BBMerge logs and ihist files
- Can extract additional statistics from ihist parsing

### 2. FUSILLI-Specific QC Metrics

#### 2.1 Library Representation Metrics
**Purpose**: Assess how well the fusion library is represented in sequencing data

**Key Metrics**:
- **Variant Coverage**:
  - `expected_variants`: Total variants in library design
  - `observed_variants`: Variants with count > 0
  - `coverage_fraction`: `observed_variants / expected_variants`
  - `missing_variants`: `expected_variants - observed_variants`
- **Breakpoint Coverage**:
  - `breakpoints_expected`: Number of breakpoint positions in design
  - `breakpoints_detected`: Breakpoint positions with at least one detection
  - `breakpoint_coverage_fraction`: `breakpoints_detected / breakpoints_expected`
- **Partner Coverage**:
  - `partners_expected`: Number of partner domains in design
  - `partners_detected`: Partners with at least one fusion detected
  - `partner_coverage_fraction`: `partners_detected / partners_expected`

**Data Sources**:
- `fusion_counts_summary.csv`: Contains all expected variants
- `partner_counts_summary.csv`: Contains partner detection counts
- Configuration: Expected variants from fusion library definition

#### 2.2 Library Diversity Metrics
**Purpose**: Quantify library complexity and evenness

**Key Metrics**:
- **Shannon Diversity Index**: `-Σ(p_i * ln(p_i))` where p_i is proportion of variant i
- **Simpson Diversity Index**: `1 - Σ(p_i²)`
- **Evenness**: `shannon_diversity / ln(observed_variants)`
- **Gini Coefficient**: Measure of inequality in variant distribution
- **Top N Fractions**: Fraction of counts in top 1, 5, 10, 25 variants
- **Median/Mean Variant Count**: Central tendency of variant counts

**Data Sources**:
- `fusion_counts_summary.csv`: Per-sample variant counts
- Calculated from count distributions

#### 2.3 Detection Performance Metrics
**Purpose**: Measure how effectively reads are being used for detection

**Key Metrics**:
- **Read Utilization**:
  - `reads_processed`: Total reads analyzed
  - `reads_with_partner_end`: Reads passing pre-filter
  - `reads_matched`: Reads with detected breakpoint
  - `prefilter_efficiency`: `reads_with_partner_end / reads_processed`
  - `detection_efficiency`: `reads_matched / reads_processed`
  - `matching_efficiency`: `reads_matched / reads_with_partner_end`
- **Detection Yield**:
  - `detections_per_read`: `total_counts / reads_processed`
  - `detections_per_merged_read`: `total_counts / merged_reads`
  - `detections_per_million_reads`: `(total_counts / reads_processed) * 1e6`
- **Sensitivity Metrics**:
  - `read_length_mean`: Mean read length (from FastQC)
  - `breakpoint_kmer_length`: `2 * breakpoint_window`
  - `fraction_reads_long_enough`: Fraction of reads ≥ breakpoint_kmer_length
  - `expected_detectable_reads`: `reads_processed * fraction_reads_long_enough`
  - `sensitivity_index`: `reads_matched / expected_detectable_reads`

**Data Sources**:
- `fusion_qc_metrics.csv`: Contains detection metrics
- `sensitivity_metrics.csv`: Contains sensitivity calculations
- `*.fusion_metrics.json`: Per-sample detection statistics
- FastQC: Read length statistics

#### 2.4 Preprocessing Efficiency Metrics
**Purpose**: Assess cumulative read loss through preprocessing pipeline

**Key Metrics**:
- **End-to-End Retention**:
  - `raw_to_merged_read_retention`: `merged_reads / raw_reads`
  - `raw_to_merged_base_retention`: `merged_bases / raw_bases`
  - `total_read_loss_fraction`: `(raw_reads - merged_reads) / raw_reads`
- **Step-by-Step Loss Breakdown**:
  - `loss_at_trim`: Reads lost during trimming
  - `loss_at_contam`: Reads lost during contaminant removal
  - `loss_at_quality`: Reads lost during quality filtering
  - `loss_at_merge`: Reads lost during merging (unmerged pairs)
  - Fraction lost at each step

**Data Sources**:
- `decay_metrics.csv`: Already tracks read/base counts through steps
- BBDuk stats/logs: Step-specific retention
- BBMerge logs: Merge statistics

#### 2.5 Unmerged Read Metrics (when enabled)
**Purpose**: Compare detection performance between merged and unmerged reads

**Key Metrics**:
- `merged_detection_rate`: `merged_detections / merged_reads`
- `unmerged_detection_rate`: `unmerged_detections / unmerged_reads`
- `merged_vs_unmerged_ratio`: `merged_detection_rate / unmerged_detection_rate`
- `total_detection_rate`: `(merged + unmerged detections) / total_reads`

**Data Sources**:
- `unmerged_qc_metrics.csv`: Unmerged-specific metrics
- `unmerged_counts_summary.csv`: Unmerged fusion counts
- BBMerge logs: Unmerged read counts

## MultiQC Configuration & Integration

### 3.1 Standard Tool Integration

#### FastQC Module
**Status**: Native MultiQC support

**Configuration**:
- Already included via FastQC ZIP files
- Can add custom sections for key statistics extraction
- Enable cross-sample comparison tables

**Enhancements**:
- Extract and tabulate critical metrics (mean quality, GC content, duplication rate)
- Flag samples with quality issues
- Compare R1 vs R2 metrics

#### BBDuk Module
**Status**: Native MultiQC support for stats files

**Current Inputs**:
- `stats/{experiment}/trim/{sample}.trim.stats.txt`
- `stats/{experiment}/contam/{sample}.contam.stats.txt`
- `stats/{experiment}/quality/{sample}.quality.stats.txt`
- Log files: `logs/{experiment}/bbduk/{sample}.*.log`

**Proposed Additions**:
- Include histogram files for visualization:
  - `stats/{experiment}/trim/{sample}.qhist` (quality histogram)
  - `stats/{experiment}/trim/{sample}.bhist` (base composition)
  - `stats/{experiment}/trim/{sample}.gchist` (GC content)
  - `stats/{experiment}/trim/{sample}.aqhist` (average quality)
  - `stats/{experiment}/trim/{sample}.lhist` (length histogram)

**Implementation**:
- MultiQC has a generic histogram parser that can handle these files
- May need custom parsing function for optimal visualization
- Create summary tables for retention rates across steps

#### BBMerge Module
**Status**: Native MultiQC support

**Current Inputs**:
- `stats/{experiment}/merge/{sample}.ihist` (insert size histogram)
- `logs/{experiment}/bbmerge/{sample}.log` (merge statistics)

**Enhancements**:
- Extract additional statistics from ihist (median, mode, std dev)
- Create summary table of merge rates across samples
- Visualize insert size distributions

### 3.2 Custom MultiQC Modules for FUSILLI Metrics

#### Module 1: Fusion Detection Metrics
**Purpose**: Visualize fusion detection performance and library representation

**Input Files**:
- `results/{experiment}/fusion_qc_metrics.csv`
- `results/{experiment}/sensitivity_metrics.csv`
- `results/{experiment}/fusion_counts_summary.csv`

**Visualizations**:
1. **Detection Efficiency Plot**:
   - X-axis: Sample
   - Y-axis: Detection efficiency, prefilter efficiency, matching efficiency
   - Multi-line plot comparing different efficiency metrics

2. **Library Coverage Plot**:
   - X-axis: Sample
   - Y-axis: Coverage fraction (variants, breakpoints, partners)
   - Stacked bar chart or multi-line plot

3. **Sensitivity Analysis**:
   - X-axis: Sample
   - Y-axis: Sensitivity index, expected detection fraction
   - Scatter plot or bar chart

4. **Detection Yield Table**:
   - Sample | Reads Processed | Detections | Detections/Read | Detections/Million

#### Module 2: Library Diversity Metrics
**Purpose**: Visualize library complexity and evenness

**Input Files**:
- `results/{experiment}/fusion_counts_summary.csv`

**Calculated Metrics**:
- Shannon diversity index
- Simpson diversity index
- Evenness
- Top N fractions

**Visualizations**:
1. **Diversity Comparison**:
   - X-axis: Sample
   - Y-axis: Shannon/Simpson diversity
   - Bar chart comparing diversity across samples

2. **Evenness Plot**:
   - X-axis: Sample
   - Y-axis: Evenness score
   - Bar chart

3. **Top N Fractions**:
   - X-axis: Sample
   - Y-axis: Fraction of counts
   - Stacked bar chart (Top 1, Top 5, Top 10, Top 25, Other)

4. **Variant Count Distribution**:
   - Histogram of variant counts (log scale)
   - Shows distribution of variant abundances

#### Module 3: Preprocessing Pipeline Metrics
**Purpose**: Visualize read retention through preprocessing steps

**Input Files**:
- `results/{experiment}/decay_metrics.csv`
- `results/{experiment}/trim_metrics.csv`
- `results/{experiment}/contam_metrics.csv`
- `results/{experiment}/quality_metrics.csv`

**Visualizations**:
1. **Read Decay Plot**:
   - X-axis: Processing step (Raw, Trimmed, Cleaned, Quality, Merged)
   - Y-axis: Read count (log scale)
   - Line plot showing read loss at each step

2. **Retention Rate Comparison**:
   - X-axis: Sample
   - Y-axis: Retention rate (%)
   - Grouped bar chart for each preprocessing step

3. **Step Loss Breakdown**:
   - Stacked bar chart showing fraction lost at each step
   - X-axis: Sample
   - Y-axis: Fraction of reads lost
   - Stacked by step (trim, contam, quality, merge)

#### Module 4: Partner Detection Metrics
**Purpose**: Visualize partner domain detection

**Input Files**:
- `results/{experiment}/partner_counts_summary.csv`
- `results/{experiment}/counts/{sample}.partner_counts.csv`

**Visualizations**:
1. **Partner Detection Heatmap**:
   - Rows: Partners
   - Columns: Samples
   - Values: Detection counts (log scale)
   - Color-coded heatmap

2. **Partner Coverage**:
   - X-axis: Partner
   - Y-axis: Number of samples with detection
   - Bar chart

3. **Partner End vs Linker Detection**:
   - Scatter plot comparing partner end counts vs linker counts
   - Per partner, per sample

### 3.3 MultiQC Configuration File

Create `config/multiqc_config.yaml` with:

```yaml
# MultiQC Configuration for FUSILLI

# Report title
report_title: "FUSILLI QC Report - {experiment}"

# Custom module order
module_order:
  - fastqc
  - bbduk
  - bbmerge
  - fusilli_detection
  - fusilli_diversity
  - fusilli_preprocessing
  - fusilli_partners

# FastQC settings
fastqc:
  # Extract key statistics
  extract_stats: true

# BBDuk settings
bbduk:
  # Include histogram files
  include_histograms: true
  # Histogram file patterns
  qhist_pattern: "*.qhist"
  bhist_pattern: "*.bhist"
  gchist_pattern: "*.gchist"
  aqhist_pattern: "*.aqhist"
  lhist_pattern: "*.lhist"

# BBMerge settings
bbmerge:
  # Extract additional statistics from ihist
  extract_ihist_stats: true

# Custom FUSILLI modules
# Note: FUSILLI custom modules are now packaged separately as 'fusilli-multiqc'
# and registered via setuptools entry points. They are automatically discovered
# by MultiQC when the package is installed. No custom_modules configuration needed.
# The package is installed via workflow/envs/qc.yaml from Git repository, or manually with:
#   pip install git+https://github.com/user/fusilli-multiqc.git@main
# Once published to PyPI, install with:
#   pip install fusilli-multiqc>=1.0.0

# Table settings
table_columns_visible:
  # Show key metrics in summary table
  - Sample
  - Total Reads
  - Merged Reads
  - Detection Efficiency
  - Coverage Fraction
  - Shannon Diversity

# Plot settings
plots:
  # Default plot settings
  default_plot_type: "bar"
  show_points: true
```

## Implementation Plan

### Phase 1: Enhanced Standard Tool Integration
**Priority**: High
**Effort**: Low-Medium

1. **Update MultiQC inputs** (`workflow/rules/qc.smk`):
   - Add BBDuk histogram files to `get_multiqc_inputs()`
   - Ensure all stats and log files are included

2. **Create MultiQC config** (`config/multiqc_config.yaml`):
   - Configure standard modules (FastQC, BBDuk, BBMerge)
   - Set up histogram file patterns
   - Configure table columns

3. **Test MultiQC aggregation**:
   - Verify all standard tool outputs are aggregated
   - Check histogram visualizations
   - Validate summary tables

### Phase 2: Custom FUSILLI Modules
**Priority**: High
**Effort**: Medium-High

1. **Create custom module structure**:
   - `fusilli-multiqc/` package in separate Git repository
   - `fusilli_multiqc/modules/` subdirectory for modules
   - `fusilli_multiqc/utils.py` for shared utilities (replaces `fusilli_base.py`)
   - Package setup files (`setup.py`, `pyproject.toml`)
   - Package installed from Git repository or PyPI (once published)

2. **Implement detection metrics module**:
   - Parse `fusion_qc_metrics.csv` and `sensitivity_metrics.csv`
   - Create detection efficiency plots
   - Create library coverage plots
   - Create sensitivity analysis plots

3. **Implement diversity metrics module**:
   - Parse `fusion_counts_summary.csv`
   - Calculate diversity indices
   - Create diversity comparison plots
   - Create evenness plots

4. **Implement preprocessing module**:
   - Parse `decay_metrics.csv` and step-specific metrics
   - Create read decay plots
   - Create retention rate comparisons

5. **Implement partner detection module**:
   - Parse `partner_counts_summary.csv`
   - Create partner detection heatmap
   - Create partner coverage plots

### Phase 3: Metric Calculation Enhancements
**Priority**: Medium
**Effort**: Medium

1. **Enhance metric calculation** (`workflow/rules/process_strings.smk`):
   - Add diversity index calculations
   - Add breakpoint coverage calculations
   - Add partner coverage calculations
   - Add detection yield metrics

2. **Update CSV outputs**:
   - Add new metrics to `fusion_qc_metrics.csv`
   - Create new summary files if needed
   - Ensure all metrics are MultiQC-compatible

### Phase 4: Documentation & Validation
**Priority**: Medium
**Effort**: Low

1. **Documentation**:
   - Update README with QC metrics explanation
   - Document MultiQC configuration
   - Create guide for interpreting QC reports

2. **Validation**:
   - Test with multiple experiments
   - Validate metric calculations
   - Check visualization quality

## File Structure

```
fusilli-multiqc/                  # Standalone pip-installable package (separate repository)
├── fusilli_multiqc/
│   ├── __init__.py
│   ├── utils.py                  # Shared utilities (replaces fusilli_base.py)
│   └── modules/
│       ├── __init__.py
│       ├── detection.py          # Detection metrics module
│       ├── diversity.py           # Diversity metrics module
│       ├── preprocessing.py       # Preprocessing module
│       └── partners.py           # Partner detection module
├── setup.py                      # Setuptools configuration
├── pyproject.toml                 # Modern Python packaging
└── README.md                     # Package documentation

workflow/
├── rules/
│   └── qc.smk                    # Updated with histogram inputs
├── envs/
│   └── qc.yaml                   # Conda environment (installs fusilli-multiqc from Git/PyPI)
config/
└── multiqc_config.yaml           # MultiQC configuration
```

## Expected Outputs

### MultiQC Report Sections

1. **FastQC**: Standard FastQC visualizations and tables
2. **BBDuk**:
   - Summary table of read/base retention across steps
   - Histogram visualizations (quality, GC, length)
   - Step-by-step comparison plots
3. **BBMerge**:
   - Merge rate summary table
   - Insert size distribution plots
   - Overlap statistics
4. **FUSILLI Detection Metrics**:
   - Detection efficiency plots
   - Library coverage plots
   - Sensitivity analysis
   - Detection yield table
5. **FUSILLI Diversity Metrics**:
   - Diversity index comparison
   - Evenness plots
   - Top N fractions
   - Variant count distribution
6. **FUSILLI Preprocessing**:
   - Read decay plot
   - Retention rate comparison
   - Step loss breakdown
7. **FUSILLI Partner Detection**:
   - Partner detection heatmap
   - Partner coverage plot
   - Partner end vs linker comparison

## Quality Thresholds & Flags

Define thresholds for key metrics to flag potential issues:

| Metric | Expected | Warning | Error |
|--------|----------|---------|-------|
| Coverage Fraction | > 0.8 | 0.6-0.8 | < 0.6 |
| Detection Efficiency | > 0.1 | 0.05-0.1 | < 0.05 |
| Merge Rate | > 0.7 | 0.5-0.7 | < 0.5 |
| Mean Quality (Q) | > 30 | 20-30 | < 20 |
| Shannon Diversity | > 5 | 3-5 | < 3 |
| Read Retention | > 0.5 | 0.3-0.5 | < 0.3 |

MultiQC can be configured to highlight samples that exceed warning/error thresholds.

## Benefits

1. **Comprehensive QC**: Single report aggregating all QC metrics
2. **Standard Tool Integration**: Leverages MultiQC's native support for FastQC, BBDuk, BBMerge
3. **Custom Visualizations**: FUSILLI-specific metrics with tailored visualizations
4. **Cross-Sample Comparison**: Easy comparison of metrics across samples
5. **Reproducibility**: Standardized QC reporting across experiments
6. **Automation**: Integrated into Snakemake workflow
7. **Extensibility**: Easy to add new metrics or visualizations

## Next Steps

1. Review and approve this proposal
2. Implement Phase 1 (enhanced standard tool integration)
3. Develop custom MultiQC modules (Phase 2)
4. Test with existing data
5. Iterate based on feedback
6. Document and deploy
