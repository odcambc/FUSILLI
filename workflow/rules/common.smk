"""
Common functions and configuration loading for FUSILLI pipeline.

This module handles:
- Configuration parsing and validation
- Sample metadata loading
- Helper functions for rule inputs/outputs
"""

import warnings
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from snakemake.utils import validate


# =============================================================================
# CONFIGURATION VALIDATION
# =============================================================================

# Validate main config against schema
validate(config, "../schemas/config.schema.yaml")


# =============================================================================
# CONFIGURATION EXTRACTION
# =============================================================================

# Experiment identification
EXPERIMENT = config["experiment"]

# Paths
DATA_DIR = config["data_dir"]
REF_DIR = config.get("ref_dir", "references")

# Fusion library configuration
FUSION_CONFIG = config["fusion_library"]
ANCHOR_NAME = FUSION_CONFIG["anchor"]["name"]
ANCHOR_POSITION = FUSION_CONFIG["anchor"].get("position", "downstream")
TRUNCATED_COMPONENT = FUSION_CONFIG["anchor"].get("truncated_component", "partner")
LINKER_SEQUENCE = FUSION_CONFIG.get("linker_sequence", "")
PARTNERS_FILE = FUSION_CONFIG["partners_file"]
SEQUENCES_FILE = FUSION_CONFIG["sequences_file"]
UNFUSED_SEQUENCES_FILE = FUSION_CONFIG.get("unfused_sequences_file", None)
EXON_PARTNERS_FILE = FUSION_CONFIG.get("exon_partners_file", None)
VARIANT_ANCHORS = FUSION_CONFIG.get("variant_anchors", [])

# Detection parameters
DETECTION_CONFIG = config.get("detection", {})
DETECTION_METHOD = DETECTION_CONFIG.get("method", "string")
BREAKPOINT_WINDOW = DETECTION_CONFIG.get("breakpoint_window", 12)
MAINTAIN_FRAME = DETECTION_CONFIG.get("maintain_frame", True)
KMER_SIZE = DETECTION_CONFIG.get("kmer_size", 15)
ORIENTATION_CHECK = DETECTION_CONFIG.get("orientation_check", False)
LINKER_FIRST = DETECTION_CONFIG.get("linker_first", False)
PREFILTER_FALLBACK = DETECTION_CONFIG.get("prefilter_fallback", False)
UNMERGED_DETECTION = DETECTION_CONFIG.get("unmerged_detection", False)

# Sequencing parameters
SEQ_CONFIG = config.get("sequencing", {})
PAIRED = SEQ_CONFIG.get("paired", True)
MIN_QUALITY = SEQ_CONFIG.get("min_quality", 30)

# Preprocessing
PREPROC_CONFIG = config.get("preprocessing", {})
ADAPTERS = PREPROC_CONFIG.get("adapters", "resources/adapters.fa")
CONTAMINANTS = PREPROC_CONFIG.get("contaminants", [])

# QC options
QC_CONFIG = config.get("qc", {})
RUN_QC = QC_CONFIG.get("run_qc", False)
BASELINE_CONDITION = QC_CONFIG.get("baseline_condition", "baseline")

# Pipeline behavior
PIPELINE_CONFIG = config.get("pipeline", {})
SHOW_PROGRESS = PIPELINE_CONFIG.get("show_progress", True)
PROGRESS_INTERVAL = PIPELINE_CONFIG.get("progress_interval", 1)

# Quick mode
QUICK_CONFIG = config.get("quick", {})
QUICK_ENABLED = QUICK_CONFIG.get("enabled", False)
QUICK_MAX_READS = QUICK_CONFIG.get("max_reads", 100000)
QUICK_FRACTION = QUICK_CONFIG.get("fraction", None)
QUICK_SEED = QUICK_CONFIG.get("seed", 1337)

# Resources
RESOURCES_CONFIG = config.get("resources", {})
DEFAULT_MEMORY = RESOURCES_CONFIG.get("memory_mb", 16000)
DEFAULT_THREADS = RESOURCES_CONFIG.get("threads", 16)


# =============================================================================
# LOAD SAMPLES
# =============================================================================

def load_samples_csv(filepath: str) -> pd.DataFrame:
    """Load and validate samples CSV file."""
    # Read CSV, skipping comment lines
    with open(filepath, 'r') as f:
        lines = [l for l in f if not l.strip().startswith('#')]

    from io import StringIO
    df = pd.read_csv(StringIO(''.join(lines)))

    # Validate required columns
    required = {'sample', 'condition', 'file'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in samples file: {missing}")

    # Set sample as index
    df = df.set_index('sample', drop=False, verify_integrity=True)

    return df


# Load samples
SAMPLES_DF = load_samples_csv(config["samples_file"])
SAMPLES = list(SAMPLES_DF["sample"])
FILES = list(SAMPLES_DF["file"])


# =============================================================================
# LOAD FUSION PARTNERS
# =============================================================================

def load_partners_csv(filepath: str) -> pd.DataFrame:
    """Load and validate fusion partners CSV file."""
    with open(filepath, 'r') as f:
        lines = [l for l in f if not l.strip().startswith('#')]

    from io import StringIO
    df = pd.read_csv(StringIO(''.join(lines)))

    # Filter to included partners only
    df = df[df['include'].astype(str).str.lower().isin(['true', 'yes', '1'])]

    return df.set_index('partner_name', drop=False)


PARTNERS_DF = load_partners_csv(PARTNERS_FILE)
PARTNERS = list(PARTNERS_DF["partner_name"])


# =============================================================================
# SAMPLE HELPER FUNCTIONS
# =============================================================================

def get_sample_info(sample_name: str) -> dict:
    """Get all info for a sample."""
    if sample_name not in SAMPLES_DF.index:
        raise ValueError(f"Unknown sample: {sample_name}")
    return SAMPLES_DF.loc[sample_name].to_dict()


def get_samples_by_condition(condition: str) -> list:
    """Get all samples with a given condition."""
    return list(SAMPLES_DF[SAMPLES_DF["condition"] == condition]["sample"])


def get_baseline_samples() -> list:
    """Get samples with baseline condition."""
    return get_samples_by_condition(BASELINE_CONDITION)


def get_experiment_samples() -> list:
    """Get non-baseline samples."""
    return list(SAMPLES_DF[SAMPLES_DF["condition"] != BASELINE_CONDITION]["sample"])


# =============================================================================
# INPUT FILE FUNCTIONS
# =============================================================================

def get_raw_fastq_r1(wildcards):
    """Get R1 FASTQ file for a sample."""
    sample = wildcards.sample
    file_prefix = SAMPLES_DF.loc[sample, "file"]
    return str(Path(DATA_DIR) / f"{file_prefix}_R1_001.fastq.gz")


def get_raw_fastq_r2(wildcards):
    """Get R2 FASTQ file for a sample."""
    sample = wildcards.sample
    file_prefix = SAMPLES_DF.loc[sample, "file"]
    return str(Path(DATA_DIR) / f"{file_prefix}_R2_001.fastq.gz")


def get_raw_fastqs(wildcards):
    """Get both FASTQ files for a sample."""
    return {
        "R1": get_raw_fastq_r1(wildcards),
        "R2": get_raw_fastq_r2(wildcards)
    }


def get_quick_fastq_r1(wildcards):
    """Get subsampled R1 FASTQ for a sample (quick mode)."""
    sample = wildcards.sample
    return f"results/{EXPERIMENT}/quick/{sample}_R1.sub.fastq.gz"


def get_quick_fastq_r2(wildcards):
    """Get subsampled R2 FASTQ for a sample (quick mode)."""
    sample = wildcards.sample
    return f"results/{EXPERIMENT}/quick/{sample}_R2.sub.fastq.gz"


def get_input_fastq_r1(wildcards):
    """Return appropriate R1 depending on quick mode."""
    if QUICK_ENABLED:
        return get_quick_fastq_r1(wildcards)
    return get_raw_fastq_r1(wildcards)


def get_input_fastq_r2(wildcards):
    """Return appropriate R2 depending on quick mode."""
    if QUICK_ENABLED:
        return get_quick_fastq_r2(wildcards)
    return get_raw_fastq_r2(wildcards)


def get_input_fastqs(wildcards):
    """Return both FASTQs depending on quick mode."""
    return {
        "R1": get_input_fastq_r1(wildcards),
        "R2": get_input_fastq_r2(wildcards)
    }


# =============================================================================
# REFERENCE FILE FUNCTIONS
# =============================================================================

def get_reference_fasta():
    """Get full path to reference FASTA."""
    return f"{REF_DIR}/{SEQUENCES_FILE}"


def get_reference_names() -> list:
    """Get list of sequence names from reference FASTA."""
    fasta_path = get_reference_fasta()
    names = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        names.append(record.id)
    return names


# =============================================================================
# HELPER FUNCTIONS FOR RULES
# =============================================================================

def pass_names(names) -> str:
    """Convert list of paths to comma-separated string for bbduk."""
    if isinstance(names, str):
        return names
    return ",".join(names)


# Computed values for rules
ADAPTERS_REF = pass_names(ADAPTERS)
CONTAMINANTS_REF = pass_names(CONTAMINANTS) if CONTAMINANTS else ""


# =============================================================================
# TARGET FILE FUNCTIONS
# =============================================================================

def get_all_targets(wildcards):
    """Generate all target files for rule all."""
    targets = []

    # Fusion counts for each sample
    targets.extend(expand(
        "results/{experiment}/counts/{sample}.fusion_counts.csv",
        experiment=EXPERIMENT,
        sample=SAMPLES
    ))

    if UNMERGED_DETECTION:
        targets.extend(expand(
            "results/{experiment}/counts/{sample}.{mate}.unmerged_fusion_counts.csv",
            experiment=EXPERIMENT,
            sample=SAMPLES,
            mate=["R1", "R2"]
        ))
        targets.append(f"results/{EXPERIMENT}/unmerged_counts_summary.csv")
        targets.append(f"results/{EXPERIMENT}/unmerged_qc_metrics.csv")
        targets.append(f"results/{EXPERIMENT}/unmerged_partner_counts_summary.csv")

    # Aggregated summary file with sample columns
    targets.append(f"results/{EXPERIMENT}/fusion_counts_summary.csv")
    targets.append(f"results/{EXPERIMENT}/sensitivity_metrics.csv")
    targets.append(f"results/{EXPERIMENT}/decay_metrics.csv")

    # QC reports if enabled
    if RUN_QC:
        targets.extend(expand(
            "stats/{experiment}/{experiment}_multiqc.html",
            experiment=EXPERIMENT
        ))

    # Reproducibility metadata
    targets.extend([
        f"results/{EXPERIMENT}/repro/metadata.json",
        f"results/{EXPERIMENT}/repro/metadata.txt",
        f"results/{EXPERIMENT}/repro/conda-env.yaml",
        f"results/{EXPERIMENT}/repro/pip-freeze.txt"
    ])

    return targets


# =============================================================================
# DEPRECATED COMPATIBILITY LAYER
# =============================================================================
# These provide backward compatibility with old config format
# Will be removed in future versions

# Legacy variable names (deprecated)
experiment = EXPERIMENT
samples = SAMPLES
files = FILES
paired = PAIRED

# Legacy function names (deprecated)
def get_input_strings(wildcards):
    """DEPRECATED: Use get_all_targets instead."""
    warnings.warn(
        "get_input_strings is deprecated, use get_all_targets",
        DeprecationWarning
    )
    return get_all_targets(wildcards)
