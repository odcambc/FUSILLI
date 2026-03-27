"""
String-matching fusion detection rules.

This module contains rules for:
1. Generating breakpoint reference sequences
2. Detecting fusions via string matching
"""


rule generate_breakpoint_sequences:
    """
    Generate all possible breakpoint sequences for the fusion library.

    This rule creates:
    - breakpoint_sequences.csv: K-mers spanning each possible breakpoint
    - domain_ends.csv: 3' end k-mers for pre-filtering
    - unfused_sequences.csv: K-mers for unfused control sequences (if configured)

    These files are used by the string matching detection step.
    """
    input:
        sequences=get_reference_fasta()
    output:
        breakpoints="results/{experiment}/references/breakpoint_sequences.csv",
        ends="results/{experiment}/references/domain_ends.csv",
        unfused="results/{experiment}/references/unfused_sequences.csv",
        variants="results/{experiment}/references/variant_catalog.csv",
        expected_counts="results/{experiment}/references/expected_counts_template.csv"
    params:
        partners_file=PARTNERS_FILE,
        exon_partners_file=EXON_PARTNERS_FILE,
        anchor_name=ANCHOR_NAME,
        anchor_position=ANCHOR_POSITION,
        truncated_component=TRUNCATED_COMPONENT,
        linker_sequence=LINKER_SEQUENCE,
        breakpoint_window=BREAKPOINT_WINDOW,
        maintain_frame=MAINTAIN_FRAME,
        kmer_size=KMER_SIZE,
        unfused_sequences_file=UNFUSED_SEQUENCES_FILE,
        variant_anchors=VARIANT_ANCHORS
    log:
        "logs/{experiment}/generate_breakpoints.log"
    script:
        "../scripts/fusion_sequences.py"


rule detect_fusions_string:
    """
    Detect fusion breakpoints in merged reads using string matching.

    This is a two-stage algorithm:
    1. Pre-filter: Check if read contains any partner domain 3' end
    2. Match: Search for specific breakpoint sequences
    3. Unfused: Detect unfused control sequences

    Input is the merged (error-corrected) reads from preprocessing.
    """
    input:
        fastq="results/{experiment}/merged/{sample}_merged.fastq.gz",
        breakpoints="results/{experiment}/references/breakpoint_sequences.csv",
        ends="results/{experiment}/references/domain_ends.csv",
        unfused="results/{experiment}/references/unfused_sequences.csv"
    output:
        counts="results/{experiment}/counts/{sample}.fusion_counts.csv",
        metrics="results/{experiment}/counts/{sample}.fusion_metrics.json",
        partner_counts="results/{experiment}/counts/{sample}.partner_counts.csv"
    params:
        show_progress=SHOW_PROGRESS,
        progress_interval=PROGRESS_INTERVAL,
        orientation_check=ORIENTATION_CHECK,
        prefilter_fallback=PREFILTER_FALLBACK,
        linker_sequence=LINKER_SEQUENCE,
        breakpoint_window=BREAKPOINT_WINDOW
    log:
        "logs/{experiment}/string_match/{sample}.log"
    threads: 1  # String matching is I/O bound, not CPU bound
    script:
        "../scripts/string_matcher.py"


rule detect_fusions_unmerged_string:
    """
    Detect fusion breakpoints in unmerged reads using string matching.

    This rule processes R1 and R2 unmerged reads separately (via {mate} wildcard).
    It uses the same two-stage algorithm as merged detection:
    1. Pre-filter: Check if read contains any partner domain 3' end
    2. Match: Search for specific breakpoint sequences
    3. Unfused: Detect unfused control sequences

    Input is the unmerged reads from bbmerge (reads that failed to merge).
    Output files use .unmerged_fusion_counts.csv suffix to keep counts distinct
    from merged read counts. Empty unmerged files are handled gracefully.
    """
    input:
        fastq="results/{experiment}/merged/{sample}_{mate}.unmerged.fastq.gz",
        breakpoints="results/{experiment}/references/breakpoint_sequences.csv",
        ends="results/{experiment}/references/domain_ends.csv",
        unfused="results/{experiment}/references/unfused_sequences.csv"
    output:
        counts="results/{experiment}/counts/{sample}.{mate}.unmerged_fusion_counts.csv",
        metrics="results/{experiment}/counts/{sample}.{mate}.unmerged_fusion_metrics.json",
        partner_counts="results/{experiment}/counts/{sample}.{mate}.unmerged_partner_counts.csv"
    params:
        show_progress=SHOW_PROGRESS,
        progress_interval=PROGRESS_INTERVAL,
        orientation_check=ORIENTATION_CHECK,
        prefilter_fallback=PREFILTER_FALLBACK,
        linker_sequence=LINKER_SEQUENCE,
        breakpoint_window=BREAKPOINT_WINDOW
    log:
        "logs/{experiment}/string_match/{sample}.{mate}.unmerged.log"
    threads: 1  # String matching is I/O bound, not CPU bound
    script:
        "../scripts/string_matcher.py"

rule aggregate_counts:
    """
    Aggregate fusion counts across all samples into a single summary file.
    """
    input:
        counts=expand(
            "results/{{experiment}}/counts/{sample}.fusion_counts.csv",
            sample=SAMPLES
        ),
        metrics=expand(
            "results/{{experiment}}/counts/{sample}.fusion_metrics.json",
            sample=SAMPLES
        ),
        partner_counts=expand(
            "results/{{experiment}}/counts/{sample}.partner_counts.csv",
            sample=SAMPLES
        ),
        variant_catalog="results/{experiment}/references/variant_catalog.csv",
        merge_logs=expand(
            "logs/{{experiment}}/bbmerge/{sample}.log",
            sample=SAMPLES
        ),
        trim_stats=expand(
            "stats/{{experiment}}/trim/{sample}.trim.stats.txt",
            sample=SAMPLES
        ),
        contam_stats=expand(
            "stats/{{experiment}}/contam/{sample}.contam.stats.txt",
            sample=SAMPLES
        ),
        quality_stats=expand(
            "stats/{{experiment}}/quality/{sample}.quality.stats.txt",
            sample=SAMPLES
        ),
        trim_logs=expand(
            "logs/{{experiment}}/bbduk/{sample}.trim.log",
            sample=SAMPLES
        ),
        contam_logs=expand(
            "logs/{{experiment}}/bbduk/{sample}.clean.log",
            sample=SAMPLES
        ),
        quality_logs=expand(
            "logs/{{experiment}}/bbduk/{sample}.quality.log",
            sample=SAMPLES
        )
    output:
        summary="results/{experiment}/fusion_counts_summary.csv",
        qc_metrics="results/{experiment}/fusion_qc_metrics.csv",
        partner_summary="results/{experiment}/partner_counts_summary.csv",
        sensitivity_metrics="results/{experiment}/sensitivity_metrics.csv",
        decay_metrics="results/{experiment}/decay_metrics.csv",
        trim_metrics="results/{experiment}/trim_metrics.csv",
        contam_metrics="results/{experiment}/contam_metrics.csv",
        quality_metrics="results/{experiment}/quality_metrics.csv"
    params:
        mode="merged",
        breakpoint_window=BREAKPOINT_WINDOW,
        ihist_base_path=f"stats/{EXPERIMENT}/merge"
    script:
        "../scripts/aggregate_counts.py"


rule aggregate_unmerged_counts:
    """
    Aggregate unmerged fusion counts across all samples into summary files.
    """
    input:
        counts=expand(
            "results/{{experiment}}/counts/{sample}.{mate}.unmerged_fusion_counts.csv",
            sample=SAMPLES,
            mate=["R1", "R2"]
        ),
        metrics=expand(
            "results/{{experiment}}/counts/{sample}.{mate}.unmerged_fusion_metrics.json",
            sample=SAMPLES,
            mate=["R1", "R2"]
        ),
        partner_counts=expand(
            "results/{{experiment}}/counts/{sample}.{mate}.unmerged_partner_counts.csv",
            sample=SAMPLES,
            mate=["R1", "R2"]
        )
    output:
        summary="results/{experiment}/unmerged_counts_summary.csv",
        qc_metrics="results/{experiment}/unmerged_qc_metrics.csv",
        partner_summary="results/{experiment}/unmerged_partner_counts_summary.csv"
    params:
        mode="unmerged"
    script:
        "../scripts/aggregate_counts.py"
