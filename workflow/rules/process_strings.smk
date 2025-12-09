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

    These files are used by the string matching detection step.
    """
    input:
        sequences=get_reference_fasta()
    output:
        breakpoints="results/{experiment}/references/breakpoint_sequences.csv",
        ends="results/{experiment}/references/domain_ends.csv"
    params:
        partners_file=PARTNERS_FILE,
        anchor_name=ANCHOR_NAME,
        anchor_position=ANCHOR_POSITION,
        linker_sequence=LINKER_SEQUENCE,
        breakpoint_window=BREAKPOINT_WINDOW,
        maintain_frame=MAINTAIN_FRAME,
        kmer_size=KMER_SIZE
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

    Input is the merged (error-corrected) reads from preprocessing.
    """
    input:
        fastq="results/{experiment}/merged/{sample}_merged.fastq.gz",
        breakpoints="results/{experiment}/references/breakpoint_sequences.csv",
        ends="results/{experiment}/references/domain_ends.csv"
    output:
        "results/{experiment}/counts/{sample}.fusion_counts.csv"
    params:
        show_progress=SHOW_PROGRESS,
        progress_interval=PROGRESS_INTERVAL
    log:
        "logs/{experiment}/string_match/{sample}.log"
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
        )
    output:
        "results/{experiment}/fusion_counts_summary.csv"
    run:
        import pandas as pd

        # Load and merge all count files
        dfs = []
        for count_file in input.counts:
            sample = Path(count_file).stem.replace('.fusion_counts', '')
            df = pd.read_csv(count_file)
            df = df.rename(columns={'count': sample})
            dfs.append(df.set_index('fusion_id'))

        # Merge on fusion_id
        merged = pd.concat(dfs, axis=1).fillna(0).astype(int)
        merged = merged.reset_index()

        # Sort by total counts
        merged['total'] = merged.drop('fusion_id', axis=1).sum(axis=1)
        merged = merged.sort_values('total', ascending=False)

        merged.to_csv(output[0], index=False)
