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
        unfused="results/{experiment}/references/unfused_sequences.csv"
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
        metrics="results/{experiment}/counts/{sample}.fusion_metrics.json"
    params:
        show_progress=SHOW_PROGRESS,
        progress_interval=PROGRESS_INTERVAL,
        orientation_check=ORIENTATION_CHECK,
        prefilter_fallback=PREFILTER_FALLBACK,
        linker_first=LINKER_FIRST,
        linker_sequence=LINKER_SEQUENCE,
        breakpoint_window=BREAKPOINT_WINDOW
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
        summary="results/{experiment}/fusion_counts_summary.csv",
        qc_metrics="results/{experiment}/fusion_qc_metrics.csv"
    run:
        import pandas as pd

        # Load and merge all count files
        dfs = []
        for count_file in input.counts:
            sample = Path(count_file).stem.replace('.fusion_counts', '')
            df = pd.read_csv(count_file)
            # Preserve type column if present
            if 'type' in df.columns:
                df = df[['fusion_id', 'type', 'count']]
            else:
                # Backward compatibility: add default type
                df['type'] = 'fusion'
                df = df[['fusion_id', 'type', 'count']]
            df = df.rename(columns={'count': sample})
            dfs.append(df.set_index('fusion_id'))

        # Merge on fusion_id
        # First, extract type column (should be same across all samples)
        type_col = None
        for df in dfs:
            if 'type' in df.columns:
                type_col = df['type']
                break

        # Merge count columns
        count_dfs = [df.drop(columns=['type'], errors='ignore') for df in dfs]
        merged = pd.concat(count_dfs, axis=1).fillna(0).astype(int)
        merged = merged.reset_index()

        # Add type column back if available
        if type_col is not None:
            merged = merged.merge(
                type_col.reset_index()[['fusion_id', 'type']],
                on='fusion_id',
                how='left'
            )
            # Reorder columns: fusion_id, type, then sample columns
            sample_cols = [c for c in merged.columns if c not in ['fusion_id', 'type']]
            merged = merged[['fusion_id', 'type'] + sample_cols]
        else:
            # Backward compatibility: add default type
            merged['type'] = 'fusion'
            sample_cols = [c for c in merged.columns if c not in ['fusion_id', 'type']]
            merged = merged[['fusion_id', 'type'] + sample_cols]

        # Sort by total counts
        merged['total'] = merged[sample_cols].sum(axis=1)
        merged = merged.sort_values('total', ascending=False)

        merged.to_csv(output.summary, index=False)

        # QC composition metrics per sample
        metrics = []
        sample_cols = [c for c in merged.columns if c not in ['fusion_id', 'type', 'total']]
        for sample in sample_cols:
            counts_series = merged[[sample]].fillna(0)[sample]
            total = counts_series.sum()
            if total == 0:
                metrics.append({
                    "sample": sample,
                    "total_counts": 0,
                    "unique_fusions": 0,
                    "top1_fraction": 0.0,
                    "top10_fraction": 0.0,
                })
                continue

            sorted_counts = counts_series.sort_values(ascending=False)
            top1 = sorted_counts.iloc[0] if len(sorted_counts) else 0
            top10 = sorted_counts.head(10).sum()

            metrics.append({
                "sample": sample,
                "total_counts": int(total),
                "unique_fusions": int((counts_series > 0).sum()),
                "top1_fraction": float(top1 / total),
                "top10_fraction": float(top10 / total),
            })

        pd.DataFrame(metrics).to_csv(output.qc_metrics, index=False)
