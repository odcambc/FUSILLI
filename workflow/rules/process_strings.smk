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
        linker_first=LINKER_FIRST,
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
        linker_first=LINKER_FIRST,
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
    run:
        import json
        import pandas as pd
        import re
        import zipfile

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
        def _parse_bbmerge_log(path):
            total_pairs = None
            joined = None
            with open(path, 'r') as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith("Pairs:"):
                        match = re.search(r"Pairs:\s*([0-9,]+)", line)
                        if match:
                            total_pairs = int(match.group(1).replace(",", ""))
                    elif line.startswith("Joined:"):
                        match = re.search(r"Joined:\s*([0-9,]+)", line)
                        if match:
                            joined = int(match.group(1).replace(",", ""))
            if total_pairs is None or joined is None:
                return None
            unmerged = max(total_pairs - joined, 0)
            unmerged_fraction = unmerged / total_pairs if total_pairs else 0.0
            return {
                "total_pairs": total_pairs,
                "merged_pairs": joined,
                "unmerged_pairs": unmerged,
                "unmerged_fraction": float(unmerged_fraction),
            }

        merge_stats = {}
        for log_path in input.merge_logs:
            sample = Path(log_path).stem
            try:
                stats = _parse_bbmerge_log(log_path)
            except FileNotFoundError:
                stats = None
            merge_stats[sample] = stats

        for sample in sample_cols:
            counts_series = merged[[sample]].fillna(0)[sample]
            total = counts_series.sum()

            fusion_mask = merged['type'] == 'fusion'
            unfused_mask = merged['type'] == 'unfused'
            fusion_series = merged.loc[fusion_mask, sample].fillna(0)
            unfused_series = merged.loc[unfused_mask, sample].fillna(0)

            expected_total = len(counts_series)
            expected_fusions = len(fusion_series)
            expected_unfused = len(unfused_series)

            zero_fraction = (
                float((counts_series == 0).sum() / expected_total)
                if expected_total else 0.0
            )
            fusion_zero_fraction = (
                float((fusion_series == 0).sum() / expected_fusions)
                if expected_fusions else 0.0
            )
            unfused_zero_fraction = (
                float((unfused_series == 0).sum() / expected_unfused)
                if expected_unfused else 0.0
            )

            sorted_fusions = fusion_series.sort_values(ascending=False)
            top1 = sorted_fusions.iloc[0] if len(sorted_fusions) else 0
            top10 = sorted_fusions.head(10).sum()
            total_fusion_counts = fusion_series.sum()

            merge_stat = merge_stats.get(sample)

            metrics.append({
                "sample": sample,
                "total_counts": int(total),
                "total_fusion_counts": int(total_fusion_counts),
                "total_unfused_counts": int(unfused_series.sum()),
                "expected_variants": int(expected_total),
                "expected_fusions": int(expected_fusions),
                "expected_unfused": int(expected_unfused),
                "unique_fusions": int((fusion_series > 0).sum()),
                "unique_unfused": int((unfused_series > 0).sum()),
                "zero_fraction": float(zero_fraction),
                "fusion_zero_fraction": float(fusion_zero_fraction),
                "unfused_zero_fraction": float(unfused_zero_fraction),
                "top1_fraction": float(top1 / total_fusion_counts) if total_fusion_counts else 0.0,
                "top10_fraction": float(top10 / total_fusion_counts) if total_fusion_counts else 0.0,
                "total_pairs": int(merge_stat["total_pairs"]) if merge_stat else 0,
                "merged_pairs": int(merge_stat["merged_pairs"]) if merge_stat else 0,
                "unmerged_pairs": int(merge_stat["unmerged_pairs"]) if merge_stat else 0,
                "unmerged_fraction": float(merge_stat["unmerged_fraction"]) if merge_stat else 0.0,
            })

        metrics_df = pd.DataFrame(metrics)

        # Merge in per-sample detection metrics from JSON
        json_metrics = []
        for metrics_file in input.metrics:
            sample = Path(metrics_file).stem.replace('.fusion_metrics', '')
            try:
                with open(metrics_file, 'r') as f:
                    data = json.load(f)
            except FileNotFoundError:
                data = {}
            json_metrics.append({
                "sample": sample,
                "partner_end_reads": int(data.get("partner_end_reads", 0)),
                "partner_linker_reads": int(data.get("partner_linker_reads", 0)),
                "unique_partners_detected": int(data.get("unique_partners_detected", 0)),
                "unique_partner_linker_detected": int(data.get("unique_partner_linker_detected", 0)),
                "matched_reads": int(data.get("matched_reads", 0)),
                "reads_processed": int(data.get("reads_processed", 0)),
            })

        json_df = pd.DataFrame(json_metrics)
        metrics_df = metrics_df.merge(json_df, on="sample", how="left")
        metrics_df.to_csv(output.qc_metrics, index=False)

        # Aggregate partner counts across samples
        partner_dfs = []
        for partner_file in input.partner_counts:
            sample = Path(partner_file).stem.replace('.partner_counts', '')
            df = pd.read_csv(partner_file)
            df = df.rename(columns={
                "partner_end_count": f"{sample}_end",
                "partner_linker_count": f"{sample}_linker"
            })
            partner_dfs.append(df.set_index("partner_name"))

        if partner_dfs:
            partner_merged = pd.concat(partner_dfs, axis=1).fillna(0).astype(int)
            partner_merged = partner_merged.reset_index()
            partner_merged.to_csv(output.partner_summary, index=False)
        else:
            pd.DataFrame(
                columns=["partner_name"]
            ).to_csv(output.partner_summary, index=False)

        def _parse_fastqc_lengths(zip_path):
            if not Path(zip_path).exists():
                return None
            with zipfile.ZipFile(zip_path) as zf:
                data_files = [n for n in zf.namelist() if n.endswith("fastqc_data.txt")]
                if not data_files:
                    return None
                with zf.open(data_files[0]) as fh:
                    lines = [l.decode("utf-8").strip() for l in fh.readlines()]
            in_section = False
            lengths = []
            for line in lines:
                if line.startswith(">>Sequence Length Distribution"):
                    in_section = True
                    continue
                if in_section:
                    if line.startswith(">>END_MODULE"):
                        break
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) < 2:
                        continue
                    length_token, count_token = parts[0], parts[1]
                    try:
                        count = float(count_token)
                    except ValueError:
                        continue
                    if "-" in length_token:
                        start, end = length_token.split("-", 1)
                        try:
                            start_v = float(start)
                            end_v = float(end)
                        except ValueError:
                            continue
                        length = (start_v + end_v) / 2.0
                        lengths.append((length, count, start_v, end_v))
                    else:
                        try:
                            length = float(length_token)
                        except ValueError:
                            continue
                        lengths.append((length, count, length, length))
            if not lengths:
                return None
            total = sum(c for _, c, _, _ in lengths)
            mean = sum(length * c for length, c, _, _ in lengths) / total if total else 0.0
            return {"lengths": lengths, "total": total, "mean": mean}

        def _fraction_reads_long_enough(lengths, threshold):
            if not lengths:
                return 0.0
            total = sum(c for _, c, _, _ in lengths)
            if total == 0:
                return 0.0
            long_count = 0.0
            for _, count, start_v, end_v in lengths:
                if end_v < threshold:
                    continue
                if start_v >= threshold:
                    long_count += count
                else:
                    # partial overlap within range
                    if end_v > start_v:
                        frac = (end_v - threshold) / (end_v - start_v)
                        long_count += count * max(0.0, min(1.0, frac))
            return long_count / total

        def _parse_ihist(path):
            if not Path(path).exists():
                return None
            rows = []
            with open(path, "r") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) < 2:
                        continue
                    try:
                        size = float(parts[0])
                        count = float(parts[1])
                    except ValueError:
                        continue
                    rows.append((size, count))
            if not rows:
                return None
            total = sum(c for _, c in rows)
            rows.sort(key=lambda x: x[0])
            cum = 0.0
            median = 0.0
            for size, count in rows:
                cum += count
                if cum >= total / 2.0:
                    median = size
                    break
            return {"median": median}

        sensitivity_rows = []
        for sample in sample_cols:
            fastqc_zip = Path(f"stats/{EXPERIMENT}/fastqc/{sample}_R1.fastqc.zip")
            lengths_info = _parse_fastqc_lengths(fastqc_zip)
            lengths = lengths_info["lengths"] if lengths_info else None
            mean_len = lengths_info["mean"] if lengths_info else 0.0

            threshold = 2 * BREAKPOINT_WINDOW
            expected_fraction = _fraction_reads_long_enough(lengths, threshold)

            matched_reads = int(json_df.loc[json_df["sample"] == sample, "matched_reads"].fillna(0).iloc[0]) if not json_df.empty else 0
            reads_processed = int(json_df.loc[json_df["sample"] == sample, "reads_processed"].fillna(0).iloc[0]) if not json_df.empty else 0
            expected_reads = reads_processed * expected_fraction
            sensitivity_index = (matched_reads / expected_reads) if expected_reads else 0.0

            ihist_path = Path(f"stats/{EXPERIMENT}/merge/{sample}.ihist")
            ihist_info = _parse_ihist(ihist_path)
            overlap_median = ihist_info["median"] if ihist_info else 0.0

            sensitivity_rows.append({
                "sample": sample,
                "read_length_mean": float(mean_len),
                "overlap_median": float(overlap_median),
                "breakpoint_kmer_length": int(threshold),
                "expected_detection_fraction": float(expected_fraction),
                "sensitivity_index": float(sensitivity_index),
            })

        pd.DataFrame(sensitivity_rows).to_csv(output.sensitivity_metrics, index=False)

        def _parse_bbduk_stats(path):
            if not Path(path).exists():
                return None
            input_reads = None
            input_bases = None
            output_reads = None
            output_bases = None
            with open(path, "r") as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith("Input:"):
                        match = re.search(r"Input:\s*([0-9,]+)\s+reads\s+([0-9,]+)\s+bases", line)
                        if match:
                            input_reads = int(match.group(1).replace(",", ""))
                            input_bases = int(match.group(2).replace(",", ""))
                    elif line.startswith("Output:"):
                        match = re.search(r"Output:\s*([0-9,]+)\s+reads\s+([0-9,]+)\s+bases", line)
                        if match:
                            output_reads = int(match.group(1).replace(",", ""))
                            output_bases = int(match.group(2).replace(",", ""))
            if input_reads is None or output_reads is None:
                return None
            return {
                "input_reads": input_reads,
                "input_bases": input_bases,
                "output_reads": output_reads,
                "output_bases": output_bases,
            }

        def _parse_bbduk_log(path):
            if not Path(path).exists():
                return None
            input_reads = None
            input_bases = None
            result_reads = None
            result_bases = None
            with open(path, "r") as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith("Input:"):
                        match = re.search(r"Input:\s*([0-9,]+)\s+reads\s+([0-9,]+)\s+bases", line)
                        if match:
                            input_reads = int(match.group(1).replace(",", ""))
                            input_bases = int(match.group(2).replace(",", ""))
                    elif line.startswith("Result:"):
                        match = re.search(r"Result:\s*([0-9,]+)\s+reads\s+([0-9,]+)\s+bases", line)
                        if match:
                            result_reads = int(match.group(1).replace(",", ""))
                            result_bases = int(match.group(2).replace(",", ""))
            if input_reads is None or result_reads is None:
                return None
            return {
                "input_reads": input_reads,
                "input_bases": input_bases,
                "output_reads": result_reads,
                "output_bases": result_bases,
            }

        def _parse_bbmerge_stats(path):
            if not Path(path).exists():
                return None
            joined = None
            avg_insert = None
            with open(path, "r") as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith("Joined:"):
                        match = re.search(r"Joined:\s*([0-9,]+)", line)
                        if match:
                            joined = int(match.group(1).replace(",", ""))
                    elif line.startswith("Avg Insert"):
                        match = re.search(r"Avg Insert:\s*([0-9.]+)", line)
                        if match:
                            avg_insert = float(match.group(1))
            if joined is None:
                return None
            merged_bases = int(round(joined * avg_insert)) if avg_insert else None
            return {
                "merged_reads": joined,
                "merged_bases": merged_bases,
            }

        trim_stats = {}
        for path in input.trim_stats:
            sample = Path(path).stem.replace(".trim", "")
            trim_stats[sample] = _parse_bbduk_stats(path)

        contam_stats = {}
        for path in input.contam_stats:
            sample = Path(path).stem.replace(".contam", "")
            contam_stats[sample] = _parse_bbduk_stats(path)

        quality_stats = {}
        for path in input.quality_stats:
            sample = Path(path).stem.replace(".quality", "")
            quality_stats[sample] = _parse_bbduk_stats(path)

        trim_logs = {}
        for path in input.trim_logs:
            sample = Path(path).stem.replace(".trim", "")
            trim_logs[sample] = _parse_bbduk_log(path)

        contam_logs = {}
        for path in input.contam_logs:
            sample = Path(path).stem.replace(".clean", "")
            contam_logs[sample] = _parse_bbduk_log(path)

        quality_logs = {}
        for path in input.quality_logs:
            sample = Path(path).stem.replace(".quality", "")
            quality_logs[sample] = _parse_bbduk_log(path)

        merge_stats_simple = {}
        for path in input.merge_logs:
            sample = Path(path).stem
            merge_stats_simple[sample] = _parse_bbmerge_stats(path)

        decay_rows = []
        for sample in sample_cols:
            trim = trim_stats.get(sample)
            contam = contam_stats.get(sample)
            quality = quality_stats.get(sample)
            merged = merge_stats_simple.get(sample)

            trim_log = trim_logs.get(sample)
            contam_log = contam_logs.get(sample)
            quality_log = quality_logs.get(sample)

            raw_reads = trim_log["input_reads"] if trim_log else (trim["input_reads"] if trim else 0)
            raw_bases = trim_log["input_bases"] if trim_log else (trim["input_bases"] if trim else None)

            steps = [
                ("raw", raw_reads, raw_bases),
                ("trimmed", trim_log["output_reads"] if trim_log else (trim["output_reads"] if trim else 0),
                 trim_log["output_bases"] if trim_log else (trim["output_bases"] if trim else None)),
                ("cleaned", contam_log["output_reads"] if contam_log else (contam["output_reads"] if contam else 0),
                 contam_log["output_bases"] if contam_log else (contam["output_bases"] if contam else None)),
                ("quality", quality_log["output_reads"] if quality_log else (quality["output_reads"] if quality else 0),
                 quality_log["output_bases"] if quality_log else (quality["output_bases"] if quality else None)),
                ("merged", merged["merged_reads"] if merged else 0, merged["merged_bases"] if merged else None),
            ]

            for step, reads, bases in steps:
                decay_rows.append({
                    "sample": sample,
                    "step": step,
                    "reads": int(reads) if reads is not None else 0,
                    "bases": int(bases) if bases is not None else None,
                    "read_fraction": float(reads / raw_reads) if raw_reads else 0.0,
                    "base_fraction": float(bases / raw_bases) if raw_bases and bases is not None else 0.0,
                })

        decay_df = pd.DataFrame(decay_rows)
        decay_df.to_csv(output.decay_metrics, index=False)

        for step_name, out_path in [
            ("trimmed", output.trim_metrics),
            ("cleaned", output.contam_metrics),
            ("quality", output.quality_metrics),
        ]:
            step_df = decay_df[decay_df["step"] == step_name].copy()
            step_df = step_df.drop(columns=["step"], errors="ignore")
            step_df.to_csv(out_path, index=False)


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
    run:
        import json
        import pandas as pd

        dfs = []
        for count_file in input.counts:
            sample = Path(count_file).stem.replace('.unmerged_fusion_counts', '')
            df = pd.read_csv(count_file)
            if 'type' in df.columns:
                df = df[['fusion_id', 'type', 'count']]
            else:
                df['type'] = 'fusion'
                df = df[['fusion_id', 'type', 'count']]
            df = df.rename(columns={'count': sample})
            dfs.append(df.set_index('fusion_id'))

        type_col = None
        for df in dfs:
            if 'type' in df.columns:
                type_col = df['type']
                break

        count_dfs = [df.drop(columns=['type'], errors='ignore') for df in dfs]
        merged = pd.concat(count_dfs, axis=1).fillna(0).astype(int)
        merged = merged.reset_index()

        if type_col is not None:
            merged = merged.merge(
                type_col.reset_index()[['fusion_id', 'type']],
                on='fusion_id',
                how='left'
            )
            sample_cols = [c for c in merged.columns if c not in ['fusion_id', 'type']]
            merged = merged[['fusion_id', 'type'] + sample_cols]
        else:
            merged['type'] = 'fusion'
            sample_cols = [c for c in merged.columns if c not in ['fusion_id', 'type']]
            merged = merged[['fusion_id', 'type'] + sample_cols]

        merged['total'] = merged[sample_cols].sum(axis=1)
        merged = merged.sort_values('total', ascending=False)
        merged.to_csv(output.summary, index=False)

        metrics = []
        sample_cols = [c for c in merged.columns if c not in ['fusion_id', 'type', 'total']]
        for sample in sample_cols:
            counts_series = merged[[sample]].fillna(0)[sample]
            total = counts_series.sum()

            fusion_mask = merged['type'] == 'fusion'
            unfused_mask = merged['type'] == 'unfused'
            fusion_series = merged.loc[fusion_mask, sample].fillna(0)
            unfused_series = merged.loc[unfused_mask, sample].fillna(0)

            expected_total = len(counts_series)
            expected_fusions = len(fusion_series)
            expected_unfused = len(unfused_series)

            zero_fraction = (
                float((counts_series == 0).sum() / expected_total)
                if expected_total else 0.0
            )
            fusion_zero_fraction = (
                float((fusion_series == 0).sum() / expected_fusions)
                if expected_fusions else 0.0
            )
            unfused_zero_fraction = (
                float((unfused_series == 0).sum() / expected_unfused)
                if expected_unfused else 0.0
            )

            sorted_fusions = fusion_series.sort_values(ascending=False)
            top1 = sorted_fusions.iloc[0] if len(sorted_fusions) else 0
            top10 = sorted_fusions.head(10).sum()
            total_fusion_counts = fusion_series.sum()

            metrics.append({
                "sample": sample,
                "total_counts": int(total),
                "total_fusion_counts": int(total_fusion_counts),
                "total_unfused_counts": int(unfused_series.sum()),
                "expected_variants": int(expected_total),
                "expected_fusions": int(expected_fusions),
                "expected_unfused": int(expected_unfused),
                "unique_fusions": int((fusion_series > 0).sum()),
                "unique_unfused": int((unfused_series > 0).sum()),
                "zero_fraction": float(zero_fraction),
                "fusion_zero_fraction": float(fusion_zero_fraction),
                "unfused_zero_fraction": float(unfused_zero_fraction),
                "top1_fraction": float(top1 / total_fusion_counts) if total_fusion_counts else 0.0,
                "top10_fraction": float(top10 / total_fusion_counts) if total_fusion_counts else 0.0,
            })

        metrics_df = pd.DataFrame(metrics)

        json_metrics = []
        for metrics_file in input.metrics:
            sample = Path(metrics_file).stem.replace('.unmerged_fusion_metrics', '')
            try:
                with open(metrics_file, 'r') as f:
                    data = json.load(f)
            except FileNotFoundError:
                data = {}
            json_metrics.append({
                "sample": sample,
                "partner_end_reads": int(data.get("partner_end_reads", 0)),
                "partner_linker_reads": int(data.get("partner_linker_reads", 0)),
                "unique_partners_detected": int(data.get("unique_partners_detected", 0)),
                "unique_partner_linker_detected": int(data.get("unique_partner_linker_detected", 0)),
            })

        json_df = pd.DataFrame(json_metrics)
        metrics_df = metrics_df.merge(json_df, on="sample", how="left")
        metrics_df.to_csv(output.qc_metrics, index=False)

        partner_dfs = []
        for partner_file in input.partner_counts:
            sample = Path(partner_file).stem.replace('.unmerged_partner_counts', '')
            df = pd.read_csv(partner_file)
            df = df.rename(columns={
                "partner_end_count": f"{sample}_end",
                "partner_linker_count": f"{sample}_linker"
            })
            partner_dfs.append(df.set_index("partner_name"))

        if partner_dfs:
            partner_merged = pd.concat(partner_dfs, axis=1).fillna(0).astype(int)
            partner_merged = partner_merged.reset_index()
            partner_merged.to_csv(output.partner_summary, index=False)
        else:
            pd.DataFrame(
                columns=["partner_name"]
            ).to_csv(output.partner_summary, index=False)
