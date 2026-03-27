#!/usr/bin/env python3
"""
Aggregate fusion counts across samples into summary files.

This script aggregates per-sample fusion counts and QC metrics into:
- fusion_counts_summary.csv: Per-fusion counts across all samples
- fusion_qc_metrics.csv: Per-sample QC metrics
- partner_counts_summary.csv: Partner-level counts across samples
- sensitivity_metrics.csv: Detection sensitivity metrics (merged mode only)
- decay_metrics.csv: Read count decay through pipeline steps (merged mode only)

USAGE:
    # Via Snakemake (automatic parameter passing)
    script: "scripts/aggregate_counts.py"

    # Standalone merged mode
    python aggregate_counts.py merged \
        --counts results/{exp}/counts/*.fusion_counts.csv \
        --metrics results/{exp}/counts/*.fusion_metrics.json \
        --partner-counts results/{exp}/counts/*.partner_counts.csv \
        --variant-catalog results/{exp}/references/variant_catalog.csv \
        --merge-logs logs/{exp}/bbmerge/*.log \
        --trim-logs logs/{exp}/bbduk/*.trim.log \
        --contam-logs logs/{exp}/bbduk/*.clean.log \
        --quality-logs logs/{exp}/bbduk/*.quality.log \
        --trim-stats stats/{exp}/trim/*.trim.stats.txt \
        --contam-stats stats/{exp}/contam/*.contam.stats.txt \
        --quality-stats stats/{exp}/quality/*.quality.stats.txt \
        --experiment {exp} \
        --breakpoint-window {breakpoint_window} \
        --output-summary results/{exp}/fusion_counts_summary.csv \
        --output-qc results/{exp}/fusion_qc_metrics.csv \
        --output-partner results/{exp}/partner_counts_summary.csv \
        --output-sensitivity results/{exp}/sensitivity_metrics.csv \
        --output-decay results/{exp}/decay_metrics.csv \
        --output-trim results/{exp}/trim_metrics.csv \
        --output-contam results/{exp}/contam_metrics.csv \
        --output-quality results/{exp}/quality_metrics.csv

    # Standalone unmerged mode
    python aggregate_counts.py unmerged \
        --counts results/{exp}/counts/*.R1.unmerged_fusion_counts.csv \
        --metrics results/{exp}/counts/*.R1.unmerged_fusion_metrics.json \
        --partner-counts results/{exp}/counts/*.R1.unmerged_partner_counts.csv \
        --output-summary results/{exp}/unmerged_counts_summary.csv \
        --output-qc results/{exp}/unmerged_qc_metrics.csv \
        --output-partner results/{exp}/unmerged_partner_counts_summary.csv
"""

import argparse
import json
from pathlib import Path
from typing import Any

try:
    from utils import (
        parse_bbmerge_log,
        parse_bbduk_log,
        parse_bbmerge_stats,
        parse_bbduk_stats,
        parse_ihist,
    )
except ImportError:
    from workflow.scripts.utils import (
        parse_bbmerge_log,
        parse_bbduk_log,
        parse_bbmerge_stats,
        parse_bbduk_stats,
        parse_ihist,
    )


def load_and_merge_counts(
    count_files: list[str],
    type_suffix: str = "fusion_counts",
) -> tuple[Any, list[str]]:
    """
    Load and merge all count files into a single DataFrame.

    Args:
        count_files: List of paths to count CSV files
        type_suffix: Suffix to strip from filename to get sample name

    Returns:
        Tuple of (merged DataFrame, list of sample column names)
    """
    import pandas as pd

    dfs = []
    for count_file in count_files:
        sample = Path(count_file).stem.replace(f".{type_suffix}", "")
        df = pd.read_csv(count_file)
        if "type" in df.columns:
            df = df[["fusion_id", "type", "count"]]
        else:
            df["type"] = "fusion"
            df = df[["fusion_id", "type", "count"]]
        df = df.rename(columns={"count": sample})
        dfs.append(df.set_index("fusion_id"))

    type_col = None
    for df in dfs:
        if "type" in df.columns:
            type_col = df["type"]
            break

    count_dfs = [df.drop(columns=["type"], errors="ignore") for df in dfs]
    merged = pd.concat(count_dfs, axis=1).fillna(0).astype(int)
    merged = merged.reset_index()

    if type_col is not None:
        merged = merged.merge(
            type_col.reset_index()[["fusion_id", "type"]], on="fusion_id", how="left"
        )
        sample_cols = [c for c in merged.columns if c not in ["fusion_id", "type"]]
        merged = merged[["fusion_id", "type"] + sample_cols]
    else:
        merged["type"] = "fusion"
        sample_cols = [c for c in merged.columns if c not in ["fusion_id", "type"]]
        merged = merged[["fusion_id", "type"] + sample_cols]

    merged["total"] = merged[sample_cols].sum(axis=1)
    merged = merged.sort_values("total", ascending=False)

    return merged, sample_cols


def calculate_qc_metrics(
    merged: Any,
    sample_cols: list[str],
    merge_stats: dict | None = None,
    json_df: Any | None = None,
) -> Any:
    """
    Calculate QC metrics per sample.

    Args:
        merged: Merged counts DataFrame
        sample_cols: List of sample column names
        merge_stats: Dict mapping sample name to merge stats from parse_bbmerge_log
        json_df: DataFrame with per-sample detection metrics from JSON files

    Returns:
        DataFrame with QC metrics per sample
    """
    import pandas as pd

    metrics = []
    for sample in sample_cols:
        counts_series = merged[[sample]].fillna(0)[sample]
        total = counts_series.sum()

        fusion_mask = merged["type"] == "fusion"
        unfused_mask = merged["type"] == "unfused"
        fusion_series = merged.loc[fusion_mask, sample].fillna(0)
        unfused_series = merged.loc[unfused_mask, sample].fillna(0)

        expected_total = len(counts_series)
        expected_fusions = len(fusion_series)
        expected_unfused = len(unfused_series)

        zero_fraction = (
            float((counts_series == 0).sum() / expected_total)
            if expected_total
            else 0.0
        )
        fusion_zero_fraction = (
            float((fusion_series == 0).sum() / expected_fusions)
            if expected_fusions
            else 0.0
        )
        unfused_zero_fraction = (
            float((unfused_series == 0).sum() / expected_unfused)
            if expected_unfused
            else 0.0
        )

        sorted_fusions = fusion_series.sort_values(ascending=False)
        top1 = sorted_fusions.iloc[0] if len(sorted_fusions) else 0
        top10 = sorted_fusions.head(10).sum()
        total_fusion_counts = fusion_series.sum()

        merge_stat = merge_stats.get(sample) if merge_stats else None

        row = {
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
            "top1_fraction": float(top1 / total_fusion_counts)
            if total_fusion_counts
            else 0.0,
            "top10_fraction": float(top10 / total_fusion_counts)
            if total_fusion_counts
            else 0.0,
        }

        if merge_stat:
            row["total_pairs"] = int(merge_stat["total_pairs"])
            row["merged_pairs"] = int(merge_stat["merged_pairs"])
            row["unmerged_pairs"] = int(merge_stat["unmerged_pairs"])
            row["unmerged_fraction"] = float(merge_stat["unmerged_fraction"])
        else:
            row["total_pairs"] = 0
            row["merged_pairs"] = 0
            row["unmerged_pairs"] = 0
            row["unmerged_fraction"] = 0.0

        metrics.append(row)

    metrics_df = pd.DataFrame(metrics)

    if json_df is not None and not json_df.empty:
        metrics_df = metrics_df.merge(json_df, on="sample", how="left")

    return metrics_df


def calculate_diversity_metrics(
    metrics_df: Any,
    merged: Any,
    fusion_catalog: Any,
    expected_breakpoints: int,
    expected_partners: int,
) -> Any:
    """
    Calculate diversity, coverage, and yield metrics.

    Args:
        metrics_df: DataFrame with base QC metrics
        merged: Merged counts DataFrame
        fusion_catalog: Variant catalog DataFrame
        expected_breakpoints: Expected number of unique breakpoints
        expected_partners: Expected number of unique partners

    Returns:
        Updated metrics DataFrame with diversity metrics added
    """
    import numpy as np

    for idx, row in metrics_df.iterrows():
        sample = row["sample"]

        fusion_data = merged[merged["type"] == "fusion"].copy()
        if sample in fusion_data.columns:
            counts = fusion_data[sample].fillna(0).values
            total_counts = counts.sum()

            non_zero_counts = counts[counts > 0]

            if len(non_zero_counts) > 0 and total_counts > 0:
                proportions = non_zero_counts / total_counts
                shannon_diversity = float(-np.sum(proportions * np.log(proportions)))
                simpson_diversity = float(1 - np.sum(proportions**2))
                observed_variants = len(non_zero_counts)
                if observed_variants > 1:
                    evenness = float(shannon_diversity / np.log(observed_variants))
                else:
                    evenness = 0.0
            else:
                shannon_diversity = 0.0
                simpson_diversity = 0.0
                evenness = 0.0
                observed_variants = 0

            unique_fusions = int(row.get("unique_fusions", 0))
            expected_fusions = int(row.get("expected_fusions", 0))
            variant_coverage = (
                float(unique_fusions / expected_fusions)
                if expected_fusions > 0
                else 0.0
            )

            detected_fusions_df = fusion_data[fusion_data[sample] > 0]
            if (
                fusion_catalog is not None
                and "breakpoint_nt" in fusion_catalog.columns
                and "fusion_id" in fusion_catalog.columns
            ):
                detected_breakpoints = (
                    detected_fusions_df.merge(
                        fusion_catalog[["fusion_id", "breakpoint_nt"]],
                        on="fusion_id",
                        how="left",
                    )["breakpoint_nt"]
                    .dropna()
                    .nunique()
                )
                if detected_breakpoints == 0 and unique_fusions > 0:
                    detected_breakpoints = unique_fusions
            else:
                detected_breakpoints = unique_fusions

            breakpoint_coverage = (
                float(detected_breakpoints / expected_breakpoints)
                if expected_breakpoints > 0
                else 0.0
            )

            unique_partners_detected = int(row.get("unique_partners_detected", 0))
            partner_coverage = (
                float(unique_partners_detected / expected_partners)
                if expected_partners > 0
                else 0.0
            )

            matched_reads = int(row.get("matched_reads", 0))
            reads_processed = int(row.get("reads_processed", 0))
            total_fusion_counts = int(row.get("total_fusion_counts", 0))

            detections_per_read = (
                float(total_fusion_counts / reads_processed)
                if reads_processed > 0
                else 0.0
            )
            detections_per_million = (
                float(detections_per_read * 1000000) if reads_processed > 0 else 0.0
            )
        else:
            shannon_diversity = 0.0
            simpson_diversity = 0.0
            evenness = 0.0
            variant_coverage = 0.0
            breakpoint_coverage = 0.0
            partner_coverage = 0.0
            detections_per_read = 0.0
            detections_per_million = 0.0
            observed_variants = 0

        metrics_df.at[idx, "shannon_diversity"] = shannon_diversity
        metrics_df.at[idx, "simpson_diversity"] = simpson_diversity
        metrics_df.at[idx, "evenness"] = evenness
        metrics_df.at[idx, "variant_coverage"] = variant_coverage
        metrics_df.at[idx, "breakpoint_coverage"] = breakpoint_coverage
        metrics_df.at[idx, "partner_coverage"] = partner_coverage
        metrics_df.at[idx, "detections_per_read"] = detections_per_read
        metrics_df.at[idx, "detections_per_million"] = detections_per_million
        if "observed_variants" not in metrics_df.columns:
            metrics_df["observed_variants"] = 0
        metrics_df.at[idx, "observed_variants"] = observed_variants

    return metrics_df


def load_partner_counts(partner_files: list[str]) -> Any:
    """
    Load and merge partner counts across samples.

    Args:
        partner_files: List of paths to partner counts CSV files

    Returns:
        Merged partner counts DataFrame
    """
    import pandas as pd

    partner_dfs = []
    for partner_file in partner_files:
        sample = Path(partner_file).stem.replace(".partner_counts", "")
        df = pd.read_csv(partner_file)
        df = df.rename(
            columns={
                "partner_end_count": f"{sample}_end",
                "partner_linker_count": f"{sample}_linker",
            }
        )
        partner_dfs.append(df.set_index("partner_name"))

    if partner_dfs:
        partner_merged = pd.concat(partner_dfs, axis=1).fillna(0).astype(int)
        return partner_merged.reset_index()
    else:
        return pd.DataFrame(columns=["partner_name"])


def load_json_metrics(metrics_files: list[str], suffix: str = "fusion_metrics") -> Any:
    """
    Load per-sample metrics from JSON files.

    Args:
        metrics_files: List of paths to metrics JSON files
        suffix: Suffix to strip from filename

    Returns:
        DataFrame with per-sample metrics
    """
    import pandas as pd

    json_metrics = []
    for metrics_file in metrics_files:
        sample = Path(metrics_file).stem.replace(f".{suffix}", "")
        try:
            with open(metrics_file, "r") as f:
                data = json.load(f)
        except FileNotFoundError:
            data = {}
        json_metrics.append(
            {
                "sample": sample,
                "partner_end_reads": int(data.get("partner_end_reads", 0)),
                "partner_linker_reads": int(data.get("partner_linker_reads", 0)),
                "unique_partners_detected": int(
                    data.get("unique_partners_detected", 0)
                ),
                "unique_partner_linker_detected": int(
                    data.get("unique_partner_linker_detected", 0)
                ),
                "matched_reads": int(data.get("matched_reads", 0)),
                "reads_processed": int(data.get("reads_processed", 0)),
            }
        )

    return pd.DataFrame(json_metrics)


def calculate_sensitivity_metrics(
    sample_cols: list[str],
    merge_stats: dict,
    trim_logs: dict,
    json_df: Any,
    variant_catalog_path: str | Path,
    breakpoint_window: int,
    ihist_base_path: str | Path | None = None,
) -> Any:
    """
    Calculate sensitivity metrics for merged detection.

    Args:
        sample_cols: List of sample column names
        merge_stats: Dict mapping sample name to merge stats
        trim_logs: Dict mapping sample name to trim log stats
        json_df: DataFrame with per-sample detection metrics
        variant_catalog_path: Path to variant catalog CSV
        breakpoint_window: Breakpoint window size
        ihist_base_path: Base path for ihist files (e.g., stats/{exp}/merge/)

    Returns:
        DataFrame with sensitivity metrics
    """
    import pandas as pd

    variant_catalog_df = pd.read_csv(variant_catalog_path)
    fusion_lengths = variant_catalog_df.loc[
        variant_catalog_df["type"] == "fusion", "full_fusion_length"
    ].dropna()
    mean_fusion_length = (
        float(fusion_lengths.mean()) if not fusion_lengths.empty else 0.0
    )

    kmer_length = 2 * breakpoint_window

    def _calculate_expected_fraction(
        avg_insert: float, lengths: pd.Series, k: int
    ) -> float:
        """Calculate expected detection fraction per fusion, then average (per-fusion calibration)."""
        if avg_insert <= 0 or lengths.empty:
            return 0.0

        probs = pd.Series(0.0, index=lengths.index)

        full_mask = lengths <= avg_insert
        probs[full_mask] = 1.0

        partial_mask = ~full_mask & (avg_insert >= k)
        probs[partial_mask] = (avg_insert - k + 1) / (
            lengths[partial_mask] - avg_insert + 1
        )

        return float(probs.mean())

    sensitivity_rows = []
    for sample in sample_cols:
        merge = merge_stats.get(sample)
        trim_log = trim_logs.get(sample)

        avg_insert = merge["avg_insert"] if merge and merge.get("avg_insert") else 0.0
        merged_reads = merge["merged_reads"] if merge else 0
        raw_reads = trim_log["input_reads"] if trim_log else 0
        raw_pairs = raw_reads // 2 if raw_reads else 0

        expected_fraction = _calculate_expected_fraction(
            avg_insert, fusion_lengths, kmer_length
        )

        matched_reads = 0
        if not json_df.empty:
            sample_rows = json_df[json_df["sample"] == sample]
            if not sample_rows.empty:
                matched_reads = int(sample_rows["matched_reads"].fillna(0).iloc[0])

        expected_from_merged = merged_reads * expected_fraction
        expected_from_raw = raw_pairs * expected_fraction

        matching_efficiency = (
            (matched_reads / expected_from_merged) if expected_from_merged else 0.0
        )
        end_to_end_efficiency = (
            (matched_reads / expected_from_raw) if expected_from_raw else 0.0
        )

        overlap_median = 0.0
        if ihist_base_path:
            ihist_path = Path(ihist_base_path) / f"{sample}.ihist"
            if ihist_path.exists():
                ihist_info = parse_ihist(ihist_path)
                overlap_median = ihist_info["median"] if ihist_info else 0.0

        sensitivity_rows.append(
            {
                "sample": sample,
                "avg_insert_size": float(avg_insert),
                "overlap_median": float(overlap_median),
                "breakpoint_kmer_length": int(kmer_length),
                "mean_fusion_length": float(mean_fusion_length),
                "expected_detection_fraction": float(expected_fraction),
                "matching_efficiency": float(matching_efficiency),
                "end_to_end_efficiency": float(end_to_end_efficiency),
            }
        )

    return pd.DataFrame(sensitivity_rows)


def calculate_decay_metrics(
    sample_cols: list[str],
    trim_stats: dict,
    contam_stats: dict,
    quality_stats: dict,
    merge_stats: dict,
    trim_logs: dict,
    contam_logs: dict,
    quality_logs: dict,
    json_df: Any,
) -> Any:
    """
    Calculate decay metrics for read counts through pipeline steps.

    Args:
        sample_cols: List of sample column names
        trim_stats: Dict mapping sample name to trim stats
        contam_stats: Dict mapping sample name to contam stats
        quality_stats: Dict mapping sample name to quality stats
        merge_stats: Dict mapping sample name to merge stats
        trim_logs: Dict mapping sample name to trim log stats
        contam_logs: Dict mapping sample name to contam log stats
        quality_logs: Dict mapping sample name to quality log stats
        json_df: DataFrame with per-sample detection metrics

    Returns:
        DataFrame with decay metrics
    """
    import pandas as pd

    decay_rows = []
    for sample in sample_cols:
        trim = trim_stats.get(sample)
        contam = contam_stats.get(sample)
        quality = quality_stats.get(sample)
        merged = merge_stats.get(sample)

        trim_log = trim_logs.get(sample)
        contam_log = contam_logs.get(sample)
        quality_log = quality_logs.get(sample)

        raw_reads = (
            trim_log["input_reads"]
            if trim_log
            else (trim["input_reads"] if trim else 0)
        )
        raw_bases = (
            trim_log["input_bases"]
            if trim_log
            else (trim["input_bases"] if trim else None)
        )

        raw_pairs = raw_reads // 2 if raw_reads else 0

        matched_reads = 0
        if not json_df.empty and sample in json_df["sample"].values:
            sample_rows = json_df[json_df["sample"] == sample]
            if not sample_rows.empty:
                matched_reads = int(sample_rows["matched_reads"].fillna(0).iloc[0])

        steps = [
            ("raw", raw_pairs, raw_bases),
            (
                "trimmed",
                (
                    trim_log["output_reads"]
                    if trim_log
                    else (trim["output_reads"] if trim else 0)
                )
                // 2,
                (
                    trim_log["output_bases"]
                    if trim_log
                    else (trim["output_bases"] if trim else None)
                ),
            ),
            (
                "cleaned",
                (
                    contam_log["output_reads"]
                    if contam_log
                    else (contam["output_reads"] if contam else 0)
                )
                // 2,
                (
                    contam_log["output_bases"]
                    if contam_log
                    else (contam["output_bases"] if contam else None)
                ),
            ),
            (
                "quality",
                (
                    quality_log["output_reads"]
                    if quality_log
                    else (quality["output_reads"] if quality else 0)
                )
                // 2,
                (
                    quality_log["output_bases"]
                    if quality_log
                    else (quality["output_bases"] if quality else None)
                ),
            ),
            (
                "merged",
                merged["merged_reads"] if merged else 0,
                merged["merged_bases"] if merged else None,
            ),
            ("matched", matched_reads, None),
        ]

        for step, pairs, bases in steps:
            read_fraction = float(pairs / raw_pairs) if raw_pairs else 0.0
            decay_rows.append(
                {
                    "sample": sample,
                    "step": step,
                    "reads": int(pairs) if pairs is not None else 0,
                    "bases": int(bases) if bases is not None else None,
                    "read_fraction": read_fraction,
                    "base_fraction": float(bases / raw_bases)
                    if raw_bases and bases is not None
                    else 0.0,
                }
            )

    return pd.DataFrame(decay_rows)


def aggregate_merged(
    counts: list[str],
    metrics: list[str],
    partner_counts: list[str],
    variant_catalog: str,
    merge_logs: list[str],
    trim_logs: list[str],
    contam_logs: list[str],
    quality_logs: list[str],
    trim_stats: list[str],
    contam_stats: list[str],
    quality_stats: list[str],
    breakpoint_window: int,
    output_summary: str,
    output_qc: str,
    output_partner: str,
    output_sensitivity: str,
    output_decay: str,
    output_trim: str,
    output_contam: str,
    output_quality: str,
    ihist_base_path: str | None = None,
) -> None:
    """
    Aggregate merged counts and generate all output files.
    """
    import pandas as pd

    merged, sample_cols = load_and_merge_counts(counts, "fusion_counts")

    merged.to_csv(output_summary, index=False)

    merge_stats = {}
    for log_path in merge_logs:
        sample = Path(log_path).stem
        merge_stats[sample] = parse_bbmerge_log(log_path)

    json_df = load_json_metrics(metrics, "fusion_metrics")

    metrics_df = calculate_qc_metrics(merged, sample_cols, merge_stats, json_df)

    try:
        variant_catalog_df = pd.read_csv(variant_catalog)
        fusion_catalog = variant_catalog_df[
            variant_catalog_df["type"] == "fusion"
        ].copy()
        if "breakpoint_nt" in fusion_catalog.columns:
            expected_breakpoints = fusion_catalog["breakpoint_nt"].nunique()
        else:
            expected_breakpoints = len(fusion_catalog)
        if "partner_name" in fusion_catalog.columns:
            expected_partners = fusion_catalog["partner_name"].nunique()
        else:
            expected_partners = 0
    except (FileNotFoundError, KeyError, pd.errors.EmptyDataError):
        fusion_data = merged[merged["type"] == "fusion"].copy()
        expected_breakpoints = len(fusion_data)
        expected_partners = 0
        fusion_catalog = None

    metrics_df = calculate_diversity_metrics(
        metrics_df, merged, fusion_catalog, expected_breakpoints, expected_partners
    )

    metrics_df.to_csv(output_qc, index=False)

    partner_merged = load_partner_counts(partner_counts)
    partner_merged.to_csv(output_partner, index=False)

    sensitivity_trim_logs = {}
    for path in trim_logs:
        sample = Path(path).stem.replace(".trim", "")
        sensitivity_trim_logs[sample] = parse_bbduk_log(path)

    sensitivity_merge_stats = {}
    for path in merge_logs:
        sample = Path(path).stem
        sensitivity_merge_stats[sample] = parse_bbmerge_stats(path)

    sensitivity_df = calculate_sensitivity_metrics(
        sample_cols=sample_cols,
        merge_stats=sensitivity_merge_stats,
        trim_logs=sensitivity_trim_logs,
        json_df=json_df,
        variant_catalog_path=variant_catalog,
        breakpoint_window=breakpoint_window,
        ihist_base_path=ihist_base_path,
    )
    sensitivity_df.to_csv(output_sensitivity, index=False)

    trim_stats_dict = {}
    for path in trim_stats:
        sample = Path(path).stem.replace(".trim", "")
        trim_stats_dict[sample] = parse_bbduk_stats(path)

    contam_stats_dict = {}
    for path in contam_stats:
        sample = Path(path).stem.replace(".contam", "")
        contam_stats_dict[sample] = parse_bbduk_stats(path)

    quality_stats_dict = {}
    for path in quality_stats:
        sample = Path(path).stem.replace(".quality", "")
        quality_stats_dict[sample] = parse_bbduk_stats(path)

    contam_logs_dict = {}
    for path in contam_logs:
        sample = Path(path).stem.replace(".clean", "")
        contam_logs_dict[sample] = parse_bbduk_log(path)

    quality_logs_dict = {}
    for path in quality_logs:
        sample = Path(path).stem.replace(".quality", "")
        quality_logs_dict[sample] = parse_bbduk_log(path)

    decay_df = calculate_decay_metrics(
        sample_cols=sample_cols,
        trim_stats=trim_stats_dict,
        contam_stats=contam_stats_dict,
        quality_stats=quality_stats_dict,
        merge_stats=sensitivity_merge_stats,
        trim_logs=sensitivity_trim_logs,
        contam_logs=contam_logs_dict,
        quality_logs=quality_logs_dict,
        json_df=json_df,
    )
    decay_df.to_csv(output_decay, index=False)

    for step_name, out_path in [
        ("trimmed", output_trim),
        ("cleaned", output_contam),
        ("quality", output_quality),
    ]:
        step_df = decay_df[decay_df["step"] == step_name].copy()
        step_df = step_df.drop(columns=["step"], errors="ignore")
        step_df.to_csv(out_path, index=False)


def aggregate_unmerged(
    counts: list[str],
    metrics: list[str],
    partner_counts: list[str],
    output_summary: str,
    output_qc: str,
    output_partner: str,
) -> None:
    """
    Aggregate unmerged counts and generate output files.
    """
    merged, sample_cols = load_and_merge_counts(counts, "unmerged_fusion_counts")

    merged.to_csv(output_summary, index=False)

    metrics_df = calculate_qc_metrics(merged, sample_cols)

    json_df = load_json_metrics(metrics, "unmerged_fusion_metrics")

    if not json_df.empty:
        metrics_df = metrics_df.merge(json_df, on="sample", how="left")

    metrics_df.to_csv(output_qc, index=False)

    partner_merged = load_partner_counts(partner_counts)
    partner_merged.to_csv(output_partner, index=False)


def main_snakemake(snakemake) -> None:
    """Entry point when called from Snakemake."""
    mode = snakemake.params.get("mode", "merged")
    breakpoint_window = snakemake.params.get("breakpoint_window", 12)
    ihist_base = snakemake.params.get("ihist_base_path", None)

    if mode == "merged":
        aggregate_merged(
            counts=snakemake.input.counts,
            metrics=snakemake.input.metrics,
            partner_counts=snakemake.input.partner_counts,
            variant_catalog=snakemake.input.variant_catalog,
            merge_logs=snakemake.input.merge_logs,
            trim_logs=snakemake.input.trim_logs,
            contam_logs=snakemake.input.contam_logs,
            quality_logs=snakemake.input.quality_logs,
            trim_stats=snakemake.input.trim_stats,
            contam_stats=snakemake.input.contam_stats,
            quality_stats=snakemake.input.quality_stats,
            breakpoint_window=breakpoint_window,
            output_summary=snakemake.output.summary,
            output_qc=snakemake.output.qc_metrics,
            output_partner=snakemake.output.partner_summary,
            output_sensitivity=snakemake.output.sensitivity_metrics,
            output_decay=snakemake.output.decay_metrics,
            output_trim=snakemake.output.trim_metrics,
            output_contam=snakemake.output.contam_metrics,
            output_quality=snakemake.output.quality_metrics,
            ihist_base_path=ihist_base,
        )
    else:
        aggregate_unmerged(
            counts=snakemake.input.counts,
            metrics=snakemake.input.metrics,
            partner_counts=snakemake.input.partner_counts,
            output_summary=snakemake.output.summary,
            output_qc=snakemake.output.qc_metrics,
            output_partner=snakemake.output.partner_summary,
        )


def main_cli() -> None:
    """Entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Aggregate fusion counts across samples",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument("mode", choices=["merged", "unmerged"], help="Aggregation mode")
    parser.add_argument("--counts", nargs="+", required=True, help="Count CSV files")
    parser.add_argument(
        "--metrics", nargs="+", required=True, help="Metrics JSON files"
    )
    parser.add_argument(
        "--partner-counts", nargs="+", required=True, help="Partner counts CSV files"
    )
    parser.add_argument("--variant-catalog", help="Variant catalog CSV (merged mode)")
    parser.add_argument(
        "--merge-logs", nargs="+", help="BBMerge log files (merged mode)"
    )
    parser.add_argument(
        "--trim-logs", nargs="+", help="BBDuk trim log files (merged mode)"
    )
    parser.add_argument(
        "--contam-logs", nargs="+", help="BBDuk contam log files (merged mode)"
    )
    parser.add_argument(
        "--quality-logs", nargs="+", help="BBDuk quality log files (merged mode)"
    )
    parser.add_argument(
        "--trim-stats", nargs="+", help="BBDuk trim stats files (merged mode)"
    )
    parser.add_argument(
        "--contam-stats", nargs="+", help="BBDuk contam stats files (merged mode)"
    )
    parser.add_argument(
        "--quality-stats", nargs="+", help="BBDuk quality stats files (merged mode)"
    )
    parser.add_argument(
        "--breakpoint-window", type=int, default=12, help="Breakpoint window size"
    )
    parser.add_argument(
        "--ihist-base-path", help="Base path for ihist files (merged mode)"
    )
    parser.add_argument("--output-summary", required=True, help="Output summary CSV")
    parser.add_argument("--output-qc", required=True, help="Output QC metrics CSV")
    parser.add_argument(
        "--output-partner", required=True, help="Output partner summary CSV"
    )
    parser.add_argument(
        "--output-sensitivity", help="Output sensitivity metrics CSV (merged mode)"
    )
    parser.add_argument("--output-decay", help="Output decay metrics CSV (merged mode)")
    parser.add_argument("--output-trim", help="Output trim metrics CSV (merged mode)")
    parser.add_argument(
        "--output-contam", help="Output contam metrics CSV (merged mode)"
    )
    parser.add_argument(
        "--output-quality", help="Output quality metrics CSV (merged mode)"
    )

    args = parser.parse_args()

    if args.mode == "merged":
        if not all(
            [
                args.variant_catalog,
                args.merge_logs,
                args.trim_logs,
                args.output_sensitivity,
                args.output_decay,
            ]
        ):
            parser.error("merged mode requires additional arguments")
        aggregate_merged(
            counts=args.counts,
            metrics=args.metrics,
            partner_counts=args.partner_counts,
            variant_catalog=args.variant_catalog,
            merge_logs=args.merge_logs,
            trim_logs=args.trim_logs,
            contam_logs=args.contam_logs or [],
            quality_logs=args.quality_logs or [],
            trim_stats=args.trim_stats or [],
            contam_stats=args.contam_stats or [],
            quality_stats=args.quality_stats or [],
            breakpoint_window=args.breakpoint_window,
            output_summary=args.output_summary,
            output_qc=args.output_qc,
            output_partner=args.output_partner,
            output_sensitivity=args.output_sensitivity,
            output_decay=args.output_decay,
            output_trim=args.output_trim,
            output_contam=args.output_contam,
            output_quality=args.output_quality,
            ihist_base_path=args.ihist_base_path,
        )
    else:
        aggregate_unmerged(
            counts=args.counts,
            metrics=args.metrics,
            partner_counts=args.partner_counts,
            output_summary=args.output_summary,
            output_qc=args.output_qc,
            output_partner=args.output_partner,
        )


if __name__ == "__main__":
    main_cli()
