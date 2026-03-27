#!/usr/bin/env python3
"""
Unit tests for aggregate_counts.py

Tests the aggregation logic for fusion counts and QC metrics.
"""

import json
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))

from aggregate_counts import (
    load_and_merge_counts,
    calculate_qc_metrics,
    calculate_diversity_metrics,
    load_partner_counts,
    load_json_metrics,
    aggregate_merged,
    aggregate_unmerged,
)


class TestLoadAndMergeCounts:
    """Tests for the count loading and merging logic."""

    def test_loads_and_merges_multiple_count_files(self, tmp_path):
        """Should load multiple count files and merge on fusion_id."""
        count1 = tmp_path / "sample1.fusion_counts.csv"
        count1.write_text(
            "fusion_id,type,count\nTPR_1_Met_WT,fusion,10\nTPR_2_Met_WT,fusion,20\n"
        )

        count2 = tmp_path / "sample2.fusion_counts.csv"
        count2.write_text(
            "fusion_id,type,count\nTPR_1_Met_WT,fusion,15\nTPR_3_Met_WT,fusion,30\n"
        )

        merged, sample_cols = load_and_merge_counts(
            [str(count1), str(count2)], "fusion_counts"
        )

        assert set(sample_cols) == {"sample1", "sample2"}
        assert len(merged) == 3
        assert (
            merged.loc[merged["fusion_id"] == "TPR_1_Met_WT", "sample1"].iloc[0] == 10
        )
        assert (
            merged.loc[merged["fusion_id"] == "TPR_1_Met_WT", "sample2"].iloc[0] == 15
        )
        assert merged.loc[merged["fusion_id"] == "TPR_3_Met_WT", "sample1"].iloc[0] == 0
        assert (
            merged.loc[merged["fusion_id"] == "TPR_3_Met_WT", "sample2"].iloc[0] == 30
        )

    def test_handles_missing_type_column(self, tmp_path):
        """Should add default type='fusion' when type column is missing."""
        count1 = tmp_path / "sample1.fusion_counts.csv"
        count1.write_text("fusion_id,count\nTPR_1_Met_WT,10\n")

        merged, sample_cols = load_and_merge_counts([str(count1)], "fusion_counts")

        assert "type" in merged.columns
        assert all(merged["type"] == "fusion")

    def test_sorts_by_total_counts(self, tmp_path):
        """Should sort merged data by total counts descending."""
        count1 = tmp_path / "sample1.fusion_counts.csv"
        count1.write_text(
            "fusion_id,type,count\nTPR_low,fusion,5\nTPR_high,fusion,100\n"
        )

        count2 = tmp_path / "sample2.fusion_counts.csv"
        count2.write_text(
            "fusion_id,type,count\nTPR_low,fusion,10\nTPR_high,fusion,50\n"
        )

        merged, _ = load_and_merge_counts([str(count1), str(count2)], "fusion_counts")

        assert merged.iloc[0]["fusion_id"] == "TPR_high"
        assert merged.iloc[1]["fusion_id"] == "TPR_low"

    def test_unmerged_suffix_handling(self, tmp_path):
        """Should correctly extract sample names from unmerged count files."""
        count1 = tmp_path / "sample1.R1.unmerged_fusion_counts.csv"
        count1.write_text("fusion_id,type,count\nTPR_1_Met_WT,fusion,10\n")

        count2 = tmp_path / "sample1.R2.unmerged_fusion_counts.csv"
        count2.write_text("fusion_id,type,count\nTPR_1_Met_WT,fusion,20\n")

        merged, sample_cols = load_and_merge_counts(
            [str(count1), str(count2)], "unmerged_fusion_counts"
        )

        assert "sample1.R1" in sample_cols
        assert "sample1.R2" in sample_cols


class TestCalculateQcMetrics:
    """Tests for QC metrics calculation."""

    def test_calculates_basic_metrics(self, tmp_path):
        """Should calculate basic count metrics per sample."""
        count1 = tmp_path / "sample1.fusion_counts.csv"
        count1.write_text(
            "fusion_id,type,count\nTPR_1,fusion,100\nTPR_2,fusion,50\nTPR_3,fusion,0\n"
        )

        count2 = tmp_path / "sample2.fusion_counts.csv"
        count2.write_text(
            "fusion_id,type,count\nTPR_1,fusion,200\nTPR_2,fusion,0\nTPR_3,fusion,0\n"
        )

        merged, sample_cols = load_and_merge_counts(
            [str(count1), str(count2)], "fusion_counts"
        )

        metrics_df = calculate_qc_metrics(merged, sample_cols)

        sample1_metrics = metrics_df[metrics_df["sample"] == "sample1"].iloc[0]
        assert sample1_metrics["total_counts"] == 150
        assert sample1_metrics["total_fusion_counts"] == 150
        assert sample1_metrics["unique_fusions"] == 2

    def test_calculates_zero_fractions(self, tmp_path):
        """Should calculate zero fraction metrics correctly."""
        count1 = tmp_path / "sample1.fusion_counts.csv"
        count1.write_text(
            "fusion_id,type,count\nTPR_1,fusion,100\nTPR_2,fusion,0\nTPR_3,fusion,0\nTPR_4,fusion,0\n"
        )

        merged, sample_cols = load_and_merge_counts([str(count1)], "fusion_counts")

        metrics_df = calculate_qc_metrics(merged, sample_cols)

        sample1_metrics = metrics_df[metrics_df["sample"] == "sample1"].iloc[0]
        assert sample1_metrics["zero_fraction"] == 0.75
        assert sample1_metrics["fusion_zero_fraction"] == 0.75

    def test_includes_merge_stats_when_provided(self, tmp_path):
        """Should include merge statistics when provided."""
        count1 = tmp_path / "sample1.fusion_counts.csv"
        count1.write_text("fusion_id,type,count\nTPR_1,fusion,100\n")

        merged, sample_cols = load_and_merge_counts([str(count1)], "fusion_counts")

        merge_stats = {
            "sample1": {
                "total_pairs": 1000,
                "merged_pairs": 800,
                "unmerged_pairs": 200,
                "unmerged_fraction": 0.2,
            }
        }

        metrics_df = calculate_qc_metrics(merged, sample_cols, merge_stats)

        sample1_metrics = metrics_df[metrics_df["sample"] == "sample1"].iloc[0]
        assert sample1_metrics["total_pairs"] == 1000
        assert sample1_metrics["merged_pairs"] == 800
        assert sample1_metrics["unmerged_fraction"] == 0.2

    def test_includes_json_metrics_when_provided(self, tmp_path):
        """Should merge in JSON detection metrics when provided."""
        count1 = tmp_path / "sample1.fusion_counts.csv"
        count1.write_text("fusion_id,type,count\nTPR_1,fusion,100\n")

        merged, sample_cols = load_and_merge_counts([str(count1)], "fusion_counts")

        metrics_json = tmp_path / "sample1.fusion_metrics.json"
        metrics_json.write_text(
            json.dumps(
                {
                    "partner_end_reads": 500,
                    "matched_reads": 150,
                    "reads_processed": 1000,
                }
            )
        )

        json_df = load_json_metrics([str(metrics_json)], "fusion_metrics")
        metrics_df = calculate_qc_metrics(merged, sample_cols, json_df=json_df)

        sample1_metrics = metrics_df[metrics_df["sample"] == "sample1"].iloc[0]
        assert sample1_metrics["partner_end_reads"] == 500
        assert sample1_metrics["matched_reads"] == 150


class TestLoadPartnerCounts:
    """Tests for partner counts loading."""

    def test_loads_and_merges_partner_counts(self, tmp_path):
        """Should load and merge partner counts across samples."""
        partner1 = tmp_path / "sample1.partner_counts.csv"
        partner1.write_text(
            "partner_name,partner_end_count,partner_linker_count\nTPR,100,50\nCCDC6,80,40\n"
        )

        partner2 = tmp_path / "sample2.partner_counts.csv"
        partner2.write_text(
            "partner_name,partner_end_count,partner_linker_count\nTPR,150,75\n"
        )

        merged = load_partner_counts([str(partner1), str(partner2)])

        assert "partner_name" in merged.columns
        assert "sample1_end" in merged.columns
        assert "sample2_end" in merged.columns
        assert len(merged) == 2

        tpr_row = merged[merged["partner_name"] == "TPR"].iloc[0]
        assert tpr_row["sample1_end"] == 100
        assert tpr_row["sample2_end"] == 150

    def test_handles_empty_partner_counts(self, tmp_path):
        """Should return empty DataFrame with correct columns when no files."""
        merged = load_partner_counts([])
        assert "partner_name" in merged.columns
        assert len(merged) == 0


class TestLoadJsonMetrics:
    """Tests for JSON metrics loading."""

    def test_loads_json_metrics(self, tmp_path):
        """Should load metrics from JSON files."""
        metrics1 = tmp_path / "sample1.fusion_metrics.json"
        metrics1.write_text(
            json.dumps(
                {
                    "partner_end_reads": 100,
                    "partner_linker_reads": 50,
                    "unique_partners_detected": 5,
                    "matched_reads": 80,
                    "reads_processed": 1000,
                }
            )
        )

        metrics2 = tmp_path / "sample2.fusion_metrics.json"
        metrics2.write_text(
            json.dumps(
                {
                    "partner_end_reads": 200,
                    "matched_reads": 160,
                    "reads_processed": 2000,
                }
            )
        )

        df = load_json_metrics([str(metrics1), str(metrics2)], "fusion_metrics")

        assert len(df) == 2
        assert "sample" in df.columns
        sample1 = df[df["sample"] == "sample1"].iloc[0]
        assert sample1["partner_end_reads"] == 100
        assert sample1["unique_partners_detected"] == 5

    def test_handles_missing_files_gracefully(self, tmp_path):
        """Should handle missing JSON files gracefully."""
        missing = tmp_path / "missing.json"

        df = load_json_metrics([str(missing)], "fusion_metrics")

        assert len(df) == 1
        assert df.iloc[0]["partner_end_reads"] == 0


class TestAggregateUnmerged:
    """Tests for unmerged aggregation."""

    def test_aggregate_unmerged_creates_outputs(self, tmp_path):
        """Should create all expected output files."""
        count1 = tmp_path / "sample1.R1.unmerged_fusion_counts.csv"
        count1.write_text("fusion_id,type,count\nTPR_1_Met_WT,fusion,10\n")

        metrics1 = tmp_path / "sample1.R1.unmerged_fusion_metrics.json"
        metrics1.write_text(json.dumps({"matched_reads": 100}))

        partner1 = tmp_path / "sample1.R1.unmerged_partner_counts.csv"
        partner1.write_text(
            "partner_name,partner_end_count,partner_linker_count\nTPR,50,25\n"
        )

        summary_out = tmp_path / "summary.csv"
        qc_out = tmp_path / "qc.csv"
        partner_out = tmp_path / "partner.csv"

        aggregate_unmerged(
            counts=[str(count1)],
            metrics=[str(metrics1)],
            partner_counts=[str(partner1)],
            output_summary=str(summary_out),
            output_qc=str(qc_out),
            output_partner=str(partner_out),
        )

        assert summary_out.exists()
        assert qc_out.exists()
        assert partner_out.exists()

        summary_df = __import__("pandas").read_csv(summary_out)
        assert len(summary_df) == 1
        assert summary_df.iloc[0]["fusion_id"] == "TPR_1_Met_WT"


class TestLogParsing:
    """Tests for log parsing functions."""

    def test_parse_bbmerge_log(self, tmp_path):
        """Should parse BBMerge log files correctly."""
        from utils import parse_bbmerge_log

        log_file = tmp_path / "bbmerge.log"
        log_file.write_text("""
Pairs: 1,000,000
Joined: 850,000
Avg Insert: 125.5
""")

        result = parse_bbmerge_log(log_file)

        assert result is not None
        assert result["total_pairs"] == 1000000
        assert result["merged_pairs"] == 850000
        assert result["unmerged_pairs"] == 150000
        assert result["unmerged_fraction"] == 0.15

    def test_parse_bbmerge_log_missing_file(self, tmp_path):
        """Should return None for missing files."""
        from utils import parse_bbmerge_log

        result = parse_bbmerge_log(tmp_path / "missing.log")
        assert result is None

    def test_parse_bbduk_log(self, tmp_path):
        """Should parse BBDuk log files correctly."""
        from utils import parse_bbduk_log

        log_file = tmp_path / "bbduk.log"
        log_file.write_text(
            "Input: 1,000,000 reads (150,000,000 bases)\nResult: 980,000 reads (147,000,000 bases)\n"
        )

        result = parse_bbduk_log(log_file)

        assert result is not None
        assert result["input_reads"] == 1000000
        assert result["output_reads"] == 980000

    def test_parse_bbduk_stats(self, tmp_path):
        """Should parse BBDuk stats files correctly."""
        from utils import parse_bbduk_stats

        stats_file = tmp_path / "stats.txt"
        stats_file.write_text(
            "Input: 500,000 reads (75,000,000 bases)\nOutput: 490,000 reads (73,500,000 bases)\n"
        )

        result = parse_bbduk_stats(stats_file)

        assert result is not None
        assert result["input_reads"] == 500000
        assert result["output_reads"] == 490000

    def test_parse_ihist(self, tmp_path):
        """Should parse ihist files correctly."""
        from utils import parse_ihist

        ihist_file = tmp_path / "merge.ihist"
        ihist_file.write_text("""# ihist version 1.0
# Created by BBMerge
# Non-overlapping: true
10\t100
20\t200
30\t300
40\t400
50\t500
""")

        result = parse_ihist(ihist_file)

        assert result is not None
        assert "median" in result
        assert result["median"] == 40.0

    def test_parse_ihist_missing_file(self, tmp_path):
        """Should return None for missing ihist files."""
        from utils import parse_ihist

        result = parse_ihist(tmp_path / "missing.ihist")
        assert result is None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
