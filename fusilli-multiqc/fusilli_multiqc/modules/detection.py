"""
MultiQC module for FUSILLI detection metrics.

This module visualizes:
- Detection efficiency (detection_efficiency, prefilter_efficiency, matching_efficiency)
- Library coverage (variant, breakpoint, partner coverage)
- Sensitivity analysis (sensitivity_index vs expected_detection_fraction)
"""

import logging
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseModule
from multiqc.plots import linegraph, bargraph, scatter

from fusilli_multiqc.utils import parse_csv_file

logger = logging.getLogger(__name__)


class MultiqcModule(BaseModule):
    """
    FUSILLI Detection Metrics Module

    Parses fusion_qc_metrics.csv, sensitivity_metrics.csv, and fusion_counts_summary.csv
    to generate detection efficiency, coverage, and sensitivity visualizations.
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="FUSILLI Detection Metrics",
            anchor="fusilli_detection",
            info="Fusion detection efficiency, library coverage, and sensitivity metrics.",
        )

        # Find and parse input files
        self.fusion_qc_data = None
        self.sensitivity_data = None
        self.fusion_counts_data = None

        # Find fusion_qc_metrics.csv
        for f in self.find_log_files("fusion_qc_metrics.csv"):
            self.fusion_qc_data = parse_csv_file(f["fn"])
            if self.fusion_qc_data is not None:
                break

        # Find sensitivity_metrics.csv
        for f in self.find_log_files("sensitivity_metrics.csv"):
            self.sensitivity_data = parse_csv_file(f["fn"])
            if self.sensitivity_data is not None:
                break

        # Find fusion_counts_summary.csv
        for f in self.find_log_files("fusion_counts_summary.csv"):
            self.fusion_counts_data = parse_csv_file(f["fn"])
            if self.fusion_counts_data is not None:
                break

        # If no data found, exit
        if self.fusion_qc_data is None and self.sensitivity_data is None:
            raise UserWarning

        # Calculate detection efficiency metrics
        if self.fusion_qc_data is not None:
            self.calculate_detection_efficiency()
            self.calculate_coverage_metrics()

        # Generate plots
        if self.fusion_qc_data is not None:
            self.detection_efficiency_plot()
            self.coverage_plot()

        if self.sensitivity_data is not None:
            self.sensitivity_plot()

        # Add summary table
        if self.fusion_qc_data is not None:
            self.add_summary_table()

    def calculate_detection_efficiency(self):
        """Calculate detection efficiency metrics from fusion_qc_metrics."""
        if self.fusion_qc_data is None:
            return

        # Calculate detection efficiency if not present
        if "detection_efficiency" not in self.fusion_qc_data.columns:
            if "matched_reads" in self.fusion_qc_data.columns and "reads_processed" in self.fusion_qc_data.columns:
                self.fusion_qc_data["detection_efficiency"] = (
                    self.fusion_qc_data["matched_reads"] / self.fusion_qc_data["reads_processed"]
                ).fillna(0.0)

        # Calculate prefilter efficiency
        if "prefilter_efficiency" not in self.fusion_qc_data.columns:
            if "partner_end_reads" in self.fusion_qc_data.columns and "reads_processed" in self.fusion_qc_data.columns:
                self.fusion_qc_data["prefilter_efficiency"] = (
                    self.fusion_qc_data["partner_end_reads"] / self.fusion_qc_data["reads_processed"]
                ).fillna(0.0)

        # Calculate matching efficiency
        if "matching_efficiency" not in self.fusion_qc_data.columns:
            if "matched_reads" in self.fusion_qc_data.columns and "partner_end_reads" in self.fusion_qc_data.columns:
                self.fusion_qc_data["matching_efficiency"] = (
                    self.fusion_qc_data["matched_reads"] / self.fusion_qc_data["partner_end_reads"]
                ).fillna(0.0)

    def calculate_coverage_metrics(self):
        """Calculate library coverage metrics."""
        if self.fusion_qc_data is None:
            return

        # Calculate variant coverage
        if "variant_coverage" not in self.fusion_qc_data.columns:
            if "unique_fusions" in self.fusion_qc_data.columns and "expected_fusions" in self.fusion_qc_data.columns:
                self.fusion_qc_data["variant_coverage"] = (
                    self.fusion_qc_data["unique_fusions"] / self.fusion_qc_data["expected_fusions"]
                ).fillna(0.0)

        # Calculate breakpoint coverage (simplified - using unique_fusions as proxy)
        if "breakpoint_coverage" not in self.fusion_qc_data.columns:
            if "unique_fusions" in self.fusion_qc_data.columns and "expected_fusions" in self.fusion_qc_data.columns:
                # For now, use variant coverage as breakpoint coverage
                # In future, could parse fusion_id to extract unique breakpoints
                self.fusion_qc_data["breakpoint_coverage"] = self.fusion_qc_data["variant_coverage"]

        # Calculate partner coverage
        if "partner_coverage" not in self.fusion_qc_data.columns:
            if "unique_partners_detected" in self.fusion_qc_data.columns:
                # Assume we know expected partners from config, or use a reasonable estimate
                # For now, use unique_partners_detected as a proxy
                self.fusion_qc_data["partner_coverage"] = (
                    self.fusion_qc_data["unique_partners_detected"] / 100.0
                ).clip(upper=1.0)

    def detection_efficiency_plot(self):
        """Create detection efficiency multi-line plot."""
        if self.fusion_qc_data is None or "sample" not in self.fusion_qc_data.columns:
            return

        plot_data = OrderedDict()
        for _, row in self.fusion_qc_data.iterrows():
            sample = row["sample"]
            plot_data[sample] = {}
            if "detection_efficiency" in row:
                plot_data[sample]["Detection Efficiency"] = row["detection_efficiency"]
            if "prefilter_efficiency" in row:
                plot_data[sample]["Prefilter Efficiency"] = row["prefilter_efficiency"]
            if "matching_efficiency" in row:
                plot_data[sample]["Matching Efficiency"] = row["matching_efficiency"]

        if not plot_data:
            return

        pconfig = {
            "id": "fusilli_detection_efficiency",
            "title": "FUSILLI: Detection Efficiency Metrics",
            "ylab": "Efficiency",
            "xlab": "Sample",
            "tt_label": "<b>{point.x}</b>: {point.y:.3f}",
            "data_labels": [
                {"name": "Detection Efficiency", "ylab": "Efficiency"},
            ],
        }

        self.add_section(
            name="Detection Efficiency",
            anchor="fusilli_detection_efficiency",
            description="Detection efficiency metrics across samples.",
            helptext="""
            This plot shows three efficiency metrics:
            - **Detection Efficiency**: Fraction of processed reads that matched fusion sequences
            - **Prefilter Efficiency**: Fraction of reads containing partner domain 3' ends
            - **Matching Efficiency**: Fraction of prefiltered reads that matched specific breakpoints
            """,
            plot=linegraph.plot(plot_data, pconfig),
        )

    def coverage_plot(self):
        """Create library coverage stacked bar chart."""
        if self.fusion_qc_data is None or "sample" not in self.fusion_qc_data.columns:
            return

        plot_data = OrderedDict()
        for _, row in self.fusion_qc_data.iterrows():
            sample = row["sample"]
            plot_data[sample] = {}
            if "variant_coverage" in row:
                plot_data[sample]["Variant Coverage"] = row["variant_coverage"]
            if "breakpoint_coverage" in row:
                plot_data[sample]["Breakpoint Coverage"] = row["breakpoint_coverage"]
            if "partner_coverage" in row:
                plot_data[sample]["Partner Coverage"] = row["partner_coverage"]

        if not plot_data:
            return

        pconfig = {
            "id": "fusilli_library_coverage",
            "title": "FUSILLI: Library Coverage Metrics",
            "ylab": "Coverage Fraction",
            "xlab": "Sample",
            "stacking": None,
            "tt_label": "<b>{point.x}</b>: {point.y:.3f}",
        }

        self.add_section(
            name="Library Coverage",
            anchor="fusilli_library_coverage",
            description="Library representation coverage across samples.",
            helptext="""
            This plot shows coverage metrics for different aspects of the fusion library:
            - **Variant Coverage**: Fraction of expected fusion variants detected
            - **Breakpoint Coverage**: Fraction of expected breakpoint positions detected
            - **Partner Coverage**: Fraction of expected fusion partners detected
            """,
            plot=bargraph.plot(plot_data, pconfig),
        )

    def sensitivity_plot(self):
        """Create sensitivity analysis scatter plot."""
        if self.sensitivity_data is None or "sample" not in self.sensitivity_data.columns:
            return

        plot_data = OrderedDict()
        for _, row in self.sensitivity_data.iterrows():
            sample = row["sample"]
            plot_data[sample] = {
                "x": row.get("expected_detection_fraction", 0.0),
                "y": row.get("sensitivity_index", 0.0),
            }

        if not plot_data:
            return

        pconfig = {
            "id": "fusilli_sensitivity",
            "title": "FUSILLI: Sensitivity Analysis",
            "ylab": "Sensitivity Index",
            "xlab": "Expected Detection Fraction",
            "tt_label": "<b>{point.sample}</b><br>Expected: {point.x:.3f}<br>Sensitivity: {point.y:.3f}",
        }

        self.add_section(
            name="Sensitivity Analysis",
            anchor="fusilli_sensitivity",
            description="Sensitivity index vs expected detection fraction.",
            helptext="""
            This plot shows the relationship between:
            - **Expected Detection Fraction**: Fraction of reads expected to be long enough to detect breakpoints
            - **Sensitivity Index**: Ratio of actual detections to expected detections

            Values closer to 1.0 indicate good sensitivity.
            """,
            plot=scatter.plot(plot_data, pconfig),
        )

    def add_summary_table(self):
        """Add detection yield summary table."""
        if self.fusion_qc_data is None or "sample" not in self.fusion_qc_data.columns:
            return

        # Prepare table data
        headers = OrderedDict()
        headers["reads_processed"] = {
            "title": "Reads Processed",
            "description": "Total reads processed for fusion detection",
            "format": "{:,.0f}",
        }
        headers["matched_reads"] = {
            "title": "Matched Reads",
            "description": "Reads that matched fusion sequences",
            "format": "{:,.0f}",
        }

        if "detections_per_read" in self.fusion_qc_data.columns:
            headers["detections_per_read"] = {
                "title": "Detections/Read",
                "description": "Number of detections per read",
                "format": "{:.4f}",
            }

        if "detections_per_million" in self.fusion_qc_data.columns:
            headers["detections_per_million"] = {
                "title": "Detections/Million",
                "description": "Number of detections per million reads",
                "format": "{:.2f}",
            }

        if "detection_efficiency" in self.fusion_qc_data.columns:
            headers["detection_efficiency"] = {
                "title": "Detection Efficiency",
                "description": "Fraction of reads that matched fusions",
                "format": "{:.3f}",
                "max": 1.0,
                "min": 0.0,
            }

        self.general_stats_addcols(self.fusion_qc_data, headers)
