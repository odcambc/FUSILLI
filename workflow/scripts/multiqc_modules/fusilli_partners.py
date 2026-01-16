"""
MultiQC module for FUSILLI partner detection metrics.

This module visualizes:
- Partner detection heatmap (samples × partners)
- Partner coverage plots
- Partner end vs linker detection comparison
"""

import logging
import numpy as np
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseModule
from multiqc.plots import heatmap, bargraph, scatter

from .fusilli_base import parse_csv_file

logger = logging.getLogger(__name__)


class MultiqcModule(BaseModule):
    """
    FUSILLI Partner Detection Module
    
    Parses partner_counts_summary.csv to visualize partner detection across samples.
    """
    
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="FUSILLI Partner Detection",
            anchor="fusilli_partners",
            info="Fusion partner detection metrics and coverage across samples.",
        )
        
        # Find and parse input files
        self.partner_data = None
        
        # Find partner_counts_summary.csv
        for f in self.find_log_files("partner_counts_summary.csv"):
            self.partner_data = parse_csv_file(f["fn"])
            if self.partner_data is not None:
                break
        
        # If no data found, exit
        if self.partner_data is None:
            raise UserWarning
        
        # Process partner data
        self.processed_data = self.process_partner_data()
        
        # Generate plots
        if self.processed_data is not None:
            self.partner_detection_heatmap()
            self.partner_coverage_plot()
            self.partner_end_vs_linker_plot()
    
    def process_partner_data(self):
        """
        Process partner counts data into a format suitable for plotting.
        
        Returns:
            Dictionary with processed data structures
        """
        if self.partner_data is None or "partner_name" not in self.partner_data.columns:
            return None
        
        import pandas as pd
        
        # Separate end and linker columns
        end_cols = [c for c in self.partner_data.columns if c.endswith("_end")]
        linker_cols = [c for c in self.partner_data.columns if c.endswith("_linker")]
        
        # Extract sample names from column names
        samples = set()
        for col in end_cols + linker_cols:
            # Remove _end or _linker suffix to get sample name
            sample = col.rsplit("_", 1)[0]
            samples.add(sample)
        samples = sorted(list(samples))
        
        # Create detection matrix (binary: detected or not)
        detection_matrix = OrderedDict()
        for _, row in self.partner_data.iterrows():
            partner = row["partner_name"]
            detection_matrix[partner] = {}
            for sample in samples:
                end_col = f"{sample}_end"
                linker_col = f"{sample}_linker"
                end_count = row.get(end_col, 0) if end_col in row else 0
                linker_count = row.get(linker_col, 0) if linker_col in row else 0
                # Consider detected if either end or linker has counts > 0
                detection_matrix[partner][sample] = 1 if (end_count > 0 or linker_count > 0) else 0
        
        # Create count matrix (total counts per partner per sample)
        count_matrix = OrderedDict()
        for _, row in self.partner_data.iterrows():
            partner = row["partner_name"]
            count_matrix[partner] = {}
            for sample in samples:
                end_col = f"{sample}_end"
                linker_col = f"{sample}_linker"
                end_count = row.get(end_col, 0) if end_col in row else 0
                linker_col = row.get(linker_col, 0) if linker_col in row else 0
                count_matrix[partner][sample] = end_count + linker_col
        
        # Create partner coverage data (number of samples with detection per partner)
        partner_coverage = []
        for partner in detection_matrix:
            num_samples = sum(detection_matrix[partner].values())
            total_counts = sum(count_matrix[partner].values())
            partner_coverage.append({
                "partner": partner,
                "samples_detected": num_samples,
                "total_counts": total_counts,
            })
        
        # Create end vs linker comparison data
        end_linker_data = []
        for sample in samples:
            end_col = f"{sample}_end"
            linker_col = f"{sample}_linker"
            if end_col in self.partner_data.columns and linker_col in self.partner_data.columns:
                total_end = self.partner_data[end_col].sum()
                total_linker = self.partner_data[linker_col].sum()
                end_linker_data.append({
                    "sample": sample,
                    "end_counts": total_end,
                    "linker_counts": total_linker,
                })
        
        return {
            "detection_matrix": detection_matrix,
            "count_matrix": count_matrix,
            "partner_coverage": pd.DataFrame(partner_coverage),
            "end_linker_data": pd.DataFrame(end_linker_data),
            "samples": samples,
        }
    
    def partner_detection_heatmap(self):
        """Create partner detection heatmap."""
        if self.processed_data is None:
            return
        
        detection_matrix = self.processed_data["detection_matrix"]
        if not detection_matrix:
            return
        
        # Convert to format expected by MultiQC heatmap
        # MultiQC expects: {sample: {partner: value}}
        heatmap_data = OrderedDict()
        samples = self.processed_data["samples"]
        
        for sample in samples:
            heatmap_data[sample] = OrderedDict()
            for partner in detection_matrix:
                heatmap_data[sample][partner] = detection_matrix[partner].get(sample, 0)
        
        pconfig = {
            "id": "fusilli_partner_heatmap",
            "title": "FUSILLI: Partner Detection Heatmap",
            "xlab": "Partner",
            "ylab": "Sample",
            "square": False,
            "colstops": [
                {"0": "#ffffff", "1": "#1f77b4"},
            ],
        }
        
        self.add_section(
            name="Partner Detection Heatmap",
            anchor="fusilli_partner_heatmap",
            description="Binary detection matrix of partners across samples.",
            helptext="""
            This heatmap shows which partners are detected in which samples.
            - **Blue**: Partner detected in sample
            - **White**: Partner not detected in sample
            
            Partners are sorted by detection frequency (most common first).
            """,
            plot=heatmap.plot(heatmap_data, pconfig),
        )
    
    def partner_coverage_plot(self):
        """Create partner coverage bar chart."""
        if self.processed_data is None:
            return
        
        partner_coverage = self.processed_data["partner_coverage"]
        if partner_coverage is None or partner_coverage.empty:
            return
        
        # Sort by number of samples detected
        partner_coverage = partner_coverage.sort_values(
            "samples_detected", ascending=False
        ).head(50)  # Limit to top 50 for readability
        
        plot_data = OrderedDict()
        for _, row in partner_coverage.iterrows():
            partner = row["partner"]
            plot_data[partner] = {"Samples Detected": row["samples_detected"]}
        
        if not plot_data:
            return
        
        pconfig = {
            "id": "fusilli_partner_coverage",
            "title": "FUSILLI: Partner Coverage (Top 50)",
            "ylab": "Number of Samples with Detection",
            "xlab": "Partner",
            "tt_label": "<b>{point.x}</b><br>{point.y} samples",
            "rotate_labels": True,
        }
        
        self.add_section(
            name="Partner Coverage",
            anchor="fusilli_partner_coverage",
            description="Number of samples with detection per partner (top 50).",
            helptext="""
            This plot shows how many samples have detected each partner.
            Partners are sorted by detection frequency (most common first).
            Only the top 50 partners are shown for readability.
            """,
            plot=bargraph.plot(plot_data, pconfig),
        )
    
    def partner_end_vs_linker_plot(self):
        """Create partner end vs linker detection scatter plot."""
        if self.processed_data is None:
            return
        
        end_linker_data = self.processed_data["end_linker_data"]
        if end_linker_data is None or end_linker_data.empty:
            return
        
        plot_data = OrderedDict()
        for _, row in end_linker_data.iterrows():
            sample = row["sample"]
            plot_data[sample] = {
                "x": row["end_counts"],
                "y": row["linker_counts"],
            }
        
        if not plot_data:
            return
        
        pconfig = {
            "id": "fusilli_partner_end_linker",
            "title": "FUSILLI: Partner End vs Linker Detection",
            "ylab": "Linker Counts",
            "xlab": "End Counts",
            "tt_label": "<b>{point.sample}</b><br>End: {point.x}<br>Linker: {point.y}",
        }
        
        self.add_section(
            name="Partner End vs Linker",
            anchor="fusilli_partner_end_linker",
            description="Comparison of partner end and linker detection counts.",
            helptext="""
            This scatter plot compares detection counts for:
            - **End counts**: Detections via partner domain 3' end sequences
            - **Linker counts**: Detections via linker sequences
            
            Points along the diagonal indicate similar detection rates for both methods.
            """,
            plot=scatter.plot(plot_data, pconfig),
        )
