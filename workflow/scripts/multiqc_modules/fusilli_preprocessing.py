"""
MultiQC module for FUSILLI preprocessing metrics.

This module visualizes:
- Read decay through preprocessing steps
- Retention rates per step
- Step loss breakdown
"""

import logging
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseModule
from multiqc.plots import linegraph, bargraph

from .fusilli_base import parse_csv_file

logger = logging.getLogger(__name__)


class MultiqcModule(BaseModule):
    """
    FUSILLI Preprocessing Metrics Module
    
    Parses decay_metrics.csv and step-specific metrics to visualize read retention
    through the preprocessing pipeline.
    """
    
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="FUSILLI Preprocessing Metrics",
            anchor="fusilli_preprocessing",
            info="Read retention and loss metrics through preprocessing steps.",
        )
        
        # Find and parse input files
        self.decay_data = None
        self.trim_data = None
        self.contam_data = None
        self.quality_data = None
        
        # Find decay_metrics.csv
        for f in self.find_log_files("decay_metrics.csv"):
            self.decay_data = parse_csv_file(f["fn"])
            if self.decay_data is not None:
                break
        
        # Find step-specific metrics
        for f in self.find_log_files("trim_metrics.csv"):
            self.trim_data = parse_csv_file(f["fn"])
            if self.trim_data is not None:
                break
        
        for f in self.find_log_files("contam_metrics.csv"):
            self.contam_data = parse_csv_file(f["fn"])
            if self.contam_data is not None:
                break
        
        for f in self.find_log_files("quality_metrics.csv"):
            self.quality_data = parse_csv_file(f["fn"])
            if self.quality_data is not None:
                break
        
        # If no data found, exit
        if self.decay_data is None:
            raise UserWarning
        
        # Calculate retention rates
        self.retention_data = self.calculate_retention_rates()
        
        # Generate plots
        if self.decay_data is not None:
            self.read_decay_plot()
        
        if self.retention_data is not None and not self.retention_data.empty:
            self.retention_rate_plot()
            self.step_loss_plot()
            self.add_summary_table()
    
    def calculate_retention_rates(self):
        """
        Calculate retention rates for each preprocessing step.
        
        Returns:
            DataFrame with columns: sample, step, retention_rate
        """
        if self.decay_data is None:
            return None
        
        import pandas as pd
        
        # Get raw reads for each sample
        raw_data = self.decay_data[self.decay_data["step"] == "raw"].copy()
        if raw_data.empty:
            return None
        
        raw_reads = dict(zip(raw_data["sample"], raw_data["reads"]))
        
        # Calculate retention rates for each step
        retention_list = []
        steps = ["trimmed", "cleaned", "quality", "merged"]
        
        for step in steps:
            step_data = self.decay_data[self.decay_data["step"] == step].copy()
            for _, row in step_data.iterrows():
                sample = row["sample"]
                reads = row["reads"]
                raw = raw_reads.get(sample, 0)
                retention = (reads / raw) if raw > 0 else 0.0
                retention_list.append({
                    "sample": sample,
                    "step": step,
                    "reads": reads,
                    "raw_reads": raw,
                    "retention_rate": retention,
                    "loss_rate": 1.0 - retention,
                })
        
        return pd.DataFrame(retention_list)
    
    def read_decay_plot(self):
        """Create read decay line plot through preprocessing steps."""
        if self.decay_data is None:
            return
        
        # Prepare data for line plot
        plot_data = OrderedDict()
        steps = ["raw", "trimmed", "cleaned", "quality", "merged"]
        
        for sample in self.decay_data["sample"].unique():
            sample_data = self.decay_data[self.decay_data["sample"] == sample]
            plot_data[sample] = {}
            for step in steps:
                step_row = sample_data[sample_data["step"] == step]
                if not step_row.empty:
                    reads = step_row.iloc[0]["reads"]
                    plot_data[sample][step] = reads
        
        if not plot_data:
            return
        
        pconfig = {
            "id": "fusilli_read_decay",
            "title": "FUSILLI: Read Decay Through Preprocessing",
            "ylab": "Number of Reads",
            "xlab": "Preprocessing Step",
            "tt_label": "<b>{point.x}</b><br>{point.y:,.0f} reads",
            "yLog": True,
        }
        
        self.add_section(
            name="Read Decay",
            anchor="fusilli_read_decay",
            description="Read counts at each preprocessing step (log scale).",
            helptext="""
            This plot shows how read counts change through the preprocessing pipeline:
            - **raw**: Initial read count
            - **trimmed**: After adapter trimming
            - **cleaned**: After contaminant removal
            - **quality**: After quality filtering
            - **merged**: After read merging
            
            The y-axis is on a log scale to better visualize the decay pattern.
            """,
            plot=linegraph.plot(plot_data, pconfig),
        )
    
    def retention_rate_plot(self):
        """Create retention rate grouped bar chart."""
        if self.retention_data is None or self.retention_data.empty:
            return
        
        plot_data = OrderedDict()
        for _, row in self.retention_data.iterrows():
            sample = row["sample"]
            step = row["step"]
            if sample not in plot_data:
                plot_data[sample] = {}
            plot_data[sample][step] = row["retention_rate"]
        
        if not plot_data:
            return
        
        pconfig = {
            "id": "fusilli_retention_rates",
            "title": "FUSILLI: Retention Rates by Step",
            "ylab": "Retention Rate",
            "xlab": "Sample",
            "stacking": None,
            "tt_label": "<b>{point.x}</b><br>{series.name}: {point.y:.3f}",
            "max": 1.0,
            "min": 0.0,
        }
        
        self.add_section(
            name="Retention Rates",
            anchor="fusilli_retention_rates",
            description="Read retention rates at each preprocessing step.",
            helptext="""
            This plot shows the fraction of reads retained after each preprocessing step:
            - **trimmed**: Retention after adapter trimming
            - **cleaned**: Retention after contaminant removal
            - **quality**: Retention after quality filtering
            - **merged**: Retention after read merging (fraction of reads that merged)
            
            Higher retention rates indicate less aggressive filtering.
            """,
            plot=bargraph.plot(plot_data, pconfig),
        )
    
    def step_loss_plot(self):
        """Create step loss stacked bar chart."""
        if self.retention_data is None or self.retention_data.empty:
            return
        
        plot_data = OrderedDict()
        for sample in self.retention_data["sample"].unique():
            sample_data = self.retention_data[
                self.retention_data["sample"] == sample
            ]
            plot_data[sample] = {}
            
            # Calculate loss at each step
            prev_retention = 1.0
            for step in ["trimmed", "cleaned", "quality", "merged"]:
                step_row = sample_data[sample_data["step"] == step]
                if not step_row.empty:
                    retention = step_row.iloc[0]["retention_rate"]
                    loss = prev_retention - retention
                    plot_data[sample][step] = loss
                    prev_retention = retention
                else:
                    plot_data[sample][step] = 0.0
        
        if not plot_data:
            return
        
        pconfig = {
            "id": "fusilli_step_loss",
            "title": "FUSILLI: Step Loss Breakdown",
            "ylab": "Fraction of Raw Reads",
            "xlab": "Sample",
            "stacking": "normal",
            "tt_label": "<b>{point.x}</b><br>{series.name}: {point.y:.3f}",
        }
        
        self.add_section(
            name="Step Loss Breakdown",
            anchor="fusilli_step_loss",
            description="Fraction of reads lost at each preprocessing step.",
            helptext="""
            This stacked bar chart shows the cumulative loss of reads through preprocessing:
            - **trimmed**: Loss during adapter trimming
            - **cleaned**: Additional loss during contaminant removal
            - **quality**: Additional loss during quality filtering
            - **merged**: Additional loss during read merging (unmerged reads)
            
            The total height shows the final retention rate.
            """,
            plot=bargraph.plot(plot_data, pconfig),
        )
    
    def add_summary_table(self):
        """Add retention rate metrics to general stats table."""
        if self.retention_data is None or self.retention_data.empty:
            return
        
        # Calculate final retention (merged step)
        merged_data = self.retention_data[
            self.retention_data["step"] == "merged"
        ].copy()
        
        if merged_data.empty:
            return
        
        headers = OrderedDict()
        headers["retention_rate"] = {
            "title": "Final Retention",
            "description": "Fraction of raw reads retained after all preprocessing",
            "format": "{:.3f}",
            "max": 1.0,
            "min": 0.0,
        }
        
        self.general_stats_addcols(merged_data, headers)
