"""
MultiQC module for FUSILLI diversity metrics.

This module visualizes:
- Shannon and Simpson diversity indices
- Evenness (Pielou's evenness)
- Top N fractions (top1, top5, top10, top25)
- Variant count distribution
"""

import logging
import numpy as np
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseModule
from multiqc.plots import linegraph, bargraph, scatter

from fusilli_multiqc.utils import parse_csv_file

logger = logging.getLogger(__name__)


class MultiqcModule(BaseModule):
    """
    FUSILLI Diversity Metrics Module

    Parses fusion_counts_summary.csv to calculate and visualize library diversity metrics.
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="FUSILLI Diversity Metrics",
            anchor="fusilli_diversity",
            info="Library diversity metrics including Shannon/Simpson diversity and evenness.",
        )

        # Find and parse input files
        self.fusion_counts_data = None

        # Find fusion_counts_summary.csv
        for f in self.find_log_files("fusion_counts_summary.csv"):
            self.fusion_counts_data = parse_csv_file(f["fn"])
            if self.fusion_counts_data is not None:
                break

        # If no data found, exit
        if self.fusion_counts_data is None:
            raise UserWarning

        # Calculate diversity metrics
        self.diversity_metrics = self.calculate_diversity_metrics()

        # Generate plots
        if self.diversity_metrics is not None and not self.diversity_metrics.empty:
            self.diversity_comparison_plot()
            self.evenness_plot()
            self.top_n_fractions_plot()
            self.variant_count_distribution_plot()
            self.add_summary_table()

    def calculate_diversity_metrics(self):
        """
        Calculate diversity metrics for each sample.

        Returns:
            DataFrame with columns: sample, shannon_diversity, simpson_diversity, evenness,
            top1_fraction, top5_fraction, top10_fraction, top25_fraction
        """
        if self.fusion_counts_data is None:
            return None

        import pandas as pd

        # Get sample columns (exclude fusion_id, type, total)
        sample_cols = [
            c for c in self.fusion_counts_data.columns
            if c not in ["fusion_id", "type", "total"]
        ]

        if not sample_cols:
            return None

        # Filter to fusion type only
        fusion_data = self.fusion_counts_data[
            self.fusion_counts_data["type"] == "fusion"
        ].copy()

        if fusion_data.empty:
            return None

        metrics_list = []

        for sample in sample_cols:
            # Get counts for this sample
            counts = fusion_data[sample].fillna(0).values
            total_counts = counts.sum()

            if total_counts == 0:
                metrics_list.append({
                    "sample": sample,
                    "shannon_diversity": 0.0,
                    "simpson_diversity": 0.0,
                    "evenness": 0.0,
                    "top1_fraction": 0.0,
                    "top5_fraction": 0.0,
                    "top10_fraction": 0.0,
                    "top25_fraction": 0.0,
                    "observed_variants": 0,
                })
                continue

            # Filter out zero counts
            non_zero_counts = counts[counts > 0]
            if len(non_zero_counts) == 0:
                metrics_list.append({
                    "sample": sample,
                    "shannon_diversity": 0.0,
                    "simpson_diversity": 0.0,
                    "evenness": 0.0,
                    "top1_fraction": 0.0,
                    "top5_fraction": 0.0,
                    "top10_fraction": 0.0,
                    "top25_fraction": 0.0,
                    "observed_variants": 0,
                })
                continue

            # Calculate proportions
            proportions = non_zero_counts / total_counts

            # Shannon diversity: H' = -Σ(p_i * ln(p_i))
            shannon = -np.sum(proportions * np.log(proportions))

            # Simpson diversity: 1 - Σ(p_i²)
            simpson = 1 - np.sum(proportions ** 2)

            # Evenness: H' / ln(S) where S is number of observed variants
            observed_variants = len(non_zero_counts)
            if observed_variants > 1:
                evenness = shannon / np.log(observed_variants)
            else:
                evenness = 0.0

            # Top N fractions
            sorted_counts = np.sort(non_zero_counts)[::-1]  # Descending order
            top1 = sorted_counts[0] if len(sorted_counts) > 0 else 0
            top5 = sorted_counts[:5].sum() if len(sorted_counts) >= 5 else sorted_counts.sum()
            top10 = sorted_counts[:10].sum() if len(sorted_counts) >= 10 else sorted_counts.sum()
            top25 = sorted_counts[:25].sum() if len(sorted_counts) >= 25 else sorted_counts.sum()

            metrics_list.append({
                "sample": sample,
                "shannon_diversity": float(shannon),
                "simpson_diversity": float(simpson),
                "evenness": float(evenness),
                "top1_fraction": float(top1 / total_counts),
                "top5_fraction": float(top5 / total_counts),
                "top10_fraction": float(top10 / total_counts),
                "top25_fraction": float(top25 / total_counts),
                "observed_variants": int(observed_variants),
            })

        return pd.DataFrame(metrics_list)

    def diversity_comparison_plot(self):
        """Create diversity comparison bar chart."""
        if self.diversity_metrics is None or self.diversity_metrics.empty:
            return

        plot_data = OrderedDict()
        for _, row in self.diversity_metrics.iterrows():
            sample = row["sample"]
            plot_data[sample] = {
                "Shannon Diversity": row["shannon_diversity"],
                "Simpson Diversity": row["simpson_diversity"],
            }

        if not plot_data:
            return

        pconfig = {
            "id": "fusilli_diversity_comparison",
            "title": "FUSILLI: Diversity Index Comparison",
            "ylab": "Diversity Index",
            "xlab": "Sample",
            "stacking": None,
            "tt_label": "<b>{point.x}</b>: {point.y:.3f}",
        }

        self.add_section(
            name="Diversity Comparison",
            anchor="fusilli_diversity_comparison",
            description="Shannon and Simpson diversity indices across samples.",
            helptext="""
            This plot compares two diversity metrics:
            - **Shannon Diversity**: Measures both richness and evenness (higher = more diverse)
            - **Simpson Diversity**: Measures dominance (higher = less dominance, more diverse)

            Higher values indicate greater library diversity.
            """,
            plot=bargraph.plot(plot_data, pconfig),
        )

    def evenness_plot(self):
        """Create evenness bar chart."""
        if self.diversity_metrics is None or self.diversity_metrics.empty:
            return

        plot_data = OrderedDict()
        for _, row in self.diversity_metrics.iterrows():
            sample = row["sample"]
            plot_data[sample] = {"Evenness": row["evenness"]}

        if not plot_data:
            return

        pconfig = {
            "id": "fusilli_evenness",
            "title": "FUSILLI: Library Evenness",
            "ylab": "Evenness (Pielou's)",
            "xlab": "Sample",
            "tt_label": "<b>{point.x}</b>: {point.y:.3f}",
            "max": 1.0,
            "min": 0.0,
        }

        self.add_section(
            name="Evenness",
            anchor="fusilli_evenness",
            description="Pielou's evenness index across samples.",
            helptext="""
            Evenness measures how evenly distributed variant counts are.
            Values range from 0 (highly uneven) to 1 (perfectly even).
            Higher evenness indicates more uniform representation of variants.
            """,
            plot=bargraph.plot(plot_data, pconfig),
        )

    def top_n_fractions_plot(self):
        """Create top N fractions stacked bar chart."""
        if self.diversity_metrics is None or self.diversity_metrics.empty:
            return

        plot_data = OrderedDict()
        for _, row in self.diversity_metrics.iterrows():
            sample = row["sample"]
            plot_data[sample] = {
                "Top 1": row["top1_fraction"],
                "Top 5": row["top5_fraction"] - row["top1_fraction"],
                "Top 10": row["top10_fraction"] - row["top5_fraction"],
                "Top 25": row["top25_fraction"] - row["top10_fraction"],
                "Other": 1.0 - row["top25_fraction"],
            }

        if not plot_data:
            return

        pconfig = {
            "id": "fusilli_top_n_fractions",
            "title": "FUSILLI: Top N Variant Fractions",
            "ylab": "Fraction of Total Counts",
            "xlab": "Sample",
            "stacking": "normal",
            "tt_label": "<b>{point.x}</b><br>{series.name}: {point.y:.3f}",
        }

        self.add_section(
            name="Top N Fractions",
            anchor="fusilli_top_n_fractions",
            description="Fraction of counts contributed by top N variants.",
            helptext="""
            This stacked bar chart shows the distribution of counts across variants:
            - **Top 1**: Fraction from the most abundant variant
            - **Top 5**: Fraction from variants ranked 2-5
            - **Top 10**: Fraction from variants ranked 6-10
            - **Top 25**: Fraction from variants ranked 11-25
            - **Other**: Fraction from all remaining variants

            More even distribution (less in Top 1) indicates better library diversity.
            """,
            plot=bargraph.plot(plot_data, pconfig),
        )

    def variant_count_distribution_plot(self):
        """Create variant count distribution histogram."""
        if self.fusion_counts_data is None:
            return

        # Get fusion counts (exclude type column and sample columns for now)
        fusion_data = self.fusion_counts_data[
            self.fusion_counts_data["type"] == "fusion"
        ].copy()

        if fusion_data.empty:
            return

        # Get sample columns
        sample_cols = [
            c for c in fusion_data.columns
            if c not in ["fusion_id", "type", "total"]
        ]

        if not sample_cols:
            return

        # Aggregate counts across all samples
        all_counts = []
        for sample in sample_cols:
            counts = fusion_data[sample].fillna(0).values
            all_counts.extend(counts[counts > 0].tolist())

        if not all_counts:
            return

        # Create histogram data (log scale bins)
        import numpy as np
        log_counts = np.log10(np.array(all_counts) + 1)
        bins = np.linspace(0, np.ceil(log_counts.max()), 20)
        hist, bin_edges = np.histogram(log_counts, bins=bins)

        # Convert to plot format
        plot_data = OrderedDict()
        for i in range(len(hist)):
            bin_label = f"10^{bin_edges[i]:.1f} - 10^{bin_edges[i+1]:.1f}"
            plot_data[bin_label] = {"Count": int(hist[i])}

        pconfig = {
            "id": "fusilli_variant_distribution",
            "title": "FUSILLI: Variant Count Distribution",
            "ylab": "Number of Variants",
            "xlab": "Count (log10 scale)",
            "tt_label": "<b>{point.x}</b>: {point.y} variants",
        }

        self.add_section(
            name="Variant Count Distribution",
            anchor="fusilli_variant_distribution",
            description="Distribution of variant counts (log scale).",
            helptext="""
            This histogram shows the distribution of variant abundances across all samples.
            The x-axis is on a log10 scale to better visualize the wide range of counts.
            A more uniform distribution indicates better library representation.
            """,
            plot=bargraph.plot(plot_data, pconfig),
        )

    def add_summary_table(self):
        """Add diversity metrics to general stats table."""
        if self.diversity_metrics is None or self.diversity_metrics.empty:
            return

        headers = OrderedDict()
        headers["shannon_diversity"] = {
            "title": "Shannon Diversity",
            "description": "Shannon diversity index (H')",
            "format": "{:.3f}",
        }
        headers["simpson_diversity"] = {
            "title": "Simpson Diversity",
            "description": "Simpson diversity index (1-D)",
            "format": "{:.3f}",
        }
        headers["evenness"] = {
            "title": "Evenness",
            "description": "Pielou's evenness index",
            "format": "{:.3f}",
            "max": 1.0,
            "min": 0.0,
        }
        headers["observed_variants"] = {
            "title": "Observed Variants",
            "description": "Number of unique variants detected",
            "format": "{:,.0f}",
        }

        self.general_stats_addcols(self.diversity_metrics, headers)
