#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chromosome Heatmap Visualization Module
Used for generating heatmap visualizations of all chromosomes for individual samples
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from typing import Dict, List, Tuple
import seaborn as sns
from snp_data_parser import SNPDataParser

# Set matplotlib to support fonts
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

class ChromosomeHeatmapVisualizer:
    """
    Chromosome Heatmap Visualizer
    """
    
    def __init__(self, parser: SNPDataParser):
        """
        Initialize the visualizer
        
        Args:
            parser: SNP data parser instance
        """
        self.parser = parser
        self.colors = plt.cm.RdYlBu_r  # Red-Yellow-Blue spectrum, suitable for scientific publication
        
    def create_single_sample_heatmap(self, sample: str, window_size: int = 50000, 
                                   figsize: Tuple[int, int] = (16, 12),
                                   save_path: str = None) -> plt.Figure:
        """
        Create heatmap for all chromosomes of a single sample
        
        Args:
            sample: Sample name
            window_size: Window size (bp)
            figsize: Figure size
            save_path: Save path
            
        Returns:
            matplotlib Figure object
        """
        print(f"Generating chromosome heatmap for sample {sample}...")
        
        # Get sample data across all chromosomes
        chrom_data = self.parser.get_sample_all_chromosomes_data(sample, window_size)
        
        if not chrom_data:
            raise ValueError(f"Data for sample {sample} is empty")
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Set chromosome spacing and height
        chrom_height = 0.8
        chrom_spacing = 1.2
        
        # Get maximum chromosome length for normalization
        max_length = max(self.parser.chromosome_lengths.values())
        
        # Draw each chromosome
        y_positions = []
        chrom_labels = []
        
        for i, chromosome in enumerate(sorted(self.parser.chromosomes)):
            if chromosome not in chrom_data:
                continue
                
            positions, rates = chrom_data[chromosome]
            
            if not positions:
                continue
            
            # Calculate y position
            y_pos = len(self.parser.chromosomes) - i - 1
            y_positions.append(y_pos)
            chrom_labels.append(chromosome)
            
            # Chromosome length
            chrom_length = self.parser.chromosome_lengths[chromosome]
            
            # Draw chromosome background
            chrom_width = chrom_length / max_length * 0.9  # Normalized width
            
            # Draw heatmap bands
            for j, (pos, rate) in enumerate(zip(positions, rates)):
                # Calculate window position and width in the figure
                window_start = (pos - window_size//2) / max_length * 0.9
                window_width = window_size / max_length * 0.9
                
                # Ensure not exceeding chromosome boundaries
                if window_start < 0:
                    window_width += window_start
                    window_start = 0
                if window_start + window_width > chrom_width:
                    window_width = chrom_width - window_start
                
                if window_width > 0:
                    # Draw heatmap rectangle
                    rect = patches.Rectangle(
                        (window_start, y_pos - chrom_height/2),
                        window_width, chrom_height,
                        facecolor=self.colors(rate),
                        edgecolor='none',
                        alpha=0.8
                    )
                    ax.add_patch(rect)
            
            # Draw chromosome border
            border = patches.Rectangle(
                (0, y_pos - chrom_height/2),
                chrom_width, chrom_height,
                facecolor='none',
                edgecolor='black',
                linewidth=1,
                alpha=0.7
            )
            ax.add_patch(border)
        
        # Set axes
        ax.set_xlim(0, 1.0)
        ax.set_ylim(-0.5, len(self.parser.chromosomes) - 0.5)
        
        # Set y-axis labels
        ax.set_yticks(y_positions)
        ax.set_yticklabels(chrom_labels, fontsize=14)
        ax.set_ylabel('Chromosomes', fontsize=15, fontweight='bold')
        
        # Set x-axis labels (genomic position) - using 5/10/20/50 as units
        max_mb = max_length / 1e6
        # Calculate appropriate interval, at least in units of 50
        if max_mb <= 50:
            tick_interval = 5
        elif max_mb <= 100:
            tick_interval = 10
        elif max_mb <= 200:
            tick_interval = 20
        else:
            tick_interval = 50
        
        x_tick_values = np.arange(0, max_mb + tick_interval, tick_interval)
        x_ticks = x_tick_values / max_mb * 0.9
        x_labels = [f"{int(val)}" for val in x_tick_values]
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_labels, fontsize=13, rotation=0)
        ax.set_xlabel('Genomic Position (Mb)', fontsize=15, fontweight='bold')
        
        # Set title
        ax.set_title(f'SNP Match Rate Heatmap for Sample: {sample}\n'
                    f'Window Size: {window_size/1000:.0f} kb', 
                    fontsize=16, fontweight='bold', pad=20)
        
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=self.colors, 
                                  norm=plt.Normalize(vmin=0, vmax=1))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.6, aspect=20, pad=0.02)
        cbar.set_label('SNP Match Rate', fontsize=14, fontweight='bold')
        cbar.ax.tick_params(labelsize=14)
        
        # Remove grid lines
        ax.grid(False)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save figure
        if save_path:
            # Save PNG format
            plt.savefig(save_path, dpi=300, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
        
        return fig
    
    def create_multiple_samples_comparison(self, samples: List[str], 
                                         window_size: int = 50000,
                                         figsize: Tuple[int, int] = (20, 15),
                                         save_path: str = None) -> plt.Figure:
        """
        Create chromosome heatmap comparison for multiple samples
        
        Args:
            samples: List of sample names
            window_size: Window size
            figsize: Figure size
            save_path: Save path
            
        Returns:
            matplotlib Figure object
        """
        print(f"Generating chromosome heatmap comparison for {len(samples)} samples...")
        
        # Create subplots
        n_samples = len(samples)
        fig, axes = plt.subplots(n_samples, 1, figsize=figsize, 
                                sharex=True, sharey=True)
        
        if n_samples == 1:
            axes = [axes]
        
        # Get maximum chromosome length
        max_length = max(self.parser.chromosome_lengths.values())
        
        for idx, sample in enumerate(samples):
            ax = axes[idx]
            
            # Get sample data
            chrom_data = self.parser.get_sample_all_chromosomes_data(sample, window_size)
            
            # Draw chromosome heatmap
            self._draw_chromosome_heatmap(ax, chrom_data, sample, window_size, max_length)
        
        # Set overall title
        fig.suptitle('SNP Match Rate Comparison Across Samples\n'
                    f'Window Size: {window_size/1000:.0f} kb', 
                    fontsize=18, fontweight='bold', y=0.98)
        
        # Add shared colorbar
        sm = plt.cm.ScalarMappable(cmap=self.colors, 
                                  norm=plt.Normalize(vmin=0, vmax=1))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes, shrink=0.6, aspect=25, pad=0.02)
        cbar.set_label('SNP Match Rate', fontsize=14, fontweight='bold')
        cbar.ax.tick_params(labelsize=14)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save figure
        if save_path:
            # Save PNG format
            plt.savefig(save_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            print(f"Comparison chart saved to: {save_path}")
            
            # Save SVG format
            svg_path = save_path.replace('.png', '.svg')
            plt.savefig(svg_path, format='svg', bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            print(f"SVG comparison chart saved to: {svg_path}")
        
        return fig
    
    def _draw_chromosome_heatmap(self, ax, chrom_data: Dict, sample: str, 
                               window_size: int, max_length: int):
        """
        Draw chromosome heatmap on specified axis
        
        Args:
            ax: matplotlib axis object
            chrom_data: Chromosome data dictionary
            sample: Sample name
            window_size: Window size
            max_length: Maximum chromosome length
        """
        chrom_height = 0.8
        
        # Draw each chromosome
        y_positions = []
        chrom_labels = []
        
        for i, chromosome in enumerate(sorted(self.parser.chromosomes)):
            if chromosome not in chrom_data:
                continue
                
            positions, rates = chrom_data[chromosome]
            
            if not positions:
                continue
            
            y_pos = len(self.parser.chromosomes) - i - 1
            y_positions.append(y_pos)
            chrom_labels.append(chromosome)
            
            chrom_length = self.parser.chromosome_lengths[chromosome]
            chrom_width = chrom_length / max_length * 0.9
            
            # Draw heatmap
            for pos, rate in zip(positions, rates):
                window_start = (pos - window_size//2) / max_length * 0.9
                window_width = window_size / max_length * 0.9
                
                if window_start < 0:
                    window_width += window_start
                    window_start = 0
                if window_start + window_width > chrom_width:
                    window_width = chrom_width - window_start
                
                if window_width > 0:
                    rect = patches.Rectangle(
                        (window_start, y_pos - chrom_height/2),
                        window_width, chrom_height,
                        facecolor=self.colors(rate),
                        edgecolor='none',
                        alpha=0.8
                    )
                    ax.add_patch(rect)
            
            # Chromosome border
            border = patches.Rectangle(
                (0, y_pos - chrom_height/2),
                chrom_width, chrom_height,
                facecolor='none',
                edgecolor='black',
                linewidth=1,
                alpha=0.7
            )
            ax.add_patch(border)
        
        # Set axes
        ax.set_xlim(-0.05, 1.0)
        ax.set_ylim(-0.5, len(self.parser.chromosomes) - 0.5)
        ax.set_yticks(y_positions)
        ax.set_yticklabels(chrom_labels, fontsize=10)
        
        # Set x-axis
        x_ticks = np.linspace(0, 0.9, 6)
        x_labels = [f"{int(x * max_length / 1e6):.0f}" for x in x_ticks]
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_labels, fontsize=10)
        ax.set_xlabel('Genomic Position (Mb)', fontsize=12)
        
        # Sample label
        ax.set_ylabel(f'Sample: {sample}', fontsize=12, fontweight='bold')
        
        # Grid
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)


if __name__ == "__main__":
    # Test code
    parser = SNPDataParser("input.txt")
    parser.load_data()
    
    visualizer = ChromosomeHeatmapVisualizer(parser)
    
    # Generate heatmap for the first sample
    sample = parser.samples[0]
    fig = visualizer.create_single_sample_heatmap(
        sample, 
        window_size=50000,
        save_path=f"chromosome_heatmap_{sample}.png"
    )
    plt.show()