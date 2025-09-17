#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sample Comparison Heatmap Visualization Module
Used for generating comparison heatmap visualizations of multiple samples on the same chromosome
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

class SampleComparisonHeatmapVisualizer:
    """
    Sample Comparison Heatmap Visualizer
    """
    
    def __init__(self, parser: SNPDataParser):
        """
        Initialize the visualizer
        
        Args:
            parser: SNP data parser instance
        """
        self.parser = parser
        self.colors = plt.cm.RdYlBu_r  # Red-Yellow-Blue spectrum
        
    def create_chromosome_comparison_heatmap(self, chromosome: str, samples: List[str] = None,
                                           window_size: int = 50000,
                                           figsize: Tuple[int, int] = (16, 12),
                                           save_path: str = None) -> plt.Figure:
        """
        Create comparison heatmap for multiple samples on the same chromosome
        
        Args:
            chromosome: Chromosome name
            samples: List of sample names, if None use all samples
            window_size: Window size (bp)
            figsize: Figure size
            save_path: Save path
            
        Returns:
            matplotlib Figure object
        """
        if samples is None:
            samples = self.parser.samples
        
        print(f"Generating sample comparison heatmap for chromosome {chromosome}...")
        print(f"Including samples: {samples}")
        
        # Get data for multiple samples on the same chromosome
        chrom_data = self.parser.get_chromosome_data_for_samples(
            chromosome, samples, window_size
        )
        
        if not chrom_data:
            raise ValueError(f"Data for chromosome {chromosome} is empty")
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Get chromosome length
        chrom_length = self.parser.chromosome_lengths[chromosome]
        
        # Calculate heatmap matrix
        heatmap_matrix, window_positions, sample_labels = self._prepare_heatmap_matrix(
            chrom_data, samples, chrom_length, window_size
        )
        
        if heatmap_matrix.size == 0:
            raise ValueError(f"Cannot generate heatmap data for chromosome {chromosome}")
        
        # Create heatmap
        im = ax.imshow(heatmap_matrix, cmap=self.colors, aspect='auto', 
                      vmin=0, vmax=1, interpolation='nearest')
        
        # Set axis labels
        self._set_axis_labels(ax, window_positions, sample_labels, chrom_length)
        
        # Set title
        chrom_display = chromosome
        ax.set_title(f'SNP Match Rate Comparison - {chrom_display}\n'
                    f'Window Size: {window_size/1000:.0f} kb | '
                    f'Chromosome Length: {chrom_length/1e6:.1f} Mb',
                    fontsize=16, fontweight='bold', pad=20)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.6, aspect=20, pad=0.02)
        cbar.set_label('SNP Match Rate', fontsize=14, fontweight='bold')
        cbar.ax.tick_params(labelsize=14)
        
        # Remove grid lines
        # self._add_grid_lines(ax, len(sample_labels), len(window_positions))
        
        # Adjust layout
        plt.tight_layout()
        
        # Save figure
        if save_path:
            # Save PNG format
            plt.savefig(save_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            print(f"Figure saved to: {save_path}")
        
        return fig
    
    def create_all_chromosomes_comparison(self, samples: List[str] = None,
                                        window_size: int = 50000,
                                        figsize: Tuple[int, int] = (20, 25),
                                        save_path: str = None) -> plt.Figure:
        """
        Create sample comparison heatmaps for all chromosomes
        
        Args:
            samples: List of sample names
            window_size: Window size
            figsize: Figure size
            save_path: Save path
            
        Returns:
            matplotlib Figure object
        """
        if samples is None:
            samples = self.parser.samples
        
        print(f"Generating sample comparison heatmaps for all chromosomes...")
        
        # Create subplots with increased chromosome spacing
        n_chroms = len(self.parser.chromosomes)
        fig, axes = plt.subplots(n_chroms, 1, figsize=figsize, 
                                sharex=False, sharey=True)
        
        # Adjust subplot spacing
        plt.subplots_adjust(hspace=0.3)
        
        if n_chroms == 1:
            axes = [axes]
        
        # Create heatmap for each chromosome
        for idx, chromosome in enumerate(sorted(self.parser.chromosomes)):
            ax = axes[idx]
            
            try:
                # Get chromosome data
                chrom_data = self.parser.get_chromosome_data_for_samples(
                    chromosome, samples, window_size
                )
                
                if not chrom_data:
                    ax.text(0.5, 0.5, f'No data for {chromosome}', 
                           ha='center', va='center', transform=ax.transAxes)
                    continue
                
                # Create heatmap matrix
                chrom_length = self.parser.chromosome_lengths[chromosome]
                heatmap_matrix, window_positions, sample_labels = self._prepare_heatmap_matrix(
                    chrom_data, samples, chrom_length, window_size
                )
                
                if heatmap_matrix.size > 0:
                    # Draw heatmap
                    im = ax.imshow(heatmap_matrix, cmap=self.colors, aspect='auto',
                                  vmin=0, vmax=1, interpolation='nearest')
                    
                    # Set labels
                    self._set_axis_labels_compact(ax, window_positions, sample_labels, 
                                                 chrom_length, chromosome)
                    
                    # Remove grid lines, add chromosome borders
                    # self._add_grid_lines(ax, len(sample_labels), len(window_positions))
                    
                    # Add chromosome borders
                    for spine in ax.spines.values():
                        spine.set_visible(True)
                        spine.set_linewidth(1.5)
                        spine.set_edgecolor('black')
                else:
                    ax.text(0.5, 0.5, f'Insufficient data for {chromosome}',
                           ha='center', va='center', transform=ax.transAxes)
                    
            except Exception as e:
                print(f"Error processing chromosome {chromosome}: {e}")
                ax.text(0.5, 0.5, f'Error processing {chromosome}',
                       ha='center', va='center', transform=ax.transAxes)
        
        # Set overall title
        fig.suptitle(f'SNP Match Rate Comparison Across All Chromosomes\n'
                    f'Window Size: {window_size/1000:.0f} kb | '
                    f'Samples: {len(samples)}',
                    fontsize=18, fontweight='bold', y=0.98)
        
        # Add shared colorbar
        sm = plt.cm.ScalarMappable(cmap=self.colors, 
                                  norm=plt.Normalize(vmin=0, vmax=1))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes, shrink=0.5, aspect=30, pad=0.02)
        cbar.set_label('SNP Match Rate', fontsize=14, fontweight='bold')
        cbar.ax.tick_params(labelsize=14)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save figure
        if save_path:
            # Save PNG format
            plt.savefig(save_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            print(f"All chromosome comparison chart saved to: {save_path}")
        
        return fig
    
    def _prepare_heatmap_matrix(self, chrom_data: Dict, samples: List[str], 
                              chrom_length: int, window_size: int) -> Tuple[np.ndarray, List[int], List[str]]:
        """
        Prepare heatmap matrix data
        
        Args:
            chrom_data: Chromosome data dictionary
            samples: Sample list
            chrom_length: Chromosome length
            window_size: Window size
            
        Returns:
            (heatmap matrix, window position list, sample label list)
        """
        # Calculate number of windows
        num_windows = int(np.ceil(chrom_length / window_size))
        
        # Create standardized window positions
        standard_positions = [(i * window_size + window_size // 2) 
                            for i in range(num_windows)]
        
        # Initialize matrix
        valid_samples = []
        matrix_rows = []
        
        for sample in samples:
            if sample in chrom_data:
                positions, rates = chrom_data[sample]
                if positions and rates:
                    # Map data to standard windows
                    standard_rates = self._map_to_standard_windows(
                        positions, rates, standard_positions, window_size
                    )
                    matrix_rows.append(standard_rates)
                    valid_samples.append(sample)
        
        if not matrix_rows:
            return np.array([]), [], []
        
        # Convert to numpy array
        heatmap_matrix = np.array(matrix_rows)
        
        return heatmap_matrix, standard_positions, valid_samples
    
    def _map_to_standard_windows(self, positions: List[int], rates: List[float],
                               standard_positions: List[int], window_size: int) -> List[float]:
        """
        Map data to standard windows
        
        Args:
            positions: Original position list
            rates: Original match rate list
            standard_positions: Standard window position list
            window_size: Window size
            
        Returns:
            Standardized match rate list
        """
        standard_rates = []
        
        for std_pos in standard_positions:
            # Find the closest window
            closest_idx = None
            min_distance = float('inf')
            
            for i, pos in enumerate(positions):
                distance = abs(pos - std_pos)
                if distance < min_distance and distance < window_size:
                    min_distance = distance
                    closest_idx = i
            
            if closest_idx is not None:
                standard_rates.append(rates[closest_idx])
            else:
                standard_rates.append(0.0)  # Set windows without data to 0
        
        return standard_rates
    
    def _set_axis_labels(self, ax, window_positions: List[int], 
                        sample_labels: List[str], chrom_length: int):
        """
        Set axis labels
        """
        # Y-axis (samples)
        ax.set_yticks(range(len(sample_labels)))
        ax.set_yticklabels(sample_labels, fontsize=13)
        ax.set_ylabel('Samples', fontsize=15, fontweight='bold')
        
        # X-axis (genomic position) - using 5/10/20/50 as units
        max_mb = chrom_length / 1e6
        # Calculate appropriate interval, at least in units of 50
        if max_mb <= 50:
            tick_interval = 5
        elif max_mb <= 100:
            tick_interval = 10
        elif max_mb <= 200:
            tick_interval = 20
        else:
            tick_interval = 50
        
        # Generate multiples of 10 ticks
        tick_values_mb = np.arange(0, max_mb + tick_interval, tick_interval)
        # Convert to window indices
        tick_indices = []
        tick_labels = []
        for val_mb in tick_values_mb:
            if val_mb <= max_mb:
                # Find the closest window position
                target_pos = val_mb * 1e6
                if len(window_positions) > 0:
                    closest_idx = np.argmin([abs(pos - target_pos) for pos in window_positions])
                    if closest_idx < len(window_positions):
                        tick_indices.append(closest_idx)
                        tick_labels.append(f"{int(val_mb)}")
        
        if tick_indices:
            ax.set_xticks(tick_indices)
            ax.set_xticklabels(tick_labels, fontsize=13, rotation=0)
        
        ax.set_xlabel('Genomic Position (Mb)', fontsize=15, fontweight='bold')
    
    def _set_axis_labels_compact(self, ax, window_positions: List[int],
                               sample_labels: List[str], chrom_length: int, chromosome: str):
        """
        Set compact axis labels
        """
        # Y-axis label (only show chromosome name)
        chrom_display = chromosome
        ax.set_ylabel(chrom_display, fontsize=13, fontweight='bold')
        
        # Hide y-axis ticks
        ax.set_yticks([])
        
        # X-axis
        n_ticks = min(5, len(window_positions))
        if n_ticks > 0:
            tick_indices = np.linspace(0, len(window_positions)-1, n_ticks, dtype=int)
            ax.set_xticks(tick_indices)
            tick_labels = [f"{window_positions[i]/1e6:.0f}" for i in tick_indices]
            ax.set_xticklabels(tick_labels, fontsize=9)
    
    def _add_grid_lines(self, ax, n_samples: int, n_windows: int):
        """
        Add grid lines
        """
        # Add separation lines between samples
        for i in range(1, n_samples):
            ax.axhline(y=i-0.5, color='white', linewidth=1, alpha=0.7)
        
        # Add separation lines between windows (sparse display)
        step = max(1, n_windows // 20)
        for i in range(step, n_windows, step):
            ax.axvline(x=i-0.5, color='white', linewidth=0.5, alpha=0.5)


if __name__ == "__main__":
    # Test code
    parser = SNPDataParser("input.txt")
    parser.load_data()
    
    visualizer = SampleComparisonHeatmapVisualizer(parser)
    
    # Generate sample comparison heatmap for the first chromosome
    chromosome = parser.chromosomes[0]
    samples = parser.samples[:8]  # Select first 8 samples
    
    fig = visualizer.create_chromosome_comparison_heatmap(
        chromosome,
        samples=samples,
        window_size=50000,
        save_path=f"sample_comparison_{chromosome}.png"
    )
    plt.show()