#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SNP Data Parser Module
Used for reading and processing SNP comparison data, calculating window match rates
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

class SNPDataParser:
    """
    SNP Data Parser
    """
    
    def __init__(self, input_file: str):
        """
        Initialize SNP data parser
        
        Args:
            input_file: Input file path
        """
        self.input_file = input_file
        self.data = None
        self.chromosomes = []
        self.samples = []
        self.chromosome_lengths = {}
        
    def load_data(self) -> pd.DataFrame:
        """
        Load SNP data
        
        Returns:
            Loaded data DataFrame
        """
        print("Loading SNP data...")
        
        try:
            # Read data, support compressed files
            if self.input_file.endswith('.gz'):
                import gzip
                # For VCF files, skip comment lines
                with gzip.open(self.input_file, 'rt') as f:
                    lines = []
                    for line in f:
                        if not line.startswith('##'):
                            lines.append(line)
                    
                # Write data to temporary string and read with pandas
                from io import StringIO
                data_str = ''.join(lines)
                self.data = pd.read_csv(StringIO(data_str), sep='\t', low_memory=False)
            else:
                # For regular files, skip comment lines
                with open(self.input_file, 'r') as f:
                    lines = []
                    for line in f:
                        if not line.startswith('##'):
                            lines.append(line)
                    
                from io import StringIO
                data_str = ''.join(lines)
                self.data = pd.read_csv(StringIO(data_str), sep='\t', low_memory=False)
                
        except Exception as e:
            raise ValueError(f"Failed to read file: {str(e)}")
        
        # Get sample names (excluding first two columns CHROM and POS)
        self.samples = list(self.data.columns[2:])
        
        # Get chromosome list
        self.chromosomes = sorted(self.data['#CHROM'].unique())
        
        # Calculate maximum position for each chromosome (chromosome length)
        for chrom in self.chromosomes:
            chrom_data = self.data[self.data['#CHROM'] == chrom]
            self.chromosome_lengths[chrom] = chrom_data['POS'].max()
        
        print(f"Data loading completed:")
        print(f"  - Number of chromosomes: {len(self.chromosomes)}")
        print(f"  - Number of samples: {len(self.samples)}")
        print(f"  - Total SNP sites: {len(self.data)}")
        
        return self.data
    
    def calculate_window_match_rate(self, chromosome: str, sample: str, 
                                  window_size: int = 50000) -> Tuple[List[int], List[float]]:
        """
        Calculate window match rates for specified chromosome and sample
        
        Args:
            chromosome: Chromosome name
            sample: Sample name
            window_size: Window size (default 500kb)
            
        Returns:
            (window position list, match rate list)
        """
        # Get data for specified chromosome
        chrom_data = self.data[self.data['#CHROM'] == chromosome].copy()
        
        if len(chrom_data) == 0:
            return [], []
        
        # Get chromosome length
        chrom_length = self.chromosome_lengths[chromosome]
        
        # Calculate number of windows
        num_windows = int(np.ceil(chrom_length / window_size))
        
        window_positions = []
        match_rates = []
        
        for i in range(num_windows):
            start_pos = i * window_size
            end_pos = min((i + 1) * window_size, chrom_length)
            
            # Get data within window
            window_data = chrom_data[
                (chrom_data['POS'] >= start_pos) & 
                (chrom_data['POS'] < end_pos)
            ]
            
            if len(window_data) == 0:
                match_rates.append(0.0)
            else:
                # Calculate match rate
                sample_values = window_data[sample].astype(str)
                
                # Count 1s and 0s (ignore '/')
                valid_values = sample_values[sample_values.isin(['0', '1'])]
                
                if len(valid_values) == 0:
                    match_rates.append(0.0)
                else:
                    match_count = (valid_values == '1').sum()
                    match_rate = match_count / len(valid_values)
                    match_rates.append(match_rate)
            
            window_positions.append(start_pos + window_size // 2)  # Window center position
        
        return window_positions, match_rates
    
    def get_chromosome_data_for_samples(self, chromosome: str, samples: List[str], 
                                      window_size: int = 50000) -> Dict[str, Tuple[List[int], List[float]]]:
        """
        Get window match rate data for multiple samples on the same chromosome
        
        Args:
            chromosome: Chromosome name
            samples: List of sample names
            window_size: Window size
            
        Returns:
            Dictionary with sample names as keys and (window positions, match rates) tuples as values
        """
        result = {}
        
        for sample in samples:
            if sample in self.samples:
                positions, rates = self.calculate_window_match_rate(
                    chromosome, sample, window_size
                )
                result[sample] = (positions, rates)
            else:
                print(f"Warning: Sample {sample} does not exist in data")
        
        return result
    
    def get_sample_all_chromosomes_data(self, sample: str, 
                                      window_size: int = 50000) -> Dict[str, Tuple[List[int], List[float]]]:
        """
        Get window match rate data for a single sample across all chromosomes
        
        Args:
            sample: Sample name
            window_size: Window size
            
        Returns:
            Dictionary with chromosome names as keys and (window positions, match rates) tuples as values
        """
        result = {}
        
        if sample not in self.samples:
            print(f"Error: Sample {sample} does not exist in data")
            return result
        
        for chromosome in self.chromosomes:
            positions, rates = self.calculate_window_match_rate(
                chromosome, sample, window_size
            )
            result[chromosome] = (positions, rates)
        
        return result
    
    def get_summary_statistics(self) -> Dict:
        """
        Get data summary statistics
        
        Returns:
            Dictionary containing statistical information
        """
        stats = {
            'chromosomes': self.chromosomes,
            'samples': self.samples,
            'chromosome_lengths': self.chromosome_lengths,
            'total_snps': len(self.data),
            'snps_per_chromosome': {}
        }
        
        for chrom in self.chromosomes:
            chrom_data = self.data[self.data['#CHROM'] == chrom]
            stats['snps_per_chromosome'][chrom] = len(chrom_data)
        
        return stats


if __name__ == "__main__":
    # Test code
    parser = SNPDataParser("input.txt")
    data = parser.load_data()
    
    # Print summary statistics
    stats = parser.get_summary_statistics()
    print("\n=== Data Summary ===")
    print(f"Chromosomes: {stats['chromosomes']}")
    print(f"Samples: {stats['samples']}")
    print(f"Chromosome lengths: {stats['chromosome_lengths']}")
    
    # Test window match rate calculation
    sample = stats['samples'][0]
    chromosome = stats['chromosomes'][0]
    positions, rates = parser.calculate_window_match_rate(chromosome, sample)
    print(f"\nFirst 5 window match rates for sample {sample} on chromosome {chromosome}:")
    for i in range(min(5, len(rates))):
        print(f"  Position {positions[i]:,}: {rates[i]:.3f}")