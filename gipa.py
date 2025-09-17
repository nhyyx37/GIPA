#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GIPA: Genomic Identity and Parentage Authenticator
A comprehensive tool for genomic identity verification and parentage analysis

Features:
1. Identity Analysis: Compare query sample with reference samples to find genetic matches
2. Parentage Analysis: Identify potential parent combinations for hybrid samples
3. Sliding Window Correction: Advanced filtering to reduce noise and improve accuracy
4. Multi-threading Support: Parallel processing for large datasets
5. Visualization: Generate heatmaps for result interpretation
6. Flexible Input: Support for VCF files with customizable parameters

Usage Examples:
# Identity analysis
python gipa.py --vcf data.vcf --sample query_sample --refs reference_list.txt --out results

# Parentage analysis
python gipa.py --vcf data.vcf --sample hybrid_sample --refs parent_list.txt --out results --find_parents

# With visualization
python gipa.py --vcf data.vcf --sample query_sample --refs reference_list.txt --out results --generate-heatmaps

# Chromosome-specific analysis
python gipa.py --vcf data.vcf --sample query_sample --refs reference_list.txt --out results --chr Chr01
"""

import argparse
import sys
import gzip
import threading
import os
import warnings
from collections import defaultdict, Counter
from itertools import combinations
from typing import List, Optional

# Optional dependencies
pd = None
pysam = None
np = None

# Try importing numpy for sample classification
try:
    import numpy as np
except ImportError:
    print("Warning: numpy not available, sample classification features will be limited")
    np = None

# Try importing visualization modules
try:
    import matplotlib.pyplot as plt
    from snp_data_parser import SNPDataParser
    from chromosome_heatmap import ChromosomeHeatmapVisualizer
    from sample_comparison_heatmap import SampleComparisonHeatmapVisualizer
    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False
    print("Warning: Visualization modules not available, heatmap generation will be skipped")

# Suppress warnings
warnings.filterwarnings('ignore')


def calculate_sample_heterozygosity(vcf_file, samples):
    """
    Calculate heterozygosity rate for each sample
    
    Args:
        vcf_file: Path to VCF file
        samples: List of sample names
    
    Returns:
        Dictionary with sample names as keys and heterozygosity rates as values
    """
    het_counts = {sample: 0 for sample in samples}
    total_counts = {sample: 0 for sample in samples}
    
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    
    with open_func(vcf_file, 'rt') as f:
        sample_indices = {}
        
        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                # Header line
                header = line.strip().split('\t')
                for sample in samples:
                    if sample in header:
                        sample_indices[sample] = header.index(sample)
                    else:
                        print(f"Warning: Sample {sample} not found in VCF header")
                continue
            
            # Data line
            fields = line.strip().split('\t')
            
            for sample in samples:
                if sample not in sample_indices:
                    continue
                    
                idx = sample_indices[sample]
                if idx < len(fields):
                    gt_field = fields[idx]
                    gt = gt_field.split(':')[0]  # Get genotype part
                    
                    # Parse genotype
                    if '/' in gt:
                        alleles = gt.split('/')
                    elif '|' in gt:
                        alleles = gt.split('|')
                    else:
                        continue
                    
                    # Skip missing genotypes
                    if '.' in alleles or len(alleles) != 2:
                        continue
                    
                    total_counts[sample] += 1
                    
                    # Check if heterozygous
                    if alleles[0] != alleles[1]:
                        het_counts[sample] += 1
    
    # Calculate heterozygosity rates
    het_rates = {}
    for sample in samples:
        if total_counts[sample] > 0:
            het_rates[sample] = het_counts[sample] / total_counts[sample]
        else:
            het_rates[sample] = 0.0
            print(f"Warning: No valid genotypes found for sample {sample}")
    
    return {
        'heterozygosity_rates': het_rates,
        'total_snps': total_counts,
        'heterozygous_snps': het_counts
    }


def classify_samples_by_heterozygosity(het_results, min_samples_for_auto=3, fallback_thresholds=(0.05, 0.10)):
    """
    Classify samples as inbred lines, hybrids, or unknown based on heterozygosity
    
    Args:
        het_results: Results from calculate_sample_heterozygosity
        min_samples_for_auto: Minimum samples needed for automatic threshold detection
        fallback_thresholds: (inbred_threshold, hybrid_threshold) for manual classification
    
    Returns:
        Dictionary with classification results
    """
    if np is None:
        print("Warning: numpy not available, using fixed thresholds for classification")
        return classify_with_fixed_thresholds(het_results, fallback_thresholds[0], fallback_thresholds[1])
    
    het_rates = het_results['heterozygosity_rates']
    rates_list = list(het_rates.values())
    
    if len(rates_list) < min_samples_for_auto:
        print(f"Warning: Only {len(rates_list)} samples available, using fixed thresholds")
        return classify_with_fixed_thresholds(het_results, fallback_thresholds[0], fallback_thresholds[1])
    
    # Automatic threshold detection using clustering
    rates_array = np.array(rates_list)
    
    # Simple threshold detection: use percentiles
    inbred_threshold = np.percentile(rates_array, 25)
    hybrid_threshold = np.percentile(rates_array, 75)
    
    # Ensure reasonable thresholds
    inbred_threshold = max(inbred_threshold, 0.02)
    hybrid_threshold = min(hybrid_threshold, 0.20)
    
    if hybrid_threshold <= inbred_threshold:
        hybrid_threshold = inbred_threshold + 0.05
    
    print(f"Auto-detected thresholds: Inbred <= {inbred_threshold:.3f}, Hybrid >= {hybrid_threshold:.3f}")
    
    # Classify samples
    classifications = {}
    for sample, rate in het_rates.items():
        if rate <= inbred_threshold:
            classifications[sample] = 'inbred_line'
        elif rate >= hybrid_threshold:
            classifications[sample] = 'hybrid'
        else:
            classifications[sample] = 'intermediate'
    
    return {
        'classifications': classifications,
        'thresholds': {'inbred': inbred_threshold, 'hybrid': hybrid_threshold},
        'heterozygosity_rates': het_rates
    }


def classify_with_fixed_thresholds(het_results, inbred_thresh=0.05, hybrid_thresh=0.10):
    """
    Classify samples using fixed thresholds
    """
    het_rates = het_results['heterozygosity_rates']
    classifications = {}
    
    for sample, rate in het_rates.items():
        if rate <= inbred_thresh:
            classifications[sample] = 'inbred_line'
        elif rate >= hybrid_thresh:
            classifications[sample] = 'hybrid'
        else:
            classifications[sample] = 'intermediate'
    
    return {
        'classifications': classifications,
        'thresholds': {'inbred': inbred_thresh, 'hybrid': hybrid_thresh},
        'heterozygosity_rates': het_rates
    }


def apply_genotype_matching_rules(p1_gt, p2_gt, hybrid_gt):
    """
    Apply genotype matching rules for parentage analysis
    
    Args:
        p1_gt: Parent 1 genotype
        p2_gt: Parent 2 genotype  
        hybrid_gt: Hybrid genotype
    
    Returns:
        Tuple (is_informative, is_match, rule_applied)
    """
    # Skip missing genotypes
    if '.' in [p1_gt, p2_gt, hybrid_gt]:
        return False, False, "missing_data"
    
    # Parse genotypes
    def parse_gt(gt):
        if '/' in gt:
            return sorted(gt.split('/'))
        elif '|' in gt:
            return sorted(gt.split('|'))
        else:
            return [gt, gt]
    
    p1_alleles = parse_gt(p1_gt)
    p2_alleles = parse_gt(p2_gt)
    hybrid_alleles = parse_gt(hybrid_gt)
    
    # Rule 1: Both parents homozygous, same allele
    if (p1_alleles[0] == p1_alleles[1] and 
        p2_alleles[0] == p2_alleles[1] and 
        p1_alleles[0] == p2_alleles[0]):
        # Hybrid should be homozygous for the same allele
        expected = [p1_alleles[0], p1_alleles[0]]
        is_match = hybrid_alleles == expected
        return True, is_match, "both_homo_same"
    
    # Rule 2: Both parents homozygous, different alleles
    if (p1_alleles[0] == p1_alleles[1] and 
        p2_alleles[0] == p2_alleles[1] and 
        p1_alleles[0] != p2_alleles[0]):
        # Hybrid should be heterozygous with both parental alleles
        expected = sorted([p1_alleles[0], p2_alleles[0]])
        is_match = hybrid_alleles == expected
        return True, is_match, "both_homo_diff"
    
    # Rule 3: One parent homozygous, other heterozygous
    if (p1_alleles[0] == p1_alleles[1] and p2_alleles[0] != p2_alleles[1]):
        # P1 homo, P2 hetero
        homo_allele = p1_alleles[0]
        hetero_alleles = p2_alleles
        
        # Hybrid should have the homozygous allele and one of the heterozygous alleles
        possible_hybrids = [
            sorted([homo_allele, hetero_alleles[0]]),
            sorted([homo_allele, hetero_alleles[1]])
        ]
        is_match = hybrid_alleles in possible_hybrids
        return True, is_match, "p1_homo_p2_hetero"
    
    if (p2_alleles[0] == p2_alleles[1] and p1_alleles[0] != p1_alleles[1]):
        # P2 homo, P1 hetero
        homo_allele = p2_alleles[0]
        hetero_alleles = p1_alleles
        
        # Hybrid should have the homozygous allele and one of the heterozygous alleles
        possible_hybrids = [
            sorted([homo_allele, hetero_alleles[0]]),
            sorted([homo_allele, hetero_alleles[1]])
        ]
        is_match = hybrid_alleles in possible_hybrids
        return True, is_match, "p2_homo_p1_hetero"
    
    # Rule 4: Both parents heterozygous
    if (p1_alleles[0] != p1_alleles[1] and p2_alleles[0] != p2_alleles[1]):
        # Collect all possible alleles from both parents
        all_parent_alleles = set(p1_alleles + p2_alleles)
        
        # Hybrid alleles should be subset of parental alleles
        hybrid_allele_set = set(hybrid_alleles)
        is_match = hybrid_allele_set.issubset(all_parent_alleles)
        return True, is_match, "both_hetero"
    
    # Rule 5: Both parents homozygous for same allele (redundant with Rule 1, but kept for clarity)
    # This case is already handled by Rule 1
    
    # If no rule applies, consider it non-informative
    return False, False, "no_rule_applied"


class SNPVisualizationTool:
    """
    Tool for generating SNP visualization heatmaps
    """
    
    def __init__(self, input_file: str):
        """
        Initialize visualization tool
        
        Args:
            input_file: Path to the matrix file (identity or parental)
        """
        self.input_file = input_file
        self.data_parser = None
        self.chromosome_visualizer = None
        self.sample_visualizer = None
        
        if not VISUALIZATION_AVAILABLE:
            raise ImportError("Visualization modules not available")
    
    def initialize(self):
        """
        Initialize visualization components
        """
        try:
            self.data_parser = SNPDataParser(self.input_file)
            self.data_parser.load_data()  # Load the data
            self.chromosome_visualizer = ChromosomeHeatmapVisualizer(self.data_parser)
            self.sample_visualizer = SampleComparisonHeatmapVisualizer(self.data_parser)
            print(f"Visualization tool initialized with {self.input_file}")
        except Exception as e:
            print(f"Failed to initialize visualization tool: {e}")
            raise
    
    def generate_single_sample_heatmap(self, sample: str, window_size: int = 50000,
                                     output_dir: str = "output", 
                                     show_plot: bool = False) -> str:
        """
        Generate heatmap for a single sample
        
        Args:
            sample: Sample name
            window_size: Window size in bp
            output_dir: Output directory
            show_plot: Whether to display the plot
        
        Returns:
            Path to generated image file
        """
        if self.chromosome_visualizer is None:
            raise RuntimeError("Visualization tool not initialized")
        
        try:
            save_path = f"{output_dir}/chromosome_heatmap_{sample}.png"
            fig = self.chromosome_visualizer.create_single_sample_heatmap(
                sample=sample,
                window_size=window_size,
                save_path=save_path
            )
            print(f"Single sample heatmap generated: {save_path}")
            return save_path
        except Exception as e:
            print(f"Failed to generate single sample heatmap: {e}")
            return None
    
    def generate_chromosome_comparison_heatmap(self, chromosome: str, 
                                             samples: Optional[List[str]] = None,
                                             window_size: int = 50000,
                                             output_dir: str = "output",
                                             show_plot: bool = False) -> str:
        """
        Generate comparison heatmap for multiple samples on a chromosome
        
        Args:
            chromosome: Chromosome name
            samples: List of sample names (if None, use all samples)
            window_size: Window size in bp
            output_dir: Output directory
            show_plot: Whether to display the plot
        
        Returns:
            Path to generated image file
        """
        if self.sample_visualizer is None:
            raise RuntimeError("Visualization tool not initialized")
        
        try:
            save_path = f"{output_dir}/sample_comparison_{chromosome}.png"
            fig = self.sample_visualizer.create_chromosome_comparison_heatmap(
                chromosome=chromosome,
                samples=samples,
                window_size=window_size,
                save_path=save_path
            )
            print(f"Chromosome comparison heatmap generated: {save_path}")
            return save_path
        except Exception as e:
            print(f"Failed to generate chromosome comparison heatmap: {e}")
            return None


def load_samples(sample_list_file):
    """
    Load sample names from file
    """
    samples = []
    with open(sample_list_file, 'r') as f:
        for line in f:
            sample = line.strip()
            if sample:
                samples.append(sample)
    return samples


def verify_samples_in_vcf(vcf_file, samples_to_check):
    """
    Verify that samples exist in VCF file
    """
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    
    with open_func(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                header = line.strip().split('\t')
                vcf_samples = header[9:]  # Samples start from column 10
                
                missing_samples = []
                for sample in samples_to_check:
                    if sample not in vcf_samples:
                        missing_samples.append(sample)
                
                if missing_samples:
                    print(f"Error: The following samples were not found in VCF: {missing_samples}")
                    print(f"Available samples in VCF: {vcf_samples[:10]}{'...' if len(vcf_samples) > 10 else ''}")
                    return False
                
                return True
        
        print("Error: No header line found in VCF file")
        return False


def get_window_representative(window, prev_repr=None):
    """
    Get representative value for a window using majority voting
    
    Args:
        window: List of values in the window
        prev_repr: Previous representative value (for tie-breaking)
    
    Returns:
        Representative value for the window
    """
    if not window:
        return prev_repr if prev_repr is not None else 0
    
    # Count occurrences of each value
    counts = Counter(window)
    
    # Find the most common value
    most_common = counts.most_common()
    
    if len(most_common) == 1:
        return most_common[0][0]
    
    # Handle ties
    max_count = most_common[0][1]
    tied_values = [value for value, count in most_common if count == max_count]
    
    # If previous representative is among tied values, prefer it
    if prev_repr is not None and prev_repr in tied_values:
        return prev_repr
    
    # Otherwise, return the first tied value
    return tied_values[0]


def apply_sliding_window_correction(matrix_data, window_size=5, filter_times=2):
    """
    Apply sliding window correction to matrix data
    
    Args:
        matrix_data: List of (chrom, pos, values_dict) tuples
        window_size: Size of sliding window
        filter_times: Number of filtering iterations
    
    Returns:
        Corrected matrix data in the same format
    """
    if not matrix_data:
        return matrix_data
    
    print(f"Applying sliding window correction (window_size={window_size}, filter_times={filter_times})")
    
    # Group data by chromosome
    chrom_data = defaultdict(list)
    for chrom, pos, values in matrix_data:
        chrom_data[chrom].append((pos, values))
    
    # Sort by position within each chromosome
    for chrom in chrom_data:
        chrom_data[chrom].sort(key=lambda x: x[0])
    
    corrected_data = []
    
    for chrom, positions_values in chrom_data.items():
        if len(positions_values) < window_size:
            # If not enough data points, keep original
            for pos, values in positions_values:
                corrected_data.append((chrom, pos, values))
            continue
        
        # Get all sample names from first entry
        sample_names = list(positions_values[0][1].keys())
        
        # Apply filtering for each sample
        for sample in sample_names:
            # Extract values for this sample
            sample_values = [values[sample] for pos, values in positions_values]
            
            # Apply multiple rounds of filtering
            filtered_values = sample_values[:]
            
            for filter_round in range(filter_times):
                new_filtered_values = []
                
                for i in range(len(filtered_values)):
                    # Define window boundaries
                    start_idx = max(0, i - window_size // 2)
                    end_idx = min(len(filtered_values), i + window_size // 2 + 1)
                    
                    # Get window values
                    window = filtered_values[start_idx:end_idx]
                    
                    # Get representative value
                    prev_repr = new_filtered_values[-1] if new_filtered_values else None
                    repr_value = get_window_representative(window, prev_repr)
                    
                    new_filtered_values.append(repr_value)
                
                filtered_values = new_filtered_values
            
            # Update values for this sample
            for i, (pos, values) in enumerate(positions_values):
                values[sample] = filtered_values[i]
        
        # Add corrected data for this chromosome
        for pos, values in positions_values:
            corrected_data.append((chrom, pos, values))
    
    # Sort final result by chromosome and position
    def sort_key(item):
        chrom, pos, values = item
        # Extract numeric part of chromosome for sorting
        chrom_num = ''.join(filter(str.isdigit, chrom))
        chrom_num = int(chrom_num) if chrom_num else 0
        return (chrom_num, pos)
    
    corrected_data.sort(key=sort_key)
    
    print(f"Sliding window correction completed. Processed {len(corrected_data)} SNP positions.")
    return corrected_data


def parse_genotype(gt_str):
    """
    Parse genotype string and return standardized format
    
    Args:
        gt_str: Genotype string (e.g., "0/1", "1|0", "./.", etc.)
    
    Returns:
        Tuple of (allele1, allele2) or None if missing
    """
    if not gt_str or gt_str == './.' or gt_str == '.|.':
        return None
    
    # Split by / or |
    if '/' in gt_str:
        alleles = gt_str.split('/')
    elif '|' in gt_str:
        alleles = gt_str.split('|')
    else:
        return None
    
    if len(alleles) != 2 or '.' in alleles:
        return None
    
    return (alleles[0], alleles[1])


def _process_vcf_chunk(vcf_file, query_sample, ref_sample_chunk, sample_indices, target_chr, thread_id):
    """
    Process a chunk of reference samples for comparison with query sample
    
    Args:
        vcf_file: Path to VCF file
        query_sample: Query sample name
        ref_sample_chunk: List of reference sample names to process
        sample_indices: Dictionary mapping sample names to column indices
        target_chr: Target chromosome (None for all chromosomes)
        thread_id: Thread identifier for logging
    
    Returns:
        Dictionary with comparison results
    """
    results = {}
    
    # Initialize results for each reference sample
    for ref_sample in ref_sample_chunk:
        results[ref_sample] = {
            'positions': [],
            'values': []
        }
    
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    
    with open_func(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            
            # Skip if target chromosome specified and doesn't match
            if target_chr and chrom != target_chr:
                continue
            
            # Get query sample genotype
            query_idx = sample_indices.get(query_sample)
            if query_idx is None or query_idx >= len(fields):
                continue
            
            query_gt_field = fields[query_idx]
            query_gt = query_gt_field.split(':')[0]
            query_parsed = parse_genotype(query_gt)
            
            if query_parsed is None:
                continue
            
            # Compare with each reference sample in this chunk
            for ref_sample in ref_sample_chunk:
                ref_idx = sample_indices.get(ref_sample)
                if ref_idx is None or ref_idx >= len(fields):
                    continue
                
                ref_gt_field = fields[ref_idx]
                ref_gt = ref_gt_field.split(':')[0]
                ref_parsed = parse_genotype(ref_gt)
                
                if ref_parsed is None:
                    continue
                
                # Calculate identity (1 if match, 0 if no match)
                identity = 1 if sorted(query_parsed) == sorted(ref_parsed) else 0
                
                results[ref_sample]['positions'].append((chrom, pos))
                results[ref_sample]['values'].append(identity)
    
    print(f"Thread {thread_id}: Processed {len(ref_sample_chunk)} reference samples")
    return results


def _generate_comparison_matrix_multithreaded(vcf_file, query_sample, ref_samples, sample_indices, target_chr, num_threads):
    """
    Generate comparison matrix using multiple threads
    
    Args:
        vcf_file: Path to VCF file
        query_sample: Query sample name
        ref_samples: List of reference sample names
        sample_indices: Dictionary mapping sample names to column indices
        target_chr: Target chromosome (None for all chromosomes)
        num_threads: Number of threads to use
    
    Returns:
        List of (chrom, pos, values_dict) tuples
    """
    # Split reference samples into chunks for threading
    chunk_size = max(1, len(ref_samples) // num_threads)
    ref_chunks = [ref_samples[i:i + chunk_size] for i in range(0, len(ref_samples), chunk_size)]
    
    print(f"Using {len(ref_chunks)} threads to process {len(ref_samples)} reference samples")
    
    # Process chunks in parallel
    threads = []
    thread_results = []
    
    for i, chunk in enumerate(ref_chunks):
        thread_result = {}
        thread_results.append(thread_result)
        
        thread = threading.Thread(
            target=lambda chunk=chunk, result=thread_result, tid=i: 
                result.update(_process_vcf_chunk(vcf_file, query_sample, chunk, sample_indices, target_chr, tid))
        )
        threads.append(thread)
        thread.start()
    
    # Wait for all threads to complete
    for thread in threads:
        thread.join()
    
    # Merge results from all threads
    print("Merging results from all threads...")
    
    # Collect all positions
    all_positions = set()
    for thread_result in thread_results:
        for ref_sample, data in thread_result.items():
            all_positions.update(data['positions'])
    
    # Sort positions
    sorted_positions = sorted(all_positions, key=lambda x: (x[0], x[1]))
    
    # Build final matrix
    matrix_data = []
    for chrom, pos in sorted_positions:
        values_dict = {}
        
        for thread_result in thread_results:
            for ref_sample, data in thread_result.items():
                if (chrom, pos) in data['positions']:
                    pos_idx = data['positions'].index((chrom, pos))
                    values_dict[ref_sample] = data['values'][pos_idx]
                else:
                    values_dict[ref_sample] = 0  # Default to no match
        
        matrix_data.append((chrom, pos, values_dict))
    
    return matrix_data


def generate_comparison_matrix(vcf_file, query_sample, ref_samples, output_prefix, target_chr=None, num_threads=1):
    """
    Generate comparison matrix between query sample and reference samples
    
    Args:
        vcf_file: Path to VCF file
        query_sample: Query sample name
        ref_samples: List of reference sample names
        output_prefix: Output file prefix
        target_chr: Target chromosome (None for all chromosomes)
        num_threads: Number of threads to use
    
    Returns:
        List of (chrom, pos, values_dict) tuples
    """
    print(f"Generating comparison matrix for {query_sample} vs {len(ref_samples)} reference samples")
    
    # Verify samples exist in VCF
    all_samples = [query_sample] + ref_samples
    if not verify_samples_in_vcf(vcf_file, all_samples):
        print("Error: Sample verification failed")
        return None
    
    # Get sample indices from VCF header
    sample_indices = {}
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    
    with open_func(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                header = line.strip().split('\t')
                for sample in all_samples:
                    if sample in header:
                        sample_indices[sample] = header.index(sample)
                break
    
    if len(sample_indices) != len(all_samples):
        print("Error: Not all samples found in VCF header")
        return None
    
    # Generate matrix using appropriate method
    if num_threads > 1:
        matrix_data = _generate_comparison_matrix_multithreaded(
            vcf_file, query_sample, ref_samples, sample_indices, target_chr, num_threads
        )
    else:
        # Single-threaded processing
        matrix_data = []
        
        open_func = gzip.open if vcf_file.endswith('.gz') else open
        
        with open_func(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                chrom = fields[0]
                pos = int(fields[1])
                
                # Skip if target chromosome specified and doesn't match
                if target_chr and chrom != target_chr:
                    continue
                
                # Get query sample genotype
                query_idx = sample_indices[query_sample]
                query_gt_field = fields[query_idx]
                query_gt = query_gt_field.split(':')[0]
                query_parsed = parse_genotype(query_gt)
                
                if query_parsed is None:
                    continue
                
                # Compare with each reference sample
                values_dict = {}
                for ref_sample in ref_samples:
                    ref_idx = sample_indices[ref_sample]
                    ref_gt_field = fields[ref_idx]
                    ref_gt = ref_gt_field.split(':')[0]
                    ref_parsed = parse_genotype(ref_gt)
                    
                    if ref_parsed is None:
                        values_dict[ref_sample] = 0
                    else:
                        # Calculate identity (1 if match, 0 if no match)
                        identity = 1 if sorted(query_parsed) == sorted(ref_parsed) else 0
                        values_dict[ref_sample] = identity
                
                matrix_data.append((chrom, pos, values_dict))
    
    # Save raw matrix
    matrix_file = f"{output_prefix}_identity_matrix_raw.txt"
    with open(matrix_file, 'w') as f:
        header = ['#CHROM', 'POS'] + ref_samples
        f.write('\t'.join(header) + '\n')
        
        for chrom, pos, values in matrix_data:
            row_data = [chrom, str(pos)]
            for ref_sample in ref_samples:
                row_data.append(str(values[ref_sample]))
            f.write('\t'.join(row_data) + '\n')
    
    print(f"Raw comparison matrix saved to: {matrix_file}")
    print(f"Matrix contains {len(matrix_data)} SNP positions")
    
    return matrix_data


def calculate_identity_from_matrix(matrix_data, ref_samples, output_prefix):
    """
    Calculate identity statistics from comparison matrix
    
    Args:
        matrix_data: List of (chrom, pos, values_dict) tuples
        ref_samples: List of reference sample names
        output_prefix: Output file prefix
    
    Returns:
        Tuple of (overall_results, chromosome_results)
    """
    print("Calculating identity statistics...")
    
    # Overall statistics
    overall_stats = {}
    for ref_sample in ref_samples:
        overall_stats[ref_sample] = {'matched': 0, 'total': 0}
    
    # Chromosome-specific statistics
    chrom_stats = defaultdict(lambda: {ref_sample: {'matched': 0, 'total': 0} for ref_sample in ref_samples})
    
    # Process matrix data
    for chrom, pos, values in matrix_data:
        for ref_sample in ref_samples:
            value = values.get(ref_sample, 0)
            
            # Update overall stats
            overall_stats[ref_sample]['total'] += 1
            if value == 1:
                overall_stats[ref_sample]['matched'] += 1
            
            # Update chromosome stats
            chrom_stats[chrom][ref_sample]['total'] += 1
            if value == 1:
                chrom_stats[chrom][ref_sample]['matched'] += 1
    
    # Calculate percentages and prepare results
    overall_results = []
    for ref_sample in ref_samples:
        stats = overall_stats[ref_sample]
        if stats['total'] > 0:
            identity_pct = (stats['matched'] / stats['total']) * 100
        else:
            identity_pct = 0.0
        
        overall_results.append({
            'Reference': ref_sample,
            'Identity (%)': identity_pct,
            'Matched_SNPs': stats['matched'],
            'Compared_SNPs': stats['total']
        })
    
    # Sort by identity percentage (descending)
    overall_results.sort(key=lambda x: x['Identity (%)'], reverse=True)
    
    # Prepare chromosome results
    chromosome_results = []
    for chrom in sorted(chrom_stats.keys()):
        for ref_sample in ref_samples:
            stats = chrom_stats[chrom][ref_sample]
            if stats['total'] > 0:
                identity_pct = (stats['matched'] / stats['total']) * 100
            else:
                identity_pct = 0.0
            
            chromosome_results.append({
                'Chromosome': chrom,
                'Reference': ref_sample,
                'Identity (%)': identity_pct,
                'Matched_SNPs': stats['matched'],
                'Compared_SNPs': stats['total']
            })
    
    # Sort chromosome results
    chromosome_results.sort(key=lambda x: (x['Chromosome'], -x['Identity (%)']))
    
    # Save results to files
    overall_file = f"{output_prefix}_identity_results.tsv"
    with open(overall_file, 'w') as f:
        f.write('Reference\tIdentity(%)\tMatched_SNPs\tCompared_SNPs\n')
        for result in overall_results:
            f.write(f"{result['Reference']}\t{result['Identity (%)']:.2f}\t{result['Matched_SNPs']}\t{result['Compared_SNPs']}\n")
    
    chrom_file = f"{output_prefix}_identity_by_chromosome.tsv"
    with open(chrom_file, 'w') as f:
        f.write('Chromosome\tReference\tIdentity(%)\tMatched_SNPs\tCompared_SNPs\n')
        for result in chromosome_results:
            f.write(f"{result['Chromosome']}\t{result['Reference']}\t{result['Identity (%)']:.2f}\t{result['Matched_SNPs']}\t{result['Compared_SNPs']}\n")
    
    print(f"Identity results saved to: {overall_file}")
    print(f"Chromosome-specific results saved to: {chrom_file}")
    
    return overall_results, chromosome_results


def generate_parental_comparison_matrix(vcf_file, hybrid_sample, ref_samples, output_prefix, filter_times=2, window_size=5, target_chr=None, num_threads=1):
    """
    Generate parental comparison matrix for hybrid sample
    
    Args:
        vcf_file: Path to VCF file
        hybrid_sample: Hybrid sample name
        ref_samples: List of potential parent sample names
        output_prefix: Output file prefix
        filter_times: Number of sliding window filter iterations
        window_size: Sliding window size
        target_chr: Target chromosome (None for all chromosomes)
        num_threads: Number of threads to use
    
    Returns:
        Tuple of (overall_results, chromosome_results)
    """
    print(f"Generating parental comparison matrix for {hybrid_sample}")
    print(f"Testing {len(ref_samples)} potential parents in {len(list(combinations(ref_samples, 2)))} combinations")
    
    # Verify samples exist in VCF
    all_samples = [hybrid_sample] + ref_samples
    if not verify_samples_in_vcf(vcf_file, all_samples):
        print("Error: Sample verification failed")
        return None, None
    
    # Get sample indices from VCF header
    sample_indices = {}
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    
    with open_func(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                header = line.strip().split('\t')
                for sample in all_samples:
                    if sample in header:
                        sample_indices[sample] = header.index(sample)
                break
    
    if len(sample_indices) != len(all_samples):
        print("Error: Not all samples found in VCF header")
        return None, None
    
    # Generate all parent combinations
    parent_combinations = list(combinations(ref_samples, 2))
    
    # Process VCF and calculate parental matches
    combination_data = {}
    for p1, p2 in parent_combinations:
        combination_name = f"{p1}_x_{p2}"
        combination_data[combination_name] = {
            'parent1': p1,
            'parent2': p2,
            'positions': [],
            'values': []
        }
    
    print("Processing VCF file for parental analysis...")
    
    with open_func(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            
            # Skip if target chromosome specified and doesn't match
            if target_chr and chrom != target_chr:
                continue
            
            # Get hybrid sample genotype
            hybrid_idx = sample_indices[hybrid_sample]
            hybrid_gt_field = fields[hybrid_idx]
            hybrid_gt = hybrid_gt_field.split(':')[0]
            
            # Process each parent combination
            for p1, p2 in parent_combinations:
                combination_name = f"{p1}_x_{p2}"
                
                # Get parent genotypes
                p1_idx = sample_indices[p1]
                p2_idx = sample_indices[p2]
                
                p1_gt_field = fields[p1_idx]
                p2_gt_field = fields[p2_idx]
                
                p1_gt = p1_gt_field.split(':')[0]
                p2_gt = p2_gt_field.split(':')[0]
                
                # Apply genotype matching rules
                is_informative, is_match, rule = apply_genotype_matching_rules(p1_gt, p2_gt, hybrid_gt)
                
                if is_informative:
                    combination_data[combination_name]['positions'].append((chrom, pos))
                    combination_data[combination_name]['values'].append(1 if is_match else 0)
    
    # Convert to matrix format for sliding window correction
    matrix_data = []
    all_positions = set()
    
    for combo_name, data in combination_data.items():
        all_positions.update(data['positions'])
    
    sorted_positions = sorted(all_positions, key=lambda x: (x[0], x[1]))
    
    for chrom, pos in sorted_positions:
        values_dict = {}
        for combo_name, data in combination_data.items():
            if (chrom, pos) in data['positions']:
                pos_idx = data['positions'].index((chrom, pos))
                values_dict[combo_name] = data['values'][pos_idx]
            else:
                values_dict[combo_name] = 0  # Default to no match
        
        matrix_data.append((chrom, pos, values_dict))
    
    # Apply sliding window correction
    print("Applying sliding window correction to parental matrix...")
    corrected_matrix_data = apply_sliding_window_correction(
        matrix_data, 
        window_size=window_size, 
        filter_times=filter_times
    )
    
    # Save corrected matrix
    matrix_file = f"{output_prefix}_parental_matrix_corrected.txt"
    with open(matrix_file, 'w') as f:
        header = ['#CHROM', 'POS'] + list(combination_data.keys())
        f.write('\t'.join(header) + '\n')
        
        for chrom, pos, values in corrected_matrix_data:
            row_data = [chrom, str(pos)]
            for combo_name in combination_data.keys():
                row_data.append(str(values.get(combo_name, 0)))
            f.write('\t'.join(row_data) + '\n')
    
    print(f"Corrected parental matrix saved to: {matrix_file}")
    
    # Calculate parental statistics
    print("Calculating parental statistics...")
    
    # Overall statistics
    overall_stats = {}
    for combo_name in combination_data.keys():
        overall_stats[combo_name] = {'matched': 0, 'total': 0}
    
    # Chromosome-specific statistics
    chrom_stats = defaultdict(lambda: {combo_name: {'matched': 0, 'total': 0} for combo_name in combination_data.keys()})
    
    # Process corrected matrix data
    for chrom, pos, values in corrected_matrix_data:
        for combo_name in combination_data.keys():
            value = values.get(combo_name, 0)
            
            # Update overall stats
            overall_stats[combo_name]['total'] += 1
            if value == 1:
                overall_stats[combo_name]['matched'] += 1
            
            # Update chromosome stats
            chrom_stats[chrom][combo_name]['total'] += 1
            if value == 1:
                chrom_stats[chrom][combo_name]['matched'] += 1
    
    # Prepare results
    overall_results = []
    for combo_name, data in combination_data.items():
        stats = overall_stats[combo_name]
        if stats['total'] > 0:
            match_pct = (stats['matched'] / stats['total']) * 100
        else:
            match_pct = 0.0
        
        overall_results.append({
            'parent1': data['parent1'],
            'parent2': data['parent2'],
            'match_percentage': match_pct,
            'matched_snps': stats['matched'],
            'informative_snps': stats['total']
        })
    
    # Sort by match percentage (descending)
    overall_results.sort(key=lambda x: x['match_percentage'], reverse=True)
    
    # Prepare chromosome results
    chromosome_results = []
    for chrom in sorted(chrom_stats.keys()):
        for combo_name, data in combination_data.items():
            stats = chrom_stats[chrom][combo_name]
            if stats['total'] > 0:
                match_pct = (stats['matched'] / stats['total']) * 100
            else:
                match_pct = 0.0
            
            chromosome_results.append({
                'chromosome': chrom,
                'parent1': data['parent1'],
                'parent2': data['parent2'],
                'match_percentage': match_pct,
                'matched_snps': stats['matched'],
                'informative_snps': stats['total']
            })
    
    # Sort chromosome results
    chromosome_results.sort(key=lambda x: (x['chromosome'], -x['match_percentage']))
    
    # Save results to files
    overall_file = f"{output_prefix}_parental_results.tsv"
    with open(overall_file, 'w') as f:
        f.write('Parent1\tParent2\tChromosome\tMatch(%)\tInformative_SNPs\tMatched_SNPs\tNote\n')
        
        # Write overall results
        for result in overall_results:
            f.write(f"{result['parent1']}\t{result['parent2']}\tOverall\t{result['match_percentage']:.2f}\t{result['informative_snps']}\t{result['matched_snps']}\tOverall_Statistics\n")
        
        # Write chromosome-specific results
        for result in chromosome_results:
            f.write(f"{result['parent1']}\t{result['parent2']}\t{result['chromosome']}\t{result['match_percentage']:.2f}\t{result['informative_snps']}\t{result['matched_snps']}\tChromosome_Specific\n")
    
    chrom_file = f"{output_prefix}_parental_by_chromosome.tsv"
    with open(chrom_file, 'w') as f:
        f.write('Chromosome\tParent1\tParent2\tMatch(%)\tInformative_SNPs\tMatched_SNPs\n')
        for result in chromosome_results:
            f.write(f"{result['chromosome']}\t{result['parent1']}\t{result['parent2']}\t{result['match_percentage']:.2f}\t{result['informative_snps']}\t{result['matched_snps']}\n")
    
    print(f"Parental results saved to: {overall_file}")
    print(f"Chromosome-specific parental results saved to: {chrom_file}")
    
    return overall_results, chromosome_results


def main():
    parser = argparse.ArgumentParser(
        description="GIPA: Genomic Identity and Parentage Authenticator\nGenomic Identity Verification and Parentage Analysis Tool",
        formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('--vcf', '-v', required=True,
                       help='Path to the VCF file')
    parser.add_argument('--sample', '-s', required=True,
                       help='The query sample name from the VCF')
    parser.add_argument('--refs', '-r', required=True,
                       help='A text file with one reference sample name per line')
    parser.add_argument('--out', '-o', default='output',
                       help='Output prefix for result files (default: output)')
    parser.add_argument('--filter-times', '-ft', type=int, default=2,
                       help='Filter times for sliding window (default: 2)')
    parser.add_argument('--filter-window', '-fw', type=int, default=5,
                       help='Sliding window size (default: 5)')
    parser.add_argument('--threads', '-t', type=int, default=1,
                       help='Number of threads to use (default: 1)')
    parser.add_argument('--chr', '-c', type=str, default=None,
                       help='Specify chromosome to analyze (e.g., Chr01, Chr02)')
    parser.add_argument('--heatmap-window', '-hw', type=int, default=50,
                       help='Heatmap window size in kb (default: 50)')
    parser.add_argument('--find_parents', action='store_true',
                       help='Run parental comparison analysis instead of identity statistics')
    parser.add_argument('--generate-heatmaps', action='store_true',
                       help='Generate heatmap visualizations (requires visualization modules)')
    
    args = parser.parse_args()
    
    # Load reference samples
    ref_samples = load_samples(args.refs)
    print(f"Loaded {len(ref_samples)} reference samples")
    
    # If parental mode is enabled, only perform parental comparison analysis
    if args.find_parents:
        if not args.refs:
            print("Error: Parental comparison analysis requires reference sample file (--refs)")
            sys.exit(1)
    else:
        # Default mode: perform Identity analysis
        # Generate comparison matrix
        matrix_data = generate_comparison_matrix(args.vcf, args.sample, ref_samples, args.out, args.chr, args.threads)
        
        # Check if matrix generation was successful
        if matrix_data is None:
            print("Error: Failed to generate comparison matrix. Please check sample names and VCF file.")
            sys.exit(1)
        
        # Apply sliding window correction
        corrected_matrix_data = apply_sliding_window_correction(
            matrix_data, 
            window_size=args.filter_window, 
            filter_times=args.filter_times
        )
        
        # Save corrected matrix
        corrected_matrix_file = f"{args.out}_identity_matrix_corrected.txt"
        with open(corrected_matrix_file, 'w') as f:
            header = ['#CHROM', 'POS'] + ref_samples
            f.write('\t'.join(header) + '\n')
            
            for chrom, pos, values in corrected_matrix_data:
                row_data = [chrom, str(pos)]
                for ref_sample in ref_samples:
                    row_data.append(str(values[ref_sample]))
                f.write('\t'.join(row_data) + '\n')
        
        # Calculate Identity statistics
        identity_results, chromosome_results = calculate_identity_from_matrix(corrected_matrix_data, ref_samples, args.out)
        print("\nIdentity Analysis - Top 5 Results:")
        display_results = identity_results if isinstance(identity_results, list) else []
        for i, result in enumerate(display_results[:5]):
            print(f"{result['Reference']}: {result['Identity (%)']:.2f}% ({result['Matched_SNPs']}/{result['Compared_SNPs']} SNPs)")
        if len(display_results) > 5:
            print(f"... and {len(display_results)-5} more results")
        print(f"Complete results saved to: {args.out}_identity_results.tsv")
        
        # Auto-generate heatmap visualizations
        if args.generate_heatmaps and VISUALIZATION_AVAILABLE:
            try:
                # Create visualization tool
                vis_tool = SNPVisualizationTool(f"{args.out}_identity_matrix_corrected.txt")
                vis_tool.initialize()
                
                window_size = args.heatmap_window * 1000  # Convert to bp
                
                if args.chr is None:
                    # No --chr parameter: generate top 5 single sample heatmaps

                    for i, result in enumerate(display_results[:5]):
                        sample_name = result['Reference']
                        
                        original_path = vis_tool.generate_single_sample_heatmap(
                            sample=sample_name,
                            window_size=window_size,
                            output_dir=".",
                            show_plot=False
                        )
                        # Rename file to include output prefix and function classification
                        if original_path and os.path.exists(original_path):
                            new_filename = f"{args.out}_identity_{sample_name}.png"
                            new_path = os.path.join(".", new_filename)
                            os.rename(original_path, new_path)
                else:
                    # With --chr parameter: generate comparison plot with top 10 samples

                    top_samples = [result['Reference'] for result in display_results[:10]]
                    original_path = vis_tool.generate_chromosome_comparison_heatmap(
                        chromosome=args.chr,
                        samples=top_samples,
                        window_size=window_size,
                        output_dir=".",
                        show_plot=False
                    )
                    # Rename file to include output prefix and function classification
                    if original_path and os.path.exists(original_path):
                        new_filename = f"{args.out}_identity_{args.chr}.png"
                        new_path = os.path.join(".", new_filename)
                        os.rename(original_path, new_path)
                    

            except Exception as e:
                print(f"Heatmap generation failed: {e}")

    
    # If parental mode is enabled, perform parental comparison analysis
    if args.find_parents:
        if not args.refs:
            print("Error: Parental comparison analysis requires reference sample file (--refs)")
            sys.exit(1)
            
        parental_results, chromosome_parental_results = generate_parental_comparison_matrix(
            args.vcf, args.sample, ref_samples, args.out, 
            args.filter_times, args.filter_window, args.chr, args.threads
        )
        if parental_results is not None:
            print("\nParental Analysis - Top 5 Combinations:")
            display_results = parental_results if isinstance(parental_results, list) else []
            for i, result in enumerate(display_results[:5]):
                print(f"{result['parent1']} x {result['parent2']}: {result['match_percentage']:.2f}% ({result['matched_snps']}/{result['informative_snps']} SNPs)")
            if len(display_results) > 5:
                print(f"... and {len(display_results)-5} more combinations")
            print(f"Complete results saved to: {args.out}_parental_results.tsv")
            
            # Auto-generate heatmap visualizations (parental mode)
            if args.generate_heatmaps and VISUALIZATION_AVAILABLE:
                try:
                    # Read parental analysis result files
                    parental_results_file = f"{args.out}_parental_results.tsv"
                    parental_matrix_file = f"{args.out}_parental_matrix_corrected.txt"
                    
                    if not os.path.exists(parental_results_file):
                        print(f"Error: Parental results file does not exist: {parental_results_file}")
                        return
                    
                    if not os.path.exists(parental_matrix_file):
                        print(f"Error: Parental matrix file does not exist: {parental_matrix_file}")
                        return
                    
                    # Read parental analysis results to get highest matching combinations
                    import pandas as pd
                    # Manually read file and handle format issues
                    with open(parental_results_file, 'r') as f:
                        lines = f.readlines()
                    
                    # Rebuild correct TSV format
                    corrected_lines = []
                    header = "Parent1\tParent2\tChromosome\tMatch(%)\tInformative_SNPs\tMatched_SNPs\tNote\n"
                    corrected_lines.append(header)
                    
                    for line in lines[1:]:  # Skip original header line
                        if line.strip():
                            corrected_lines.append(line)
                    
                    # Write to temporary file
                    temp_file = parental_results_file + ".temp"
                    with open(temp_file, 'w') as f:
                        f.writelines(corrected_lines)
                    
                    results_df = pd.read_csv(temp_file, sep='\t')
                    
                    # Filter Overall results and sort by Match(%)
                    overall_results = results_df[results_df['Chromosome'] == 'Overall'].copy()
                    overall_results['Match(%)'] = pd.to_numeric(overall_results['Match(%)'], errors='coerce')
                    overall_results = overall_results.sort_values('Match(%)', ascending=False)
                    
                    # Clean up temporary file
                    os.remove(temp_file)
                    
                    # Create visualization tool (using parental matrix file)
                    vis_tool = SNPVisualizationTool(parental_matrix_file)
                    vis_tool.initialize()
                    
                    window_size = args.heatmap_window * 1000  # Convert to bp
                    
                    if args.chr is None:
                        # No --chr parameter: generate top 5 parental combination heatmaps

                        top_combinations = overall_results.head(5)
                        
                        for i, (_, row) in enumerate(top_combinations.iterrows()):
                            parent1 = row['Parent1']
                            parent2 = row['Parent2']
                            combination_name = f"{parent1}_x_{parent2}"
                            match_rate = row['Match(%)']
                            
    
                            
                            # Generate single sample heatmap using combination name
                            output_file = vis_tool.generate_single_sample_heatmap(
                                sample=combination_name,
                                window_size=window_size,
                                output_dir=".",
                                show_plot=False
                            )
                            
                            # Rename file to include output prefix and function classification
                            if output_file and os.path.exists(output_file):
                                old_name = output_file
                                new_name = os.path.join(".", f"{args.out}_parental_{combination_name}.png")
                                os.rename(old_name, new_name)
    
                    else:
                        # With --chr parameter: generate top 10 parental combination comparison heatmap

                        top_combinations = overall_results.head(10)
                        
                        # Build combination name list
                        combination_samples = []
                        for _, row in top_combinations.iterrows():
                            parent1 = row['Parent1']
                            parent2 = row['Parent2']
                            combination_name = f"{parent1}_x_{parent2}"
                            combination_samples.append(combination_name)
                        

                        
                        output_file = vis_tool.generate_chromosome_comparison_heatmap(
                            chromosome=args.chr,
                            samples=combination_samples,
                            window_size=window_size,
                            output_dir=".",
                            show_plot=False
                        )
                        
                        # Rename file to include output prefix and function classification
                        if output_file and os.path.exists(output_file):
                            old_name = output_file
                            new_name = os.path.join(".", f"{args.out}_parental_{args.chr}_comparison.png")
                            os.rename(old_name, new_name)

                        

                except Exception as e:
                    print(f"Parental analysis heatmap generation failed: {e}")
                    import traceback
                    traceback.print_exc()

    
    print("\n--- GIPA Analysis Complete ---")

if __name__ == '__main__':
    main()