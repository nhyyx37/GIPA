# GIPA: Genomic Identity and Parentage Analysis

**Genomic Identity and Parentage Analysis Tool**

GIPA is a professional genomic analysis tool specifically designed for sample identity verification and parentage analysis based on SNP (Single Nucleotide Polymorphism) data. This tool provides high-precision genomic comparison analysis and scientific-grade visualization capabilities.

## 🌟 Key Features

### 1. Identity Statistical Analysis
- **Sample Identity Verification**: Calculate similarity between query samples and reference samples based on SNP data
- **High-Precision Alignment**: Use sliding window correction algorithms to improve analysis accuracy
- **Batch Processing**: Support simultaneous comparison with multiple reference samples
- **Detailed Statistics**: Provide comprehensive comparison statistics and result ranking

### 2. Parental Comparison Analysis
- **Hybrid Identification**: Automatically identify hybrid and inbred line samples
- **Parental Combination Analysis**: Calculate matching rates for all possible parental combinations
- **Intelligent Thresholds**: Automatically determine classification thresholds based on sample heterozygosity
- **Genetic Law Validation**: Apply Mendelian inheritance laws to validate parental relationships

### 3. Scientific-Grade Visualization
- **Chromosome Heatmaps**: Display SNP distribution of individual samples across chromosomes
- **Sample Comparison Heatmaps**: Genomic comparison visualization between multiple samples
- **High-Quality Output**: 300 DPI resolution, suitable for academic publication
- **Customizable Windows**: Support user-defined analysis window sizes

### 4. Advanced Features
- **Multi-threading Processing**: Support parallel computing to accelerate analysis
- **Chromosome-Specific Analysis**: Analyze specific chromosomes as needed
- **Sliding Window Correction**: Reduce noise and improve result reliability
- **Flexible Output Formats**: Support multiple data output formats

## 📋 System Requirements

### Python Version
- Python 3.7 or higher

### Core Dependencies
```bash
# Required dependencies
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
pysam>=0.19.0

# Visualization dependencies
matplotlib>=3.5.0
seaborn>=0.11.0

# Optional dependencies
scikit-learn>=1.0.0
tqdm>=4.62.0
colorlog>=6.6.0
```

## 🚀 Installation Guide

### 1. Clone the Project
```bash
git clone https://github.com/nhyyx/GIPA.git
```

### 2. Install Dependencies
```bash
pip install -r requirements.txt
```

### 3. Verify Installation
```bash
python gipa.py --help
```

## 📖 Usage Guide

### Basic Usage

#### Identity Analysis
```bash
# Basic Identity analysis (using default output prefix "output")
python gipa.py --vcf input.vcf.gz --sample query_sample --refs samples.txt

# Specify output prefix
python gipa.py --vcf input.vcf.gz --sample query_sample --refs samples.txt --out my_analysis

# Specify chromosome analysis
python gipa.py --vcf input.vcf.gz --sample query_sample --refs samples.txt --chr Chr01

# Generate visualization heatmaps
python gipa.py --vcf input.vcf.gz --sample query_sample --refs samples.txt --generate-heatmaps
```

#### Parental Analysis
```bash
# Parental comparison analysis
python gipa.py --vcf input.vcf.gz --sample hybrid_sample --refs samples.txt --find_parents

# Parental analysis + visualization
python gipa.py --vcf input.vcf.gz --sample hybrid_sample --refs samples.txt --find_parents --generate-heatmaps

# Parental analysis with custom output prefix
python gipa.py --vcf input.vcf.gz --sample hybrid_sample --refs samples.txt --out parental_analysis --find_parents
```

### Advanced Parameters

```bash
# Custom heatmap window size (default 50kb)
python gipa.py --vcf input.vcf.gz --sample query_sample --refs samples.txt --heatmap-window 100

# Multi-threading processing (accelerate large dataset analysis)
python gipa.py --vcf input.vcf.gz --sample query_sample --refs samples.txt --threads 4

# Custom sliding window parameters
python gipa.py --vcf input.vcf.gz --sample query_sample --refs samples.txt --filter-times 3 --filter-window 10

# Combine multiple advanced parameters
python gipa.py --vcf input.vcf.gz --sample query_sample --refs samples.txt --out advanced_analysis --threads 8 --chr Chr01 --heatmap-window 200 --generate-heatmaps
```

### Complete Parameter List

| Parameter | Short | Type | Default | Description |
|-----------|-------|------|---------|-------------|
| `--vcf` | `-v` | Required | - | Path to VCF format SNP data file |
| `--sample` | `-s` | Required | - | Query sample name (from VCF file) |
| `--refs` | `-r` | Required | - | Reference sample list file (one sample name per line) |
| `--out` | `-o` | Optional | output | Output file prefix |
| `--filter-times` | `-ft` | Optional | 2 | Number of sliding window filter iterations |
| `--filter-window` | `-fw` | Optional | 5 | Sliding window size |
| `--threads` | `-t` | Optional | 1 | Number of threads to use |
| `--chr` | `-c` | Optional | All | Specify chromosome to analyze (e.g., Chr01, Chr02) |
| `--heatmap-window` | `-hw` | Optional | 50 | Heatmap window size (kb) |
| `--find_parents` | - | Optional | False | Enable parental comparison analysis mode |
| `--generate-heatmaps` | - | Optional | False | Generate heatmap visualizations (requires visualization modules) |

## 📁 Input File Formats

### VCF File
Standard VCF format file, supports compressed format (.vcf.gz):
```
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3
Chr01	1000	.	A	T	60	PASS	.	GT	0/0	0/1	1/1
Chr01	2000	.	G	C	60	PASS	.	GT	1/1	0/0	0/1
```

### Reference Sample File
Text file with one sample name per line:
```
sample1
sample2
sample3
```

## 📊 Output File Description

### Identity Analysis Output

#### 1. Comparison Matrix File (`*_identity_matrix.txt`)
```
Position	sample1	sample2	sample3
Chr01:1000	1	0	1
Chr01:2000	0	1	0
```

#### 2. Corrected Matrix (`*_identity_matrix_corrected.txt`)
Comparison results after sliding window correction

#### 3. Statistical Results File (`*_identity_results.tsv`)
```
Sample	Identity(%)	Matched_SNPs	Total_SNPs
query_sample	100.00	50000	50000
sample1	85.50	42750	50000
sample2	72.30	36150	50000
```

### Parental Analysis Output

#### 1. Parental Comparison Matrix (`*_parental_matrix.txt`)
Comparison results for all parental combinations

#### 2. Parental Analysis Results (`*_parental_results.tsv`)
```
Parental_Combination	Match(%)	Matched_SNPs	Total_SNPs
parent1_x_parent2	95.50	47750	50000
parent1_x_parent3	78.20	39100	50000
```

### Visualization Output

#### Heatmap Files (saved in output directory by default)
- `output/sample_comparison_{chromosome}.png`: Chromosome sample comparison heatmap
- `output/{sample}_identity_{chromosome}.png`: Single sample chromosome heatmap
- `output/{combination}_parental_{chromosome}.png`: Parental combination heatmap

**Note**: If using the `--out` parameter to specify a custom prefix, output file names will change accordingly. For example, when using `--out my_analysis`, files will be named `my_analysis_identity_matrix.txt`, etc.

