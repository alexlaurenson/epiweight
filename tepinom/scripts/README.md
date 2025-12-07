# TepiNom 1.6 - ILP-Based Epitope Optimization Tool

TepiNom 1.6 is a unified, refactored epitope optimization tool that uses Integer Linear Programming (ILP) to maximize population coverage for vaccine design.

## Features

- **Dual optimization modes**: Epitope combination selection or protein region identification
- **Phenotypic frequency-based coverage**: Calculates population coverage using intra-locus and inter-locus HLA frequency analysis
- **MHC Class I & II support**: Optimizes for CD8+ (HLA-A/B/C) or CD4+ (DRB1/DPA1/DQA1) epitopes
- **High-risk HLA prioritization**: Apply custom weights to specific HLA alleles
- **Protein diversity constraints**: Enforce selection of epitopes from different proteins
- **Adaptive region selection**: Returns all 100% coverage regions or top N by coverage
- **Sliding window analysis**: Identify optimal regions using customizable window sizes

## Installation

### Requirements
- Python 3.7+
- Dependencies: `pandas`, `numpy`, `pulp`
- CBC solver (COIN-OR Branch and Cut)

### Install Dependencies
```bash
pip install pandas numpy pulp
```

### Install CBC Solver
- **Ubuntu/Debian**: `sudo apt-get install coinor-cbc`
- **macOS**: `brew install coin-or-tools/coinor/cbc`
- **Windows**: Download from [COIN-OR](https://www.coin-or.org/download/binary/Cbc/)

## Usage

### Basic Syntax
```bash
python3 tepinom_1.6.py <epitope_file> <hla_file> [options]
```

### Required Arguments
- `epitope_file`: Path to epitope binding affinity matrix (TSV format)
- `hla_file`: Path to HLA allele frequency file (TSV format)

### Key Options

#### General Options
- `--mhc_class {i,ii}`: MHC class (default: i)
- `--postfilter`: Use postfilter input format (required for most use cases)
- `--ba_rank_hit FLOAT`: Binding affinity rank threshold (default: 0.5 for MHC I, 1.0 for MHC II)
- `--conservation_score FLOAT`: Minimum conservation score (default: 0.8)
- `--out_dir PATH`: Output directory (default: output/)
- `--out_suffix STRING`: Suffix for output filenames

#### Optimization Mode
- `--goal {combo,region}`: Optimization goal (default: combo)
  - `combo`: Select optimal epitope combinations
  - `region`: Identify high-coverage protein regions

#### Region-Specific Options
- `--window_size INT`: Sliding window size in amino acids (required for region mode)
- `--best_region_num INT`: Number of top regions to output (default: 10, or adaptive)

#### Advanced Options
- `--max_epitopes INT`: Maximum epitopes to select
- `--prioritize_protein_diversity`: Enforce ≤1 epitope per protein
- `--high_risk_hlas STRING`: Comma-separated HLA patterns to prioritize (e.g., "HLA-A*29:02,HLA-B*57:01")
- `--hla_weight FLOAT`: Weight multiplier for high-risk HLAs (default: 2.0)
- `--hard_constraint`: Require coverage of high-risk HLAs

## Usage Examples

### Example 1: Basic Epitope Combination (MHC I)
```bash
python3 tepinom_1.6.py \
  mhci_epitopes.txt \
  mhci_frequencies.txt \
  --postfilter \
  --mhc_class i \
  --goal combo \
  --ba_rank_hit 0.5 \
  --max_epitopes 15 \
  --out_dir results/ \
  --out_suffix liver_mhci
```

**Output:**
- `results/optimized_epitopes_liver_mhci.tsv`
- `results/population_coverage_liver_mhci.tsv`

### Example 2: Region Search with Sliding Window (MHC II)
```bash
python3 tepinom_1.6.py \
  mhcii_epitopes.txt \
  mhcii_frequencies.txt \
  --postfilter \
  --mhc_class ii \
  --goal region \
  --window_size 100 \
  --ba_rank_hit 1.0 \
  --out_dir results/ \
  --out_suffix blood_mhcii_100aa
```

**Output:**
- `results/region_coverage_blood_mhcii_100aa.tsv`
- `results/optimized_epitopes_region_blood_mhcii_100aa.tsv`

### Example 3: High-Risk HLA Prioritization
```bash
python3 tepinom_1.6.py \
  epitopes.txt \
  frequencies.txt \
  --postfilter \
  --mhc_class i \
  --goal combo \
  --ba_rank_hit 0.5 \
  --high_risk_hlas "HLA-A*29:02,HLA-A*30:01,HLA-A*33:01" \
  --hla_weight 5.0 \
  --hard_constraint \
  --out_dir results/
```

### Example 4: Multiple Window Sizes (Batch Processing)
```bash
for window in 50 100 150 200; do
  python3 tepinom_1.6.py \
    epitopes.txt \
    frequencies.txt \
    --postfilter \
    --mhc_class i \
    --goal region \
    --window_size $window \
    --ba_rank_hit 0.5 \
    --out_dir results/ \
    --out_suffix window_${window}aa
done
```

## Input File Formats

### Epitope File (Postfilter Format)
Tab-separated file with required columns:
```
Epitope	Position	ConservationScore	Source	Stage	HLA-A*01:01	HLA-A*02:01	HLA-B*07:02	...
SHNGFWSVP	145	0.95	CSP	Sporozoite	0.3	1.2	2.5	...
KGVNINSPK	289	0.88	TRAP	Liver	2.1	0.4	0.6	...
MARTYGSAL	412	0.92	MSP1	Blood	0.8	0.9	1.1	...
```

**Required columns:**
1. `Epitope`: Peptide sequence
2. `Position`: Starting position in protein
3. `ConservationScore`: Conservation metric (0-1)
4. `Source`: Protein identifier
5. `Stage`: Life cycle stage or category
6. HLA columns: Binding affinity ranks for each allele

### HLA Frequency File
Tab-separated file with allele frequencies:
```
HLA-A*01:01	0.1234
HLA-A*02:01	0.2456
HLA-A*03:01	0.0876
HLA-B*07:02	0.0987
HLA-B*08:01	0.1123
```

**Format:**
- Column 1: HLA allele name (must match epitope file headers)
- Column 2: Population frequency (0-1 scale)

## Output Files

### Combo Mode Outputs

#### `optimized_epitopes_*.tsv`
Selected epitopes with individual coverage metrics:
```
Epitope	Position	ConservationScore	Source	Stage	HLA-A_Cov_pct	HLA-B_Cov_pct	HLA-C_Cov_pct	Total_Cov_pct	HLA-A*01:01	...
SHNGFWSVP	145	0.95	CSP	Sporozoite	45.23	38.67	22.45	72.35	0.3	...
KGVNINSPK	289	0.88	TRAP	Liver	32.45	28.92	18.23	58.47	2.1	...
```

#### `population_coverage_*.tsv`
Coverage progression for increasing epitope counts:
```
Num_Epitopes	Best_Combination	HLA-A_Cov_pct	HLA-B_Cov_pct	HLA-C_Cov_pct	Total_Cov_pct
1	SHNGFWSVP	45.23	38.67	22.45	72.35
2	SHNGFWSVP,KGVNINSPK	67.89	58.92	41.23	88.47
3	SHNGFWSVP,KGVNINSPK,MARTYGSAL	78.45	72.34	55.67	94.28
```

### Region Mode Outputs

#### `region_coverage_*.tsv`
Top protein regions ranked by coverage:
```
Protein	Region_Start	Region_End	Region_Type	Num_Epitopes	Epitopes	HLA-A_Cov_pct	HLA-B_Cov_pct	HLA-C_Cov_pct	Total_Cov_pct
CSP	145	245	sliding	5	SHNGFWSVP,KPKDELAQ,...	78.45	72.34	55.67	94.28
TRAP	280	380	sliding	4	KGVNINSPK,GIAGGLALL,...	65.23	68.92	48.23	87.15
```

#### `optimized_epitopes_region_*.tsv`
All epitopes from top regions with individual metrics (same format as combo mode epitope output)

## Population Coverage Calculation

TepiNom calculates population coverage using **phenotypic frequencies** for accurate prevalence estimation:

### Methodology

1. **Phenotypic Frequency Conversion**: Converts allele frequencies to phenotypic presence
   - P(allele present in individual) = 2f - f² 
   - Based on Hardy-Weinberg equilibrium assumptions
   - Accounts for both homozygous and heterozygous carriers

2. **Intra-locus Coverage**: Coverage within each HLA locus (A, B, C or DRB1, DPA1, DQA1)
   - Locus coverage = 1 - ∏(1 - phenofreq) for all covered alleles in that locus
   - Represents probability that at least one allele in the locus is covered

3. **Inter-locus Coverage**: Combined coverage across all loci
   - Total coverage = 1 - ∏(1 - locus_coverage) for all loci
   - Represents probability that at least one allele across all loci is covered

### Example Calculation

For an epitope binding HLA-A\*02:01 (f=0.25) and HLA-B\*07:02 (f=0.15):

**Step 1: Phenotypic frequencies**
- HLA-A\*02:01: 2(0.25) - (0.25)² = 0.4375 (43.75%)
- HLA-B\*07:02: 2(0.15) - (0.15)² = 0.2775 (27.75%)

**Step 2: Intra-locus coverage**
- HLA-A locus: 43.75%
- HLA-B locus: 27.75%
- HLA-C locus: 0% (no binding)

**Step 3: Inter-locus coverage**
- Total = 1 - (1-0.4375)(1-0.2775)(1-0) = 1 - (0.5625)(0.7225)(1) = 59.36%

This means 59.36% of the population carries at least one of these HLA alleles.

## Optimization Algorithm

### Combo Mode
1. Filter epitopes by conservation score
2. Build binary coverage matrix (epitopes × HLA alleles)
3. Calculate allele weights based on phenotypic frequency and risk priority
4. For k = 1 to max_epitopes:
   - Solve ILP to maximize weighted allele coverage
   - Calculate exact population coverage using phenotypic frequencies
   - Stop if no improvement or 100% coverage reached

### Region Mode
1. Filter epitopes by conservation score
2. For each protein:
   - Slide window across protein positions
   - For each window:
     - Run ILP optimization on epitopes within window
     - Calculate coverage for selected epitopes
3. Rank all windows by coverage
4. Return top N regions (or all with 100% coverage)

## Adaptive Region Selection

When `--best_region_num` is **not specified**, TepiNom automatically:
1. Returns **all regions** achieving ≥100% coverage, OR
2. Returns **top 10 regions** if none reach 100%

When `--best_region_num X` **is specified**:
- Returns exactly top X regions by coverage

## Performance Considerations

- **Runtime**: Depends on dataset size
  - Combo mode: 10 min - 4 hours typical
  
- **Memory**: 8-16GB recommended for large datasets

- **Parallelization**: Use array jobs for multiple window sizes or tissues

## Troubleshooting

### "Could not locate CBC solver"
- Ensure CBC is installed and in system PATH
- Check installation: `which cbc` (Unix) or `where cbc` (Windows)

### "No epitopes pass conservation filter"
- Lower `--conservation_score` threshold
- Check that conservation scores are in correct format (0-1 scale)

### ILP solver timeouts
- Reduce dataset size by pre-filtering
- Increase time limit in code (default: 300s per ILP solve)
- Consider using `--max_epitopes` to limit search space

### Memory errors
- Process proteins separately
- Reduce window size
- Filter to specific HLA alleles of interest

## Advanced: High-Performance Computing

### SLURM Example
```bash
#!/bin/bash
#SBATCH --job-name=tepinom
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

module load python
python3 tepinom_1.6.py \
  epitopes.txt \
  frequencies.txt \
  --postfilter \
  --mhc_class i \
  --goal region \
  --window_size 100 \
  --ba_rank_hit 0.5 \
  --out_dir results/
```

## Citation

If you use TepiNom in your research, please cite:
```
[Citation details to be added - Patent Pending]
```

## Version History

- **v1.6**: Unified refactored version with improved code organization
- **v1.5**: Added region search capabilities
- **v1.0-1.4**: Original MHC I/II separate scripts

## License

Apache 2.0

## Support

For questions, bug reports, or feature requests:
- Open an issue on GitHub
- Contact: Alex Laurenson

## Acknowledgments

TepiNom uses:
- **PuLP** for linear programming
- **CBC** (COIN-OR) for ILP solving
- **Phenotypic frequency calculations** based on Hardy-Weinberg equilibrium for accurate population coverage estimation