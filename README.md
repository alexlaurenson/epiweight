# EpiWeight
A Toolkit for Epitope Selection and Population Coverage Optimization

Welcome to EpiWeight! This repository contains tools for epitope downselection and optimization, particularly focusing on T-cell epitope selection based on binding affinity, population coverage, and conservation.

## Overview

EpiWeight provides methods for:
- Filtering epitopes based on binding affinity thresholds for specific HLA alleles
- Optimizing epitope selection to maximize population coverage while minimizing redundancy
- Ensuring diversity in selected epitopes by considering protein sources and avoiding overlapping sequences
- Identifying optimal protein regions for vaccine design
- Calculating population coverage metrics using phenotypic frequency

## Tools in This Repository

### TepiNom 1.6 (Patent Pending) - **Unified ILP Optimization Tool**
The latest iteration of TepiNom provides a single, refactored script (`tepinom_1.6.py`) that handles both MHC Class I and MHC Class II epitope optimization using Integer Linear Programming (ILP).

**Key Features:**
- **Two optimization modes:**
  - **Combo mode**: Select optimal epitope combinations for maximum population coverage
  - **Region mode**: Identify protein regions with highest immunogenic potential using sliding windows
- **Computed phenotypic frequency** for population coverage calculations
- **High-risk HLA prioritization** with configurable weights
- **Protein diversity constraints** to avoid redundant epitope selection
- **Adaptive region selection**: Automatically returns all regions achieving 100% coverage, or top N if none reach 100%
- Support for both MHC Class I (CD8+ T-cell epitopes: HLA-A, HLA-B, HLA-C) and MHC Class II (CD4+ T-cell epitopes: DRB1, DPA1, DQA1)

See [tepinom/README.md](tepinom/README.md) for detailed documentation.

### Legacy Tools
- `tepinom_mhci_arg.py` – Original MHC Class I epitope selection script
- `tepinom_mhcii_arg.py` – Original MHC Class II epitope selection script

*Note: These are maintained for backward compatibility. New projects should use `tepinom_1.6.py`.*

## How TepiNom Works

TepiNom uses Integer Linear Programming to optimize epitope selection by:

1. **Filtering** epitopes based on:
   - Binding affinity to HLA alleles (rank-based thresholds)
   - Conservation scores across protein variants
   
2. **Optimizing** selections to:
   - Maximize population coverage
   - Consider intra-locus and inter-locus HLA coverage
   - Apply weights to high-risk or priority HLA alleles
   - Enforce protein diversity (optional)
   
3. **Region Analysis** (optional):
   - Identify contiguous protein regions with high immunogenic potential
   - Use sliding windows of customizable sizes (e.g., 50aa, 100aa, 200aa)
   - Rank regions by population coverage

## Quick Start

### Installation
```bash
git clone https://github.com/alexlaurenson/epiweight
cd epiweight
```

**Dependencies:**
```bash
pip install pandas numpy pulp
```

**Solver Required:** TepiNom requires the CBC (COIN-OR Branch and Cut) solver for ILP optimization.
- Ubuntu/Debian: `sudo apt-get install coinor-cbc`
- macOS: `brew install coin-or-tools/coinor/cbc`
- Windows: Download from [COIN-OR website](https://www.coin-or.org/download/binary/Cbc/)

### Basic Usage Examples

#### Epitope Combination Optimization (MHC I)
```bash
python3 tepinom_1.6.py \
  epitopes.txt \
  hla_frequencies.txt \
  --postfilter \
  --mhc_class i \
  --goal combo \
  --ba_rank_hit 0.5 \
  --max_epitopes 20 \
  --out_dir output/
```

#### Region-Based Search (MHC II, 100aa windows)
```bash
python3 tepinom_1.6.py \
  epitopes.txt \
  hla_frequencies.txt \
  --postfilter \
  --mhc_class ii \
  --goal region \
  --window_size 100 \
  --ba_rank_hit 1.0 \
  --out_dir output/
```

#### With High-Risk HLA Prioritization
```bash
python3 tepinom_1.6.py \
  epitopes.txt \
  hla_frequencies.txt \
  --postfilter \
  --mhc_class i \
  --goal combo \
  --ba_rank_hit 0.5 \
  --high_risk_hlas HLA-A*29:02,HLA-A*30:01 \
  --hla_weight 5.0 \
  --out_dir output/
```

## Input File Formats

### Epitope File (TSV, postfilter format)
```
Epitope	Position	ConservationScore	Source	Stage	HLA-A*01:01	HLA-A*02:01	...
SHNGFWSVP	145	0.95	Protein1	Sporozoite	0.3	1.2	...
KGVNINSPK	289	0.88	Protein2	Liver	2.1	0.4	...
```

### HLA Frequency File (TSV)
```
HLA-A*01:01	0.1234
HLA-A*02:01	0.2456
HLA-B*07:02	0.0987
```

## Output Files

### Combo Mode
- `optimized_epitopes_*.tsv`: Selected epitopes with individual coverage metrics
- `population_coverage_*.tsv`: Progression showing coverage for increasing epitope counts

### Region Mode
- `region_coverage_*.tsv`: Top protein regions ranked by coverage
- `optimized_epitopes_region_*.tsv`: All epitopes from top regions with individual metrics

## Example Output (Population Coverage)
```
Num_Epitopes	Best_Combination	HLA-A_Cov_pct	HLA-B_Cov_pct	HLA-C_Cov_pct	Total_Cov_pct
1	SHNGFWSVP	45.23	38.67	22.45	72.35
2	SHNGFWSVP,KGVNINSPK	67.89	58.92	41.23	88.47
3	SHNGFWSVP,KGVNINSPK,MARTYGSA	78.45	72.34	55.67	94.28
```

## Advanced Features

- **Protein diversity enforcement**: Use `--prioritize_protein_diversity` to select at most one epitope per protein
- **Hard constraints**: Use `--hard_constraint` to require coverage of specified high-risk HLAs
- **Adaptive region selection**: Automatically identifies all regions with 100% coverage (region mode only)
- **Conservation filtering**: Set minimum conservation scores with `--conservation_score`

## Contributing

Want to improve EpiWeight or add new features? Feel free to fork this repo, submit pull requests, or open issues for discussions.

## Citation

If you use TepiNom in your research, please cite:
```
[Citation details to be added]
```

## License

EpiWeight is licensed under the Apache 2.0 License.

## Contact & Support

Have questions? Found a bug? Feel free to open an issue or reach out!

**Maintainer:** Alex Laurenson  
**Repository:** https://github.com/alexlaurenson/epiweight