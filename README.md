# EpiWeight

### A Toolkit for Epitope Selection and Population Coverage Optimization

Welcome to **EpiWeight**! This repository contains tools for **epitope downselection and optimization**, particularly focusing on **T-cell epitope selection** based on binding affinity, population coverage, and conservation.

## Overview
EpiWeight provides methods for:
- **Filtering epitopes** based on binding affinity thresholds for specific **HLA alleles**.
- **Optimizing epitope selection** to **maximize population coverage** while minimizing redundancy.
- **Ensuring diversity** in selected epitopes by considering protein sources and avoiding overlapping sequences.
- **Normalizing HLA allele frequencies** within each HLA class for fair optimization.

The primary tools in this repository are **TepiNom (Patent Pending)** for general epitope selection and optimization, with two separate scripts for MHC I and MHC II alleles:
- **`tepinom_mhci_arg.py`** – for **MHC Class I** epitope selection targeting **CD8+ T-cell epitopes** (HLA-A, HLA-B, HLA-C coverage).
- **`tepinom_mhcii_arg.py`** – for **MHC Class II** epitope selection targeting **CD4+ T-cell epitopes** (HLA-DRB1, HLA-DQ, HLA-DP coverage).

## TepiNom: Epitope Downselection & Population Optimization

**TepiNom** is designed to help researchers select **optimal epitope combinations** by filtering candidates based on:
- **Binding Affinity** – Only considering epitopes with strong predicted binding to specific **HLA alleles** (based on user-defined rank cutoffs).
- **Conservation** – Filtering out poorly conserved epitopes across protein variants.
- **Population Coverage** – Prioritizing epitopes that **cover the highest proportion** of a target population.
- **Avoiding Redundancy** – Ensuring non-overlapping sequences within the same protein and prioritizing diverse HLA coverage.

### How It Works
TepiNom takes as input:
- **Protein sequences** (FASTA file)
- **Epitope predictions** (CSV file with binding affinity data)
- **HLA allele frequencies** (CSV file with population frequency data)

It processes these inputs to **select the best epitope combination** and outputs:
- **Filtered epitopes** (based on binding affinity & conservation)
- **Optimized epitope list** (maximizing population coverage)
- **Population coverage statistics** (per HLA class and overall)

### Input Files
1. **Protein Sequences:** `proteins.fasta`
2. **Epitope Predictions:** `epitopes.csv`
3. **HLA Frequencies:** `hla_frequencies.csv`

### Output Files
1. `tepinom_output_unfiltered.txt`: **All input epitopes** (before filtering).
2. `tepinom_output_filtered.txt`: **Epitopes that pass filtering** (based on affinity and conservation).
3. `tepinom_output_popcov.txt`: **Optimized epitope selection with population coverage details**.

## Installation & Usage
Clone the repository and navigate to the project directory:

`git clone https://github.com/alexlaurenson/epiweight`<br>
`cd epiweight`

### Running TepiNom for MHC Class I Epitope Selection (CD8+ T-cells)
`python tepinom_mhcii_arg.py proteins.fasta epitopes.csv hla_frequencies.csv <rank_cutoff> <conservation_cutoff> <max_epitopes> <prioritize_different_proteins> <hla_rank_cutoff>`

Example:
`python tepinom_mhci_arg.py proteins.fasta epitopes.csv hla_frequencies.csv 2.0 0.7 10 1 0.5`

Where:
- `rank_cutoff` = Max binding affinity rank for filtering (e.g., `2.0`)
- `conservation_cutoff` = Min conservation score required (e.g., `0.7`)
- `max_epitopes` = Max number of selected epitopes (e.g., `10`)
- `prioritize_different_proteins` = Whether to enforce epitope diversity across proteins (`1` for True, `0` for False)
- `hla_rank_cutoff` = Threshold for an epitope to be considered for an HLA allele

## Example Output (Population Coverage Report)
| epitope_sequence | HLA_A | HLA_B | HLA_C | population_coverage |
|------------------|-------|-------|-------|---------------------|
| SHNGFWSVP        | 0.75  | 0.85  | 0.50  | 87.0                |
| KGVNINSPKPV      | 0.60  | 0.80  | 0.40  | 80.5                |
| MARTYGSA         | 0.55  | 0.30  | 0.90  | 67.8                |
| combination      | 0.95  | 0.92  | 0.85  | 98.2                |


- **Each epitope's coverage** is reported for **HLA-A, HLA-B, and HLA-C** categories for MHC Class I or **HLA-DRB1, HLA-DQ, and HLA-DP** categories for MHC Class II.
- **Final combination row** shows the **unique** HLA coverage across the entire selected epitope set.

## Contributing
Want to improve **TepiNom** or add new features? Feel free to **fork this repo**, submit **pull requests**, or open **issues** for discussions.

## License
EpiWeight is licensed under the **Apache 2.0**.

## Contact & Support
Have questions? Found a bug? Feel free to open an **issue** or reach out!
