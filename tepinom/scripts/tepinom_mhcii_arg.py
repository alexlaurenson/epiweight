import os
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
import sys

# Function to load protein sequences from a FASTA file
def load_protein_sequences(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    print(f"Loaded {len(sequences)} protein sequences.")
    return sequences

# Function to load epitope data from a CSV file
def load_epitopes(epitope_file):
    try:
        epitopes = pd.read_csv(epitope_file)
        print(f"Loaded {len(epitopes)} epitopes with {len(epitopes.columns)} columns.")
    except Exception as e:
        print(f"Error loading epitope file: {e}")
        sys.exit(1)

    hla_DQA1ols = [col for col in epitopes.columns if 'DRB1' in col or 'HLA-DPA1' in col or 'HLA-DQA1' in col]
    return epitopes, hla_DQA1ols

# Function to load HLA frequencies from a CSV file
def load_hla_frequencies(hla_file):
    try:
        hla_freq = pd.read_csv(hla_file)
        print(f"Loaded HLA frequencies for {len(hla_freq)} alleles.")
    except Exception as e:
        print(f"Error loading HLA frequency file: {e}")
        sys.exit(1)
    freq_dict = dict(zip(hla_freq['HLA_allele'], hla_freq['freq']))
    return freq_dict

# Normalize the HLA frequencies within each class (HLA-A, HLA-B, HLA-C)
def normalize_hla_frequencies(hla_frequencies):
    # Separate HLA frequencies by class
    DRB1 = {key: val for key, val in hla_frequencies.items() if key.startswith('DRB1')}
    hla_DPA1 = {key: val for key, val in hla_frequencies.items() if key.startswith('HLA-DPA1')}
    hla_DQA1 = {key: val for key, val in hla_frequencies.items() if key.startswith('HLA-DQA1')}

    # Normalize frequencies for each class
    def normalize_class(frequencies):
        total_freq = sum(frequencies.values())
        return {key: val / total_freq for key, val in frequencies.items()} if total_freq > 0 else frequencies

    normalized_DRB1 = normalize_class(DRB1)
    normalized_hla_DPA1 = normalize_class(hla_DPA1)
    normalized_hla_DQA1 = normalize_class(hla_DQA1)

    # Combine normalized frequencies
    normalized_hla_frequencies = {**normalized_DRB1, **normalized_hla_DPA1, **normalized_hla_DQA1}

    return normalized_hla_frequencies

# Function to filter epitopes by median binding affinity (based on HLA columns)
def filter_by_median_binding_affinity(epitopes, hla_DQA1ols, rank_cutoff):
    epitopes['median_binding_rank'] = epitopes[hla_DQA1ols].median(axis=1)
    filtered_epitopes = epitopes[epitopes['median_binding_rank'] <= rank_cutoff]
    print(f"Remaining after median binding affinity filtering: {len(filtered_epitopes)} epitopes.")
    return filtered_epitopes

# Function to calculate conservation score for each epitope
def calculate_conservation_score(epitope_sequence, protein_prefix, protein_sequences):
    matching_proteins = {key: seq for key, seq in protein_sequences.items() if key.startswith(protein_prefix)}
    matches = sum(1 for seq in matching_proteins.values() if epitope_sequence in seq)
    conservation_score = matches / len(matching_proteins) if matching_proteins else 0
    return conservation_score

# Function to check if two epitopes overlap based on amino acid sequence
def check_overlap(epitope1, epitope2):
    return epitope1[-4:] == epitope2[:4]

# Function to optimize population coverage by selecting epitopes
def optimize_population_coverage(epitopes, hla_frequencies, max_epitopes, prioritize_different_proteins, hla_rank_cutoff):
    selected_epitopes = []
    used_proteins = set()

    epitopes = epitopes.sort_values(by='median_binding_rank')

    # Function to calculate coverage based on selected epitopes
    def calculate_current_coverage():
        total_coverage = 0
        for epitope in selected_epitopes:
            for hla in hla_frequencies:
                if epitope.get(hla, 100) <= hla_rank_cutoff:
                    total_coverage += hla_frequencies.get(hla, 0)
        return total_coverage

    for _, epitope in epitopes.iterrows():
        hla_DQA1overage = 0
        # Calculate how much coverage this epitope can contribute
        for hla in hla_frequencies:
            if epitope.get(hla, 100) <= hla_rank_cutoff:
                hla_DQA1overage += hla_frequencies.get(hla, 0)

        if hla_DQA1overage > 0:
            if epitope['protein'] in used_proteins:
                continue

            overlapping = any(
                selected['protein'] == epitope['protein'] and 
                check_overlap(selected['epitope_sequence'], epitope['epitope_sequence'])
                for selected in selected_epitopes
            )
            if overlapping:
                continue

            current_coverage_before = calculate_current_coverage()
            selected_epitopes.append(epitope)
            used_proteins.add(epitope['protein'])

            current_coverage_after = calculate_current_coverage()

            if current_coverage_after > current_coverage_before:
                pass  # Coverage increased, continue
            else:
                break

            if len(selected_epitopes) >= max_epitopes:
                break

    print(f"Final selected epitopes: {len(selected_epitopes)}")
    return pd.DataFrame(selected_epitopes)

# Function to save the output to a CSV file
def save_output_file(epitopes, output_file):
    epitopes.to_csv(output_file, sep='\t', index=False)
    print(f"Output file saved as {output_file}")

# Function to calculate population coverage for each epitope
def calculate_population_coverage(epitopes, hla_frequencies, hla_rank_cutoff):
    coverage_data = []

    # Track total coverage across all selected epitopes
    total_DRB1_covered = set()
    total_hla_DPA1_covered = set()
    total_hla_DQA1_covered = set()

    for _, epitope in epitopes.iterrows():
        epitope_data = {'epitope_sequence': epitope['epitope_sequence']}
        total_coverage = 0

        covered_DRB1 = set()
        covered_hla_DPA1 = set()
        covered_hla_DQA1 = set()

        for hla in hla_frequencies:
            if epitope.get(hla, 100) <= hla_rank_cutoff:
                if hla.startswith('DRB1'):
                    covered_DRB1.add(hla)
                    total_DRB1_covered.add(hla)
                elif hla.startswith('HLA-DPA1'):
                    covered_hla_DPA1.add(hla)
                    total_hla_DPA1_covered.add(hla)
                elif hla.startswith('HLA-DQA1'):
                    covered_hla_DQA1.add(hla)
                    total_hla_DQA1_covered.add(hla)

        # Calculate the epitope coverage within each HLA class
        epitope_data['DRB1'] = sum(hla_frequencies[h] for h in covered_DRB1)
        epitope_data['HLA-DPA1'] = sum(hla_frequencies[h] for h in covered_hla_DPA1)
        epitope_data['HLA-DQA1'] = sum(hla_frequencies[h] for h in covered_hla_DQA1)

        # Calculate population coverage based on the frequencies of covered alleles
        total_coverage = epitope_data['DRB1'] + epitope_data['HLA-DPA1'] + epitope_data['HLA-DQA1']
        epitope_data['population_coverage'] = (total_coverage / sum(hla_frequencies.values())) * 100
        
        coverage_data.append(epitope_data)

    # Calculate the overall combination coverage based on unique alleles covered
    total_unique_DRB1 = sum(1 for h in hla_frequencies if h.startswith('DRB1'))
    total_unique_hla_DPA1 = sum(1 for h in hla_frequencies if h.startswith('HLA-DPA1'))
    total_unique_hla_DQA1 = sum(1 for h in hla_frequencies if h.startswith('HLA-DQA1'))

    # Compute percentages of unique alleles covered
    combo_coverage_a = (len(total_DRB1_covered) / total_unique_DRB1) * 100 if total_unique_DRB1 > 0 else 0
    combo_coverage_b = (len(total_hla_DPA1_covered) / total_unique_hla_DPA1) * 100 if total_unique_hla_DPA1 > 0 else 0
    combo_coverage_c = (len(total_hla_DQA1_covered) / total_unique_hla_DQA1) * 100 if total_unique_hla_DQA1 > 0 else 0

    # Overall population coverage based on the combination
    total_combination_coverage = (
        len(total_DRB1_covered | total_hla_DPA1_covered | total_hla_DQA1_covered) / len(hla_frequencies)
    ) * 100

    # Add the combination row to the output
    combination_row = {
        'epitope_sequence': 'combination',
        'DRB1': combo_coverage_a,
        'HLA-DPA1': combo_coverage_b,
        'HLA-DQA1': combo_coverage_c,
        'population_coverage': total_combination_coverage
    }
    coverage_data.append(combination_row)

    return pd.DataFrame(coverage_data)

# Define the argument parser
def parse_args():
    parser = argparse.ArgumentParser(
        description="This script processes and optimizes epitope data based on protein sequences, epitope binding affinities, and HLA allele frequencies."
    )
    parser.add_argument('protein_fasta', type=str, help="Path to the protein FASTA file containing protein sequences.")
    parser.add_argument('epitope_file', type=str, help="Path to the CSV file containing epitopes data.")
    parser.add_argument('hla_file', type=str, help="Path to the CSV file containing HLA allele frequency data.")
    
    # New argument for output file suffix
    parser.add_argument('--out_suffix', type=str, default='', 
                        help="Suffix for output file names (e.g., 'test' will produce output files like 'tepinom_output_filtered_test.txt').")

    # New argument for output directory
    parser.add_argument('--out_dir', type=str, default='.', 
                        help="Directory to save output files (default: current directory). The directory will be created if it doesn't exist.")
    
    # Optional arguments with descriptions
    parser.add_argument('--median_ba_rank', type=float, default=1.0, 
                        help="Rank cutoff for binding affinity (default: 1.0). Only epitopes with a median binding rank lower than this value are selected.")
    parser.add_argument('--conservation_score', type=float, default=0.0, 
                        help="Minimum conservation score for an epitope (default: 0.0). Epitopes with a score below this threshold are excluded.")
    parser.add_argument('--max_epitopes', type=int, default=10, 
                        help="Maximum number of epitopes to select (default: 10). The script will select up to this many epitopes based on the criteria.")
    parser.add_argument('--prioritize_protein_diversity', action='store_true', 
                        help="Flag to prioritize protein diversity in epitope selection (default: False). If set, the script aims to maximize the diversity of the selected epitopes.")
    parser.add_argument('--ba_rank_hit', type=float, default=1.0, 
                        help="Binding affinity rank hit threshold (default: 1.0). Only epitopes with binding ranks below this value are considered for population coverage.")

    return parser.parse_args()

# Main function
def main():
    args = parse_args()

    # Print out the argument values if the script runs without errors
    print(f"Protein FASTA file: {args.protein_fasta}")
    print(f"Epitopes file: {args.epitope_file}")
    print(f"HLA frequencies file: {args.hla_file}")
    print(f"Median binding affinity rank cutoff: {args.median_ba_rank}")
    print(f"Minimum conservation score: {args.conservation_score}")
    print(f"Maximum number of epitopes: {args.max_epitopes}")
    print(f"Prioritize protein diversity: {args.prioritize_protein_diversity}")
    print(f"Binding affinity rank hit threshold: {args.ba_rank_hit}")
    print(f"Output directory: {args.out_dir}")
    print(f"Output file suffix: {args.out_suffix}")

    # Create the output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading inputs...")
    proteins = load_protein_sequences(args.protein_fasta)
    epitopes, hla_DQA1ols = load_epitopes(args.epitope_file)
    hla_frequencies = load_hla_frequencies(args.hla_file)

    # Normalize HLA frequencies
    print("Normalizing HLA frequencies...")
    hla_frequencies = normalize_hla_frequencies(hla_frequencies)

    # Calculate conservation scores for each epitope
    print("Calculating conservation scores...")
    epitopes['conservation_score'] = epitopes.apply(
        lambda row: calculate_conservation_score(row['epitope_sequence'], row['protein'], proteins), 
        axis=1
    )

    # Filtering by median binding affinity
    print("Filtering by median binding affinity...")
    filtered_epitopes = filter_by_median_binding_affinity(epitopes, hla_DQA1ols, args.median_ba_rank)

    # Filtering by conservation score
    print("Filtering by conservation score...")
    filtered_epitopes = filtered_epitopes[filtered_epitopes['conservation_score'] >= args.conservation_score]
    print(f"Remaining after conservation filtering: {len(filtered_epitopes)} epitopes.")

    # Optimizing population coverage
    print("Optimizing population coverage...")
    optimized_epitopes = optimize_population_coverage(
        filtered_epitopes, hla_frequencies, args.max_epitopes, 
        args.prioritize_protein_diversity, args.ba_rank_hit
    )

    # Calculating population coverage for optimized epitopes
    print("Calculating population coverage for selected epitopes...")
    population_coverage_df = calculate_population_coverage(
        optimized_epitopes, hla_frequencies, args.ba_rank_hit
    )

    # Define suffix and output paths
    suffix = f"_{args.out_suffix}" if args.out_suffix else ""
    out_dir = args.out_dir

    # Saving unfiltered epitopes (all the original epitopes)
    print("Saving unfiltered output file...")
    unfiltered_epitopes = epitopes[['epitope_sequence', 'protein', 'conservation_score', 'median_binding_rank'] + hla_DQA1ols]
    save_output_file(unfiltered_epitopes, os.path.join(out_dir, f"tepinom_output_unfiltered{suffix}.txt"))

    # Saving filtered epitopes (those that passed the filters)
    print("Saving filtered output file...")
    filtered_epitopes = filtered_epitopes[['epitope_sequence', 'protein', 'conservation_score', 'median_binding_rank'] + hla_DQA1ols]
    save_output_file(filtered_epitopes, os.path.join(out_dir, f"tepinom_output_filtered{suffix}.txt"))

    # Saving population coverage output
    print("Saving population coverage output...")
    save_output_file(population_coverage_df, os.path.join(out_dir, f"tepinom_output_popcov{suffix}.txt"))

if __name__ == "__main__":
    main()
