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

    hla_cols = [col for col in epitopes.columns if col.startswith('HLA_')]
    return epitopes, hla_cols

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
    hla_A = {key: val for key, val in hla_frequencies.items() if key.startswith('HLA_A')}
    hla_B = {key: val for key, val in hla_frequencies.items() if key.startswith('HLA_B')}
    hla_C = {key: val for key, val in hla_frequencies.items() if key.startswith('HLA_C')}

    # Normalize frequencies for each class
    def normalize_class(frequencies):
        total_freq = sum(frequencies.values())
        return {key: val / total_freq for key, val in frequencies.items()} if total_freq > 0 else frequencies

    normalized_hla_A = normalize_class(hla_A)
    normalized_hla_B = normalize_class(hla_B)
    normalized_hla_C = normalize_class(hla_C)

    # Combine normalized frequencies
    normalized_hla_frequencies = {**normalized_hla_A, **normalized_hla_B, **normalized_hla_C}

    return normalized_hla_frequencies

# Function to filter epitopes by median binding affinity (based on HLA columns)
def filter_by_median_binding_affinity(epitopes, hla_cols, rank_cutoff):
    epitopes['median_binding_rank'] = epitopes[hla_cols].median(axis=1)
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
        hla_coverage = 0
        # Calculate how much coverage this epitope can contribute
        for hla in hla_frequencies:
            if epitope.get(hla, 100) <= hla_rank_cutoff:
                hla_coverage += hla_frequencies.get(hla, 0)

        if hla_coverage > 0:
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
    total_hla_a_covered = set()
    total_hla_b_covered = set()
    total_hla_c_covered = set()

    for _, epitope in epitopes.iterrows():
        epitope_data = {'epitope_sequence': epitope['epitope_sequence']}
        total_coverage = 0

        covered_hla_a = set()
        covered_hla_b = set()
        covered_hla_c = set()

        for hla in hla_frequencies:
            if epitope.get(hla, 100) <= hla_rank_cutoff:
                if hla.startswith('HLA_A'):
                    covered_hla_a.add(hla)
                    total_hla_a_covered.add(hla)
                elif hla.startswith('HLA_B'):
                    covered_hla_b.add(hla)
                    total_hla_b_covered.add(hla)
                elif hla.startswith('HLA_C'):
                    covered_hla_c.add(hla)
                    total_hla_c_covered.add(hla)

        # Calculate the epitope coverage within each HLA class
        epitope_data['HLA_A'] = sum(hla_frequencies[h] for h in covered_hla_a)
        epitope_data['HLA_B'] = sum(hla_frequencies[h] for h in covered_hla_b)
        epitope_data['HLA_C'] = sum(hla_frequencies[h] for h in covered_hla_c)

        # Calculate population coverage based on the frequencies of covered alleles
        total_coverage = epitope_data['HLA_A'] + epitope_data['HLA_B'] + epitope_data['HLA_C']
        epitope_data['population_coverage'] = (total_coverage / sum(hla_frequencies.values())) * 100
        
        coverage_data.append(epitope_data)

    # Calculate the overall combination coverage based on unique alleles covered
    total_unique_hla_a = sum(1 for h in hla_frequencies if h.startswith('HLA_A'))
    total_unique_hla_b = sum(1 for h in hla_frequencies if h.startswith('HLA_B'))
    total_unique_hla_c = sum(1 for h in hla_frequencies if h.startswith('HLA_C'))

    # Compute percentages of unique alleles covered
    combo_coverage_a = (len(total_hla_a_covered) / total_unique_hla_a) * 100 if total_unique_hla_a > 0 else 0
    combo_coverage_b = (len(total_hla_b_covered) / total_unique_hla_b) * 100 if total_unique_hla_b > 0 else 0
    combo_coverage_c = (len(total_hla_c_covered) / total_unique_hla_c) * 100 if total_unique_hla_c > 0 else 0

    # Overall population coverage based on the combination
    total_combination_coverage = (
        len(total_hla_a_covered | total_hla_b_covered | total_hla_c_covered) / len(hla_frequencies)
    ) * 100

    # Add the combination row to the output
    combination_row = {
        'epitope_sequence': 'combination',
        'HLA_A': combo_coverage_a,
        'HLA_B': combo_coverage_b,
        'HLA_C': combo_coverage_c,
        'population_coverage': total_combination_coverage
    }
    coverage_data.append(combination_row)

    return pd.DataFrame(coverage_data)

# Main function to load input data, filter and optimize epitopes, and save results
def main(protein_fasta, epitope_file, hla_file, rank_cutoff, conservation_cutoff, max_epitopes, prioritize_different_proteins, hla_rank_cutoff):
    print("Loading inputs...")
    proteins = load_protein_sequences(protein_fasta)
    epitopes, hla_cols = load_epitopes(epitope_file)
    hla_frequencies = load_hla_frequencies(hla_file)

    # Normalize HLA frequencies
    print("Normalizing HLA frequencies...")
    hla_frequencies = normalize_hla_frequencies(hla_frequencies)

    # Calculate conservation scores for each epitope
    print("Calculating conservation scores...")
    epitopes['conservation_score'] = epitopes.apply(lambda row: calculate_conservation_score(row['epitope_sequence'], row['protein'], proteins), axis=1)

    # Filtering by median binding affinity
    print("Filtering by median binding affinity...")
    filtered_epitopes = filter_by_median_binding_affinity(epitopes, hla_cols, rank_cutoff)

    # Filtering by conservation score
    print("Filtering by conservation score...")
    filtered_epitopes = filtered_epitopes[filtered_epitopes['conservation_score'] >= conservation_cutoff]
    print(f"Remaining after conservation filtering: {len(filtered_epitopes)} epitopes.")

    # Optimizing population coverage
    print("Optimizing population coverage...")
    optimized_epitopes = optimize_population_coverage(filtered_epitopes, hla_frequencies, max_epitopes, prioritize_different_proteins, hla_rank_cutoff)

    # Calculating population coverage for optimized epitopes
    print("Calculating population coverage for selected epitopes...")
    population_coverage_df = calculate_population_coverage(optimized_epitopes, hla_frequencies, hla_rank_cutoff)

    # Saving unfiltered epitopes (all the original epitopes)
    print("Saving unfiltered output file...")
    unfiltered_epitopes = epitopes[['epitope_sequence', 'protein', 'conservation_score', 'median_binding_rank'] + hla_cols]
    save_output_file(unfiltered_epitopes, "tepinom_output_unfiltered.txt")

    # Saving filtered epitopes (those that passed the filters)
    print("Saving filtered output file...")
    filtered_epitopes = filtered_epitopes[['epitope_sequence', 'protein', 'conservation_score', 'median_binding_rank'] + hla_cols]
    save_output_file(filtered_epitopes, "tepinom_output_filtered.txt")

    # Saving population coverage output
    print("Saving population coverage output...")
    save_output_file(population_coverage_df, "tepinom_output_popcov.txt")

if __name__ == "__main__":
    protein_fasta = sys.argv[1]
    epitope_file = sys.argv[2]
    hla_file = sys.argv[3]
    rank_cutoff = float(sys.argv[4])  # Input rank_cutoff as argument
    conservation_cutoff = float(sys.argv[5])
    max_epitopes = int(sys.argv[6])
    prioritize_different_proteins = bool(int(sys.argv[7]))
    hla_rank_cutoff = float(sys.argv[8])

    main(protein_fasta, epitope_file, hla_file, rank_cutoff, conservation_cutoff, max_epitopes, prioritize_different_proteins, hla_rank_cutoff)

