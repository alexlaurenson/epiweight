import pandas as pd
import numpy as np
from Bio import SeqIO
import sys

def load_protein_sequences(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    print(f"Loaded {len(sequences)} protein sequences.")
    return sequences

def load_epitopes(epitope_file):
    try:
        epitopes = pd.read_csv(epitope_file)
        print(f"Loaded {len(epitopes)} epitopes with {len(epitopes.columns)} columns.")
    except Exception as e:
        print(f"Error loading epitope file: {e}")
        sys.exit(1)
    hla_cols = [col for col in epitopes.columns if col.startswith('HLA_')]
    return epitopes, hla_cols

def load_hla_frequencies(hla_file):
    try:
        hla_freq = pd.read_csv(hla_file)
        print(f"Loaded HLA frequencies for {len(hla_freq)} alleles.")
    except Exception as e:
        print(f"Error loading HLA frequency file: {e}")
        sys.exit(1)
    freq_dict = dict(zip(hla_freq['HLA_allele'], hla_freq['freq']))
    return freq_dict

def filter_by_median_binding_affinity(epitopes, hla_cols, rank_cutoff):
    epitopes['median_binding_rank'] = epitopes[hla_cols].median(axis=1)
    filtered_epitopes = epitopes[epitopes['median_binding_rank'] <= rank_cutoff]
    print(f"Remaining after median binding affinity filtering: {len(filtered_epitopes)} epitopes.")
    return filtered_epitopes

def filter_by_conservation(epitopes, conservation_cutoff):
    if 'conservation_score' in epitopes.columns:
        filtered_epitopes = epitopes[epitopes['conservation_score'] >= conservation_cutoff]
        print(f"Remaining after conservation filtering: {len(filtered_epitopes)} epitopes.")
        return filtered_epitopes
    else:
        print("No conservation score column found. Skipping conservation filtering.")
        return epitopes

def optimize_population_coverage(epitopes, hla_frequencies, max_epitopes, prioritize_different_proteins):
    selected_epitopes = []
    covered_hla = set()

    epitopes = epitopes.sort_values(by='median_binding_rank')

    for _, epitope in epitopes.iterrows():
        hla_coverage = sum([hla_frequencies.get(hla, 0) for hla in hla_frequencies if epitope.get(hla, 100) <= 100])
        if hla_coverage > 0:
            if prioritize_different_proteins:
                if epitope['protein'] not in [e['protein'] for e in selected_epitopes]:
                    selected_epitopes.append(epitope)
                    covered_hla.update([hla for hla in hla_frequencies if epitope.get(hla, 100) <= 100])
            else:
                selected_epitopes.append(epitope)
                covered_hla.update([hla for hla in hla_frequencies if epitope.get(hla, 100) <= 100])
            if len(selected_epitopes) >= max_epitopes:
                break

    print(f"Final selected epitopes: {len(selected_epitopes)}")
    return pd.DataFrame(selected_epitopes)

def calculate_conservation_score(epitope_sequence, protein_prefix, protein_sequences):
    # Find sequences that match the protein prefix in their header
    matching_proteins = {key: seq for key, seq in protein_sequences.items() if key.startswith(protein_prefix)}

    # Count how many sequences contain the epitope sequence
    matches = sum(1 for seq in matching_proteins.values() if epitope_sequence in seq)

    # Calculate conservation score as the ratio of matching sequences to total sequences
    conservation_score = matches / len(matching_proteins) if matching_proteins else 0
    return conservation_score

def save_output_file(epitopes, output_file):
    epitopes.to_csv(output_file, sep='\t', index=False)
    print(f"Output file saved as {output_file}")

def main(protein_fasta, epitope_file, hla_file, rank_cutoff, conservation_cutoff, max_epitopes, prioritize_different_proteins):
    print("Loading inputs...")
    proteins = load_protein_sequences(protein_fasta)
    epitopes, hla_cols = load_epitopes(epitope_file)
    hla_frequencies = load_hla_frequencies(hla_file)

    # Calculate conservation scores for each epitope
    print("Calculating conservation scores...")
    epitopes['conservation_score'] = epitopes.apply(lambda row: calculate_conservation_score(row['epitope_sequence'], row['protein'], proteins), axis=1)

    # Filtering by median binding affinity
    print("Filtering by median binding affinity...")
    filtered_epitopes = filter_by_median_binding_affinity(epitopes, hla_cols, rank_cutoff)

    # Filtering by conservation score
    print("Filtering by conservation score...")
    filtered_epitopes = filter_by_conservation(filtered_epitopes, conservation_cutoff)

    # Optimizing population coverage
    print("Optimizing population coverage...")
    optimized_epitopes = optimize_population_coverage(filtered_epitopes, hla_frequencies, max_epitopes, prioritize_different_proteins)

    # Saving unfiltered epitopes (all the original epitopes)
    print("Saving unfiltered output file...")
    unfiltered_epitopes = epitopes[['epitope_sequence', 'protein', 'conservation_score', 'median_binding_rank'] + hla_cols]
    save_output_file(unfiltered_epitopes, "tepinom_output_unfiltered.txt")

    # Saving filtered epitopes (those that passed the filters)
    print("Saving filtered output file...")
    filtered_epitopes = filtered_epitopes[['epitope_sequence', 'protein', 'conservation_score', 'median_binding_rank'] + hla_cols]
    save_output_file(filtered_epitopes, "tepinom_output_filtered.txt")

if __name__ == "__main__":
    protein_fasta = sys.argv[1]
    epitope_file = sys.argv[2]
    hla_file = sys.argv[3]
    rank_cutoff = float(sys.argv[4])
    conservation_cutoff = float(sys.argv[5])
    max_epitopes = int(sys.argv[6])
    prioritize_different_proteins = bool(int(sys.argv[7]))

    main(protein_fasta, epitope_file, hla_file, rank_cutoff, conservation_cutoff, max_epitopes, prioritize_different_proteins)

