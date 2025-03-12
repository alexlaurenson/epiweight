import pandas as pd
from Bio import SeqIO

def load_protein_sequences(fasta_file):
    # Load protein sequences from a multifasta file
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.id
        sequences[header] = str(record.seq)
    print(f"Loaded {len(sequences)} protein sequences.")
    return sequences

def load_epitope_hla_matrix(epitope_file):
    # Load epitope-HLA binding data as a matrix
    df = pd.read_csv(epitope_file)
    missing = df.isnull().sum().sum()
    if missing > 0:
        print(f"Warning: {missing} missing values detected. Filling with -1.")
        df.fillna(-1, inplace=True)
    print(f"Loaded {df.shape[0]} epitopes with {df.shape[1] - 3} HLA alleles.")
    return df

def load_hla_frequencies(hla_file):
    # Load HLA allele frequencies from CSV
    hla_freq = pd.read_csv(hla_file)
    freq_dict = dict(zip(hla_freq['HLA Allele'], hla_freq['Frequency']))
    print(f"Loaded {len(freq_dict)} HLA allele frequencies.")
    return freq_dict

def main(protein_fasta, epitope_file, hla_file):
    print("Loading inputs...")
    
    sequences = load_protein_sequences(protein_fasta)
    epitope_hla = load_epitope_hla_matrix(epitope_file)
    hla_frequencies = load_hla_frequencies(hla_file)
    
    missing_hla = set(epitope_hla.columns[3:]) - set(hla_frequencies.keys())
    if missing_hla:
        print(f"Warning: {len(missing_hla)} HLA alleles in epitope matrix not found in frequency file:")
        print(missing_hla)
    
    print("All input files processed successfully.")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python script.py <protein_fasta> <epitope_file> <hla_freq_file>")
        sys.exit(1)
    
    protein_fasta = sys.argv[1]
    epitope_file = sys.argv[2]
    hla_file = sys.argv[3]
    
    main(protein_fasta, epitope_file, hla_file)

