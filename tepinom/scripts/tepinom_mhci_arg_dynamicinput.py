import os
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Epitope optimization script.")
    parser.add_argument("epitope_file", help="Path to the epitope data file.")
    parser.add_argument("hla_file", help="Path to the HLA frequencies file.")
    parser.add_argument("--postfilter", action="store_true", help="Indicates the input file has already been filtered and is ready for population coverage.")
    parser.add_argument("--median_ba_rank", type=float, default=1.0, help="Median binding affinity rank cutoff.")
    parser.add_argument("--conservation_score", type=float, default=0.8, help="Minimum conservation score.")
    parser.add_argument("--max_epitopes", type=int, default=10, help="Maximum number of epitopes to select.")
    parser.add_argument("--prioritize_protein_diversity", action="store_true", help="Prioritize protein diversity in epitope selection.")
    parser.add_argument("--ba_rank_hit", type=float, default=0.5, help="Binding affinity rank hit threshold.")
    parser.add_argument("--out_dir", default="output", help="Output directory.")
    parser.add_argument("--out_suffix", default="", help="Suffix for output files.")
    return parser.parse_args()

def load_epitopes(epitope_file, delimiter=None, postfilter=False):
    try:
        if delimiter is None:
            delimiter = '\t' if postfilter else (',' if epitope_file.endswith('.csv') else '\t')
        epitopes = pd.read_csv(epitope_file, delimiter=delimiter, low_memory=False)
        print(f"Loaded {len(epitopes)} epitopes with {len(epitopes.columns)} columns.")
    except Exception as e:
        print(f"Error loading epitope file: {e}")
        sys.exit(1)

    if postfilter:
        required_columns = ['Epitope', 'Position', 'ConservationScore', 'Source', 'Stage']
        if not all(col in epitopes.columns[:5] for col in required_columns):
            print("Error: Input missing required postfilter columns.")
            sys.exit(1)
        hla_cols = epitopes.columns[5:]
    else:
        required_columns = ['Epitope', 'Position', 'ConservationScore']
        if not all(col in epitopes.columns for col in required_columns):
            print("Error: Input missing required columns.")
            sys.exit(1)
        hla_cols = [col for col in epitopes.columns if col not in required_columns]

    if len(hla_cols) == 0:
        print("Error: No HLA allele columns found.")
        sys.exit(1)

    return epitopes, hla_cols

def load_hla_frequencies(hla_file):
    try:
        hla_freq = pd.read_csv(hla_file, delimiter='\t', header=None)
        print(f"Loaded HLA frequencies for {len(hla_freq)} alleles.")
    except Exception as e:
        print(f"Error loading HLA frequency file: {e}")
        sys.exit(1)

    freq_dict = dict(zip(hla_freq.iloc[:, 0], hla_freq.iloc[:, 1]))
    return freq_dict

def normalize_hla_frequencies(hla_frequencies):
    hla_A = {k: v for k, v in hla_frequencies.items() if k.startswith("HLA-A")}
    hla_B = {k: v for k, v in hla_frequencies.items() if k.startswith("HLA-B")}
    hla_C = {k: v for k, v in hla_frequencies.items() if k.startswith("HLA-C")}

    def normalize(freqs):
        total = sum(freqs.values())
        return {k: v / total for k, v in freqs.items()} if total > 0 else freqs

    return {**normalize(hla_A), **normalize(hla_B), **normalize(hla_C)}

def optimize_population_coverage(epitopes, hla_frequencies, max_epitopes, prioritize_different_proteins, hla_rank_cutoff):
    selected = []
    used_proteins = set()

    epitopes = epitopes.sort_values(by="median_binding_rank")

    for _, ep in epitopes.iterrows():
        coverage = 0
        for hla in hla_frequencies:
            if ep.get(hla, 100) <= hla_rank_cutoff:
                coverage += hla_frequencies[hla]

        if coverage > 0:
            if prioritize_different_proteins and ep["Position"] in used_proteins:
                continue

            selected.append(ep)
            used_proteins.add(ep["Position"])

            if len(selected) >= max_epitopes:
                break

    print(f"Final selected epitopes: {len(selected)}")
    df = pd.DataFrame(selected)

    df = add_hla_class_coverage(df, hla_frequencies, hla_rank_cutoff)

    return df

def add_hla_class_coverage(df, hla_frequencies, hla_rank_cutoff):
    covA = []
    covB = []
    covC = []

    for _, ep in df.iterrows():
        A, B, C = 0, 0, 0
        for hla, freq in hla_frequencies.items():
            if ep.get(hla, 100) <= hla_rank_cutoff:
                if hla.startswith("HLA A") or hla.startswith("HLA-A"):
                    A += freq
                elif hla.startswith("HLA B") or hla.startswith("HLA-B"):
                    B += freq
                elif hla.startswith("HLA C") or hla.startswith("HLA-C"):
                    C += freq
        covA.append(A)
        covB.append(B)
        covC.append(C)

    df["Coverage_HLA_A"] = covA
    df["Coverage_HLA_B"] = covB
    df["Coverage_HLA_C"] = covC

    base_cols = ["Epitope", "Position", "ConservationScore", "Coverage_HLA_A", "Coverage_HLA_B", "Coverage_HLA_C"]
    other_cols = [c for c in df.columns if c not in base_cols]
    df = df[base_cols + other_cols]

    return df

def calculate_population_coverage(epitopes, hla_frequencies, hla_rank_cutoff, max_epitopes):
    epitopes = epitopes.copy()
    all_hla = list(hla_frequencies.keys())

    individual = []
    for _, ep in epitopes.iterrows():
        A = B = C = 0
        for hla, freq in hla_frequencies.items():
            if ep.get(hla, 100) <= hla_rank_cutoff:
                if hla.startswith("HLA-A"):
                    A += freq
                elif hla.startswith("HLA-B"):
                    B += freq
                elif hla.startswith("HLA-C"):
                    C += freq
        total = (A + B + C) / 3
        individual.append((ep["Epitope"], A, B, C, total))

    individual.sort(key=lambda x: x[4], reverse=True)

    selected = []
    covered = set()
    results = []

    for step in range(max_epitopes):
        best = None
        best_inc = -1

        for name, A, B, C, total in individual:
            if name in selected:
                continue

            ep = epitopes[epitopes["Epitope"] == name].iloc[0]

            incA = incB = incC = 0
            for hla, freq in hla_frequencies.items():
                if hla in covered:
                    continue
                if ep.get(hla, 100) <= hla_rank_cutoff:
                    if hla.startswith("HLA-A"):
                        incA += freq
                    elif hla.startswith("HLA-B"):
                        incB += freq
                    elif hla.startswith("HLA-C"):
                        incC += freq

            inc_total = (incA + incB + incC) / 3
            if inc_total > best_inc:
                best = name
                best_inc = inc_total

        if best is None:
            break

        selected.append(best)

        ep_row = epitopes[epitopes["Epitope"] == best].iloc[0]
        for hla in hla_frequencies:
            if ep_row.get(hla, 100) <= hla_rank_cutoff:
                covered.add(hla)

        cumA = sum(hla_frequencies[h] for h in covered if h.startswith("HLA-A"))
        cumB = sum(hla_frequencies[h] for h in covered if h.startswith("HLA-B"))
        cumC = sum(hla_frequencies[h] for h in covered if h.startswith("HLA-C"))
        cum_total = (cumA + cumB + cumC) / 3

        results.append({
            "Num_Epitopes": len(selected),
            "Best_Combination": ", ".join(selected),
            "HLA_A_Cov_pct": round(cumA * 100, 2),
            "HLA_B_Cov_pct": round(cumB * 100, 2),
            "HLA_C_Cov_pct": round(cumC * 100, 2),
            "Total_Cov_pct": round(cum_total * 100, 2),
        })

    return pd.DataFrame(results)

def main():
    args = parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    epitopes, hla_cols = load_epitopes(args.epitope_file, postfilter=args.postfilter)
    hla_frequencies = load_hla_frequencies(args.hla_file)

    hla_frequencies = {allele: freq for allele, freq in hla_frequencies.items() if allele in hla_cols}

    if len(hla_frequencies) == 0:
        print("Error: None of the HLA frequencies match epitope matrix columns.")
        sys.exit(1)

    hla_frequencies = normalize_hla_frequencies(hla_frequencies)


    if not args.postfilter:
        print("Conservation processing not implemented outside postfilter.")
        sys.exit(1)

    if "median_binding_rank" not in epitopes.columns:
        epitopes["median_binding_rank"] = epitopes[hla_cols].median(axis=1)

    filtered = epitopes[epitopes["ConservationScore"] >= args.conservation_score]

    optimized = optimize_population_coverage(
        filtered,
        hla_frequencies,
        args.max_epitopes,
        args.prioritize_protein_diversity,
        args.ba_rank_hit
    )

    coverage = calculate_population_coverage(
        optimized,
        hla_frequencies,
        args.ba_rank_hit,
        args.max_epitopes
    )

    suffix = f"_{args.out_suffix}" if args.out_suffix else ""

    optimized.to_csv(os.path.join(args.out_dir, f"optimized_epitopes{suffix}.csv"), index=False)
    coverage.to_csv(os.path.join(args.out_dir, f"population_coverage{suffix}.tsv"), sep="\t", index=False)

if __name__ == "__main__":
    main()
