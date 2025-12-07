import os
import argparse
import pandas as pd
import numpy as np
import sys

try:
    import pulp
except Exception:
    pulp = None

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
    parser.add_argument("--method", choices=["greedy", "beam", "ilp"], default="greedy", help="Selection method to use.")
    parser.add_argument("--beam_width", type=int, default=100, help="Beam width when method is beam.")
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
        hla_cols = list(epitopes.columns[5:])
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
    DRB1 = {k: v for k, v in hla_frequencies.items() if k.startswith("DRB1")}
    hla_DPA1 = {k: v for k, v in hla_frequencies.items() if k.startswith("HLA-DPA1")}
    hla_DQA1 = {k: v for k, v in hla_frequencies.items() if k.startswith("HLA-DQA1")}

    def normalize(freqs):
        total = sum(freqs.values())
        return {k: v / total for k, v in freqs.items()} if total > 0 else freqs

    return {**normalize(DRB1), **normalize(hla_DPA1), **normalize(hla_DQA1)}

def optimize_population_coverage(epitopes, hla_frequencies, max_epitopes, prioritize_different_proteins, hla_rank_cutoff):
    selected = []
    used_proteins = set()
    epitopes = epitopes.sort_values(by="median_binding_rank")
    for _, ep in epitopes.iterrows():
        coverage = 0.0
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
    print(f"Final selected epitopes (greedy prefilter): {len(selected)}")
    df = pd.DataFrame(selected)
    df = add_hla_class_coverage(df, hla_frequencies, hla_rank_cutoff)
    return df

def add_hla_class_coverage(df, hla_frequencies, hla_rank_cutoff):
    covDR, covDP, covDQ = [], [], []
    for _, ep in df.iterrows():
        A = B = C = 0.0
        for hla, freq in hla_frequencies.items():
            if ep.get(hla, 100) <= hla_rank_cutoff:
                if hla.startswith("DRB1"):
                    A += freq
                elif hla.startswith("HLA-DPA1"):
                    B += freq
                elif hla.startswith("HLA-DQA1"):
                    C += freq
        covDR.append(A)
        covDP.append(B)
        covDQ.append(C)
    df["Coverage_HLA_DRB1"] = covDR
    df["Coverage_HLA_DPA1"] = covDP
    df["Coverage_HLA_DQA1"] = covDQ
    base_cols = ["Epitope", "Position", "ConservationScore", "Coverage_HLA_DRB1", "Coverage_HLA_DPA1", "Coverage_HLA_DQA1"]
    other_cols = [c for c in df.columns if c not in base_cols]
    df = df[base_cols + other_cols]
    return df

def calculate_population_coverage_greedy(epitopes, hla_frequencies, hla_rank_cutoff, max_epitopes):
    epitopes = epitopes.copy()
    individual = []
    for _, ep in epitopes.iterrows():
        A = B = C = 0.0
        for hla, freq in hla_frequencies.items():
            if ep.get(hla, 100) <= hla_rank_cutoff:
                if hla.startswith("DRB1"):
                    A += freq
                elif hla.startswith("HLA-DPA1"):
                    B += freq
                elif hla.startswith("HLA-DQA1"):
                    C += freq
        total = (A + B + C) / 3.0
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
            incA = incB = incC = 0.0
            for hla, freq in hla_frequencies.items():
                if hla in covered:
                    continue
                if ep.get(hla, 100) <= hla_rank_cutoff:
                    if hla.startswith("DRB1"):
                        incA += freq
                    elif hla.startswith("HLA-DPA1"):
                        incB += freq
                    elif hla.startswith("HLA-DQA1"):
                        incC += freq
            inc_total = (incA + incB + incC) / 3.0
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
        cumA = sum(hla_frequencies[h] for h in covered if h.startswith("DRB1"))
        cumB = sum(hla_frequencies[h] for h in covered if h.startswith("HLA-DPA1"))
        cumC = sum(hla_frequencies[h] for h in covered if h.startswith("HLA-DQA1"))
        cum_total = (cumA + cumB + cumC) / 3.0
        results.append({
            "Num_Epitopes": len(selected),
            "Best_Combination": ", ".join(selected),
            "HLA_DRB1_Cov_pct": round(cumA * 100, 2),
            "HLA_DPA1_Cov_pct": round(cumB * 100, 2),
            "HLA_DQA1_Cov_pct": round(cumC * 100, 2),
            "Total_Cov_pct": round(cum_total * 100, 2),
        })
    return pd.DataFrame(results), selected

def run_ilp_progression(epitopes_df, hla_frequencies, hla_rank_cutoff, max_k, prioritize_protein_diversity):
    if pulp is None:
        print("Error: pulp is required for ILP method. Install with `pip install pulp` and ensure CBC solver is available.")
        sys.exit(1)

    epitopes_df = epitopes_df.copy().reset_index(drop=True)
    epi_names = list(epitopes_df["Epitope"])
    num_epi = len(epi_names)
    alleles = list(hla_frequencies.keys())
    num_all = len(alleles)
    a = np.zeros((num_all, num_epi), dtype=int)
    for j, ep_row in epitopes_df.iterrows():
        for i, hla in enumerate(alleles):
            if ep_row.get(hla, 100) <= hla_rank_cutoff:
                a[i, j] = 1

    protein_to_indices = {}
    if "Position" in epitopes_df.columns:
        for idx, prot in enumerate(epitopes_df["Position"]):
            protein_to_indices.setdefault(prot, []).append(idx)

    progression_rows = []
    best_selected_for_k = None

    for s in range(1, max_k + 1):
        prob = pulp.LpProblem(f"maxcov_k{s}", pulp.LpMaximize)
        x_vars = [pulp.LpVariable(f"x_{j}", cat="Binary") for j in range(num_epi)]
        y_vars = [pulp.LpVariable(f"y_{i}", cat="Binary") for i in range(num_all)]
        prob += pulp.lpSum([hla_frequencies[alleles[i]] * y_vars[i] for i in range(num_all)])
        for i in range(num_all):
            prob += y_vars[i] <= pulp.lpSum([a[i, j] * x_vars[j] for j in range(num_epi)])
        prob += pulp.lpSum(x_vars) <= s
        if prioritize_protein_diversity and len(protein_to_indices) > 0:
            for prot, idxs in protein_to_indices.items():
                prob += pulp.lpSum([x_vars[j] for j in idxs]) <= 1

        solver = pulp.COIN_CMD(msg=False, timeLimit=300, path="/opt/homebrew/Cellar/cbc/2.10.12/bin/cbc")
        prob.solve(solver)

        chosen_indices = [j for j in range(num_epi) if pulp.value(x_vars[j]) is not None and pulp.value(x_vars[j]) > 0.5]
        chosen_names = [epi_names[j] for j in chosen_indices]

        covered_alleles = set()
        for j in chosen_indices:
            for i in range(num_all):
                if a[i, j] == 1:
                    covered_alleles.add(alleles[i])

        cumA = sum(hla_frequencies[h] for h in covered_alleles if h.startswith("DRB1"))
        cumB = sum(hla_frequencies[h] for h in covered_alleles if h.startswith("HLA-DPA1"))
        cumC = sum(hla_frequencies[h] for h in covered_alleles if h.startswith("HLA-DQA1"))
        cum_total = (cumA + cumB + cumC) / 3.0

        progression_rows.append({
            "Num_Epitopes": len(chosen_indices),
            "Best_Combination": ", ".join(chosen_names),
            "HLA_DRB1_Cov_pct": round(cumA * 100, 2),
            "HLA_DPA1_Cov_pct": round(cumB * 100, 2),
            "HLA_DQA1_Cov_pct": round(cumC * 100, 2),
            "Total_Cov_pct": round(cum_total * 100, 2),
        })

        if s == max_k:
            best_selected_for_k = chosen_indices

    prog_df = pd.DataFrame(progression_rows)
    return prog_df, best_selected_for_k

def calculate_population_coverage(epitopes, hla_frequencies, hla_rank_cutoff, max_epitopes, method="greedy", beam_width=100, prioritize_protein_diversity=False):
    if method == "greedy":
        prog_df, selected_names = calculate_population_coverage_greedy(epitopes, hla_frequencies, hla_rank_cutoff, max_epitopes)
        return prog_df, selected_names
    elif method == "ilp":
        prog_df, best_indices = run_ilp_progression(epitopes, hla_frequencies, hla_rank_cutoff, max_epitopes, prioritize_protein_diversity)
        selected_names = [epitopes.iloc[i]["Epitope"] for i in best_indices] if best_indices is not None else []
        return prog_df, selected_names
    else:
        print("Beam search not implemented in this script. Use greedy or ilp.")
        sys.exit(1)

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

    filtered = epitopes[epitopes["ConservationScore"] >= args.conservation_score].reset_index(drop=True)
    if filtered.shape[0] == 0:
        print("No epitopes pass conservation filter.")
        sys.exit(1)

    if args.method == "ilp":
        print("Running ILP progression for sizes 1..", args.max_epitopes)
        population_coverage_df, selected_names = calculate_population_coverage(
            filtered, hla_frequencies, args.ba_rank_hit, args.max_epitopes, method="ilp", prioritize_protein_diversity=args.prioritize_protein_diversity
        )
        selected_df = filtered[filtered["Epitope"].isin(selected_names)].copy()
        if not selected_df.empty:
            selected_df = add_hla_class_coverage(selected_df, hla_frequencies, args.ba_rank_hit)
    else:
        print("Running greedy prefilter optimization")
        selected_df = optimize_population_coverage(filtered, hla_frequencies, args.max_epitopes, args.prioritize_protein_diversity, args.ba_rank_hit)
        population_coverage_df, _ = calculate_population_coverage(selected_df, hla_frequencies, args.ba_rank_hit, args.max_epitopes, method="greedy")

    suffix = f"_{args.out_suffix}" if args.out_suffix else ""
    out_dir = args.out_dir

    selected_df.to_csv(os.path.join(out_dir, f"optimized_epitopes{suffix}.csv"), index=False)
    population_coverage_df.to_csv(os.path.join(out_dir, f"population_coverage{suffix}.tsv"), sep="\t", index=False)

    print("Saved optimized epitopes and population coverage progression.")
    print(f"Selected epitopes count: {selected_df.shape[0]}")

if __name__ == "__main__":
    main()
