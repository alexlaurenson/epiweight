import os
import argparse
import pandas as pd
import numpy as np
import sys
import platform
import shutil

try:
    import pulp
except Exception:
    pulp = None

# ----------------- Helper to find CBC solver -----------------
def get_solver_path():
    cbc_path = shutil.which("cbc")
    if cbc_path:
        return cbc_path
    if platform.system().lower() == "darwin":
        for path in ["/opt/homebrew/bin/cbc", "/usr/local/bin/cbc"]:
            if os.path.exists(path):
                return path
    if platform.system().lower() == "windows":
        for path in [r"C:\Program Files\cbc\bin\cbc.exe", r"C:\cbc\bin\cbc.exe"]:
            if os.path.exists(path):
                return path
    return None

# ----------------- Argument parsing -----------------
def parse_args():
    parser = argparse.ArgumentParser(description="Epitope ILP optimization script with HLA clinical weighting.")
    parser.add_argument("epitope_file", help="Path to the epitope data file.")
    parser.add_argument("hla_file", help="Path to the HLA frequencies file.")
    parser.add_argument("--postfilter", action="store_true", help="Indicates the input file has already been filtered and is ready for population coverage.")
    parser.add_argument("--conservation_score", type=float, default=0.8, help="Minimum conservation score.")
    parser.add_argument("--max_epitopes", type=int, default=None, help="Maximum number of epitopes to select.")
    parser.add_argument("--prioritize_protein_diversity", action="store_true", help="Prioritize protein diversity in epitope selection.")
    parser.add_argument("--ba_rank_hit", type=float, default=0.5, help="Binding affinity rank hit threshold.")
    parser.add_argument("--out_dir", default="output", help="Output directory.")
    parser.add_argument("--out_suffix", default="", help="Suffix for output files.")
    parser.add_argument("--high_risk_hlas", default="", help="Comma-separated list of HLA alleles associated with severe malaria outcomes.")
    parser.add_argument("--hla_weight", type=float, default=2.0, help="Weight for high-risk HLAs when using soft prioritization.")
    parser.add_argument("--hard_constraint", action="store_true", help="If set, high-risk HLAs must be covered (hard constraint).")
    return parser.parse_args()

# ----------------- File loading -----------------
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

    return dict(zip(hla_freq.iloc[:, 0], hla_freq.iloc[:, 1]))

# ----------------- Phenotypic HLA coverage -----------------
def phenotypic_frequency(f):
    return 2 * f - f**2

def calculate_locus_coverage(alleles, hla_frequencies):
    # Compute locus coverage using phenotypic frequencies
    q_list = [phenotypic_frequency(hla_frequencies[h]) for h in alleles if h in hla_frequencies]
    if not q_list:
        return 0.0
    return 1 - np.prod([1 - q for q in q_list])

def add_hla_class_coverage(df, hla_frequencies, hla_rank_cutoff):
    coverage_A, coverage_B, coverage_C = [], [], []

    for _, ep in df.iterrows():
        alleles_A = [h for h in hla_frequencies if h.startswith("DRB1") and ep.get(h, 100) <= hla_rank_cutoff]
        alleles_B = [h for h in hla_frequencies if h.startswith("HLA-DPA1") and ep.get(h, 100) <= hla_rank_cutoff]
        alleles_C = [h for h in hla_frequencies if h.startswith("HLA-DQA1") and ep.get(h, 100) <= hla_rank_cutoff]

        coverage_A.append(calculate_locus_coverage(alleles_A, hla_frequencies))
        coverage_B.append(calculate_locus_coverage(alleles_B, hla_frequencies))
        coverage_C.append(calculate_locus_coverage(alleles_C, hla_frequencies))

    df["Coverage_HLA_DRB1"] = coverage_A
    df["Coverage_HLA_DPA1"] = coverage_B
    df["Coverage_HLA_DQA1"] = coverage_C
    base_cols = ["Epitope", "Position", "ConservationScore", "Coverage_HLA_DRB1", "Coverage_HLA_DPA1", "Coverage_HLA_DQA1"]
    other_cols = [c for c in df.columns if c not in base_cols]
    return df[base_cols + other_cols]

# ----------------- ILP progression -----------------
def run_ilp_progression(epitopes_df, hla_frequencies, hla_rank_cutoff, max_k, prioritize_protein_diversity, hla_weights=None, high_risk_hlas=None, hard_constraint=False):
    if pulp is None:
        print("Error: pulp is required for ILP method.")
        sys.exit(1)
    solver_path = get_solver_path()
    if solver_path is None:
        print("Error: Could not locate CBC solver.")
        sys.exit(1)
    solver = pulp.COIN_CMD(msg=False, timeLimit=300, path=solver_path)

    epitopes_df = epitopes_df.reset_index(drop=True)
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
    best_combination_indices = {}
    if max_k is None:
        max_k = num_epi
    prev_total_cov = 0.0

    for s in range(1, max_k + 1):
        prob = pulp.LpProblem(f"maxcov_k{s}", pulp.LpMaximize)
        x_vars = [pulp.LpVariable(f"x_{j}", cat="Binary") for j in range(num_epi)]
        y_vars = [pulp.LpVariable(f"y_{i}", cat="Binary") for i in range(num_all)]

        # Objective: weighted sum of allele coverage
        prob += pulp.lpSum([hla_weights.get(alleles[i], 1.0) * y_vars[i] for i in range(num_all)])

        # Coverage constraints
        for i in range(num_all):
            prob += y_vars[i] <= pulp.lpSum([a[i, j] * x_vars[j] for j in range(num_epi)])
        prob += pulp.lpSum(x_vars) <= s

        # Protein diversity
        if prioritize_protein_diversity and len(protein_to_indices) > 0:
            for prot, idxs in protein_to_indices.items():
                prob += pulp.lpSum([x_vars[j] for j in idxs]) <= 1

        # Hard constraints for high-risk HLAs
        if hard_constraint and high_risk_hlas:
            for hla in high_risk_hlas:
                if hla in alleles:
                    i = alleles.index(hla)
                    prob += y_vars[i] == 1

        prob.solve(solver)

        chosen_indices = [j for j in range(num_epi) if pulp.value(x_vars[j]) is not None and pulp.value(x_vars[j]) > 0.5]
        if not chosen_indices:
            continue

        # Compute locus-specific phenotypic coverage
        covered_alleles = set()
        for j in chosen_indices:
            for i in range(num_all):
                if a[i, j] == 1:
                    covered_alleles.add(alleles[i])

        alleles_A = [h for h in covered_alleles if h.startswith("DRB1")]
        alleles_B = [h for h in covered_alleles if h.startswith("HLA-DPA1")]
        alleles_C = [h for h in covered_alleles if h.startswith("HLA-DQA1")]

        covA = calculate_locus_coverage(alleles_A, hla_frequencies)
        covB = calculate_locus_coverage(alleles_B, hla_frequencies)
        covC = calculate_locus_coverage(alleles_C, hla_frequencies)

        # Combine loci coverage iteratively
        total_cov = covA
        total_cov += (1 - total_cov) * covB
        total_cov += (1 - total_cov) * covC
        pct = round(total_cov * 100, 2)

        best_combination_indices[s] = chosen_indices

        progression_rows.append({
            "Num_Epitopes": s,
            "Best_Combination": ",".join([epi_names[j] for j in chosen_indices]),
            "HLA_DRB1_Cov_pct": round(covA * 100, 2),
            "HLA_DPA1_Cov_pct": round(covB * 100, 2),
            "HLA_DQA1_Cov_pct": round(covC * 100, 2),
            "Total_Cov_pct": pct,
        })

        if pct <= prev_total_cov or pct >= 100.0:
            break
        prev_total_cov = pct

    final_best_indices = best_combination_indices[max(best_combination_indices.keys())]
    return pd.DataFrame(progression_rows), final_best_indices

# ----------------- Main -----------------
def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    epitopes, hla_cols = load_epitopes(args.epitope_file, postfilter=args.postfilter)
    hla_frequencies = load_hla_frequencies(args.hla_file)
    hla_frequencies = {allele: freq for allele, freq in hla_frequencies.items() if allele in hla_cols}
    if len(hla_frequencies) == 0:
        print("Error: None of the HLA frequencies match epitope matrix columns.")
        sys.exit(1)

    high_risk_hlas = [x.strip() for x in args.high_risk_hlas.split(",")] if args.high_risk_hlas else []
    hla_weights = {}
    for hla in hla_frequencies:
        weight = 1.0
        for high_risk in high_risk_hlas:
            if high_risk in hla:
                weight = args.hla_weight
                break
        hla_weights[hla] = weight

    if not args.postfilter:
        print("Conservation processing not implemented outside postfilter.")
        sys.exit(1)

    if "median_binding_rank" not in epitopes.columns:
        epitopes["median_binding_rank"] = epitopes[hla_cols].median(axis=1)

    filtered = epitopes[epitopes["ConservationScore"] >= args.conservation_score].reset_index(drop=True)
    if filtered.shape[0] == 0:
        print("No epitopes pass conservation filter.")
        sys.exit(1)

    print("Running ILP progression")
    population_coverage_df, selected_indices = run_ilp_progression(
        filtered, hla_frequencies, args.ba_rank_hit, args.max_epitopes,
        args.prioritize_protein_diversity, hla_weights, high_risk_hlas, args.hard_constraint
    )

    selected_names = [filtered.iloc[i]["Epitope"] for i in selected_indices]
    selected_df = filtered[filtered["Epitope"].isin(selected_names)].copy()
    selected_df = add_hla_class_coverage(selected_df, hla_frequencies, args.ba_rank_hit)

    suffix = f"_{args.out_suffix}" if args.out_suffix else ""
    out_dir = args.out_dir
    selected_df.to_csv(os.path.join(out_dir, f"optimized_epitopes{suffix}.tsv"), index=False, sep="\t")
    population_coverage_df.to_csv(os.path.join(out_dir, f"population_coverage{suffix}.tsv"), index=False, sep="\t")
    print(f"Saved optimized epitopes and population coverage. Selected count: {selected_df.shape[0]}")

if __name__ == "__main__":
    main()
