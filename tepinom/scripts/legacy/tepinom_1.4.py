#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np
import sys
import platform
import shutil
import re

try:
    import pulp
except Exception:
    pulp = None

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

def parse_args():
    parser = argparse.ArgumentParser(description="Epitope ILP optimization with exact PopCover coverage reporting.")
    parser.add_argument("epitope_file", help="Path to the epitope data file.")
    parser.add_argument("hla_file", help="Path to the HLA frequencies file.")
    parser.add_argument("--mhc_class", choices=["i","ii"], default="i", help="MHC class for coverage calculation (i or ii).")
    parser.add_argument("--postfilter", action="store_true")
    parser.add_argument("--conservation_score", type=float, default=0.8)
    parser.add_argument("--max_epitopes", type=int, default=None)
    parser.add_argument("--prioritize_protein_diversity", action="store_true")
    parser.add_argument("--ba_rank_hit", type=float, default=0.5)
    parser.add_argument("--out_dir", default="output")
    parser.add_argument("--out_suffix", default="")
    parser.add_argument("--high_risk_hlas", default="")
    parser.add_argument("--hla_weight", type=float, default=2.0)
    parser.add_argument("--hard_constraint", action="store_true")
    return parser.parse_args()

def load_epitopes(epitope_file, postfilter=False):
    try:
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
    freq_dict = {}
    for i in range(hla_freq.shape[0]):
        allele = str(hla_freq.iloc[i, 0]).strip()
        try:
            freq = float(hla_freq.iloc[i, 1])
        except Exception:
            print(f"Warning: frequency for allele {allele} is not numeric. Setting to 0.")
            freq = 0.0
        freq_dict[allele] = freq
    return freq_dict

def phenotypic_frequency(f):
    # P(allele present) = 1 - (1 - f)^2 = 2f - f^2
    return 2.0 * f - f * f

def alleles_to_locus(allele):
    s = str(allele).strip().upper()
    if s.startswith("HLA-A"): return "HLA-A"
    if s.startswith("HLA-B"): return "HLA-B"
    if s.startswith("HLA-C"): return "HLA-C"
    if "DRB1" in s: return "DRB1"
    if s.startswith("HLA-DPA1"): return "HLA-DPA1"
    if s.startswith("HLA-DQA1"): return "HLA-DQA1"
    m = re.search(r'D[A-Z]+[0-9]', s)
    return m.group(0) if m else None

def calculate_population_coverage_exact(selected_indices, epitopes_df, hla_frequencies, hla_cols, rank_cutoff, mhc_class):
    """
    PopCover-style exact intra-locus and inter-locus coverage calculation.
    selected_indices: list of integer row indices in epitopes_df
    """
    allele_list = []
    for idx in selected_indices:
        row = epitopes_df.iloc[idx]
        for hla in hla_cols:
            val = row.get(hla, np.nan)
            if not pd.isna(val):
                try:
                    if float(val) <= rank_cutoff:
                        allele_list.append(hla)
                except Exception:
                    # non-numeric cell -> skip
                    pass

    if mhc_class.lower() == "i":
        loci_list = ["HLA-A", "HLA-B", "HLA-C"]
    else:
        loci_list = ["DRB1", "HLA-DPA1", "HLA-DQA1"]

    loci_map = {locus: [] for locus in loci_list}
    for allele in allele_list:
        locus = alleles_to_locus(allele)
        if locus and locus in loci_map:
            loci_map[locus].append(allele)

    intra_covg = {}
    for locus, alleles_in_locus in loci_map.items():
        freqs = []
        for a in alleles_in_locus:
            if a in hla_frequencies:
                f = float(hla_frequencies[a])
                freqs.append(phenotypic_frequency(f))
        if freqs:
            intra_covg[locus] = 1.0 - np.prod([1.0 - q for q in freqs])
        else:
            intra_covg[locus] = 0.0

    inter_covg = 1.0 - np.prod([1.0 - c for c in intra_covg.values()]) if intra_covg else 0.0
    return intra_covg, inter_covg

def run_ilp_progression(epitopes_df, hla_frequencies, hla_cols, rank_cutoff, max_k,
                        prioritize_protein_diversity, hla_weights=None, high_risk_hlas=None,
                        hard_constraint=False, mhc_class="i"):
    """
    ILP that selects up to s epitopes to maximize weighted allele coverage.
    We already bake in:
      - hla_weights: multiplicative factor applied to allele contribution if it matches a high-risk HLA
      - prioritize_protein_diversity: enforced as a constraint (at most one epitope per protein group)
      - hard_constraint: requires coverage of specified HLAs if present
    """

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

    # Only consider alleles for which we have frequencies
    alleles = list(hla_frequencies.keys())
    num_all = len(alleles)

    # build binary coverage matrix a[i,j] = 1 if allele i is covered by epitope j
    a = np.zeros((num_all, num_epi), dtype=int)
    for j, ep_row in epitopes_df.iterrows():
        for i, hla in enumerate(alleles):
            if hla in epitopes_df.columns:
                val = ep_row.get(hla, np.nan)
                if not pd.isna(val):
                    try:
                        if float(val) <= rank_cutoff:
                            a[i, j] = 1
                    except Exception:
                        pass

    # allele linear weights: use -log(1 - phenofreq) as linear approx for multiplicative merging
    allele_linear_weights = {}
    for i, hla in enumerate(alleles):
        f = float(hla_frequencies.get(hla, 0.0))
        q = phenotypic_frequency(f)
        # small epsilon avoid log(0)
        weight = -np.log(1.0 - q + 1e-12)
        # apply HLA-specific multiplicative weight if it matches any high_risk pattern
        if high_risk_hlas:
            for hr in high_risk_hlas:
                if hr and hr in hla:
                    weight *= (hla_weights.get(hla, 1.0) if hla_weights else 1.0)
                    break
        allele_linear_weights[hla] = weight

    # --- Protein diversity grouping detection ---
    protein_to_indices = {}
    protein_column_candidates = ["Source", "Protein", "Gene", "ProteinID"]
    prot_col = None
    for c in protein_column_candidates:
        if c in epitopes_df.columns:
            prot_col = c
            break
    if prioritize_protein_diversity:
        if prot_col is None:
            print("Error: --prioritize_protein_diversity requested but no protein identifier column found (checked Source/Protein/Gene/ProteinID).")
            sys.exit(1)
        for idx, prot in enumerate(epitopes_df[prot_col]):
            protein_to_indices.setdefault(str(prot), []).append(idx)

    progression_rows = []
    best_combination_indices = {}
    if max_k is None:
        max_k = num_epi
    prev_total_cov = 0.0

    # We will iterate s from 1..max_k
    for s in range(1, max_k + 1):
        prob = pulp.LpProblem(f"maxcov_k{s}", pulp.LpMaximize)

        # decision vars
        x_vars = [pulp.LpVariable(f"x_{j}", cat="Binary") for j in range(num_epi)]  # pick epitope j
        y_vars = [pulp.LpVariable(f"y_{i}", cat="Binary") for i in range(num_all)]  # allele i covered

        # objective: sum over alleles weight * y_i
        prob += pulp.lpSum([allele_linear_weights[alleles[i]] * y_vars[i] for i in range(num_all)])

        # link allele coverage to selected epitopes
        for i in range(num_all):
            prob += y_vars[i] <= pulp.lpSum([a[i, j] * x_vars[j] for j in range(num_epi)])

        # limit number of epitopes
        prob += pulp.lpSum(x_vars) <= s

        # enforce at-most-one per protein if requested (hard constraint)
        if prioritize_protein_diversity and protein_to_indices:
            for prot, idxs in protein_to_indices.items():
                if len(idxs) > 0:
                    prob += pulp.lpSum([x_vars[j] for j in idxs]) <= 1

        # hard constraint: require coverage for particular HLAs (if requested)
        if hard_constraint and high_risk_hlas:
            for hr in high_risk_hlas:
                # find indices of alleles that match hr substring exactly or as prefix
                matched = [i for i, a_name in enumerate(alleles) if hr in a_name]
                for i in matched:
                    prob += y_vars[i] == 1

        # solve
        prob.solve(solver)

        chosen_indices = [j for j in range(num_epi) if pulp.value(x_vars[j]) is not None and pulp.value(x_vars[j]) > 0.5]
        if not chosen_indices:
            # no feasible or nothing chosen
            continue

        # compute PopCover-style exact coverage for chosen set
        intra_cov, total_cov = calculate_population_coverage_exact(
            chosen_indices, epitopes_df, hla_frequencies, hla_cols, rank_cutoff, mhc_class
        )
        pct = round(total_cov * 100.0, 2)
        best_combination_indices[s] = chosen_indices

        # Prepare row: ensure locus columns ordered (HLA-A, HLA-B, HLA-C) / (DRB1,HLA-DPA1,HLA-DQA1)
        row = {"Num_Epitopes": s, "Best_Combination": ",".join([epi_names[j] for j in chosen_indices])}
        loci_order = ["HLA-A", "HLA-B", "HLA-C"] if mhc_class.lower() == "i" else ["DRB1", "HLA-DPA1", "HLA-DQA1"]
        for locus in loci_order:
            row[f"{locus}_Cov_pct"] = round(intra_cov.get(locus, 0.0) * 100.0, 4)
        row["Total_Cov_pct"] = pct

        progression_rows.append(row)

        if pct <= prev_total_cov or pct >= 100.0:
            break
        prev_total_cov = pct

    # Build dataframe and ensure columns order with Total_Cov_pct last
    if progression_rows:
        df = pd.DataFrame(progression_rows)
        # desired column order
        base_cols = ["Num_Epitopes", "Best_Combination"]
        loci_cols = [c for c in df.columns if c.endswith("_Cov_pct") and c != "Total_Cov_pct"]
        # try to sort loci_cols by preferred order
        preferred = ["HLA-A_Cov_pct", "HLA-B_Cov_pct", "HLA-C_Cov_pct", "DRB1_Cov_pct", "HLA-DPA1_Cov_pct", "HLA-DQA1_Cov_pct"]
        sorted_loci = [c for c in preferred if c in loci_cols] + [c for c in loci_cols if c not in preferred]
        final_cols = base_cols + sorted_loci + ["Total_Cov_pct"]
        df = df[final_cols]
    else:
        # empty
        df = pd.DataFrame(columns=["Num_Epitopes", "Best_Combination", "HLA-A_Cov_pct", "HLA-B_Cov_pct", "HLA-C_Cov_pct", "Total_Cov_pct"])

    final_best_indices = best_combination_indices[max(best_combination_indices.keys())] if best_combination_indices else []
    return df, final_best_indices

def prepend_command_header(filepath, command):
    try:
        with open(filepath, "r") as f:
            original = f.read()
    except FileNotFoundError:
        original = ""
    header = f"# Command used to run this script:\n# {command}\n\n"
    with open(filepath, "w") as f:
        f.write(header + original)

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    command_used = " ".join(sys.argv)

    epitopes, hla_cols = load_epitopes(args.epitope_file, postfilter=args.postfilter)
    hla_frequencies = load_hla_frequencies(args.hla_file)
    # keep only frequencies that appear in epitopes matrix columns
    hla_frequencies = {allele: float(freq) for allele, freq in hla_frequencies.items() if allele in hla_cols}
    if len(hla_frequencies) == 0:
        print("Error: None of the HLA frequencies match epitope matrix columns.")
        sys.exit(1)

    # parse high risk list
    high_risk_hlas = [x.strip() for x in args.high_risk_hlas.split(",") if x.strip()] if args.high_risk_hlas else []
    # hla_weights mapping - use value only for alleles matching any high_risk pattern
    hla_weights = {}
    for allele in hla_frequencies:
        w = 1.0
        for hr in high_risk_hlas:
            if hr and hr in allele:
                w = args.hla_weight
                break
        hla_weights[allele] = w

    if not args.postfilter:
        print("Conservation processing not implemented outside postfilter.")
        sys.exit(1)

    filtered = epitopes[epitopes["ConservationScore"] >= args.conservation_score].reset_index(drop=True)
    if filtered.shape[0] == 0:
        print("No epitopes pass conservation filter.")
        sys.exit(1)

    print("Running ILP progression (priorities baked into objective/constraints)...")
    population_coverage_df, selected_indices = run_ilp_progression(
        filtered, hla_frequencies, hla_cols, args.ba_rank_hit, args.max_epitopes,
        args.prioritize_protein_diversity, hla_weights, high_risk_hlas,
        args.hard_constraint, args.mhc_class
    )

    # --- Build optimized_epitopes output: union of epitopes mentioned in Best_Combination rows ---
    best_epitopes = set()
    if "Best_Combination" in population_coverage_df.columns:
        for comb in population_coverage_df["Best_Combination"].fillna(""):
            for e in comb.split(","):
                e = e.strip()
                if e:
                    best_epitopes.add(e)

    epitope_subset = filtered[filtered["Epitope"].isin(best_epitopes)].copy()
    # ensure requested metadata columns exist (Position, Source, Stage)
    for col in ["Position", "Source", "Stage", "ConservationScore"]:
        if col not in epitope_subset.columns:
            epitope_subset[col] = np.nan

    # add per-epitope intra-locus and total cov columns
    if not epitope_subset.empty:
        loci_order = ["HLA-A", "HLA-B", "HLA-C"] if args.mhc_class.lower() == "i" else ["DRB1", "HLA-DPA1", "HLA-DQA1"]
        # compute coverage for each epitope individually
        for idx in epitope_subset.index:
            # need the original integer index in filtered: idx is already aligned because epitope_subset came from filtered with copy()
            intra_cov, total_cov = calculate_population_coverage_exact([idx], filtered, hla_frequencies, hla_cols, args.ba_rank_hit, args.mhc_class)
            for locus in loci_order:
                epitope_subset.loc[idx, f"{locus}_Cov_pct"] = round(intra_cov.get(locus, 0.0) * 100.0, 4)
            epitope_subset.loc[idx, "Total_Cov_pct"] = round(total_cov * 100.0, 2)

    # reorder columns to: Epitope | Position | ConservationScore | Source | Stage | <per-locus coverage> | Total_Cov_pct | <all HLA columns>
    locus_cols = []
    if args.mhc_class.lower() == "i":
        locus_cols = ["HLA-A_Cov_pct", "HLA-B_Cov_pct", "HLA-C_Cov_pct"]
    else:
        locus_cols = ["DRB1_Cov_pct", "HLA-DPA1_Cov_pct", "HLA-DQA1_Cov_pct"]

    # ensure locus cols present (if absent create zero)
    for lc in locus_cols:
        if lc not in epitope_subset.columns:
            epitope_subset[lc] = 0.0

    # ensure all HLA columns appear (order maintained as in hla_cols)
    hla_cols_in_df = [c for c in hla_cols if c in epitope_subset.columns]

    ordered_cols = ["Epitope", "Position", "ConservationScore", "Source", "Stage"] + locus_cols + ["Total_Cov_pct"] + hla_cols_in_df

    # Add any missing requested columns at end if they are not present in the dataframe
    for c in ordered_cols:
        if c not in epitope_subset.columns:
            epitope_subset[c] = np.nan

    epitope_subset = epitope_subset[ordered_cols]

    # file paths
    suffix = f"_{args.out_suffix}" if args.out_suffix else ""
    out_dir = args.out_dir
    optimized_path = os.path.join(out_dir, f"optimized_epitopes{suffix}.tsv")
    popcov_path = os.path.join(out_dir, f"population_coverage{suffix}.tsv")

    # write outputs
    epitope_subset.to_csv(optimized_path, index=False, sep="\t")
    population_coverage_df.to_csv(popcov_path, index=False, sep="\t")

    # prepend header with command used
    prepend_command_header(optimized_path, command_used)
    prepend_command_header(popcov_path, command_used)

    print(f"Saved optimized epitopes and population coverage. Selected count: {epitope_subset.shape[0]}")

if __name__ == "__main__":
    main()
