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
    parser.add_argument("--goal", choices=["combo", "region"], default="combo",
                        help="Optimization goal: 'combo' = epitope combinations; 'region' = find regions (by protein/position) that maximize coverage.")
    parser.add_argument("--window_size", type=int, default=None,
                        help="When --goal region: optional sliding window size (aa). If omitted, minimal contiguous region spanning chosen epitopes is used.")
    parser.add_argument("--best_region_num", type=int, default=10,
                    help="Number of top regions to include in output (default 10). Sorted with best at bottom.")
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
        if idx is None:
            continue
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
    Works on the provided epitopes_df (which can be a subset of the original filtered dataframe).
    Returns (population_coverage_dataframe, indices_of_best_combination_relative_to_epitopes_df)
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
    if "Epitope" not in epitopes_df.columns:
        print("Error: epitopes_df must contain 'Epitope' column.")
        sys.exit(1)
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

def make_region_row(protein_name, start_pos, end_pos, region_type, epi_names, intra_cov, total_cov):
    row = {
        "Protein": protein_name,
        "Region_Start": int(start_pos) if start_pos is not None else None,
        "Region_End": int(end_pos) if end_pos is not None else None,
        "Region_Type": region_type,
        "Num_Epitopes": len(epi_names),
        "Epitopes": ",".join(epi_names),
        "Total_Cov_pct": round(total_cov * 100.0, 2)
    }
    # attach per-locus if present in intra_cov
    for locus, val in intra_cov.items():
        row[f"{locus}_Cov_pct"] = round(val * 100.0, 4)
    return row

def run_region_search_by_window(filtered_df, hla_frequencies, hla_cols, args):
    """
    Modified sliding-window search:
    - For each protein, evaluate windows from size 1 to window_size
    - Compute coverage with ILP for epitopes inside each window
    - For each protein:
        • Keep all windows reaching 100% total coverage, OR
        • If none reach 100%, keep the top 10 windows by total coverage
    - Return region_df containing these best regions per protein
    - Return best_region_epitopes_df (best region globally)
    """

    if args.window_size is None:
        raise ValueError("window_size is required for sliding window mode.")

    protein_column_candidates = ["Source", "Protein", "Gene", "ProteinID"]
    prot_col = next((c for c in protein_column_candidates if c in filtered_df.columns), None)
    if prot_col is None:
        print("Error: No protein identifier column found for region analysis.")
        sys.exit(1)

    all_region_rows = []           # all regions from all proteins
    best_overall = {"total_cov": -1.0}   # globally best region
    proteins = filtered_df[prot_col].unique()

    # ============================================================
    # MAIN LOOP: PROCESS EACH PROTEIN
    # ============================================================
    for prot in proteins:
        subset = filtered_df[filtered_df[prot_col] == prot].copy()
        if subset.empty or "Position" not in subset.columns:
            continue

        # ensure positions are integers
        try:
            subset["Position"] = subset["Position"].astype(int)
        except Exception:
            subset["Position"] = subset["Position"].apply(
                lambda x: int(re.sub(r'[^0-9]', '', str(x)))
                if pd.notna(x) and re.search(r'\d', str(x))
                else np.nan
            )
        subset = subset.dropna(subset=["Position"])
        if subset.empty:
            continue

        min_pos = int(subset["Position"].min())
        max_pos = int(subset["Position"].max())

        # Store all windows from this protein
        prot_rows = []

        # ============================================================
        # SLIDING WINDOWS
        # ============================================================
        for win_size in range(1, args.window_size + 1):
            step_end = max_pos - win_size + 1
            window_starts = [min_pos] if step_end < min_pos else list(range(min_pos, step_end + 1))

            for start in window_starts:
                end = start + win_size - 1

                in_window = subset[
                    (subset["Position"] >= start) &
                    (subset["Position"] <= end)
                ].reset_index(drop=True)

                if in_window.empty:
                    continue

                # Run ILP
                try:
                    popcov_df, chosen_idx_relative = run_ilp_progression(
                        in_window, hla_frequencies, hla_cols, args.ba_rank_hit, args.max_epitopes,
                        args.prioritize_protein_diversity, None,
                        [x.strip() for x in args.high_risk_hlas.split(",") if x.strip()],
                        args.hard_constraint, args.mhc_class
                    )
                except Exception as e:
                    print(f"Warning: ILP failed for protein {prot} window {start}-{end}: {e}")
                    continue

                if not chosen_idx_relative:
                    continue

                # Coverage
                intra_cov, total_cov = calculate_population_coverage_exact(
                    chosen_idx_relative, in_window,
                    hla_frequencies, hla_cols,
                    args.ba_rank_hit, args.mhc_class
                )

                epi_names = [in_window.iloc[i]["Epitope"] for i in chosen_idx_relative]

                row = make_region_row(
                    prot, start, end, f"sliding_{win_size}",
                    epi_names, intra_cov, total_cov
                )
                prot_rows.append(row)

                # track global best region
                if total_cov > best_overall["total_cov"]:
                    best_overall = {
                        "total_cov": total_cov,
                        "protein": prot,
                        "start": start,
                        "end": end,
                        "chosen_indices_relative": chosen_idx_relative,
                        "in_window_df": in_window.copy(),
                        "intra_cov": intra_cov
                    }

        # ============================================================
        # SELECT PER-PROTEIN BEST REGIONS
        # ============================================================
        if not prot_rows:
            continue

        dfp = pd.DataFrame(prot_rows)

        # 1. Try to get all regions reaching 100%
        perfect = dfp[dfp["Total_Cov_pct"] >= 100.0]

        if len(perfect) > 0:
            # Keep ALL that hit 100 coverage
            best_for_protein = perfect.sort_values("Total_Cov_pct", ascending=False)
        else:
            # 2. Otherwise keep top 10 for this protein
            best_for_protein = dfp.sort_values("Total_Cov_pct", ascending=False).head(10)

        # Add to global list
        all_region_rows.append(best_for_protein)

    # ============================================================
    # FINAL region_df: union of per-protein best sets
    # ============================================================
    if len(all_region_rows) == 0:
        region_df = pd.DataFrame(columns=[
            "Protein","Region_Start","Region_End","Region_Type",
            "Num_Epitopes","Epitopes","Total_Cov_pct"
        ])
    else:
        region_df = pd.concat(all_region_rows, ignore_index=True)

    # ============================================================
    # best_region_epitopes_df: same logic as before (global best)
    # ============================================================
    if best_overall["total_cov"] < 0:
        best_epitope_df = pd.DataFrame(columns=filtered_df.columns)
    else:
        chosen = best_overall["chosen_indices_relative"]
        best_epitope_df = best_overall["in_window_df"].iloc[chosen].copy().reset_index(drop=True)

        loci_order = (
            ["HLA-A", "HLA-B", "HLA-C"]
            if args.mhc_class.lower() == "i"
            else ["DRB1", "HLA-DPA1", "HLA-DQA1"]
        )

        # Compute coverage per epitope
        for idx in best_epitope_df.index:
            intra_cov_each, total_cov_each = calculate_population_coverage_exact(
                [idx], best_epitope_df, hla_frequencies, hla_cols,
                args.ba_rank_hit, args.mhc_class
            )
            for locus in loci_order:
                best_epitope_df.loc[idx, f"{locus}_Cov_pct"] = round(
                    intra_cov_each.get(locus, 0.0) * 100.0, 4
                )
            best_epitope_df.loc[idx, "Total_Cov_pct"] = round(total_cov_each * 100.0, 2)

    return region_df, best_epitope_df

def run_region_search_minimal_span(filtered_df, hla_frequencies, hla_cols, args):
    """
    Minimal contiguous region approach (no window_size provided).
    For each protein, run ILP restricted to epitopes in that protein for s=1..max_k.
    For each chosen set, compute the minimal contiguous span (min start to max end across chosen epitopes).
    Record the region and coverage; pick best across proteins.
    """
    protein_column_candidates = ["Source", "Protein", "Gene", "ProteinID"]
    prot_col = None
    for c in protein_column_candidates:
        if c in filtered_df.columns:
            prot_col = c
            break
    if prot_col is None:
        print("Error: No protein identifier column found for region analysis.")
        sys.exit(1)

    region_rows = []
    best_overall = {"total_cov": -1.0}
    proteins = filtered_df[prot_col].unique()

    for prot in proteins:
        subset = filtered_df[filtered_df[prot_col] == prot].copy()
        if subset.empty:
            continue
        # ensure Position numeric
        if "Position" not in subset.columns:
            continue
        try:
            subset["Position"] = subset["Position"].astype(int)
        except Exception:
            subset["Position"] = subset["Position"].apply(lambda x: int(re.sub(r'[^0-9]', '', str(x))) if pd.notna(x) and re.search(r'\d', str(x)) else np.nan)
        subset = subset.dropna(subset=["Position"])
        if subset.empty:
            continue

        # run ILP on subset (this will itself iterate s=1..max_k)
        try:
            popcov_df, chosen_indices_rel = run_ilp_progression(
                subset, hla_frequencies, hla_cols, args.ba_rank_hit, args.max_epitopes,
                args.prioritize_protein_diversity, None, [x.strip() for x in args.high_risk_hlas.split(",") if x.strip()],
                args.hard_constraint, args.mhc_class
            )
        except Exception as e:
            print(f"Warning: ILP failed for protein {prot} in minimal-span mode: {e}")
            continue

        if not chosen_indices_rel:
            continue

        # compute the minimal contiguous region covering chosen epitopes
        chosen_positions = []
        chosen_epitopes = []
        chosen_rows = subset.reset_index(drop=True).iloc[chosen_indices_rel]
        for idx_rel, row in chosen_rows.iterrows():
            start = int(row["Position"])
            # estimate end as start + len(epitope) - 1 if Epitope sequence present
            seq = str(row.get("Epitope", ""))
            end = start + max(len(seq) - 1, 0)
            chosen_positions.append((start, end))
            chosen_epitopes.append(row["Epitope"])

        region_start = min(s for s, e in chosen_positions)
        region_end = max(e for s, e in chosen_positions)

        intra_cov, total_cov = calculate_population_coverage_exact(chosen_indices_rel, subset.reset_index(drop=True), hla_frequencies, hla_cols, args.ba_rank_hit, args.mhc_class)

        row = make_region_row(prot, region_start, region_end, "minimal_span", chosen_epitopes, intra_cov, total_cov)
        region_rows.append(row)

        if total_cov > best_overall["total_cov"]:
            best_overall = {
                "total_cov": total_cov,
                "protein": prot,
                "start": region_start,
                "end": region_end,
                "chosen_indices_relative": chosen_indices_rel,
                "subset_df": subset.reset_index(drop=True).copy(),
                "intra_cov": intra_cov
            }

    if not region_rows:
        region_df = pd.DataFrame(columns=["Protein","Region_Start","Region_End","Region_Type","Num_Epitopes","Epitopes","Total_Cov_pct"])
    else:
        region_df = pd.DataFrame(region_rows)

    if best_overall.get("total_cov", -1.0) < 0:
        best_epitope_df = pd.DataFrame(columns=filtered_df.columns)
    else:
        chosen = best_overall["chosen_indices_relative"]
        best_epitope_df = best_overall["subset_df"].iloc[chosen].copy().reset_index(drop=True)
        loci_order = ["HLA-A", "HLA-B", "HLA-C"] if args.mhc_class.lower() == "i" else ["DRB1", "HLA-DPA1", "HLA-DQA1"]
        for idx in best_epitope_df.index:
            intra_cov_each, total_cov_each = calculate_population_coverage_exact([idx], best_epitope_df, hla_frequencies, hla_cols, args.ba_rank_hit, args.mhc_class)
            for locus in loci_order:
                best_epitope_df.loc[idx, f"{locus}_Cov_pct"] = round(intra_cov_each.get(locus, 0.0)*100.0, 4)
            best_epitope_df.loc[idx, "Total_Cov_pct"] = round(total_cov_each*100.0, 2)

    return region_df, best_epitope_df

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    command_used = " ".join(sys.argv)

    epitopes, hla_cols = load_epitopes(args.epitope_file, postfilter=args.postfilter)
    hla_frequencies = load_hla_frequencies(args.hla_file)
    hla_frequencies = {allele: float(freq) for allele, freq in hla_frequencies.items() if allele in hla_cols}
    if len(hla_frequencies) == 0:
        print("Error: None of the HLA frequencies match epitope matrix columns.")
        sys.exit(1)

    high_risk_hlas = [x.strip() for x in args.high_risk_hlas.split(",") if x.strip()] if args.high_risk_hlas else []
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

    if args.goal == "combo":
        print("Running ILP progression (combo goal)...")
        population_coverage_df, selected_indices = run_ilp_progression(
            filtered, hla_frequencies, hla_cols, args.ba_rank_hit, args.max_epitopes,
            args.prioritize_protein_diversity, hla_weights, high_risk_hlas,
            args.hard_constraint, args.mhc_class
        )

        best_epitopes = set()
        if "Best_Combination" in population_coverage_df.columns:
            for comb in population_coverage_df["Best_Combination"].fillna(""):
                for e in comb.split(","):
                    e = e.strip()
                    if e:
                        best_epitopes.add(e)

        epitope_subset = filtered[filtered["Epitope"].isin(best_epitopes)].copy()
        for col in ["Position", "Source", "Stage", "ConservationScore"]:
            if col not in epitope_subset.columns:
                epitope_subset[col] = np.nan

        if not epitope_subset.empty:
            loci_order = (["HLA-A", "HLA-B", "HLA-C"] if args.mhc_class.lower() == "i"
                          else ["DRB1", "HLA-DPA1", "HLA-DQA1"])
            for idx in epitope_subset.index:
                intra_cov, total_cov = calculate_population_coverage_exact(
                    [idx], epitope_subset.reset_index(drop=True),
                    hla_frequencies, hla_cols, args.ba_rank_hit, args.mhc_class
                )
                for locus in loci_order:
                    epitope_subset.loc[idx, f"{locus}_Cov_pct"] = round(intra_cov.get(locus, 0.0) * 100.0, 4)
                epitope_subset.loc[idx, "Total_Cov_pct"] = round(total_cov * 100.0, 2)

        locus_cols = (["HLA-A_Cov_pct", "HLA-B_Cov_pct", "HLA-C_Cov_pct"]
                      if args.mhc_class.lower() == "i"
                      else ["DRB1_Cov_pct", "HLA-DPA1_Cov_pct", "HLA-DQA1_Cov_pct"])
        for lc in locus_cols:
            if lc not in epitope_subset.columns:
                epitope_subset[lc] = 0.0

        hla_cols_in_df = [c for c in hla_cols if c not in locus_cols and c in epitope_subset.columns]
        ordered_cols = ["Epitope", "Position", "ConservationScore", "Source", "Stage"] + locus_cols + ["Total_Cov_pct"] + hla_cols_in_df
        for c in ordered_cols:
            if c not in epitope_subset.columns:
                epitope_subset[c] = np.nan
        epitope_subset = epitope_subset[ordered_cols]

        suffix = f"_{args.out_suffix}" if args.out_suffix else ""
        optimized_path = os.path.join(args.out_dir, f"optimized_epitopes{suffix}.tsv")
        popcov_path = os.path.join(args.out_dir, f"population_coverage{suffix}.tsv")
        epitope_subset.to_csv(optimized_path, index=False, sep="\t")
        population_coverage_df.to_csv(popcov_path, index=False, sep="\t")
        prepend_command_header(optimized_path, command_used)
        prepend_command_header(popcov_path, command_used)
        print(f"Saved optimized epitopes and population coverage. Selected count: {epitope_subset.shape[0]}")

    else:
        # region goal
        if args.window_size:
            region_df, best_epitopes_df = run_region_search_by_window(filtered, hla_frequencies, hla_cols, args)
        else:
            # minimal span now ranked globally across all proteins
            all_regions_df, _ = run_region_search_minimal_span(filtered, hla_frequencies, hla_cols, args)
            best_n = getattr(args, "best_region_num", 10)
            # sort by Total_Cov_pct descending to get best regions globally
            region_df = all_regions_df.sort_values(by="Total_Cov_pct", ascending=False).head(best_n).reset_index(drop=True)
            
            # collect all unique epitopes in the top regions
            all_top_epitopes = set()
            for eplist in region_df["Epitopes"].fillna(""):
                for e in eplist.split(","):
                    e = e.strip()
                    if e:
                        all_top_epitopes.add(e)
            best_epitopes_df = filtered[filtered["Epitope"].isin(all_top_epitopes)].copy()

        if args.mhc_class.lower() == "i":
            locus_cols = ["HLA-A_Cov_pct", "HLA-B_Cov_pct", "HLA-C_Cov_pct"]
        else:
            locus_cols = ["DRB1_Cov_pct", "HLA-DPA1_Cov_pct", "HLA-DQA1_Cov_pct"]

        for lc in locus_cols:
            if lc not in region_df.columns:
                region_df[lc] = 0.0

        if "Epitopes" in region_df.columns:
            region_df["Num_Epitopes"] = region_df["Epitopes"].apply(lambda x: len(x.split(",")) if x else 0)
        else:
            region_df["Epitopes"] = ""
            region_df["Num_Epitopes"] = 0

        if "Total_Cov_pct" not in region_df.columns:
            region_df["Total_Cov_pct"] = 0.0

        # reorder columns with Total_Cov_pct last
        region_cols = ["Protein", "Region_Start", "Region_End", "Region_Type", "Num_Epitopes", "Epitopes"] + locus_cols + ["Total_Cov_pct"]
        for c in region_cols:
            if c not in region_df.columns:
                region_df[c] = np.nan
        region_df = region_df[region_cols]

        suffix = f"_{args.out_suffix}" if args.out_suffix else ""
        region_path = os.path.join(args.out_dir, f"region_coverage{suffix}.tsv")
        optimized_path = os.path.join(args.out_dir, f"optimized_epitopes_region{suffix}.tsv")
        region_df.to_csv(region_path, index=False, sep="\t")
        best_epitopes_df.to_csv(optimized_path, index=False, sep="\t")
        prepend_command_header(region_path, command_used)
        prepend_command_header(optimized_path, command_used)
        print(f"Saved region search results. Best regions reported: {region_df.shape[0]}. Best region epitopes: {best_epitopes_df.shape[0]}")

if __name__ == "__main__":
    main()
