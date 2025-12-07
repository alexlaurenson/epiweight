#!/usr/bin/env python3
"""
Epitope ILP Optimization Tool

This tool performs optimization of epitope selection for vaccine design using
Integer Linear Programming (ILP) to maximize population coverage based on
HLA allele frequencies.
"""

import os
import sys
import argparse
import platform
import shutil
import re
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
import logging

import pandas as pd
import numpy as np

try:
    import pulp
except ImportError:
    pulp = None


# =============================================================================
# Configuration and Data Classes
# =============================================================================

@dataclass
class OptimizationConfig:
    """Configuration parameters for epitope optimization."""
    epitope_file: str
    hla_file: str
    mhc_class: str = "i"
    postfilter: bool = False
    conservation_score: float = 0.8
    max_epitopes: Optional[int] = None
    prioritize_protein_diversity: bool = False
    ba_rank_hit: float = 0.5
    out_dir: str = "output"
    out_suffix: str = ""
    high_risk_hlas: str = ""
    hla_weight: float = 2.0
    hard_constraint: bool = False
    goal: str = "combo"
    window_size: Optional[int] = None
    best_region_num: int = 10

    @property
    def high_risk_hla_list(self) -> List[str]:
        """Parse high risk HLAs into a list."""
        return [x.strip() for x in self.high_risk_hlas.split(",") if x.strip()]

    @property
    def loci_order(self) -> List[str]:
        """Get the ordered list of loci based on MHC class."""
        if self.mhc_class.lower() == "i":
            return ["HLA-A", "HLA-B", "HLA-C"]
        return ["DRB1", "HLA-DPA1", "HLA-DQA1"]

    @property
    def locus_coverage_columns(self) -> List[str]:
        """Get locus coverage column names."""
        return [f"{locus}_Cov_pct" for locus in self.loci_order]


@dataclass
class CoverageResult:
    """Results from population coverage calculation."""
    intra_locus_coverage: Dict[str, float]
    total_coverage: float

    @property
    def total_coverage_pct(self) -> float:
        """Total coverage as percentage."""
        return round(self.total_coverage * 100.0, 2)

    def get_locus_coverage_pct(self, locus: str) -> float:
        """Get coverage percentage for a specific locus."""
        return round(self.intra_locus_coverage.get(locus, 0.0) * 100.0, 4)


# =============================================================================
# Logger Setup
# =============================================================================

def setup_logger(name: str = "epitope_optimizer") -> logging.Logger:
    """Configure and return a logger instance."""
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(levelname)s: %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    
    return logger


logger = setup_logger()


# =============================================================================
# Solver Utilities
# =============================================================================

class SolverManager:
    """Manages ILP solver detection and configuration."""
    
    @staticmethod
    def get_solver_path() -> Optional[str]:
        """Detect and return path to CBC solver."""
        # Try system PATH first
        cbc_path = shutil.which("cbc")
        if cbc_path:
            return cbc_path
        
        # Platform-specific paths
        system = platform.system().lower()
        
        if system == "darwin":
            paths = ["/opt/homebrew/bin/cbc", "/usr/local/bin/cbc"]
        elif system == "windows":
            paths = [r"C:\Program Files\cbc\bin\cbc.exe", r"C:\cbc\bin\cbc.exe"]
        else:
            paths = []
        
        for path in paths:
            if os.path.exists(path):
                return path
        
        return None
    
    @staticmethod
    def get_solver() -> pulp.COIN_CMD:
        """Get configured CBC solver instance."""
        if pulp is None:
            raise ImportError("pulp package is required for ILP optimization")
        
        solver_path = SolverManager.get_solver_path()
        if solver_path is None:
            raise RuntimeError("Could not locate CBC solver")
        
        return pulp.COIN_CMD(msg=False, timeLimit=300, path=solver_path)


# =============================================================================
# Data Loading and Validation
# =============================================================================

class DataLoader:
    """Handles loading and validation of input data files."""
    
    PROTEIN_COLUMN_CANDIDATES = ["Source", "Protein", "Gene", "ProteinID"]
    
    @staticmethod
    def load_epitopes(epitope_file: str, postfilter: bool = False) -> Tuple[pd.DataFrame, List[str]]:
        """
        Load epitope data from file.
        
        Returns:
            Tuple of (epitopes_df, hla_column_names)
        """
        try:
            delimiter = '\t' if postfilter else (',' if epitope_file.endswith('.csv') else '\t')
            epitopes = pd.read_csv(epitope_file, delimiter=delimiter, low_memory=False)
            logger.info(f"Loaded {len(epitopes)} epitopes with {len(epitopes.columns)} columns")
        except Exception as e:
            logger.error(f"Error loading epitope file: {e}")
            sys.exit(1)
        
        # Validate required columns
        if postfilter:
            required_columns = ['Epitope', 'Position', 'ConservationScore', 'Source', 'Stage']
            if not all(col in epitopes.columns[:5] for col in required_columns):
                logger.error("Input missing required postfilter columns")
                sys.exit(1)
            hla_cols = list(epitopes.columns[5:])
        else:
            required_columns = ['Epitope', 'Position', 'ConservationScore']
            if not all(col in epitopes.columns for col in required_columns):
                logger.error("Input missing required columns")
                sys.exit(1)
            hla_cols = [col for col in epitopes.columns if col not in required_columns]
        
        if len(hla_cols) == 0:
            logger.error("No HLA allele columns found")
            sys.exit(1)
        
        return epitopes, hla_cols
    
    @staticmethod
    def load_hla_frequencies(hla_file: str) -> Dict[str, float]:
        """Load HLA allele frequencies from file."""
        try:
            hla_freq = pd.read_csv(hla_file, delimiter='\t', header=None)
            logger.info(f"Loaded HLA frequencies for {len(hla_freq)} alleles")
        except Exception as e:
            logger.error(f"Error loading HLA frequency file: {e}")
            sys.exit(1)
        
        freq_dict = {}
        for _, row in hla_freq.iterrows():
            allele = str(row.iloc[0]).strip()
            try:
                freq = float(row.iloc[1])
            except (ValueError, TypeError):
                logger.warning(f"Frequency for allele {allele} is not numeric. Setting to 0")
                freq = 0.0
            freq_dict[allele] = freq
        
        return freq_dict
    
    @staticmethod
    def filter_frequencies_by_columns(
        hla_frequencies: Dict[str, float],
        hla_cols: List[str]
    ) -> Dict[str, float]:
        """Filter HLA frequencies to only include those present in epitope matrix."""
        filtered = {
            allele: float(freq)
            for allele, freq in hla_frequencies.items()
            if allele in hla_cols
        }
        
        if len(filtered) == 0:
            logger.error("None of the HLA frequencies match epitope matrix columns")
            sys.exit(1)
        
        return filtered
    
    @staticmethod
    def find_protein_column(df: pd.DataFrame) -> Optional[str]:
        """Find the protein identifier column in a dataframe."""
        for col in DataLoader.PROTEIN_COLUMN_CANDIDATES:
            if col in df.columns:
                return col
        return None


# =============================================================================
# HLA and Coverage Calculations
# =============================================================================

class HLAAnalyzer:
    """Utilities for HLA allele analysis and coverage calculations."""
    
    @staticmethod
    def phenotypic_frequency(allele_freq: float) -> float:
        """
        Calculate phenotypic frequency from allele frequency.
        P(allele present) = 1 - (1 - f)^2 = 2f - f^2
        """
        return 2.0 * allele_freq - allele_freq * allele_freq
    
    @staticmethod
    def allele_to_locus(allele: str) -> Optional[str]:
        """Map an allele name to its locus."""
        s = str(allele).strip().upper()
        
        if s.startswith("HLA-A"):
            return "HLA-A"
        if s.startswith("HLA-B"):
            return "HLA-B"
        if s.startswith("HLA-C"):
            return "HLA-C"
        if "DRB1" in s:
            return "DRB1"
        if s.startswith("HLA-DPA1"):
            return "HLA-DPA1"
        if s.startswith("HLA-DQA1"):
            return "HLA-DQA1"
        
        # Generic pattern match
        match = re.search(r'D[A-Z]+[0-9]', s)
        return match.group(0) if match else None
    
    @staticmethod
    def calculate_allele_weights(
        hla_frequencies: Dict[str, float],
        high_risk_hlas: List[str],
        hla_weight: float = 2.0
    ) -> Dict[str, float]:
        """
        Calculate linear weights for alleles based on frequencies and risk.
        
        Uses -log(1 - phenofreq) as linear approximation for multiplicative merging.
        Applies multiplicative weight for high-risk HLAs.
        """
        weights = {}
        for allele, freq in hla_frequencies.items():
            pheno_freq = HLAAnalyzer.phenotypic_frequency(freq)
            # Small epsilon to avoid log(0)
            weight = -np.log(1.0 - pheno_freq + 1e-12)
            
            # Apply high-risk weight if applicable
            for hr in high_risk_hlas:
                if hr and hr in allele:
                    weight *= hla_weight
                    break
            
            weights[allele] = weight
        
        return weights


class CoverageCalculator:
    """Calculates population coverage using phenotypic frequency methodology."""
    
    @staticmethod
    def calculate_exact_coverage(
        selected_indices: List[int],
        epitopes_df: pd.DataFrame,
        hla_frequencies: Dict[str, float],
        hla_cols: List[str],
        rank_cutoff: float,
        mhc_class: str
    ) -> CoverageResult:
        """
        Calculate exact phenotypic frequency-style population coverage.
        
        Performs intra-locus and inter-locus coverage calculation for selected epitopes.
        """
        # Collect all alleles covered by selected epitopes
        covered_alleles = []
        for idx in selected_indices:
            if idx is None:
                continue
            row = epitopes_df.iloc[idx]
            for hla in hla_cols:
                val = row.get(hla, np.nan)
                if not pd.isna(val):
                    try:
                        if float(val) <= rank_cutoff:
                            covered_alleles.append(hla)
                    except (ValueError, TypeError):
                        pass
        
        # Define loci based on MHC class
        if mhc_class.lower() == "i":
            loci_list = ["HLA-A", "HLA-B", "HLA-C"]
        else:
            loci_list = ["DRB1", "HLA-DPA1", "HLA-DQA1"]
        
        # Group alleles by locus
        loci_map = {locus: [] for locus in loci_list}
        for allele in covered_alleles:
            locus = HLAAnalyzer.allele_to_locus(allele)
            if locus and locus in loci_map:
                loci_map[locus].append(allele)
        
        # Calculate intra-locus coverage
        intra_coverage = {}
        for locus, alleles_in_locus in loci_map.items():
            pheno_freqs = []
            for allele in alleles_in_locus:
                if allele in hla_frequencies:
                    freq = float(hla_frequencies[allele])
                    pheno_freqs.append(HLAAnalyzer.phenotypic_frequency(freq))
            
            if pheno_freqs:
                # Coverage = 1 - product(1 - phenofreq)
                intra_coverage[locus] = 1.0 - np.prod([1.0 - pf for pf in pheno_freqs])
            else:
                intra_coverage[locus] = 0.0
        
        # Calculate inter-locus coverage
        if intra_coverage:
            total_coverage = 1.0 - np.prod([1.0 - c for c in intra_coverage.values()])
        else:
            total_coverage = 0.0
        
        return CoverageResult(
            intra_locus_coverage=intra_coverage,
            total_coverage=total_coverage
        )


# =============================================================================
# ILP Optimization
# =============================================================================

class ILPOptimizer:
    """Performs ILP-based epitope selection optimization."""
    
    def __init__(self, config: OptimizationConfig):
        self.config = config
        self.solver = SolverManager.get_solver()
    
    def build_coverage_matrix(
        self,
        epitopes_df: pd.DataFrame,
        alleles: List[str],
        rank_cutoff: float
    ) -> np.ndarray:
        """
        Build binary coverage matrix a[i,j].
        
        a[i,j] = 1 if allele i is covered by epitope j, 0 otherwise.
        """
        num_alleles = len(alleles)
        num_epitopes = len(epitopes_df)
        coverage_matrix = np.zeros((num_alleles, num_epitopes), dtype=int)
        
        for j, (_, epitope_row) in enumerate(epitopes_df.iterrows()):
            for i, allele in enumerate(alleles):
                if allele in epitopes_df.columns:
                    val = epitope_row.get(allele, np.nan)
                    if not pd.isna(val):
                        try:
                            if float(val) <= rank_cutoff:
                                coverage_matrix[i, j] = 1
                        except (ValueError, TypeError):
                            pass
        
        return coverage_matrix
    
    def get_protein_groups(self, epitopes_df: pd.DataFrame) -> Dict[str, List[int]]:
        """Group epitope indices by protein for diversity constraints."""
        if not self.config.prioritize_protein_diversity:
            return {}
        
        protein_col = DataLoader.find_protein_column(epitopes_df)
        if protein_col is None:
            logger.error(
                "Protein diversity requested but no protein identifier column found "
                "(checked Source/Protein/Gene/ProteinID)"
            )
            sys.exit(1)
        
        protein_to_indices = {}
        for idx, protein in enumerate(epitopes_df[protein_col]):
            protein_key = str(protein)
            protein_to_indices.setdefault(protein_key, []).append(idx)
        
        return protein_to_indices
    
    def solve_for_k_epitopes(
        self,
        k: int,
        epitopes_df: pd.DataFrame,
        hla_frequencies: Dict[str, float],
        allele_weights: Dict[str, float],
        coverage_matrix: np.ndarray,
        alleles: List[str],
        protein_groups: Dict[str, List[int]]
    ) -> List[int]:
        """
        Solve ILP to select up to k epitopes maximizing weighted coverage.
        
        Returns list of selected epitope indices.
        """
        num_epitopes = len(epitopes_df)
        num_alleles = len(alleles)
        
        # Create problem
        prob = pulp.LpProblem(f"maxcov_k{k}", pulp.LpMaximize)
        
        # Decision variables
        x_vars = [pulp.LpVariable(f"x_{j}", cat="Binary") for j in range(num_epitopes)]
        y_vars = [pulp.LpVariable(f"y_{i}", cat="Binary") for i in range(num_alleles)]
        
        # Objective: maximize sum of weighted allele coverage
        prob += pulp.lpSum([
            allele_weights[alleles[i]] * y_vars[i]
            for i in range(num_alleles)
        ])
        
        # Constraint: link allele coverage to selected epitopes
        for i in range(num_alleles):
            prob += y_vars[i] <= pulp.lpSum([
                coverage_matrix[i, j] * x_vars[j]
                for j in range(num_epitopes)
            ])
        
        # Constraint: limit number of epitopes
        prob += pulp.lpSum(x_vars) <= k
        
        # Constraint: at most one epitope per protein (if diversity enabled)
        if protein_groups:
            for protein, indices in protein_groups.items():
                if len(indices) > 0:
                    prob += pulp.lpSum([x_vars[j] for j in indices]) <= 1
        
        # Hard constraint: require coverage for high-risk HLAs
        if self.config.hard_constraint and self.config.high_risk_hla_list:
            for hr_pattern in self.config.high_risk_hla_list:
                matched_indices = [
                    i for i, allele_name in enumerate(alleles)
                    if hr_pattern in allele_name
                ]
                for i in matched_indices:
                    prob += y_vars[i] == 1
        
        # Solve
        prob.solve(self.solver)
        
        # Extract selected epitopes
        selected = [
            j for j in range(num_epitopes)
            if pulp.value(x_vars[j]) is not None and pulp.value(x_vars[j]) > 0.5
        ]
        
        return selected
    
    def run_progression(
        self,
        epitopes_df: pd.DataFrame,
        hla_frequencies: Dict[str, float],
        hla_cols: List[str]
    ) -> Tuple[pd.DataFrame, List[int]]:
        """
        Run ILP progression from k=1 to max_k epitopes.
        
        Returns:
            Tuple of (coverage_progression_df, best_combination_indices)
        """
        epitopes_df = epitopes_df.reset_index(drop=True)
        
        if "Epitope" not in epitopes_df.columns:
            logger.error("Epitopes dataframe must contain 'Epitope' column")
            sys.exit(1)
        
        epitope_names = list(epitopes_df["Epitope"])
        num_epitopes = len(epitope_names)
        
        # Filter to alleles with known frequencies
        alleles = list(hla_frequencies.keys())
        
        # Build coverage matrix
        coverage_matrix = self.build_coverage_matrix(
            epitopes_df, alleles, self.config.ba_rank_hit
        )
        
        # Calculate allele weights
        allele_weights = HLAAnalyzer.calculate_allele_weights(
            hla_frequencies,
            self.config.high_risk_hla_list,
            self.config.hla_weight
        )
        
        # Get protein groups for diversity constraints
        protein_groups = self.get_protein_groups(epitopes_df)
        
        # Run progression
        max_k = self.config.max_epitopes if self.config.max_epitopes else num_epitopes
        progression_rows = []
        best_combinations = {}
        prev_total_cov = 0.0
        
        for k in range(1, max_k + 1):
            selected_indices = self.solve_for_k_epitopes(
                k, epitopes_df, hla_frequencies, allele_weights,
                coverage_matrix, alleles, protein_groups
            )
            
            if not selected_indices:
                continue
            
            # Calculate exact coverage
            coverage = CoverageCalculator.calculate_exact_coverage(
                selected_indices, epitopes_df, hla_frequencies,
                hla_cols, self.config.ba_rank_hit, self.config.mhc_class
            )
            
            best_combinations[k] = selected_indices
            
            # Build row for progression dataframe
            row = {
                "Num_Epitopes": k,
                "Best_Combination": ",".join([epitope_names[i] for i in selected_indices])
            }
            
            for locus in self.config.loci_order:
                row[f"{locus}_Cov_pct"] = coverage.get_locus_coverage_pct(locus)
            
            row["Total_Cov_pct"] = coverage.total_coverage_pct
            progression_rows.append(row)
            
            # Stop if coverage not improving or at 100%
            if coverage.total_coverage_pct <= prev_total_cov or coverage.total_coverage_pct >= 100.0:
                break
            prev_total_cov = coverage.total_coverage_pct
        
        # Build progression dataframe
        if progression_rows:
            df = pd.DataFrame(progression_rows)
            # Reorder columns
            base_cols = ["Num_Epitopes", "Best_Combination"]
            locus_cols = self.config.locus_coverage_columns
            final_cols = base_cols + locus_cols + ["Total_Cov_pct"]
            df = df[final_cols]
        else:
            df = pd.DataFrame(columns=[
                "Num_Epitopes", "Best_Combination"
            ] + self.config.locus_coverage_columns + ["Total_Cov_pct"])
        
        # Get best combination from last iteration
        final_best = best_combinations[max(best_combinations.keys())] if best_combinations else []
        
        return df, final_best


# =============================================================================
# Region Search
# =============================================================================

class RegionSearcher:
    """Performs region-based epitope search."""
    
    def __init__(self, config: OptimizationConfig):
        self.config = config
        self.optimizer = ILPOptimizer(config)
    
    def prepare_positions(self, df: pd.DataFrame) -> pd.DataFrame:
        """Convert Position column to numeric integers."""
        if "Position" not in df.columns:
            return df
        
        df = df.copy()
        try:
            df["Position"] = df["Position"].astype(int)
        except (ValueError, TypeError):
            # Try extracting digits
            df["Position"] = df["Position"].apply(
                lambda x: int(re.sub(r'[^0-9]', '', str(x)))
                if pd.notna(x) and re.search(r'\d', str(x))
                else np.nan
            )
        
        return df.dropna(subset=["Position"])
    
    def make_region_row(
        self,
        protein: str,
        start: int,
        end: int,
        region_type: str,
        epitope_names: List[str],
        coverage: CoverageResult
    ) -> Dict:
        """Create a row dictionary for region results."""
        row = {
            "Protein": protein,
            "Region_Start": int(start) if start is not None else None,
            "Region_End": int(end) if end is not None else None,
            "Region_Type": region_type,
            "Num_Epitopes": len(epitope_names),
            "Epitopes": ",".join(epitope_names),
            "Total_Cov_pct": coverage.total_coverage_pct
        }
        
        # Add per-locus coverage
        for locus in self.config.loci_order:
            row[f"{locus}_Cov_pct"] = coverage.get_locus_coverage_pct(locus)
        
        return row
    
    def search_sliding_window(
        self,
        filtered_df: pd.DataFrame,
        hla_frequencies: Dict[str, float],
        hla_cols: List[str]
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Search for best regions using sliding window approach.
        
        Returns:
            Tuple of (region_results_df, best_epitopes_df)
        """
        if self.config.window_size is None:
            raise ValueError("window_size is required for sliding window mode")
        
        protein_col = DataLoader.find_protein_column(filtered_df)
        if protein_col is None:
            logger.error("No protein identifier column found for region analysis")
            sys.exit(1)
        
        region_rows = []
        best_overall = {"total_cov": -1.0}
        proteins = filtered_df[protein_col].unique()
        
        for protein in proteins:
            subset = filtered_df[filtered_df[protein_col] == protein].copy()
            if subset.empty:
                continue
            
            subset = self.prepare_positions(subset)
            if subset.empty:
                continue
            
            min_pos = int(subset["Position"].min())
            max_pos = int(subset["Position"].max())
            
            # Calculate window positions
            step_end = max_pos - self.config.window_size + 1
            if step_end < min_pos:
                window_starts = [min_pos]
            else:
                window_starts = list(range(min_pos, step_end + 1))
            
            for start in window_starts:
                end = start + self.config.window_size - 1
                in_window = subset[
                    (subset["Position"] >= start) & (subset["Position"] <= end)
                ].reset_index(drop=True)
                
                if in_window.empty:
                    continue
                
                # Run ILP on window
                try:
                    _, chosen_indices = self.optimizer.run_progression(
                        in_window, hla_frequencies, hla_cols
                    )
                except Exception as e:
                    logger.warning(
                        f"ILP failed for protein {protein} window {start}-{end}: {e}"
                    )
                    continue
                
                if not chosen_indices:
                    continue
                
                # Calculate coverage
                coverage = CoverageCalculator.calculate_exact_coverage(
                    chosen_indices, in_window, hla_frequencies,
                    hla_cols, self.config.ba_rank_hit, self.config.mhc_class
                )
                
                epitope_names = [in_window.iloc[i]["Epitope"] for i in chosen_indices]
                row = self.make_region_row(
                    protein, start, end, "sliding", epitope_names, coverage
                )
                region_rows.append(row)
                
                # Track best overall
                if coverage.total_coverage > best_overall["total_cov"]:
                    best_overall = {
                        "total_cov": coverage.total_coverage,
                        "protein": protein,
                        "start": start,
                        "end": end,
                        "chosen_indices": chosen_indices,
                        "in_window_df": in_window.copy(),
                        "coverage": coverage
                    }
        
        # Build results dataframes with adaptive selection
        all_regions_df = self._build_region_dataframe(region_rows)
        
        if not all_regions_df.empty:
            region_df = self._select_best_regions(all_regions_df, filtered_df)
            
            # Collect epitopes from selected regions
            all_top_epitopes = set()
            for eplist in region_df["Epitopes"].fillna(""):
                for e in eplist.split(","):
                    e = e.strip()
                    if e:
                        all_top_epitopes.add(e)
            
            best_epitopes_df = filtered_df[
                filtered_df["Epitope"].isin(all_top_epitopes)
            ].copy()
        else:
            region_df = all_regions_df
            best_epitopes_df = pd.DataFrame(columns=filtered_df.columns)
        
        return region_df, best_epitopes_df
    
    def search_minimal_span(
        self,
        filtered_df: pd.DataFrame,
        hla_frequencies: Dict[str, float],
        hla_cols: List[str]
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Search for best regions using minimal contiguous span.
        
        Returns:
            Tuple of (region_results_df, best_epitopes_df)
        """
        protein_col = DataLoader.find_protein_column(filtered_df)
        if protein_col is None:
            logger.error("No protein identifier column found for region analysis")
            sys.exit(1)
        
        region_rows = []
        proteins = filtered_df[protein_col].unique()
        
        for protein in proteins:
            subset = filtered_df[filtered_df[protein_col] == protein].copy()
            if subset.empty:
                continue
            
            subset = self.prepare_positions(subset)
            if subset.empty:
                continue
            
            # Run ILP on entire protein
            try:
                _, chosen_indices = self.optimizer.run_progression(
                    subset, hla_frequencies, hla_cols
                )
            except Exception as e:
                logger.warning(f"ILP failed for protein {protein}: {e}")
                continue
            
            if not chosen_indices:
                continue
            
            # Calculate minimal contiguous span
            chosen_rows = subset.reset_index(drop=True).iloc[chosen_indices]
            positions = []
            epitope_names = []
            
            for _, row in chosen_rows.iterrows():
                start = int(row["Position"])
                seq = str(row.get("Epitope", ""))
                end = start + max(len(seq) - 1, 0)
                positions.append((start, end))
                epitope_names.append(row["Epitope"])
            
            region_start = min(s for s, e in positions)
            region_end = max(e for s, e in positions)
            
            # Calculate coverage
            coverage = CoverageCalculator.calculate_exact_coverage(
                chosen_indices, subset.reset_index(drop=True),
                hla_frequencies, hla_cols, self.config.ba_rank_hit,
                self.config.mhc_class
            )
            
            row = self.make_region_row(
                protein, region_start, region_end, "minimal_span",
                epitope_names, coverage
            )
            region_rows.append(row)
        
        # Build results
        all_regions_df = self._build_region_dataframe(region_rows)
        
        # Select regions based on adaptive strategy
        if not all_regions_df.empty:
            region_df = self._select_best_regions(all_regions_df, filtered_df)
            
            # Collect all epitopes from selected regions
            all_top_epitopes = set()
            for eplist in region_df["Epitopes"].fillna(""):
                for e in eplist.split(","):
                    e = e.strip()
                    if e:
                        all_top_epitopes.add(e)
            
            best_epitopes_df = filtered_df[
                filtered_df["Epitope"].isin(all_top_epitopes)
            ].copy()
        else:
            region_df = all_regions_df
            best_epitopes_df = pd.DataFrame(columns=filtered_df.columns)
        
        return region_df, best_epitopes_df
    
    def _build_region_dataframe(self, region_rows: List[Dict]) -> pd.DataFrame:
        """Build standardized region results dataframe."""
        if not region_rows:
            return pd.DataFrame(columns=[
                "Protein", "Region_Start", "Region_End", "Region_Type",
                "Num_Epitopes", "Epitopes"
            ] + self.config.locus_coverage_columns + ["Total_Cov_pct"])
        
        df = pd.DataFrame(region_rows)
        
        # Ensure all expected columns exist
        for col in self.config.locus_coverage_columns:
            if col not in df.columns:
                df[col] = 0.0
        
        # Standardize column order
        cols = [
            "Protein", "Region_Start", "Region_End", "Region_Type",
            "Num_Epitopes", "Epitopes"
        ] + self.config.locus_coverage_columns + ["Total_Cov_pct"]
        
        for col in cols:
            if col not in df.columns:
                df[col] = np.nan
        
        return df[cols]
    
    def _select_best_regions(
        self,
        all_regions_df: pd.DataFrame,
        filtered_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Adaptively select best regions based on coverage.
        
        Strategy:
        1. If best_region_num was explicitly provided by user, use that value
        2. Otherwise, select all regions with 100% coverage
        3. If no regions achieve 100% coverage, default to top 10
        """
        if all_regions_df.empty:
            return all_regions_df
        
        # Sort by coverage descending
        sorted_df = all_regions_df.sort_values(
            by="Total_Cov_pct", ascending=False
        ).reset_index(drop=True)
        
        # Check if best_region_num was explicitly set by user
        # We detect this by checking if it differs from the default value
        user_specified_num = (
            self.config.best_region_num != 10 or 
            '--best_region_num' in sys.argv
        )
        
        if user_specified_num:
            # User explicitly specified a number, use it
            logger.info(
                f"Using user-specified region count: {self.config.best_region_num}"
            )
            return sorted_df.head(self.config.best_region_num)
        
        # Adaptive selection: prefer 100% coverage regions
        perfect_coverage = sorted_df[sorted_df["Total_Cov_pct"] >= 100.0]
        
        if not perfect_coverage.empty:
            logger.info(
                f"Found {len(perfect_coverage)} region(s) with 100% coverage"
            )
            return perfect_coverage
        
        # No perfect coverage, use default top 10
        logger.info(
            f"No regions achieved 100% coverage. "
            f"Returning top {self.config.best_region_num} regions "
            f"(max coverage: {sorted_df['Total_Cov_pct'].iloc[0]:.2f}%)"
        )
        return sorted_df.head(self.config.best_region_num)
    
    def _build_best_epitopes_dataframe(
        self,
        best_overall: Dict,
        hla_frequencies: Dict[str, float],
        hla_cols: List[str],
        filtered_df: pd.DataFrame
    ) -> pd.DataFrame:
        """Build epitopes dataframe for best region."""
        if best_overall.get("total_cov", -1.0) < 0:
            return pd.DataFrame(columns=filtered_df.columns)
        
        chosen = best_overall["chosen_indices"]
        best_df = best_overall["in_window_df"].iloc[chosen].copy().reset_index(drop=True)
        
        # Add individual coverage for each epitope
        for idx in best_df.index:
            coverage = CoverageCalculator.calculate_exact_coverage(
                [idx], best_df, hla_frequencies, hla_cols,
                self.config.ba_rank_hit, self.config.mhc_class
            )
            
            for locus in self.config.loci_order:
                best_df.loc[idx, f"{locus}_Cov_pct"] = coverage.get_locus_coverage_pct(locus)
            best_df.loc[idx, "Total_Cov_pct"] = coverage.total_coverage_pct
        
        return best_df


# =============================================================================
# Output Management
# =============================================================================

class OutputManager:
    """Handles file output and formatting."""
    
    def __init__(self, config: OptimizationConfig):
        self.config = config
        os.makedirs(config.out_dir, exist_ok=True)
    
    def get_output_path(self, filename: str) -> str:
        """Get full output path with suffix."""
        suffix = f"_{self.config.out_suffix}" if self.config.out_suffix else ""
        name, ext = os.path.splitext(filename)
        return os.path.join(self.config.out_dir, f"{name}{suffix}{ext}")
    
    @staticmethod
    def prepend_command_header(filepath: str, command: str) -> None:
        """Prepend command used to the beginning of file."""
        try:
            with open(filepath, "r") as f:
                original = f.read()
        except FileNotFoundError:
            original = ""
        
        header = f"# Command used to run this script:\n# {command}\n\n"
        with open(filepath, "w") as f:
            f.write(header + original)
    
    def save_combo_results(
        self,
        epitope_subset: pd.DataFrame,
        coverage_df: pd.DataFrame,
        command: str
    ) -> None:
        """Save results from combo optimization."""
        optimized_path = self.get_output_path("optimized_epitopes.tsv")
        popcov_path = self.get_output_path("population_coverage.tsv")
        
        epitope_subset.to_csv(optimized_path, index=False, sep="\t")
        coverage_df.to_csv(popcov_path, index=False, sep="\t")
        
        self.prepend_command_header(optimized_path, command)
        self.prepend_command_header(popcov_path, command)
        
        logger.info(
            f"Saved optimized epitopes and population coverage. "
            f"Selected count: {epitope_subset.shape[0]}"
        )
    
    def save_region_results(
        self,
        region_df: pd.DataFrame,
        epitopes_df: pd.DataFrame,
        command: str
    ) -> None:
        """Save results from region optimization."""
        region_path = self.get_output_path("region_coverage.tsv")
        optimized_path = self.get_output_path("optimized_epitopes_region.tsv")
        
        region_df.to_csv(region_path, index=False, sep="\t")
        epitopes_df.to_csv(optimized_path, index=False, sep="\t")
        
        self.prepend_command_header(region_path, command)
        self.prepend_command_header(optimized_path, command)
        
        logger.info(
            f"Saved region search results. "
            f"Best regions reported: {region_df.shape[0]}. "
            f"Best region epitopes: {epitopes_df.shape[0]}"
        )


# =============================================================================
# Combo Optimization Pipeline
# =============================================================================

class ComboOptimizationPipeline:
    """Pipeline for combination-based epitope optimization."""
    
    def __init__(self, config: OptimizationConfig):
        self.config = config
        self.optimizer = ILPOptimizer(config)
    
    def run(
        self,
        filtered_df: pd.DataFrame,
        hla_frequencies: Dict[str, float],
        hla_cols: List[str]
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Run combo optimization pipeline.
        
        Returns:
            Tuple of (optimized_epitopes_df, population_coverage_df)
        """
        logger.info("Running ILP progression (combo goal)...")
        
        # Run ILP progression
        coverage_df, selected_indices = self.optimizer.run_progression(
            filtered_df, hla_frequencies, hla_cols
        )
        
        # Extract best epitopes
        best_epitopes = self._extract_best_epitopes(coverage_df)
        epitope_subset = filtered_df[
            filtered_df["Epitope"].isin(best_epitopes)
        ].copy()
        
        # Add required columns if missing
        for col in ["Position", "Source", "Stage", "ConservationScore"]:
            if col not in epitope_subset.columns:
                epitope_subset[col] = np.nan
        
        # Calculate individual coverage for each epitope
        if not epitope_subset.empty:
            epitope_subset = self._add_individual_coverage(
                epitope_subset, hla_frequencies, hla_cols
            )
        
        # Format columns
        epitope_subset = self._format_output_columns(epitope_subset, hla_cols)
        
        return epitope_subset, coverage_df
    
    def _extract_best_epitopes(self, coverage_df: pd.DataFrame) -> Set[str]:
        """Extract unique epitopes from Best_Combination column."""
        best_epitopes = set()
        if "Best_Combination" in coverage_df.columns:
            for comb in coverage_df["Best_Combination"].fillna(""):
                for epitope in comb.split(","):
                    epitope = epitope.strip()
                    if epitope:
                        best_epitopes.add(epitope)
        return best_epitopes
    
    def _add_individual_coverage(
        self,
        epitope_df: pd.DataFrame,
        hla_frequencies: Dict[str, float],
        hla_cols: List[str]
    ) -> pd.DataFrame:
        """Add individual coverage columns for each epitope."""
        for idx in epitope_df.index:
            coverage = CoverageCalculator.calculate_exact_coverage(
                [idx], epitope_df.reset_index(drop=True),
                hla_frequencies, hla_cols, self.config.ba_rank_hit,
                self.config.mhc_class
            )
            
            for locus in self.config.loci_order:
                epitope_df.loc[idx, f"{locus}_Cov_pct"] = coverage.get_locus_coverage_pct(locus)
            epitope_df.loc[idx, "Total_Cov_pct"] = coverage.total_coverage_pct
        
        return epitope_df
    
    def _format_output_columns(
        self,
        epitope_df: pd.DataFrame,
        hla_cols: List[str]
    ) -> pd.DataFrame:
        """Format and order output columns."""
        # Ensure locus coverage columns exist
        for col in self.config.locus_coverage_columns:
            if col not in epitope_df.columns:
                epitope_df[col] = 0.0
        
        # Define column order
        hla_cols_in_df = [
            c for c in hla_cols
            if c not in self.config.locus_coverage_columns
            and c in epitope_df.columns
        ]
        
        ordered_cols = [
            "Epitope", "Position", "ConservationScore", "Source", "Stage"
        ] + self.config.locus_coverage_columns + ["Total_Cov_pct"] + hla_cols_in_df
        
        # Add missing columns
        for col in ordered_cols:
            if col not in epitope_df.columns:
                epitope_df[col] = np.nan
        
        return epitope_df[ordered_cols]


# =============================================================================
# Main Application
# =============================================================================

class EpitopeOptimizer:
    """Main application controller."""
    
    def __init__(self, config: OptimizationConfig):
        self.config = config
        self.output_manager = OutputManager(config)
    
    def run(self) -> None:
        """Execute the complete optimization pipeline."""
        command_used = " ".join(sys.argv)
        
        # Load data
        logger.info("Loading data files...")
        epitopes, hla_cols = DataLoader.load_epitopes(
            self.config.epitope_file,
            self.config.postfilter
        )
        hla_frequencies = DataLoader.load_hla_frequencies(self.config.hla_file)
        hla_frequencies = DataLoader.filter_frequencies_by_columns(
            hla_frequencies, hla_cols
        )
        
        # Apply conservation filter
        if not self.config.postfilter:
            logger.error("Conservation processing not implemented outside postfilter")
            sys.exit(1)
        
        filtered = epitopes[
            epitopes["ConservationScore"] >= self.config.conservation_score
        ].reset_index(drop=True)
        
        if filtered.shape[0] == 0:
            logger.error("No epitopes pass conservation filter")
            sys.exit(1)
        
        logger.info(f"Filtered to {filtered.shape[0]} epitopes")
        
        # Run optimization based on goal
        if self.config.goal == "combo":
            pipeline = ComboOptimizationPipeline(self.config)
            epitope_subset, coverage_df = pipeline.run(
                filtered, hla_frequencies, hla_cols
            )
            self.output_manager.save_combo_results(
                epitope_subset, coverage_df, command_used
            )
        
        else:  # region goal
            logger.info("Running region search...")
            searcher = RegionSearcher(self.config)
            
            if self.config.window_size:
                region_df, best_epitopes_df = searcher.search_sliding_window(
                    filtered, hla_frequencies, hla_cols
                )
            else:
                region_df, best_epitopes_df = searcher.search_minimal_span(
                    filtered, hla_frequencies, hla_cols
                )
            
            self.output_manager.save_region_results(
                region_df, best_epitopes_df, command_used
            )


# =============================================================================
# CLI Entry Point
# =============================================================================

def parse_arguments() -> OptimizationConfig:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Epitope ILP optimization with exact phenotypic frequency coverage reporting."
    )
    
    # Required arguments
    parser.add_argument(
        "epitope_file",
        help="Path to the epitope data file"
    )
    parser.add_argument(
        "hla_file",
        help="Path to the HLA frequencies file"
    )
    
    # Optional arguments
    parser.add_argument(
        "--mhc_class",
        choices=["i", "ii"],
        default="i",
        help="MHC class for coverage calculation (i or ii)"
    )
    parser.add_argument(
        "--postfilter",
        action="store_true",
        help="Use postfilter mode for input file parsing"
    )
    parser.add_argument(
        "--conservation_score",
        type=float,
        default=0.8,
        help="Minimum conservation score threshold (default: 0.8)"
    )
    parser.add_argument(
        "--max_epitopes",
        type=int,
        default=None,
        help="Maximum number of epitopes to select"
    )
    parser.add_argument(
        "--prioritize_protein_diversity",
        action="store_true",
        help="Enforce at most one epitope per protein"
    )
    parser.add_argument(
        "--ba_rank_hit",
        type=float,
        default=0.5,
        help="Binding affinity rank threshold (default: 0.5)"
    )
    parser.add_argument(
        "--out_dir",
        default="output",
        help="Output directory (default: output)"
    )
    parser.add_argument(
        "--out_suffix",
        default="",
        help="Suffix to append to output filenames"
    )
    parser.add_argument(
        "--high_risk_hlas",
        default="",
        help="Comma-separated list of high-risk HLA patterns"
    )
    parser.add_argument(
        "--hla_weight",
        type=float,
        default=2.0,
        help="Weight multiplier for high-risk HLAs (default: 2.0)"
    )
    parser.add_argument(
        "--hard_constraint",
        action="store_true",
        help="Enforce hard constraint requiring coverage of high-risk HLAs"
    )
    parser.add_argument(
        "--goal",
        choices=["combo", "region"],
        default="combo",
        help="Optimization goal: 'combo' = epitope combinations; 'region' = find regions"
    )
    parser.add_argument(
        "--window_size",
        type=int,
        default=None,
        help="Sliding window size in amino acids (for --goal region)"
    )
    parser.add_argument(
        "--best_region_num",
        type=int,
        default=10,
        help="Number of top regions to include in output (default: 10)"
    )
    
    args = parser.parse_args()
    
    # Convert to config object
    return OptimizationConfig(
        epitope_file=args.epitope_file,
        hla_file=args.hla_file,
        mhc_class=args.mhc_class,
        postfilter=args.postfilter,
        conservation_score=args.conservation_score,
        max_epitopes=args.max_epitopes,
        prioritize_protein_diversity=args.prioritize_protein_diversity,
        ba_rank_hit=args.ba_rank_hit,
        out_dir=args.out_dir,
        out_suffix=args.out_suffix,
        high_risk_hlas=args.high_risk_hlas,
        hla_weight=args.hla_weight,
        hard_constraint=args.hard_constraint,
        goal=args.goal,
        window_size=args.window_size,
        best_region_num=args.best_region_num
    )


def main():
    """Main entry point."""
    try:
        config = parse_arguments()
        optimizer = EpitopeOptimizer(config)
        optimizer.run()
    except KeyboardInterrupt:
        logger.info("Interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
