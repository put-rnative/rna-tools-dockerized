#! /usr/bin/env python
import argparse
import os.path
from typing import List

SCORING_METHODS = [
    "3dRNAscore",
    "DFIRE",
    "RASP",
    "RNA-BRiQ",
    "RNA3DCNN",
    "cgRNASP",
    "cgRNASP-C",
    "cgRNASP-CN",
    "cgRNASP-PC",
    "lociPARSE",
    "rsRNASP",
]


def score_3drnascore(pdb_path: str) -> float:
    """Score RNA structure using 3dRNAscore method"""
    return 0.0

def score_dfire(pdb_path: str) -> float:
    """Score RNA structure using DFIRE method"""
    return 0.0

def score_rasp(pdb_path: str) -> float:
    """Score RNA structure using RASP method"""
    return 0.0

def score_rna_briq(pdb_path: str) -> float:
    """Score RNA structure using RNA-BRiQ method"""
    return 0.0

def score_rna3dcnn(pdb_path: str) -> float:
    """Score RNA structure using RNA3DCNN method"""
    return 0.0

def score_cgrnasp(pdb_path: str) -> float:
    """Score RNA structure using cgRNASP method"""
    return 0.0

def score_cgrnasp_c(pdb_path: str) -> float:
    """Score RNA structure using cgRNASP-C method"""
    return 0.0

def score_cgrnasp_cn(pdb_path: str) -> float:
    """Score RNA structure using cgRNASP-CN method"""
    return 0.0

def score_cgrnasp_pc(pdb_path: str) -> float:
    """Score RNA structure using cgRNASP-PC method"""
    return 0.0

def score_lociparse(pdb_path: str) -> float:
    """Score RNA structure using lociPARSE method"""
    return 0.0

def score_rsrnasp(pdb_path: str) -> float:
    """Score RNA structure using rsRNASP method"""
    return 0.0

# Mapping of method names to scoring functions
SCORING_FUNCTIONS = {
    "3dRNAscore": score_3drnascore,
    "DFIRE": score_dfire,
    "RASP": score_rasp,
    "RNA-BRiQ": score_rna_briq,
    "RNA3DCNN": score_rna3dcnn,
    "cgRNASP": score_cgrnasp,
    "cgRNASP-C": score_cgrnasp_c,
    "cgRNASP-CN": score_cgrnasp_cn,
    "cgRNASP-PC": score_cgrnasp_pc,
    "lociPARSE": score_lociparse,
    "rsRNASP": score_rsrnasp,
}

def main():
    parser = argparse.ArgumentParser(
        description="RNA structure scoring wrapper supporting multiple scoring methods"
    )
    parser.add_argument(
        "--scoring-method",
        choices=SCORING_METHODS,
        action="append",
        help="Scoring method to use (can be specified multiple times). If none given, all methods will be used.",
    )
    parser.add_argument(
        "pdb_files",
        nargs="+",
        help="One or more PDB files to score",
    )

    args = parser.parse_args()

    # Use all methods if none specified
    methods = args.scoring_method if args.scoring_method else SCORING_METHODS
    print(f"Using scoring methods: {methods}")

    # Validate PDB files exist
    for pdb_file in args.pdb_files:
        if not os.path.isfile(pdb_file):
            parser.error(f"PDB file not found: {pdb_file}")

    print(f"Processing PDB files: {args.pdb_files}")

    # Score each PDB file with selected methods
    for pdb_file in args.pdb_files:
        print(f"\nScoring {pdb_file}:")
        for method in methods:
            score = SCORING_FUNCTIONS[method](pdb_file)
            print(f"{method}: {score:.3f}")


if __name__ == "__main__":
    main()
