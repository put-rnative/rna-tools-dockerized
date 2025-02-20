#! /usr/bin/env python3
import argparse
import os.path
import subprocess
import tempfile
from multiprocessing.pool import ThreadPool
from typing import Dict, List, Tuple

import pandas as pd

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
    with tempfile.NamedTemporaryFile(suffix=".pdb") as formatted:
        # Format the PDB file
        subprocess.run(
            ["perl", "/opt/3dRNAscore/lib/format.pl", pdb_path],
            stdout=formatted,
            check=True,
            text=True,
        )

        # Run 3dRNAscore
        result = subprocess.run(
            ["/opt/3dRNAscore/bin/3dRNAscore", "-s", formatted.name],
            capture_output=True,
            text=True,
            check=True,
        )

        # Parse the score from output
        try:
            # Assuming the score is the last number in the output
            score = float(result.stdout.strip().split()[-1])
            return score
        except (ValueError, IndexError):
            raise RuntimeError(f"Failed to parse 3dRNAscore output: {result.stdout}")


def score_dfire(pdb_path: str) -> float:
    """Score RNA structure using DFIRE method"""
    result = subprocess.run(
        ["/opt/dfire/bin/DFIRE_RNA", pdb_path],
        capture_output=True,
        text=True,
        check=True,
    )

    try:
        # Parse the score from the second column of output
        score = float(result.stdout.strip().split()[1])
        return score
    except (ValueError, IndexError):
        raise RuntimeError(f"Failed to parse DFIRE output: {result.stdout}")


def score_rasp(pdb_path: str) -> float:
    """Score RNA structure using RASP method"""
    result = subprocess.run(
        ["/opt/rasp-fd-1.0/bin/rasp_fd", "-p", pdb_path],
        capture_output=True,
        text=True,
        check=True,
    )

    try:
        # Parse the score from output
        score = float(result.stdout.strip())
        return score
    except (ValueError, IndexError):
        raise RuntimeError(f"Failed to parse RASP output: {result.stdout}")


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

    # Prepare all scoring tasks
    tasks: List[Tuple[str, str]] = [
        (pdb_file, method) for pdb_file in args.pdb_files for method in methods
    ]

    # Process files in parallel while preserving order
    print("\nScoring files in parallel...")
    with ThreadPool() as pool:
        scores = pool.starmap(lambda f, m: (f, m, SCORING_FUNCTIONS[m](f)), tasks)

    # Collect results preserving file order
    results: Dict[str, Dict[str, float]] = {pdb_file: {} for pdb_file in args.pdb_files}
    for pdb_file, method, score in scores:
        results[pdb_file][method] = score

    # Create DataFrame and display results
    df = pd.DataFrame.from_dict(results, orient="index")
    print("\nScoring Results:")
    print(df.round(3))


if __name__ == "__main__":
    main()
