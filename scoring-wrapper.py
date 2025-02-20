#! /usr/bin/env python
import argparse

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
    import os.path

    for pdb_file in args.pdb_files:
        if not os.path.isfile(pdb_file):
            parser.error(f"PDB file not found: {pdb_file}")

    print(f"Processing PDB files: {args.pdb_files}")


if __name__ == "__main__":
    main()
