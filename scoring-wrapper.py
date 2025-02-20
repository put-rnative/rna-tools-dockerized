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
    "rsRNASP"
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--scoring-method",
        choices=SCORING_METHODS,
        action="append",
        help="Scoring method to use (can be specified multiple times). If none given, all methods will be used."
    )
