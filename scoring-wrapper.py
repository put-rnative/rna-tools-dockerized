#! /usr/bin/env python3
import argparse
import os.path
import shutil
import subprocess
import tempfile
from multiprocessing.pool import ThreadPool
from pathlib import Path
from typing import Dict, List, Tuple, Optional

from lociPARSE import lociparse
import pandas as pd
from tqdm import tqdm


def run_command(
    cmd: List[str],
    *,
    stdout: Optional[int] = subprocess.PIPE,
    stderr: Optional[int] = subprocess.PIPE,
    expected_returncode: int = 0,
    **kwargs,
) -> subprocess.CompletedProcess:
    """Run command with better error handling"""
    try:
        result = subprocess.run(
            cmd,
            stdout=stdout,
            stderr=stderr,
            text=True,
            **kwargs,
        )
        if result.returncode != expected_returncode:
            raise subprocess.CalledProcessError(
                result.returncode, cmd, result.stdout, result.stderr
            )
        return result
    except subprocess.CalledProcessError as e:
        msg = f"Command failed with exit code {e.returncode}:\n"
        msg += f"Command: {' '.join(cmd)}\n"
        if e.stdout:
            msg += f"stdout:\n{e.stdout}\n"
        if e.stderr:
            msg += f"stderr:\n{e.stderr}\n"
        raise RuntimeError(msg) from e


SCORING_METHODS = [
    "3dRNAscore",
    "DFIRE",
    "RASP",
    "RNA-BRiQ",
    "RNA3DCNN_MD",
    "RNA3DCNN_MDMC",
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
        run_command(
            ["perl", "/opt/3dRNAscore/lib/format.pl", pdb_path],
            stdout=formatted,
        )

        # Run 3dRNAscore
        result = run_command(
            ["/opt/3dRNAscore/bin/3dRNAscore", "-s", formatted.name],
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
    result = run_command(["/opt/dfire/bin/DFIRE_RNA", pdb_path])

    try:
        # Parse the score from the second column of output
        score = float(result.stdout.strip().split()[1])
        return score
    except (ValueError, IndexError):
        raise RuntimeError(f"Failed to parse DFIRE output: {result.stdout}")


def score_rasp(pdb_path: str) -> float:
    """Score RNA structure using RASP method"""
    result = run_command(
        ["/opt/rasp-fd-1.0/bin/rasp_fd", "-p", pdb_path],
        stderr=subprocess.DEVNULL,
    )

    try:
        # Parse the score from first column of output
        score = float(result.stdout.strip().split()[0])
        return score
    except (ValueError, IndexError):
        raise RuntimeError(f"Failed to parse RASP output: {result.stdout}")


def score_rna_briq(pdb_path: str) -> float:
    """Score RNA structure using RNA-BRiQ method"""
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".txt") as ss_file:
        # Run BRiQ_AssignSS
        run_command(["/opt/RNA-BRiQ/build/bin/BRiQ_AssignSS", pdb_path, ss_file.name])

        # Add pdb path at the beginning of the file
        ss_file.seek(0)
        content = ss_file.read()
        ss_file.seek(0)
        ss_file.write(f"pdb {pdb_path}\n{content}")
        ss_file.flush()

        # Run BRiQ_Energy and capture stderr for score
        result = run_command(
            ["/opt/RNA-BRiQ/build/bin/BRiQ_Energy", ss_file.name],
            stdout=subprocess.DEVNULL,
        )

        try:
            # Parse score from second column
            score = float(result.stderr.strip().split()[1])
            return score
        except (ValueError, IndexError):
            raise RuntimeError(f"Failed to parse RNA-BRiQ output: {result.stderr}")


def _run_rna3dcnn(pdb_path: str, model_path: str) -> float:
    """Helper function to run RNA3DCNN with specified model"""
    result = run_command(
        [
            "/opt/RNA3DCNN/venv/bin/python",
            "/opt/RNA3DCNN/Main.py",
            "-pn",
            pdb_path,
            "-model",
            model_path,
            "-local",
            "0",
        ]
    )

    try:
        # Find the last line starting with "Total score for"
        for line in result.stdout.splitlines()[::-1]:
            if line.startswith("Total score for"):
                score = float(line.split()[-1])
                return score
        raise ValueError("No score line found")
    except (ValueError, IndexError):
        raise RuntimeError(f"Failed to parse RNA3DCNN output: {result.stdout}")


def score_rna3dcnn_md(pdb_path: str) -> float:
    """Score RNA structure using RNA3DCNN_MD method"""
    return _run_rna3dcnn(pdb_path, "/opt/RNA3DCNN/RNA3DCNN_MD.hdf5")


def score_rna3dcnn_mdmc(pdb_path: str) -> float:
    """Score RNA structure using RNA3DCNN_MDMC method"""
    return _run_rna3dcnn(pdb_path, "/opt/RNA3DCNN/RNA3DCNN_MDMC.hdf5")


def _run_cgrnasp(pdb_path: str, executable: str) -> float:
    """Helper function to run cgRNASP variants"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Copy PDB file to temp directory
        tmp_pdb = os.path.join(tmpdir, "input.pdb")
        shutil.copy2(pdb_path, tmp_pdb)

        # Create temp file for output
        tmp_out = os.path.join(tmpdir, "output.txt")

        # Get executable directory to use as cwd
        exe_dir = os.path.dirname(executable)

        # Run scoring from executable directory
        result = run_command(
            [executable, tmpdir, "1", tmp_out],
            cwd=exe_dir,
            expected_returncode=6,
        )

        # Read score from second column
        try:
            with open(tmp_out) as f:
                score = float(f.read().strip().split()[1])
                return score
        except (ValueError, IOError, IndexError) as e:
            raise RuntimeError(f"Failed to read score from {tmp_out}: {e}")


def score_cgrnasp(pdb_path: str) -> float:
    """Score RNA structure using cgRNASP method"""
    return _run_cgrnasp(pdb_path, "/opt/cgRNASP/cgRNASP/cgRNASP")


def score_cgrnasp_c(pdb_path: str) -> float:
    """Score RNA structure using cgRNASP-C method"""
    return _run_cgrnasp(pdb_path, "/opt/cgRNASP/cgRNASP-C/cgRNASP-C")


def score_cgrnasp_cn(pdb_path: str) -> float:
    """Score RNA structure using cgRNASP-CN method"""
    return _run_cgrnasp(pdb_path, "/opt/cgRNASP-CN/cgRNASP-CN")


def score_cgrnasp_pc(pdb_path: str) -> float:
    """Score RNA structure using cgRNASP-PC method"""
    return _run_cgrnasp(pdb_path, "/opt/cgRNASP/cgRNASP-PC/cgRNASP-PC")


def score_lociparse(pdb_path: str) -> float:
    """Score RNA structure using lociPARSE method"""
    lp = lociparse()
    return lp.score(pdb_path).pMoL.value


def score_rsrnasp(pdb_path: str) -> float:
    """Score RNA structure using rsRNASP method"""
    executable = "/opt/rsRNASP/rsRNASP"
    exe_dir = os.path.dirname(executable)

    with tempfile.NamedTemporaryFile(suffix=".txt") as tmp_out:
        result = run_command(
            [executable, pdb_path, tmp_out.name],
            expected_returncode=6,
            cwd=exe_dir,
        )

        # Read score from second column
        try:
            with open(tmp_out.name) as f:
                score = float(f.read().strip().split()[1])
                return score
        except (ValueError, IOError, IndexError) as e:
            raise RuntimeError(f"Failed to read score from {tmp_out.name}: {e}")


# Mapping of method names to scoring functions
SCORING_FUNCTIONS = {
    "3dRNAscore": score_3drnascore,
    "DFIRE": score_dfire,
    "RASP": score_rasp,
    "RNA-BRiQ": score_rna_briq,
    "RNA3DCNN_MD": score_rna3dcnn_md,
    "RNA3DCNN_MDMC": score_rna3dcnn_mdmc,
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
    parser.add_argument(
        "--output",
        help="Save results to CSV file",
        type=str,
        metavar="FILE.csv",
    )
    parser.add_argument(
        "--checkpoint",
        help="Load/save partial results from/to this file",
        type=str,
        metavar="FILE.csv",
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

    # Setup default checkpoint in user's home directory
    default_checkpoint = Path.home() / ".rna_scoring_checkpoint.csv"
    checkpoint_file = Path(args.checkpoint) if args.checkpoint else default_checkpoint

    # Load existing results if checkpoint exists
    results: Dict[str, Dict[str, float]] = {pdb_file: {} for pdb_file in args.pdb_files}
    if checkpoint_file.exists():
        checkpoint_df = pd.read_csv(checkpoint_file, index_col=0)
        for pdb_file in args.pdb_files:
            if pdb_file in checkpoint_df.index:
                for method in methods:
                    if method in checkpoint_df.columns and not pd.isna(
                        checkpoint_df.loc[pdb_file, method]
                    ):
                        results[pdb_file][method] = checkpoint_df.loc[pdb_file, method]

    # Count total tasks and already completed ones
    total_tasks = len(args.pdb_files) * len(methods)
    completed_tasks = sum(
        1 for pdb_file in args.pdb_files
        for method in methods
        if method in results.get(pdb_file, {})
    )

    # Filter out already computed tasks
    tasks = [
        (pdb_file, method)
        for pdb_file in args.pdb_files
        for method in methods
        if method not in results.get(pdb_file, {})
    ]

    if not tasks:
        print("\nAll requested scores already computed!")
    else:
        # Process remaining files in parallel, handling results as they complete
        with ThreadPool() as pool:
            with tqdm(
                total=total_tasks,
                initial=completed_tasks,
                desc="Scoring files"
            ) as pbar:
                for pdb_file, method, score in pool.imap_unordered(
                    lambda t: (t[0], t[1], SCORING_FUNCTIONS[t[1]](t[0])), tasks
                ):
                    pbar.update(1)
                # Update results and checkpoint as soon as each score is ready
                results[pdb_file][method] = score
                pd.DataFrame.from_dict(results, orient="index").to_csv(checkpoint_file)

    # Create DataFrame and display/save results
    df = pd.DataFrame.from_dict(results, orient="index")
    print("\nScoring Results:")
    print(df.round(3))

    if args.output:
        df.to_csv(args.output)
        print(f"\nResults saved to: {args.output}")


if __name__ == "__main__":
    main()
