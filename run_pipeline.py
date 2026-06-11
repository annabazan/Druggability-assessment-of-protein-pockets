#!/usr/bin/env python3
"""
Pipeline runner for the project.
Run from the repository root.

Default behavior: no PyMOL visualizations, quiet mode. Use `--loud` for verbose logs
and `--visualize` to enable PyMOL steps where available.
"""
import argparse
import subprocess
import os
import sys

def run_cmd(cmd, cwd=None, loud=False):
    if loud:
        print(f"Running in {cwd or os.getcwd()}: {' '.join(cmd)}")
    res = subprocess.run(cmd, cwd=cwd)
    if res.returncode != 0:
        raise RuntimeError(f"Command failed (rc={res.returncode}): {' '.join(cmd)}")
    return res.returncode

def run_all(visualize=False, loud=False, skip_fpocket=False, skip_analysis=False):
    repo_root = os.path.dirname(os.path.abspath(__file__))

    # 1) Download structures
    print("[1/8] Downloading experimental PDB structures...")
    cmd = [sys.executable, "pdb_download.py"]
    if loud:
        cmd.append("--loud")
    run_cmd(cmd, cwd=os.path.join(repo_root, "targets"), loud=loud)

    print("[2/8] Downloading AlphaFold models...")
    cmd = [sys.executable, "af_download.py"]
    if loud:
        cmd.append("--loud")
    run_cmd(cmd, cwd=os.path.join(repo_root, "targets"), loud=loud)

    # 2) Preprocess structures
    print("[3/8] Cleaning PDB structures (clear_pdb.py)...")
    cmd = [sys.executable, "clear_pdb.py", "--pdb-dir", "pdb", "--output-dir", "filtered_pdb"]
    if loud:
        cmd.append("--loud")
    run_cmd(cmd, cwd=os.path.join(repo_root, "targets"), loud=loud)

    print("[4/8] Filtering AlphaFold models (filter_af.py)...")
    cmd = [sys.executable, "filter_af.py", "--mode", "range"]
    if visualize:
        cmd.append("--visualize")
    if loud:
        cmd.append("--loud")
    run_cmd(cmd, cwd=os.path.join(repo_root, "targets"), loud=loud)

    # 3) Sequence comparison
    print("[5/8] Comparing sequences (sequence_compare.py)...")
    cmd = [sys.executable, "sequence_compare.py", "--pdb_dir", "filtered_pdb", "--af_dir", "cut_alpha_fold"]
    if loud:
        cmd.append("--loud")
    run_cmd(cmd, cwd=os.path.join(repo_root, "targets"), loud=loud)

    # 4) 3D alignment
    print("[6/8] Running 3D alignment (3D_align_superimposer.py)...")
    cmd = [sys.executable, "3D_align_superimposer.py"]
    if visualize:
        cmd.append("--visualize")
    if loud:
        cmd.append("--loud")
    run_cmd(cmd, cwd=os.path.join(repo_root, "targets"), loud=loud)

    # 5) fpocket on PDB and AlphaFold
    if not skip_fpocket:
        print("[7/8] Running fpocket on PDB structures...")
        run_cmd(["bash", "run_fpocket_pdb.sh"], cwd=os.path.join(repo_root, "pocket_detection/fpocket"), loud=loud)

        print("\n[8/8] Running fpocket on aligned AlphaFold structures...")
        run_cmd(["bash", "run_fpocket_alpha_fold.sh"], cwd=os.path.join(repo_root, "pocket_detection/fpocket"), loud=loud)

    # 6) Analysis
    if not skip_analysis:
        print("\n[final] Running analysis (run_analysis.py)...")
        analysis_cmd = [sys.executable, os.path.join("analysis", "pocket_comparison", "run_analysis.py"),
                        "--targets", os.path.join("targets", "targets_list.csv"),
                        "--pdb-fpocket", os.path.join("pocket_detection", "fpocket", "pdb_out"),
                        "--af-fpocket", os.path.join("pocket_detection", "fpocket", "alpha_fold_out"),
                        "--pdb-structures", os.path.join("targets", "filtered_pdb"),
                        "--af-structures", os.path.join("targets", "3D_aligned_alpha_fold"),
                        "--out", os.path.join("analysis", "pocket_comparison", "outputs")]
        run_cmd(analysis_cmd, cwd=repo_root, loud=loud)

    print("Pipeline finished successfully.")


def parse_arguments():
    parser = argparse.ArgumentParser(description="Run full project pipeline from repository root.")
    parser.add_argument("--visualize", action="store_true", help="Enable PyMOL visualizations where supported.")
    parser.add_argument("--loud", action="store_true", help="Print verbose progress for each step.")
    parser.add_argument("--skip-fpocket", action="store_true", help="Skip fpocket runs.")
    parser.add_argument("--skip-analysis", action="store_true", help="Skip final analysis step.")
    return parser.parse_args()


def main():
    args = parse_arguments()
    run_all(visualize=args.visualize, loud=args.loud, skip_fpocket=args.skip_fpocket, skip_analysis=args.skip_analysis)


if __name__ == "__main__":
    main()
