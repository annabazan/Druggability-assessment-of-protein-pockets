#!/usr/bin/env python3
import argparse
import os
import subprocess
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select

CSV_PATH="targets_list.csv"
AF_DIR = "alpha_fold"
CUT_DIR = "cut_alpha_fold"
FILTERED_DIR = "filtered_alpha_fold"
COMPLETE_DIR = "filtered_cut_alpha_fold"
PYMOL_DIR = "pymol"
PYMOL_SCRIPT_PATH = "temp_script.pml"

class RangeSelect(Select):
    def __init__(self, start=None, end=None):
        self.start = start
        self.end = end

    def accept_residue(self, residue):
        if residue.id[0] != " ":
            return False
        resseq = residue.id[1]
        if self.start is not None and resseq < self.start:
            return False
        if self.end is not None and resseq > self.end:
            return False
        return True
    

class PLDDTSelect(Select):
    def __init__(self, scores, threshold):
        self.scores = scores
        self.threshold = threshold

    def accept_residue(self, residue):
        if residue.id[0] != " ":
            return False
        key = (residue.get_parent().id, residue.id)
        if key not in self.scores:
            return False
        return self.scores[key]>=self.threshold


class CompleteSelect(Select):
    def __init__(self, selectors):
        self.selectors = selectors

    def accept_residue(self, residue):
        return all(
            sel.accept_residue(residue) 
            for sel in self.selectors
            )             


def compute_plddt(structure):
    scores = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                b = [atom.get_bfactor() for atom in residue]
                if b:
                    scores[(chain.id, residue.id)] = sum(b)/len(b)
    return scores


def build_seqres(structure, select):
    seqres = []
    for model in structure:
        for chain in model:
            residues = []
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                if not select.accept_residue(residue):
                    continue
                residues.append(residue.get_resname())
            if not residues: continue

            num_res = len(residues)
            serial = 1
            for i in range(0, num_res, 13):
                chunk = residues[i : i + 13]
                line = f"SEQRES {serial:>3} {chain.id} {num_res:>4}  {' '.join(chunk)}\n"
                seqres.append(line)
                serial += 1
    return seqres


def process_structure(structure, output_path, select):
    io = PDBIO()
    io.set_structure(structure)
    tmp_atom_file = output_path+".tmp"
    io.save(tmp_atom_file, select)

    # --- SEQRES correction ---
    seqres = build_seqres(structure, select)

    # --- saving new PDB ---
    with open(output_path, "w") as out:
        for line in seqres:
            out.write(line)
        with open(tmp_atom_file) as tmp:
            for line in tmp:
                if line.startswith("ATOM"):
                    out.write(line)
        out.write("END\n")
    os.remove(tmp_atom_file)

# load alpha_fold/O60885.pdb, full
# load filtered_alpha_fold_70_test/O60885.pdb, filtered

def generate_pymol_script(input_pdb, filtered_pdb, output_png):
    script = f"""
load {input_pdb}, full 
load {filtered_pdb}, filtered

hide everything
show cartoon, full
show cartoon, filtered

spectrum b, blue_cyan_green_yellow_orange_red, minimum=0, maximum=100
set cartoon_smooth_loops, 1
set cartoon_fancy_helices, 1
bg_color white

align filtered, full

python
from pymol import cmd

min_full, max_full = cmd.get_extent("full")
width_full = max_full[0] - min_full[0]

min_filt, max_filt = cmd.get_extent("filtered")
width_filt = max_filt[0] - min_filt[0]

shift = width_full + width_filt + 10

cmd.translate([shift, 0, 0], "filtered")
python end

zoom all
ray 1600,800
png {output_png}
quit
"""
    return script


def run_pymol(input_pdb, filtered_pdb, output_png):
    script_content = generate_pymol_script(input_pdb, filtered_pdb, output_png)
    with open(PYMOL_SCRIPT_PATH, "w") as f:
        f.write(script_content)
    subprocess.run(["pymol", "-cq", PYMOL_SCRIPT_PATH])
    os.remove(PYMOL_SCRIPT_PATH)


def process_target(target, output_dir, args):
    id = target["AF_ID"]
    input_pdb = f"{AF_DIR}/{id}.pdb"
    output_pdb = f"{output_dir}/{id}.pdb"
    print(f"Processing {id}...")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)

    selectors = []
    if args.mode in ["range", "complete"]:
        start = int(target["START"]) if pd.notna(target["START"]) else None
        end = int(target["END"]) if pd.notna(target["END"]) else None
        selectors.append(RangeSelect(start, end))
    if args.mode in ["plddt", "complete"]:
        scores = compute_plddt(structure)
        selectors.append(PLDDTSelect(scores, args.plddt))
    if len(selectors)==1:
        select = selectors[0]
    else:
        select = CompleteSelect(selectors)

    process_structure(structure, output_pdb, select)
    if args.visualize:
        print(f"Visualizing results for {id}...")
        output_png = f"{output_dir}/{PYMOL_DIR}/{id}.png"
        run_pymol(input_pdb, output_pdb, output_png)
        print(f"Visualization saved in {output_png} file.")


def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter AlphaFold PDB files.")
    parser.add_argument(
        "--mode",
        choices=["range", "plddt", "complete"],
        default="range",
        help="Filtering mode"
    )
    parser.add_argument(
        "--plddt",
        type=float,
        default=70.0,
        help="pLDDT threshold"
    )
    parser.add_argument(
        "--visualize",
        action="store_true",
        help="Enable visualization"
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.mode=="range":
        output_dir = CUT_DIR
    elif args.mode=="plddt":
        output_dir = FILTERED_DIR
    elif args.mode=="complete":
        output_dir = COMPLETE_DIR

    targets = pd.read_csv(CSV_PATH)
    os.makedirs(output_dir, exist_ok=True)

    if args.visualize:
        pymol_dir = f"{output_dir}/{PYMOL_DIR}"
        os.makedirs(pymol_dir, exist_ok=True)

    for _, target in targets.iterrows():
        process_target(target, output_dir, args)


if __name__ == "__main__":
    main()

