import pandas as pd
import os
from Bio.PDB import PDBParser, PDBIO, Select


class ProteinSelect(Select):
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


def build_seqres_lines(structure, start=None, end=None):
    seqres_lines = []
    for model in structure:
        for chain in model:
            residues = []
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                resseq = residue.id[1]
                if start is not None and resseq < start:
                    continue
                if end is not None and resseq > end:
                    continue
                residues.append(residue.get_resname())

            if not residues:
                continue

            num_res = len(residues)
            serial = 1
            for i in range(0, num_res, 13):
                chunk = residues[i : i + 13]
                line = f"SEQRES {serial:>3} {chain.id} {num_res:>4}  {' '.join(chunk)}\n"
                seqres_lines.append(line)
                serial += 1

    return seqres_lines


def process_pdb(pdb_path, output_path, start=None, end=None):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    # --- ATOM correction ---
    io = PDBIO()
    io.set_structure(structure)
    tmp_atom_file = output_path + ".tmp"
    io.save(tmp_atom_file, ProteinSelect(start, end))

    # --- SEQRES correction ---
    seqres_lines = build_seqres_lines(structure, start, end)

    # --- saving new PDB ---
    with open(output_path, "w") as out:
        for line in seqres_lines:
            out.write(line)
        with open(tmp_atom_file) as tmp:
            for line in tmp:
                if line.startswith("ATOM"):
                    out.write(line)
        out.write("END\n")

    os.remove(tmp_atom_file)


def run(csv_path, af_dir, output_dir):
    df = pd.read_csv(csv_path)
    os.makedirs(output_dir, exist_ok=True)

    for _, row in df.iterrows():
        af_id = row["AF_ID"]
        start = int(row["START"]) if pd.notna(row["START"]) else None
        end = int(row["END"]) if pd.notna(row["END"]) else None

        input_pdb = f"{af_dir}/{af_id}.pdb"
        output_pdb = f"{output_dir}/{af_id}.pdb"

        print(f"Processing {af_id} (start={start}, end={end})...")
        process_pdb(input_pdb, output_pdb, start, end)

run(
    csv_path="targets_list.csv",
    af_dir="alpha_fold",
    output_dir="cut_alpha_fold",
)