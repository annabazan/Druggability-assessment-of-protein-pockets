import pandas as pd
import os
from Bio.PDB import PDBParser, PDBIO, Select
import io as python_io


class ProteinSelect(Select):
    def __init__(self, chain_id=None):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        if self.chain_id is None:
            return True
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        # only standard amino acids (not heteroatoms or water)
        return residue.id[0] == " "


def extract_seqres(lines, selected_chain=None):
    """
    Filter SEQRES lines only for the selected chain 
    (or all chains if selected_chain is None).
    """
    seqres_lines = []
    for line in lines:
        if line.startswith("SEQRES"):
            chain = line[11]
            if selected_chain is None or chain == selected_chain:
                seqres_lines.append(line)
    return seqres_lines


def process_pdb(pdb_path, output_path, chain_id):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    selected_chain = None if chain_id == "-" else chain_id

    # save selected chain to string buffer
    string_io = python_io.StringIO()
    io_pdb = PDBIO()
    io_pdb.set_structure(structure)
    io_pdb.save(string_io, ProteinSelect(selected_chain))
    
    # get ATOM lines from buffer
    atom_lines = string_io.getvalue().splitlines()

    # extract SEQRES lines from original file
    with open(pdb_path) as f:
        original_lines = f.readlines()
    seqres_lines = extract_seqres(original_lines, selected_chain)

    # final file
    with open(output_path, "w") as out:
        # Zapisujemy SEQRES
        for line in seqres_lines:
            out.write(line.rstrip() + "\n")
        # Zapisujemy tylko linie ATOM z pamięci
        for line in atom_lines:
            if line.startswith("ATOM"):
                out.write(line + "\n")
        out.write("END\n")


def run(csv_path, pdb_dir, output_dir):
    df = pd.read_csv(csv_path)
    os.makedirs(output_dir, exist_ok=True)

    for _, row in df.iterrows():
        pdb_id = row["PDB_ID"]
        chain = row["CHAIN"]

        input_pdb = f"{pdb_dir}/{pdb_id}.pdb"
        output_pdb = f"{output_dir}/{pdb_id}.pdb"

        print(f"Processing {pdb_id} (chain={chain})...")
        process_pdb(input_pdb, output_pdb, chain)

run(
    csv_path="targets_list.csv",
    pdb_dir="pdb",
    output_dir="filtered_pdb"
)