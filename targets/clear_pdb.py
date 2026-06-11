import argparse
import pandas as pd
import os
from Bio.PDB import PDBParser, PDBIO, Select
import io as python_io

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Clean and standardize downloaded PDB structures using chain definitions from targets_list.csv."
    )
    parser.add_argument(
        "--targets-file",
        default="targets_list.csv",
        help="CSV file containing target mapping and chain information.",
    )
    parser.add_argument(
        "--pdb-dir",
        default="pdb",
        help="Directory containing downloaded experimental PDB files.",
    )
    parser.add_argument(
        "--output-dir",
        default="filtered_pdb",
        help="Directory where cleaned PDB files will be written.",
    )
    parser.add_argument(
        "--loud",
        action="store_true",
        help="Print detailed progress for each processed structure.",
    )
    return parser.parse_args()

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

def main():
    args = parse_arguments()
    df = pd.read_csv(args.targets_file)
    os.makedirs(args.output_dir, exist_ok=True)

    if not args.loud:
        print(f"Processing {len(df)} structures from {args.targets_file}...\n")

    processed = 0
    for _, row in df.iterrows():
        pdb_id = row["PDB_ID"]
        chain = row["CHAIN"]

        input_pdb = f"{args.pdb_dir}/{pdb_id}.pdb"
        output_pdb = f"{args.output_dir}/{pdb_id}.pdb"

        if args.loud:
            print(f"Processing {pdb_id} (chain={chain})...")
        process_pdb(input_pdb, output_pdb, chain)
        processed += 1

    print("=" * 40)
    print("Clean PDB summary")
    print(f"Processed: {processed}")
    print(f"Saved cleaned files in: {args.output_dir}")
    print("=" * 40)
    print()

if __name__ == "__main__":
    main()