import os
import pandas as pd
import requests
from Bio.PDB import MMCIFParser, PDBIO

df = pd.read_csv("targets_list.csv")
pdb_ids = df["PDB_ID"].dropna().unique()

output_dir = "pdb"
os.makedirs(output_dir, exist_ok=True)

parser = MMCIFParser(QUIET=True)

for pdb_id in pdb_ids:
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    pdb_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    cif_path = os.path.join(output_dir, f"{pdb_id}.cif")
    
    response = requests.get(pdb_url)
    
    if response.status_code == 200:
        with open(pdb_path, "w") as f:
            f.write(response.text)
        print(f"{pdb_id} -> downloaded PDB")
    else:
        print(f"{pdb_id} -> PDB not found, trying CIF")
        response = requests.get(cif_url)
        if response.status_code == 200:
            with open(cif_path, "w") as f:
                f.write(response.text)
            print(f"{pdb_id} -> downloaded CIF")
            try:
                structure = parser.get_structure(pdb_id, cif_path)
                io = PDBIO()
                io.set_structure(structure)
                io.save(pdb_path)
                print(f"{pdb_id} -> converted CIF to PDB")
            except Exception as e:
                print(f"{pdb_id} -> failed to convert CIF to PDB: {e}")
        else:
            print(f"{pdb_id} failed")