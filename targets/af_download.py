import os
import pandas as pd
import requests

# === 1. Wczytaj CSV ===
df = pd.read_csv("targets_list.csv")

# === 2. UniProt IDs ===
af_ids = df["AF ID"].dropna().unique()

# === 3. Folder ===
output_dir = "alpha_fold"
os.makedirs(output_dir, exist_ok=True)

# === 4. Pobieranie ===
for af_id in af_ids:
    url = f"https://alphafold.ebi.ac.uk/files/AF-{af_id}-F1-model_v6.pdb"
    file_path = os.path.join(output_dir, f"{af_id}.pdb")
    
    response = requests.get(url)
    
    if response.status_code == 200:
        with open(file_path, "w") as f:
            f.write(response.text)
        print(f"{af_id} downloaded")
    else:
        print(f"{af_id} failed")