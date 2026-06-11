import os
import argparse
import pandas as pd
import requests
from Bio.PDB import MMCIFParser, PDBIO

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Download experimental PDB structures listed in targets_list.csv."
    )
    parser.add_argument(
        "--targets-file",
        default="targets_list.csv",
        help="CSV file containing PDB_ID values.",
    )
    parser.add_argument(
        "--output-dir",
        default="pdb",
        help="Directory in which to save downloaded PDB files.",
    )
    parser.add_argument(
        "--loud",
        action="store_true",
        help="Print detailed progress for each download.",
    )
    return parser.parse_args()

def download_pdb_structure(pdb_id, output_dir, parser, loud=False):
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    pdb_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    cif_path = os.path.join(output_dir, f"{pdb_id}.cif")

    try:
        response = requests.get(pdb_url, timeout=30)
    except requests.RequestException as err:
        if loud:
            print(f"{pdb_id} -> PDB request failed: {err}")
        response = None

    if response is not None and response.status_code == 200:
        with open(pdb_path, "w") as f:
            f.write(response.text)
        return "downloaded_pdb"

    if loud:
        print(f"{pdb_id} -> PDB not found, trying CIF")

    try:
        response = requests.get(cif_url, timeout=30)
    except requests.RequestException as err:
        if loud:
            print(f"{pdb_id} -> CIF request failed: {err}")
        return "failed"

    if response.status_code != 200:
        if loud:
            print(f"{pdb_id} -> CIF not found")
        return "failed"

    with open(cif_path, "w") as f:
        f.write(response.text)

    try:
        structure = parser.get_structure(pdb_id, cif_path)
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_path)
        if loud:
            print(f"{pdb_id} -> converted CIF to PDB")
        return "converted_cif"
    except Exception as err:
        if loud:
            print(f"{pdb_id} -> failed to convert CIF to PDB: {err}")
        return "failed"

def main():
    args = parse_arguments()
    df = pd.read_csv(args.targets_file)
    pdb_ids = df["PDB_ID"].dropna().unique()

    os.makedirs(args.output_dir, exist_ok=True)
    parser = MMCIFParser(QUIET=True)

    if not args.loud:
        print(f"Downloading {len(pdb_ids)} experimental structures from {args.targets_file}...\n")

    summary = {"downloaded_pdb": 0, "converted_cif": 0, "failed": 0}
    for pdb_id in pdb_ids:
        status = download_pdb_structure(pdb_id, args.output_dir, parser, loud=args.loud)
        if status in summary:
            summary[status] += 1
        else:
            summary["failed"] += 1
        if args.loud:
            print(f"{pdb_id}: {status}")

    print("=" * 40)
    print("Download summary:")
    print(f"Total requested: {len(pdb_ids)}")
    print(f"Downloaded PDB: {summary['downloaded_pdb']}")
    print(f"Failed: {summary['failed']}")
    print(f"Saved structures in: {args.output_dir}")
    print("=" * 40)
    print()

if __name__ == "__main__":
    main()