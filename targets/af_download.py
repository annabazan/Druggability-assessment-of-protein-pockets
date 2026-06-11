import os
import argparse
import pandas as pd
import requests

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Download AlphaFold models listed in targets_list.csv."
    )
    parser.add_argument(
        "--targets-file",
        default="targets_list.csv",
        help="CSV file containing AF_ID values.",
    )
    parser.add_argument(
        "--output-dir",
        default="alpha_fold",
        help="Directory in which to save downloaded AlphaFold PDB files.",
    )
    parser.add_argument(
        "--loud",
        action="store_true",
        help="Print detailed progress for each download.",
    )
    return parser.parse_args()

def download_af_model(af_id, output_dir, loud=False):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{af_id}-F1-model_v6.pdb"
    file_path = os.path.join(output_dir, f"{af_id}.pdb")

    try:
        response = requests.get(url, timeout=30)
    except requests.RequestException as err:
        if loud:
            print(f"{af_id} -> request failed: {err}")
        return False

    if response.status_code == 200:
        with open(file_path, "w") as f:
            f.write(response.text)
        return True

    if loud:
        print(f"{af_id} -> download failed with status code {response.status_code}")
    return False

def main():
    args = parse_arguments()
    df = pd.read_csv(args.targets_file)
    af_ids = df["AF_ID"].dropna().unique()

    os.makedirs(args.output_dir, exist_ok=True)

    if not args.loud:
        print(f"Downloading {len(af_ids)} AlphaFold models from {args.targets_file}...\n")

    success = 0
    failed = 0
    for af_id in af_ids:
        ok = download_af_model(af_id, args.output_dir, loud=args.loud)
        if ok:
            success += 1
            if args.loud:
                print(f"{af_id} downloaded")
        else:
            failed += 1
            if args.loud:
                print(f"{af_id} failed")

    print("=" * 40)
    print("Download summary")
    print(f"Total requested: {len(af_ids)}")
    print(f"Downloaded: {success}")
    print(f"Failed: {failed}")
    print(f"Saved models in: {args.output_dir}")
    print("=" * 40)
    print()
    
if __name__ == "__main__":
    main()