import os
import argparse
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner, substitution_matrices

def extract_sequence(pdb_path, chain_id=None):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    sequence = ""

    for model in structure:
        for chain in model:
            if chain_id and chain.id != chain_id:
                continue
            for residue in chain:
                if residue.id[0] == " ":
                    try:
                        sequence += seq1(residue.resname)
                    except:
                        sequence += "X"
    return sequence

def is_subsequence(pdb_seq, af_seq):
    return pdb_seq in af_seq

def check_alignment(pdb_seq, af_seq):
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    best = aligner.align(af_seq, pdb_seq)[0]

    matches = 0
    aligned_length = 0

    for (a_start, a_end), (b_start, b_end) in zip(best.aligned[0], best.aligned[1]):
        af_fragment = af_seq[a_start:a_end]
        pdb_fragment = pdb_seq[b_start:b_end]
        aligned_length += len(pdb_fragment)

        for a, b in zip(af_fragment, pdb_fragment):
            if a == b:
                matches += 1

    identity = matches / aligned_length if aligned_length > 0 else 0
    coverage = aligned_length / len(pdb_seq) if len(pdb_seq) > 0 else 0

    return {
        "score": best.score,
        "identity": identity,
        "coverage": coverage,
        "alignment": best,
    }

def save_alignment(pdb_id, af_id, dir_path, score, identity, coverage, alignment):
    with open(os.path.join(dir_path, f"{pdb_id}_{af_id}.txt"), "w") as f:
        f.write(f"Alignment score: {score:.2f}\n")
        f.write(f"Identity: {identity:.2f}\n")
        f.write(f"Coverage: {coverage:.2f}\n")
        f.write(str(alignment))

def process_alignment(pdb_id, af_id, pdb_seq, af_seq, loud=False):
    is_substring = 1 if is_subsequence(pdb_seq, af_seq) else 0
    alignment_result = check_alignment(pdb_seq, af_seq)
    score = alignment_result["score"]
    identity = alignment_result["identity"]
    coverage = alignment_result["coverage"]
    alignment = alignment_result["alignment"]

    if loud:
        if is_substring:
            print("PDB sequence is a subsequence of AF sequence.")
        print(f"Alignment score: {score:.2f}")
        print(f"Identity: {identity:.2f}")
        print(f"Coverage: {coverage:.2f}")
        if len(pdb_seq) > len(af_seq):
            print("Warning: PDB sequence is longer than AF sequence.")

    low_path = "alignment_results/low_identity"
    med_path = "alignment_results/med_identity"
    high_path = "alignment_results/high_identity"
    os.makedirs(low_path, exist_ok=True)
    os.makedirs(med_path, exist_ok=True)
    os.makedirs(high_path, exist_ok=True)

    if identity < 0.8:
        save_alignment(pdb_id, af_id, low_path, score, identity, coverage, alignment)
        classification = "low"
    elif identity < 0.9:
        save_alignment(pdb_id, af_id, med_path, score, identity, coverage, alignment)
        classification = "medium_80_90"
    elif identity < 0.95:
        save_alignment(pdb_id, af_id, med_path, score, identity, coverage, alignment)
        classification = "medium_90_95"
    else:
        save_alignment(pdb_id, af_id, high_path, score, identity, coverage, alignment)
        classification = "high"

    return {
        "PDB_ID": pdb_id,
        "AF_ID": af_id,
        "pdb_length": len(pdb_seq),
        "af_length": len(af_seq),
        "is_subsequence": is_substring,
        "score": score,
        "identity": identity,
        "coverage": coverage,
        "classification": classification,
        "full_identity": 1 if identity == 1.0 else 0,
        "pdb_longer_than_af": 1 if len(pdb_seq) > len(af_seq) else 0,
    }

def main():
    parser = argparse.ArgumentParser(description="Compare PDB and AF sequences.")
    parser.add_argument("--pdb_dir", default="pdb", help="Directory with PDB files")
    parser.add_argument("--af_dir", default="alpha_fold", help="Directory with AlphaFold files")
    parser.add_argument("--loud", action="store_true", help="Print detailed comparison output for each pair.")
    args = parser.parse_args()

    df = pd.read_csv("targets_list.csv")
    if not args.loud:
        print(f"Starting sequence comparison for {len(df)} targets...\n")

    results = []
    subseq_count = 0
    identity_values = []
    coverage_values = []

    for index, row in df.iterrows():
        pdb_id = row["PDB_ID"]
        af_id = row["AF_ID"]

        pdb_file = os.path.join(args.pdb_dir, f"{pdb_id}.pdb")
        if not os.path.exists(pdb_file):
            if args.loud:
                print(f"File {pdb_file} does not exist. Skipping.")
            continue

        af_file = os.path.join(args.af_dir, f"{af_id}.pdb")
        if not os.path.exists(af_file):
            if args.loud:
                print(f"File {af_file} does not exist. Skipping.")
            continue

        pdb_seq = extract_sequence(pdb_file)
        af_seq = extract_sequence(af_file)

        if args.loud:
            print(f"\n{index}: {pdb_id} - {af_id}")
            print(f"PDB sequence length: {len(pdb_seq)}")
            print(f"AF sequence length: {len(af_seq)}")

        row_result = process_alignment(pdb_id, af_id, pdb_seq, af_seq, loud=args.loud)
        results.append(row_result)
        subseq_count += row_result["is_subsequence"]
        identity_values.append(row_result["identity"])
        coverage_values.append(row_result["coverage"])

    results_df = pd.DataFrame(results)
    results_df.to_csv("sequence_compare_results.csv", index=False)

    count_0_8 = results_df[results_df["identity"] < 0.8].shape[0]
    count_0_8_0_9 = results_df[(results_df["identity"] >= 0.8) & (results_df["identity"] < 0.9)].shape[0]
    count_0_9_0_95 = results_df[(results_df["identity"] >= 0.9) & (results_df["identity"] < 0.95)].shape[0]
    count_0_95_1_0 = results_df[results_df["identity"] >= 0.95].shape[0]
    full_identities = results_df[results_df["identity"] == 1.0].shape[0]

    print("=" * 40)
    print("Sequence comparison summary")
    print(f"Total PDB sequences that are subsequences of AF sequences: {subseq_count} out of {len(results_df)}")
    print(f"Average identity: {sum(identity_values) / len(identity_values):.2f}" if identity_values else "Average identity: 0.00")
    print(f"Full identities: {full_identities}")
    print(f"Average coverage: {sum(coverage_values) / len(coverage_values):.2f}" if coverage_values else "Average coverage: 0.00")
    print("\nIdentity ranges:")
    print(f"  0.0 <= identity < 0.8: {count_0_8}")
    print(f"  0.8 <= identity < 0.9: {count_0_8_0_9}")
    print(f"  0.9 <= identity < 0.95: {count_0_9_0_95}")
    print(f"  0.95 <= identity <= 1.0: {count_0_95_1_0}")
    print("Saved sequence comparison report: sequence_compare_results.csv")
    print("=" * 40)
    print()

if __name__ == "__main__":
    main()
