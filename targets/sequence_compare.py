import os
import argparse
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner, substitution_matrices

sub_count = 0
identity_results = []
coverage_results = []

low_results = []
med_results = []
high_results = []

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
        "alignment": best
    }

def save_alignment(pdb_id, af_id, dir_path, score, identity, coverage, alignment):
    with open(os.path.join(dir_path, f"{pdb_id}_{af_id}.txt"), "w") as f:
        f.write(f"Alignment score: {score:.2f}\n")
        f.write(f"Identity: {identity:.2f}\n")
        f.write(f"Coverage: {coverage:.2f}\n")
        f.write(str(alignment))

def process_alignment(pdb_seq, af_seq):
    is_substring = 0
    # Basic compare
    if is_subsequence(pdb_seq, af_seq):
        print("PDB sequence is a subsequence of AF sequence.")
        is_substring = 1

    # Alignment-based compare
    align_result = check_alignment(pdb_seq, af_seq)
    score = align_result['score']
    identity = align_result['identity']
    coverage = align_result['coverage']
    alignment = align_result['alignment']
    print(f"Alignment score: {score:.2f}")
    print(f"Identity: {identity:.2f}")
    print(f"Coverage: {coverage:.2f}")

    # Saving results
    low_path = "alignment_results/low_identity"
    med_path = "alignment_results/med_identity"
    high_path = "alignment_results/high_identity"
    os.makedirs(low_path, exist_ok=True)
    os.makedirs(med_path, exist_ok=True)
    os.makedirs(high_path, exist_ok=True)

    identity_results.append(identity)
    coverage_results.append(coverage)
    if identity < 0.8:
        low_results.append((pdb_id, af_id, identity))
        save_alignment(pdb_id, af_id, low_path, score, identity, coverage, alignment)
    elif identity < 0.95:
        med_results.append((pdb_id, af_id, identity))
        save_alignment(pdb_id, af_id, med_path, score, identity, coverage, alignment)
    else:
        high_results.append((pdb_id, af_id, identity))
        save_alignment(pdb_id, af_id, high_path, score, identity, coverage, alignment)

    return is_substring

parser = argparse.ArgumentParser(description="Compare PDB and AF sequences.")
parser.add_argument("--pdb_dir", default="pdb", help="Directory with PDB files")
parser.add_argument("--af_dir", default="alpha_fold", help="Directory with AlphaFold files")
args = parser.parse_args()

df = pd.read_csv("targets_list.csv")

for id, row in df.iterrows():
    pdb_id = row["PDB_ID"]
    af_id = row["AF_ID"]
    print(f"\n{id}: {pdb_id} - {af_id}")

    # PDB sequence
    pdb_file = os.path.join(args.pdb_dir, f"{pdb_id}.pdb")
    if not os.path.exists(pdb_file):
        print(f"File {pdb_file} does not exist. Skipping.")
        continue
    pdb_seq = extract_sequence(pdb_file)
    print(f"PDB sequence length: {len(pdb_seq)}")

    # AF sequence
    af_file = os.path.join(args.af_dir, f"{af_id}.pdb")
    if not os.path.exists(af_file):
        print(f"File {af_file} does not exist. Skipping.")
        continue
    af_seq = extract_sequence(af_file)
    print(f"AF sequence length: {len(af_seq)}")

    sub_count +=process_alignment(pdb_seq, af_seq)
    if len(pdb_seq)>len(af_seq):
        print("Warning: PDB sequence is longer than AF sequence.")

print("\n", "*" * 50)
print(f"Total PDB sequences that are subsequences of AF sequences: {sub_count} out of {len(df)}")
print(f"Average identity: {sum(identity_results)/len(identity_results):.2f}")
print(f"Full identities: {len([score for score in identity_results if score == 1.0])}")
print(f"Average coverage: {sum(coverage_results)/len(coverage_results):.2f}")

print("\n", "*" * 50)
print(f"Identities under 0.8 (LOW): {len(low_results)}")
low_results.sort(key=lambda x: x[2])
for pdb_id, af_id, identity in low_results:
    print(f"{pdb_id} - {af_id}: {identity:.2f}")

print("\n", "*" * 50)
print(f"Identities between 0.8 and 0.95 (MEDIUM): {len(med_results)}")
med_results.sort(key=lambda x: x[2])
for pdb_id, af_id, identity in med_results:
    print(f"{pdb_id} - {af_id}: {identity:.2f}")    

print("\n", "*" * 50)
print(f"Identities above 0.95 (HIGH): {len(high_results)}")
high_results.sort(key=lambda x: x[2])
for pdb_id, af_id, identity in high_results:
    print(f"{pdb_id} - {af_id}: {identity:.2f}")