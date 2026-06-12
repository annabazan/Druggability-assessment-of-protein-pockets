import argparse
from tabnanny import verbose
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio.SeqUtils import seq1
import os
import pandas as pd

try:
    import pymol
    from pymol import cmd
    HAVE_PYMOL = True
except ImportError:
    pymol = None
    cmd = None
    HAVE_PYMOL = False

def get_ca_sequence(structure):
    residues = []
    for residue in structure.get_residues():
        if "CA" not in residue:
            continue
        try:
            aa = seq1(residue.get_resname())
        except Exception:
            aa = "X"
        residues.append((
            residue.get_parent().id,
            residue.get_id()[1],
            residue.get_id()[2].strip(),
            aa,
            residue["CA"],
        ))
    return residues

def map_atoms_by_alignment(ref_list, sample_list, alignment):
    ref_atoms = []
    sample_atoms = []

    ref_seq = "".join([item[3] for item in ref_list])
    sample_seq = "".join([item[3] for item in sample_list])

    for (r0, r1), (s0, s1) in zip(alignment.aligned[0], alignment.aligned[1]):
        for ri, si in zip(range(r0, r1), range(s0, s1)):
            if ref_seq[ri] != sample_seq[si]:
                continue
            ref_atoms.append(ref_list[ri][4])
            sample_atoms.append(sample_list[si][4])

    return ref_atoms, sample_atoms

def superimpose_robust(pdb_path, af_path, output_path, direct_numbering=False, verbose=False):
    parser = PDBParser(QUIET=True)
    ref_struct = parser.get_structure("ref", pdb_path)
    sample_struct = parser.get_structure("sample", af_path)

    ref_list = get_ca_sequence(ref_struct)
    sample_list = get_ca_sequence(sample_struct)
    if verbose:
        print(f"Reference residues: {len(ref_list)}, Sample residues: {len(sample_list)}")

    ref_by_number = {item[1]: item[4] for item in ref_list}
    sample_by_number = {item[1]: item[4] for item in sample_list}
    common_res_ids = sorted(set(ref_by_number.keys()) & set(sample_by_number.keys()))

    if direct_numbering:
        if common_res_ids:
            ref_atoms = [ref_by_number[res_id] for res_id in common_res_ids]
            sample_atoms = [sample_by_number[res_id] for res_id in common_res_ids]
            method = "direct numbering"
            if verbose:
                print(f"Using direct numbering for {len(ref_atoms)} matched CA atoms.")
        else:
            if verbose:
                print("Direct numbering requested but no common residue numbers were found; falling back to sequence alignment.")
            ref_atoms, sample_atoms, method = None, None, None
    else:
        ref_atoms, sample_atoms, method = None, None, None

    if ref_atoms is None:
        ref_seq = "".join([item[3] for item in ref_list])
        sample_seq = "".join([item[3] for item in sample_list])

        aligner = PairwiseAligner()
        aligner.mode = "local"
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5

        alignment = aligner.align(ref_seq, sample_seq)[0]
        ref_atoms, sample_atoms = map_atoms_by_alignment(ref_list, sample_list, alignment)
        method = "sequence alignment"
        if verbose:
            print(f"Using sequence alignment for {len(ref_atoms)} matched CA atoms.")

    if len(ref_atoms) < 3:
        if verbose:
            print("Not enough matched CA atoms for superposition.")
        return None, method, len(ref_list), len(sample_list), len(ref_atoms)

    if verbose:
        print(f"Superimposing using {len(ref_atoms)} CA atoms ({method}).")

    si = Superimposer()
    si.set_atoms(ref_atoms, sample_atoms)
    si.apply(sample_struct.get_atoms())

    io = PDBIO()
    io.set_structure(sample_struct)
    io.save(output_path)
    
    return si.rms, method, len(ref_list), len(sample_list), len(ref_atoms)

def visualize_alignment(pdb_file, alphafold_file, output_name="comparison.png", verbose=False):
    # PyMOL without GUI
    pymol.finish_launching(['pymol', '-cq']) 

    cmd.reinitialize()
    cmd.load(pdb_file, "experimental")
    cmd.load(alphafold_file, "alphafold")
    cmd.show_as("cartoon")
    cmd.color("slate", "experimental")
    cmd.color("orange", "alphafold")
    
    # structure next to each other for better visualization
    cmd.translate([80, 0, 0], "alphafold")
    
    # centering the view
    cmd.zoom("all", buffer=5)
    
    # render settings for better quality
    cmd.set("ray_opaque_background", "on")
    cmd.set("antialias", 2)
    
    cmd.png(output_name, width=1200, height=800, dpi=300, ray=1)
    if verbose:
        print(f"Saved: {output_name}")

def run_pipeline(args):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(script_dir, "targets_list.csv")
    pdb_dir = os.path.join(script_dir, "filtered_pdb")
    af_dir = os.path.join(script_dir, "cut_alpha_fold")
    output_dir = os.path.join(script_dir, "3D_aligned_alpha_fold")
    out_png_dir = os.path.join(script_dir, "3D_alignment_visualisations")
    
    # creating output directories if they don't exist
    os.makedirs(output_dir, exist_ok=True)
    if args.visualize:
        os.makedirs(out_png_dir, exist_ok=True)

    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"targets_list.csv not found at: {csv_path}")

    df = pd.read_csv(csv_path)
    results = []
    success_count = 0

    print(f"Beginning alignment of {len(df)} proteins...")

    for index, row in df.iterrows():
        pdb_id = row['PDB_ID']
        af_id = row['AF_ID']
        
        pdb_file = os.path.join(pdb_dir, f"{pdb_id}.pdb")
        af_file = os.path.join(af_dir, f"{af_id}.pdb")
        out_file = os.path.join(output_dir, f"{af_id}.pdb")

        if not os.path.exists(pdb_file):
            print(f"Warning: missing PDB file {pdb_file}; skipping {pdb_id}.")
            results.append({
                'PDB_ID': pdb_id,
                'AF_ID': af_id,
                'ref_ca_count': 0,
                'sample_ca_count': 0,
                'matched_ca_count': 0,
                'method': None,
                'RMSD': None,
            })
            continue

        if not os.path.exists(af_file):
            if args.loud:
                print(f"Warning: missing AlphaFold file {af_file}; skipping {af_id}.")
            results.append({
                'PDB_ID': pdb_id,
                'AF_ID': af_id,
                'ref_ca_count': 0,
                'sample_ca_count': 0,
                'matched_ca_count': 0,
                'method': None,
                'RMSD': None,
            })
            continue

        if args.loud:
            print(f"[{index+1}/{len(df)}] {pdb_id} vs {af_id}...")

        rmsd, method, ref_count, sample_count, matched_count = superimpose_robust(
            pdb_file,
            af_file,
            out_file,
            direct_numbering=args.direct_numbering,
            verbose=args.loud,
        )

        results.append({
            'PDB_ID': pdb_id,
            'AF_ID': af_id,
            'ref_ca_count': ref_count,
            'sample_ca_count': sample_count,
            'matched_ca_count': matched_count,
            'method': method,
            'RMSD': rmsd,
        })

        if rmsd is not None:
            success_count += 1

        status = f"RMSD: {rmsd:.2f}" if rmsd is not None else "Error"
        if args.loud:
            print(f"[{index+1}/{len(df)}] {pdb_id} vs {af_id} -> {status}")

        if rmsd is not None and args.visualize:
            if HAVE_PYMOL:
                try:
                    visualize_alignment(pdb_file, out_file, output_name=os.path.join(out_png_dir, f"{pdb_id}_vs_{af_id}.png"), verbose=args.loud)
                except Exception as err:
                    print(f"Warning: visualization failed for {pdb_id} vs {af_id}: {err}")
                    print("Continuing without visualization.")
            else:
                print("PyMOL is not available: skipping visualization.")

    # save results to CSV
    results_df = pd.DataFrame(results)
    results_df.to_csv("rmsd_results.csv", index=False)

    count_below_1 = results_df[results_df['RMSD'] <= 1.0].shape[0]
    count_1_2 = results_df[(results_df['RMSD'] > 1.0) & (results_df['RMSD'] <= 2.0)].shape[0]
    count_2_5 = results_df[(results_df['RMSD'] > 2.0) & (results_df['RMSD'] <= 5.0)].shape[0]

    print("\n" + "="*30)
    print("3D alignment summary")
    print(f"Successful alignments: {success_count} / {len(df)}")
    print(f"Number of pairs with RMSD <= 1.0: {count_below_1}")
    print(f"Number of pairs with 1.0 < RMSD <= 2.0: {count_1_2}")
    print(f"Number of pairs with 2.0 < RMSD <= 5.0: {count_2_5}")
    print(f"Rotated AlphaFold models saved to: {output_dir}")
    print(f"Summary CSV: rmsd_results.csv")
    print("="*30)
    print()

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Align trimmed AlphaFold models to their experimental PDB structures."
    )
    parser.add_argument(
        "--visualize",
        action="store_true",
        help="Generate optional PyMOL visualizations of aligned structure pairs.",
    )
    parser.add_argument(
        "--loud",
        action="store_true",
        help="Print detailed progress and alignment diagnostics.",
    )
    parser.add_argument(
        "--direct-numbering",
        action="store_true",
        help="Attempt direct residue-number-based C-alpha matching when possible.",
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    run_pipeline(args)
