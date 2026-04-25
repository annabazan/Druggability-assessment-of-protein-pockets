import sys
from Bio.PDB import PDBParser, PDBIO, Superimposer
import os
import pandas as pd 
import pymol
from pymol import cmd

def superimpose_robust(pdb_path, af_path, output_path):
    parser = PDBParser(QUIET=True)
    ref_struct = parser.get_structure("ref", pdb_path)
    sample_struct = parser.get_structure("sample", af_path)

    ref_ca = {res.get_id()[1]: res['CA'] for res in ref_struct.get_residues() if 'CA' in res}
    sample_ca = {res.get_id()[1]: res['CA'] for res in sample_struct.get_residues() if 'CA' in res}
    print(f"Reference residues: {len(ref_ca)}, Sample residues: {len(sample_ca)}")
    common_res_ids = sorted(set(ref_ca.keys()) & set(sample_ca.keys()))

    if not common_res_ids:
        return None  # None if no common residues found

    ref_atoms = [ref_ca[res_id] for res_id in common_res_ids]
    sample_atoms = [sample_ca[res_id] for res_id in common_res_ids]

    si = Superimposer()
    si.set_atoms(ref_atoms, sample_atoms)
    si.apply(sample_struct.get_atoms())

    io = PDBIO()
    io.set_structure(sample_struct)
    io.save(output_path)
    
    return si.rms # Return the RMSD value

def visualize_alignment(pdb_file, alphafold_file, output_name="comparison.png"):
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
    print(f"Saved: {output_name}")


def run_pipeline():
    csv_path = "../targets_list.csv"
    pdb_dir = "../filtered_pdb"
    af_dir = "../cut_alpha_fold"
    output_dir = "./3D_aligned_alpha_fold"
    out_png_dir = "./3D_alignment_visualisations"
    
    # creating output directories if they don't exist
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(out_png_dir, exist_ok=True)

    df = pd.read_csv(csv_path)
    results = []

    print(f"Beginning alignment of {len(df)} proteins...\n")

    for index, row in df.iterrows():
        pdb_id = row['PDB_ID']
        af_id = row['AF_ID']
        
        pdb_file = os.path.join(pdb_dir, f"{pdb_id}.pdb")
        af_file = os.path.join(af_dir, f"{af_id}.pdb")
        out_file = os.path.join(output_dir, f"{af_id}.pdb")

        rmsd = superimpose_robust(pdb_file, af_file, out_file)
        
        results.append({
            'PDB_ID': pdb_id,
            'AF_ID': af_id,
            'RMSD': rmsd
        })
        
        status = f"RMSD: {rmsd:.2f}" if rmsd is not None else "Error"
        print(f"[{index+1}/{len(df)}] {pdb_id} vs {af_id} -> {status}")

        if rmsd is not None:
            visualize_alignment(pdb_file, out_file, output_name=os.path.join(out_png_dir, f"{pdb_id}_vs_{af_id}.png"))

    # save results to CSV
    results_df = pd.DataFrame(results)
    results_df.to_csv("rmsd_results.csv", index=False)

    count_below_1 = results_df[results_df['RMSD'] <= 1.0].shape[0]
    count_below_2 = results_df[results_df['RMSD'] <= 2.0].shape[0]
    
    print("\n" + "="*30)
    print(f"Results saved to: rmsd_results.csv")
    print(f"Number of pairs with RMSD <= 1.0: {count_below_1}")
    print(f"Number of pairs with RMSD <= 2.0: {count_below_2}")
    print("="*30)

if __name__ == "__main__":
    run_pipeline()