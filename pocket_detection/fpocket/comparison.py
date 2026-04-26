import os
import pandas as pd
import glob
import numpy as np
from Bio.PDB import PDBParser
import pymol
from pymol import cmd

def get_pocket_residues(pdb_file):
    """Gets the set of residue numbers that are part of the pocket from a PDB file."""
    residues = set()
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    res_num = int(line[22:26].strip()) 
                    residues.add(res_num)
    except (FileNotFoundError, ValueError):
        return None
    return residues

def calculate_jaccard(set1, set2):
    if not set1 or not set2: return 0.0
    return len(set1 & set2) / len(set1 | set2)

def calculate_rmsd(struct_pdb, struct_af, common_res_ids):
    """Calculates RMSD between two structures based on common residue IDs."""
    coords_pdb = []
    coords_af = []
    
    # we assume both structures have the same chain ID (e.g., 'A') and that residue numbering is consistent
    try:
        chain_pdb = struct_pdb[0]['A']
        chain_af = struct_af[0]['A']
    except KeyError:
        return None

    for res_id in common_res_ids:
        try:
            # coordinates of CA atoms for the common residues
            p_coord = chain_pdb[res_id]['CA'].get_coord()
            a_coord = chain_af[res_id]['CA'].get_coord()
            coords_pdb.append(p_coord)
            coords_af.append(a_coord)
        except KeyError:
            continue

    if len(coords_pdb) < 3:  # RMSD is not meaningful with fewer than 3 points
        return None

    p = np.array(coords_pdb)
    q = np.array(coords_af)
    
    # RMSD calculation (without superimposition, just direct distance)
    diff = p - q
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))

def compare_pockets(pdb_out_dir, af_out_dir, full_pdb_path, full_af_path):
    # structure loading
    parser = PDBParser(QUIET=True)
    try:
        struct_pdb = parser.get_structure('PDB', full_pdb_path)
        struct_af = parser.get_structure('AF', full_af_path)
    except Exception as e:
        print(f"Error loading structures: {e}")
        return pd.DataFrame()

    pdb_path = os.path.join(pdb_out_dir, "pockets")
    af_path = os.path.join(af_out_dir, "pockets")
    
    pdb_files = sorted(glob.glob(os.path.join(pdb_path, "pocket*_atm.pdb")))
    af_files = sorted(glob.glob(os.path.join(af_path, "pocket*_atm.pdb")))
    
    pdb_data = {os.path.basename(f): get_pocket_residues(f) for f in pdb_files}
    af_data = {os.path.basename(f): get_pocket_residues(f) for f in af_files}

    results = []

    for p_name, p_res in pdb_data.items():
        best_jaccard = 0
        best_match = "None"
        
        for a_name, a_res in af_data.items():
            score = calculate_jaccard(p_res, a_res)
            if score > best_jaccard:
                best_jaccard = score
                best_match = a_name
        
        if best_jaccard > 0:
            common_res = p_res & af_data[best_match]
            rmsd_val = calculate_rmsd(struct_pdb, struct_af, common_res)
            
            results.append({
                "PDB_Pocket": p_name.replace("_atm.pdb", ""),
                "AF_Best_Match": best_match.replace("_atm.pdb", ""),
                "Jaccard_Index": round(best_jaccard, 3),
                "Shared_Residues": len(common_res),
                "RMSD": round(rmsd_val, 3) if rmsd_val is not None else "N/A"
            })
    
    return pd.DataFrame(results)


def visualize_pockets(pdb_file, af_file, common_res, output_png):
    """Generates a PNG visualization of the pocket comparison using PyMOL."""
    cmd.reinitialize()
    
    cmd.load(pdb_file, "PDB_struct")
    cmd.load(af_file, "AF_struct")
    
    cmd.hide("everything")
    cmd.show("cartoon", "all")
    cmd.set("cartoon_transparency", 0.7)
    cmd.color("gray80", "all")

    # selecting pocket residues for both structures
    res_str = "+".join(map(str, common_res))
    cmd.select("pocket_pdb", f"PDB_struct and resi {res_str}")
    cmd.select("pocket_af", f"AF_struct and resi {res_str}")

    cmd.show("sticks", "pocket_pdb")
    cmd.show("sticks", "pocket_af")
    cmd.color("marine", "pocket_pdb")
    cmd.color("orange", "pocket_af")

    
    # surfaces and transparency for better visualization
    cmd.show("surface", "pocket_pdb")
    cmd.set("transparency", 0.5, "pocket_pdb")
    cmd.show("surface", "pocket_af")
    cmd.set("transparency", 0.5, "pocket_af")
    
    cmd.zoom("pocket_pdb", buffer=5)
    cmd.png(output_png, width=1200, height=900, ray=1)
    print(f"Generated visualization: {output_png}")

def visualize_pockets_side_by_side(pdb_file, af_file, common_res, output_png):
    """Generates a PNG visualization of the pocket comparison side-by-side using translation."""
    cmd.reinitialize()
    
    # Loading structures
    cmd.load(pdb_file, "PDB_struct")
    cmd.load(af_file, "AF_struct")
    
    # Translating AF structure to the right for side-by-side view
    cmd.translate([50, 0, 0], "AF_struct")

    cmd.hide("everything")
    cmd.show("cartoon", "all")
    cmd.set("cartoon_transparency", 0.7)
    cmd.color("gray80", "all")

    res_str = "+".join(map(str, common_res))
    cmd.select("pocket_pdb", f"PDB_struct and resi {res_str}")
    cmd.select("pocket_af", f"AF_struct and resi {res_str}")

    cmd.show("sticks", "pocket_pdb")
    cmd.show("sticks", "pocket_af")
    cmd.color("marine", "pocket_pdb")
    cmd.color("orange", "pocket_af")
    
    cmd.show("surface", "pocket_pdb")
    cmd.show("surface", "pocket_af")
    cmd.set("transparency", 0.5, "pocket_pdb")
    cmd.set("transparency", 0.5, "pocket_af")

    # Centering the view on both structures
    cmd.zoom("all", buffer=1) 
    
    cmd.png(output_png, width=1600, height=900, ray=1)
    print(f"Generated side-by-side visualization: {output_png}")

def compare_pockets(pdb_out_dir, af_out_dir, full_pdb_path, full_af_path, viz_dir, pair_id):
    parser = PDBParser(QUIET=True)
    try:
        struct_pdb = parser.get_structure('PDB', full_pdb_path)
        struct_af = parser.get_structure('AF', full_af_path)
    except Exception: return pd.DataFrame()

    pdb_path = os.path.join(pdb_out_dir, "pockets")
    af_path = os.path.join(af_out_dir, "pockets")
    
    pdb_files = sorted(glob.glob(os.path.join(pdb_path, "pocket*_atm.pdb")))
    af_files = sorted(glob.glob(os.path.join(af_path, "pocket*_atm.pdb")))
    
    pdb_data = {os.path.basename(f): get_pocket_residues(f) for f in pdb_files}
    af_data = {os.path.basename(f): get_pocket_residues(f) for f in af_files}

    results = []

    for p_name, p_res in pdb_data.items():
        best_jaccard, best_match = 0, "None"
        for a_name, a_res in af_data.items():
            score = calculate_jaccard(p_res, a_res)
            if score > best_jaccard:
                best_jaccard, best_match = score, a_name
        
        if best_jaccard > 0:
            common_res = p_res & af_data[best_match]
            rmsd_val = calculate_rmsd(struct_pdb, struct_af, common_res)
            
            p_short = p_name.replace("_atm.pdb", "")
            a_short = best_match.replace("_atm.pdb", "")
            
            # generating visualizations for the best match
            viz_name = f"{pair_id}_{p_short}_vs_{a_short}.png"
            visualize_pockets(full_pdb_path, full_af_path, common_res, os.path.join(viz_dir, viz_name))
            visualize_pockets_side_by_side(full_pdb_path, full_af_path, common_res, os.path.join(viz_dir, f"side_by_side_{viz_name}"))
            
            results.append({
                "PDB_Pocket": p_short,
                "AF_Best_Match": a_short,
                "Jaccard_Index": round(best_jaccard, 3),
                "RMSD": round(rmsd_val, 3) if rmsd_val is not None else "N/A",
                "Visualization": viz_name
            })
    
    return pd.DataFrame(results)


def run_pipeline():
    csv_path = "../../targets/targets_list.csv"
    
    pdb_raw_dir = "../../targets/filtered_pdb" 
    af_raw_dir = "../../targets/3D_alignment/3D_aligned_alpha_fold"

    pdb_out_root = "pdb_out"
    af_out_root = "alpha_fold_out"
    output_dir = "comparison_results"
    out_png_dir = "pocket_visualisations"
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(out_png_dir, exist_ok=True)
    df = pd.read_csv(csv_path)

    for index, row in df.iterrows():
        pdb_id = row['PDB_ID']
        af_id = row['AF_ID']
        

        pdb_pocket_dir = os.path.join(pdb_out_root, f"{pdb_id}_out")
        af_pocket_dir = os.path.join(af_out_root, f"{af_id}_out")
        
        full_pdb = os.path.join(pdb_raw_dir, f"{pdb_id}.pdb")
        full_af = os.path.join(af_raw_dir, f"{af_id}.pdb")

        out_file = os.path.join(output_dir, f"{pdb_id}_vs_{af_id}.csv")

        if os.path.exists(pdb_pocket_dir) and os.path.exists(af_pocket_dir):
            df_comparison = compare_pockets(pdb_pocket_dir, af_pocket_dir, full_pdb, full_af, out_png_dir, f"{pdb_id}_{af_id}")
            df_comparison.to_csv(out_file, index=False)
            print(f"Saved comparison with RMSD for {pdb_id} vs {af_id}")
        else:
            print(f"Skipping {pdb_id} - missing fpocket directories")

if __name__ == "__main__":
    run_pipeline()