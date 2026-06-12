import pandas as pd
from pathlib import Path
import py3Dmol

def get_pocket_files(pockets_dir):
    pockets_dir = Path(pockets_dir)
    pocket_files = {}

    for vert_file in sorted(pockets_dir.glob("pocket*_vert.pqr")):
        pocket_id = int(
            vert_file.stem.replace("pocket", "").replace("_vert", "")
        )
        pocket_files[pocket_id] = {
            "vert": vert_file,
            "atm": pockets_dir / f"pocket{pocket_id}_atm.pdb",
        }
    return pocket_files

def get_pdb_files(
        PDB_id, AF_id,
        PDB_pockets = "pdb_out",
        AF_pockets = "alpha_fold_out",
        PDB_dir = "filtered_pdb",
        AF_dir = "3D_aligned_alpha_fold"
        ):
    structure_pdb = Path("targets")/PDB_dir/f"{PDB_id}.pdb"
    structure_af = Path("targets")/AF_dir/f"{AF_id}.pdb"

    fpocket_dir_pdb = Path("pocket_detection/fpocket")/PDB_pockets/f"{PDB_id}_out"
    fpocket_dir_af = Path("pocket_detection/fpocket")/AF_pockets/f"{AF_id}_out"

    out_pdb = fpocket_dir_pdb/f"{PDB_id}_out.pdb"
    out_af = fpocket_dir_af/f"{AF_id}_out.pdb"

    pockets_dir_pdb = fpocket_dir_pdb/"pockets"
    pockets_dir_af = fpocket_dir_af/"pockets"

    pocket_files_pdb = get_pocket_files(pockets_dir_pdb)
    pocket_files_af = get_pocket_files(pockets_dir_af)

    files = {
        "pdb" : {
            "structure" : structure_pdb,
            "out" : out_pdb,
            "pocket_files" : pocket_files_pdb
        },
        "af" : {
            "structure" : structure_af,
            "out" : out_af,
            "pocket_files" : pocket_files_af
        },
    }
    return files    

def add_pocket(
    view, viewer,
    structure_pdb,
    pocket_files,
    pocket_id
):
    pocket_vert = pocket_files[pocket_id]["vert"]
    pocket_atm = pocket_files[pocket_id]["atm"]

    with open(structure_pdb) as f:
        view.addModel(f.read(), "pdb", viewer=viewer)
    with open(pocket_vert) as f:
        view.addModel(f.read(), "pqr", viewer=viewer)
    with open(pocket_atm) as f:
        view.addModel(f.read(), "pdb", viewer=viewer)

    # proteine structure (cartoon)
    view.setStyle(
        {"model": 0},
        {"cartoon": {"color": "white"}},
        viewer=viewer
    )
    # proteine surface
    view.addSurface(
        py3Dmol.MS,
        {
            "color": "white",
            "opacity": 0.7
        },
        {"model": 0},
        viewer=viewer
    )

    # pocket spheres
    view.setStyle(
        {"model": 1},
        {"sphere": {"radius": 1.0, "color": "red"}},
        viewer=viewer
    )
    # residue sticks
    view.setStyle(
        {"model": 2},
        {"stick": {"colorscheme": "amino"}},
        viewer=viewer
    )
    # pocket surface
    view.addSurface(
        py3Dmol.VDW,
        {
            "color": "red",
            "opacity": 0.7
        },
        {"model": 1},
        viewer=viewer
    )

def show_pocket(
    structure_pdb,
    pocket_files,
    pocket_id,
    width=500,
    height=400
):
    view = py3Dmol.view(
        width=width,
        height=height
    )

    add_pocket(
        view=view,
        viewer=None,
        structure_pdb=structure_pdb,
        pocket_files=pocket_files,
        pocket_id=pocket_id
    )

    view.setBackgroundColor("black")
    view.zoomTo({"model": 1})

    return view

def show_pocket_comparison(
    left_structure_pdb,
    left_pocket_files,
    left_pocket_id,
    right_structure_pdb,
    right_pocket_files,
    right_pocket_id,
    width=500,
    height=400
):
    view = py3Dmol.view(
        width=2*width,
        height=height,
        viewergrid=(1, 2)
    )

    add_pocket(
        view=view,
        viewer=(0, 0),
        structure_pdb=left_structure_pdb,
        pocket_files=left_pocket_files,
        pocket_id=left_pocket_id
    )

    add_pocket(
        view=view,
        viewer=(0, 1),
        structure_pdb=right_structure_pdb,
        pocket_files=right_pocket_files,
        pocket_id=right_pocket_id
    )

    view.setBackgroundColor("black")
    view.zoomTo()

    return view

def add_fpocket_structure(
    view, viewer,
    structure_pdb,
    out_pdb,
    pocket_files,
    pocket_ids=None
):
    colors = [
        "red", "blue", "green", "orange", "purple",
        "cyan", "magenta", "yellow", "lime", "pink"
    ]

    with open(structure_pdb) as f:
        view.addModel(f.read(), "pdb", viewer=viewer)
    with open(out_pdb) as f:
        view.addModel(f.read(), "pdb", viewer=viewer)

    # proteine structure (cartoon)
    view.setStyle(
        {"model": 0},
        {"cartoon": {"color": "white"}},
        viewer=viewer
    )
    # protein surface
    view.addSurface(
        py3Dmol.MS,
        {
            "color": "white",
            "opacity": 0.7
        },
        {"model": 0},
        viewer=viewer
    )
    view.setStyle(
        {"model": 1},
        {},
        viewer=viewer
    )

    if pocket_ids is None:
        pocket_ids = sorted(pocket_files.keys())
    
    for i, pocket_id in enumerate(pocket_ids):
        view.setStyle(
            {
                "model": 1,
                "resn": "STP",
                "resi": str(pocket_id)
            },
            {
                "sphere": {
                    "radius": 0.8,
                    "color": colors[i % len(colors)]
                }
            },
            viewer=viewer
        )

def show_fpocket_structure(
    structure_pdb,
    out_pdb,
    pocket_files,
    pocket_ids=None,
    width=500,
    height=400,
):
    view = py3Dmol.view(width=width, height=height)
    add_fpocket_structure(
        view, viewer=None, 
        structure_pdb=structure_pdb, 
        out_pdb=out_pdb, 
        pocket_files=pocket_files, 
        pocket_ids=pocket_ids
        )
    
    view.setBackgroundColor("black")
    view.zoomTo()
    
    return view

def show_fpocket_comparison(
    left_structure, left_out, left_pockets, 
    right_structure, right_out, right_pockets,
    left_pocket_ids=None, right_pocket_ids=None,
    width=500, height=400,
):
    view = py3Dmol.view(
        width=2*width,
        height=height,
        viewergrid=(1, 2)
    )
    add_fpocket_structure(
        view, viewer=(0,0), 
        structure_pdb=left_structure, 
        out_pdb = left_out,
        pocket_files = left_pockets, 
        pocket_ids = left_pocket_ids
    )
    add_fpocket_structure(
        view, viewer=(0,1), 
        structure_pdb=right_structure, 
        out_pdb = right_out,
        pocket_files = right_pockets, 
        pocket_ids = right_pocket_ids
    )
    view.setBackgroundColor("black")
    view.zoomTo()
    
    return view

def filter_pockets_by_jaccard(pockets_df, jaccard_threshold=0.75):
    filtered_df = pockets_df.loc[
        pockets_df["jaccard"] >= jaccard_threshold
    ].copy()

    pdb_pockets = (
        filtered_df["pdb_pocket_num"]
        .dropna()
        .astype(int)
        .tolist()
    )

    af_pockets = (
        filtered_df["af_pocket_num"]
        .dropna()
        .astype(int)
        .tolist()
    )

    return filtered_df, pdb_pockets, af_pockets

def best_matched(PDB_id, AF_id, jaccard_threshold=0.75):
    ANALYSIS_DIR = Path("analysis/pocket_comparison/outputs/local")/f"{PDB_id}_vs_{AF_id}"
    pockets_info = ANALYSIS_DIR/"pocket_pairs_detailed.csv"
    df = pd.read_csv(pockets_info)

    df, pdb_pockets, af_pockets = filter_pockets_by_jaccard(df, jaccard_threshold)

    report_df = df[
        [
            "pdb_pocket_num",
            "af_pocket_num",
            "jaccard",
            "rmsd_common_ca_local_aligned",
            "af_pocket_mean_plddt",
            "pdb_fpocket_pocket_score",
            "af_fpocket_pocket_score",
            "pdb_fpocket_drug_score",
            "af_fpocket_drug_score",
            "pdb_fpocket_volume",
            "af_fpocket_volume",
        ]
    ].copy()

    report_df = report_df.rename(
        columns={
            "pdb_pocket_num": "PDB",
            "af_pocket_num": "AF",
            "jaccard": "Jaccard",
            "rmsd_common_ca_local_aligned": "RMSD",
            "af_pocket_mean_plddt": "pLDDT",
            "pdb_fpocket_pocket_score": "PDB score",
            "af_fpocket_pocket_score": "AF score",
            "pdb_fpocket_drug_score": "PDB drug",
            "af_fpocket_drug_score": "AF drug",
            "pdb_fpocket_volume": "PDB vol",
            "af_fpocket_volume": "AF vol",
        }
    )

    report_df["PDB"] = report_df["PDB"].astype("Int64")
    report_df["AF"] = report_df["AF"].astype("Int64")

    report_df["Jaccard"] = report_df["Jaccard"].round(2)
    report_df["RMSD"] = report_df["RMSD"].round(2)
    report_df["pLDDT"] = report_df["pLDDT"].round(2)

    report_df["PDB score"] = report_df["PDB score"].round(2)
    report_df["AF score"] = report_df["AF score"].round(2)

    report_df["PDB drug"] = report_df["PDB drug"].round(2)
    report_df["AF drug"] = report_df["AF drug"].round(2)

    report_df["PDB vol"] = report_df["PDB vol"].round(1)
    report_df["AF vol"] = report_df["AF vol"].round(1)

    return report_df, pdb_pockets, af_pockets

def prepare_df(PDB_id, AF_id):
    ANALYSIS_DIR = Path("analysis/pocket_comparison/outputs/local")/f"{PDB_id}_vs_{AF_id}"
    pockets_info = ANALYSIS_DIR/"pocket_pairs_detailed.csv"
    df = pd.read_csv(pockets_info)

    df = df[
        [
        "pdb_pocket_num",
        "af_pocket_num",
        "pair_status",

        "jaccard",
        "shared_residues",
        "pdb_residue_count",
        "af_residue_count",
        "rmsd_common_ca_local_aligned",
        "af_pocket_mean_plddt",

        "pdb_fpocket_pocket_score",
        "af_fpocket_pocket_score",

        "pdb_fpocket_drug_score",
        "af_fpocket_drug_score",

        "pdb_fpocket_volume",
        "af_fpocket_volume",

        "pdb_fpocket_hydrophobicity_score",
        "af_fpocket_hydrophobicity_score",

        "pdb_fpocket_polarity_score",
        "af_fpocket_polarity_score",

        "pdb_fpocket_charge_score",
        "af_fpocket_charge_score",

        "pdb_fpocket_flexibility",
        "af_fpocket_flexibility",
        ]
    ].copy()

    df = df.rename(
        columns={
            "pdb_pocket_num": "PDB",
            "af_pocket_num": "AF",

            "jaccard": "Jaccard",

            "shared_residues": "Shared",
            "pdb_residue_count": "PDB res",
            "af_residue_count": "AF res",

            "rmsd_common_ca_local_aligned": "RMSD",
            "af_pocket_mean_plddt": "pLDDT",

            "pdb_fpocket_pocket_score": "PDB score",
            "af_fpocket_pocket_score": "AF score",

            "pdb_fpocket_drug_score": "PDB drug",
            "af_fpocket_drug_score": "AF drug",

            "pdb_fpocket_volume": "PDB vol",
            "af_fpocket_volume": "AF vol",

            "pdb_fpocket_hydrophobicity_score": "PDB hydro",
            "af_fpocket_hydrophobicity_score": "AF hydro",

            "pdb_fpocket_polarity_score": "PDB polar",
            "af_fpocket_polarity_score": "AF polar",

            "pdb_fpocket_charge_score": "PDB charge",
            "af_fpocket_charge_score": "AF charge",

            "pdb_fpocket_flexibility": "PDB flex",
            "af_fpocket_flexibility": "AF flex",
        }
    )

    df["PDB"] = df["PDB"].astype("Int64")
    df["AF"] = df["AF"].astype("Int64")

    df["Shared"] = df["Shared"].astype("Int64")
    df["PDB res"] = df["PDB res"].astype("Int64")
    df["AF res"] = df["AF res"].astype("Int64")

    df["Jaccard"] = df["Jaccard"].round(2)

    df["RMSD"] = df["RMSD"].round(2)
    df["pLDDT"] = df["pLDDT"].round(2)

    df["PDB score"] = df["PDB score"].round(2)
    df["AF score"] = df["AF score"].round(2)

    df["PDB drug"] = df["PDB drug"].round(2)
    df["AF drug"] = df["AF drug"].round(2)

    df["PDB vol"] = df["PDB vol"].round(1)
    df["AF vol"] = df["AF vol"].round(1)

    df["PDB hydro"] = df["PDB hydro"].round(2)
    df["AF hydro"] = df["AF hydro"].round(2)

    df["PDB polar"] = df["PDB polar"].round(2)
    df["AF polar"] = df["AF polar"].round(2)

    df["PDB charge"] = df["PDB charge"].round(2)
    df["AF charge"] = df["AF charge"].round(2)

    df["PDB flex"] = df["PDB flex"].round(2)
    df["AF flex"] = df["AF flex"].round(2)

    return df

def prepare_status_summary(df):
    status_order = [
        "matched",
        "weak_match",
        "pdb_only",
        "af_only"
    ]

    status_summary = (
        df.assign(
            pair_status=pd.Categorical(
                df["pair_status"],
                categories=status_order,
                ordered=True
            )
        )
        ["pair_status"]
        .value_counts(sort=False)
        .rename_axis("Status")
        .reset_index(name="Count")
    )

    return status_summary

def prepare_score_summary(df):
    score_summary = (
        df[
            [
                "PDB score",
                "AF score",
                "PDB drug",
                "AF drug"
            ]
        ]
        .agg(["mean", "median", "max"])
        .round(2)
    )

    return score_summary
