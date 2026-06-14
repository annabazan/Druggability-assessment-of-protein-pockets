import pandas as pd
from pathlib import Path
import py3Dmol
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

# analysis utils

def prepare_pockets_df(df):
    rows = []
    for _, row in df.iterrows():
        status = row["pair_status"]
        # PDB pockets
        if status in ["matched", "weak_match", "pdb_only"]:
            rows.append({
                "Origin": "PDB",
                "ID": row["pdb_id"],
                "Pocket ID": int(row["pdb_pocket_num"]),
                "Status": status,
                "Match": (
                    int(row["af_pocket_num"])
                    if status in ["matched", "weak_match"]
                    else "-"
                ),
                "Jaccard": (
                    row["jaccard"] 
                    if status in ["matched", "weak_match"]
                    else "-"                
                ),
                "Res": row["pdb_residue_count"],
                "Mean pLDDT": 0,
                "Median pLDDT": 0,
                "Score": row["pdb_fpocket_pocket_score"],
                "Drug": row["pdb_fpocket_drug_score"],
                "Vol": row["pdb_fpocket_volume"],
                "Hydro": row["pdb_fpocket_hydrophobicity_score"],
                "Polar": row["pdb_fpocket_polarity_score"],
                "Charge": row["pdb_fpocket_charge_score"],
                "Flex": row["pdb_fpocket_flexibility"],
            })
        # AF pockets
        if status in ["matched", "weak_match", "af_only"]:
            rows.append({
                "Origin": "AF",
                "ID": row["af_id"],
                "Pocket ID": int(row["af_pocket_num"]),
                "Status": status,
                "Match": (
                    int(row["pdb_pocket_num"])
                    if status in ["matched", "weak_match"]
                    else "-"
                ),
                "Jaccard": (
                    row["jaccard"] 
                    if status in ["matched", "weak_match"]
                    else "-"                
                ),
                "Res": row["af_residue_count"],
                "Mean pLDDT": row["af_pocket_mean_plddt"],
                "Median pLDDT": row["af_pocket_median_plddt"],
                "Score": row["af_fpocket_pocket_score"],
                "Drug": row["af_fpocket_drug_score"],
                "Vol": row["af_fpocket_volume"],
                "Hydro": row["af_fpocket_hydrophobicity_score"],
                "Polar": row["af_fpocket_polarity_score"],
                "Charge": row["af_fpocket_charge_score"],
                "Flex": row["af_fpocket_flexibility"],
            })

    pockets_df = pd.DataFrame(rows)

    pockets_pdb = (
        pockets_df[pockets_df["Origin"] == "PDB"]
        .drop(columns="Origin")
        .reset_index(drop=True)
    ).copy()
    pockets_af = (
        pockets_df[pockets_df["Origin"] == "AF"]
        .drop(columns="Origin")
        .reset_index(drop=True)
    ).copy()

    return pockets_df, pockets_pdb, pockets_af

def viz_status(pockets_pdb, pockets_af):
    status_order = ["matched", "weak_match", "pdb_only", "af_only"]
    label_map = {
        "matched": "Matched",
        "weak_match": "Weak match",
        "pdb_only": "PDB only",
        "af_only": "AF only"
    }

    pdb_counts = (
        pockets_pdb["Status"]
        .value_counts()
        .reindex(status_order, fill_value=0)
    )
    af_counts = (
        pockets_af["Status"]
        .value_counts()
        .reindex(status_order, fill_value=0)
    )
    pdb_counts = pdb_counts[pdb_counts > 0]
    af_counts = af_counts[af_counts > 0]    

    _, axes = plt.subplots(1, 2, figsize=(8, 4))

    axes[0].pie(
        pdb_counts,
        labels=[label_map[status] for status in pdb_counts.index],
        autopct="%1.1f%%",
        startangle=90,
        colors=["tab:blue", "tab:green", "tab:red"]
    )
    axes[1].pie(
        af_counts,
        labels=[label_map[status] for status in af_counts.index],
        autopct="%1.1f%%",
        startangle=90,
        colors=["tab:blue", "tab:green", "tab:orange"]
    )

    axes[0].set_title("PDB pockets")
    axes[1].set_title("AF pockets")
    plt.tight_layout()
    plt.show()

def show_score_stats(pockets_pdb, pockets_af):
    score_stats = pd.concat(
        [
            pockets_pdb[["Score", "Drug"]]
            .agg(["mean", "median", "std", "min", "max"])
            .T
            .assign(Origin="PDB"),

            pockets_af[["Score", "Drug"]]
            .agg(["mean", "median", "std", "min", "max"])
            .T
            .assign(Origin="AF"),
        ]
    )

    score_stats = (
        score_stats
        .reset_index(names="Metric")
        .loc[:, ["Metric", "Origin", "mean", "median", "max", "min", "std"]]
    )
    score_stats["Metric"] = pd.Categorical(
        score_stats["Metric"],
        categories=["Score", "Drug"],
        ordered=True
    )
    score_stats["Origin"] = pd.Categorical(
        score_stats["Origin"],
        categories=["PDB", "AF"],
        ordered=True
    )
    score_stats = (
        score_stats
        .sort_values(["Metric", "Origin"])
        .reset_index(drop=True)
    )

    numeric_cols = ["mean", "median", "max", "min", "std"]
    score_stats[numeric_cols] = score_stats[numeric_cols].round(3)

    return score_stats

def show_characteristics_stats(pockets_pdb, pockets_af):
    score_stats = pd.concat(
        [
            pockets_pdb[["Vol", "Hydro", "Polar", "Charge", "Flex"]]
            .agg(["mean", "median", "std", "min", "max"])
            .T
            .assign(Origin="PDB"),

            pockets_af[["Vol", "Hydro", "Polar", "Charge"]]
            .agg(["mean", "median", "std", "min", "max"])
            .T
            .assign(Origin="AF"),
        ]
    )

    score_stats = (
        score_stats
        .reset_index(names="Metric")
        .loc[:, ["Metric", "Origin", "mean", "median", "max", "min", "std"]]
    )
    score_stats["Metric"] = pd.Categorical(
        score_stats["Metric"],
        categories=["Vol", "Hydro", "Polar", "Charge", "Flex"],
        ordered=True
    )
    score_stats["Origin"] = pd.Categorical(
        score_stats["Origin"],
        categories=["PDB", "AF"],
        ordered=True
    )
    score_stats = (
        score_stats
        .sort_values(["Metric", "Origin"])
        .reset_index(drop=True)
    )

    numeric_cols = ["mean", "median", "max", "min", "std"]
    score_stats[numeric_cols] = score_stats[numeric_cols].round(1)

    return score_stats

def score_filter(pockets_df, min_score=0.1, min_drug=0.3):
    return pockets_df[
        (pockets_df["Score"] >= min_score)
        & (pockets_df["Drug"] >= min_drug)
    ].copy()

def show_score_quantiles(pockets_pdb, pockets_af, quantiles=None):
    if quantiles is None:
        quantiles = [0.25, 0.5, 0.75, 0.8, 0.9, 0.95]

    score_quantiles = pd.concat(
        {
            "Score PDB": pockets_pdb["Score"].quantile(quantiles),
            "Score AF": pockets_af["Score"].quantile(quantiles),
            "Drug PDB": pockets_pdb["Drug"].quantile(quantiles),
            "Drug AF": pockets_af["Drug"].quantile(quantiles),
        },
        axis=1
    ).round(2)

    score_quantiles["Count PDB"] = [
        (
            (pockets_pdb["Score"] >= pockets_pdb["Score"].quantile(q))
            &
            (pockets_pdb["Drug"] >= pockets_pdb["Drug"].quantile(q))
        ).sum()
        for q in quantiles
    ]
    score_quantiles["Count AF"] = [
        (
            (pockets_af["Score"] >= pockets_af["Score"].quantile(q))
            &
            (pockets_af["Drug"] >= pockets_af["Drug"].quantile(q))
        ).sum()
        for q in quantiles
    ]

    return score_quantiles

def plot_distribution(pockets_pdb, pockets_af, bins=60, metric="Score"):
    _, ax = plt.subplots(figsize=(6, 4)) 

    ax.hist(
        pockets_pdb[metric],
        bins=bins,
        alpha=0.6,
        label="PDB",
        density=False
    )
    ax.hist(
        pockets_af[metric],
        bins=bins,
        alpha=0.6,
        label="AF",
        density=False
    )

    pdb_max = pockets_pdb[metric].max()
    af_max = pockets_af[metric].max()
    ax.axvline(
        pdb_max,
        color="tab:blue",
        linestyle="--",
        linewidth=1.5,
        label=f"PDB max ({pdb_max:.2f})"
    )
    ax.axvline(
        af_max,
        color="tab:orange",
        linestyle="--",
        linewidth=1.5,
        label=f"AF max ({af_max:.2f})"
    )

    ax.set_xlabel(f"{metric}")
    ax.set_ylabel("Density")
    ax.set_title(f"Distribution of {metric}")
    ax.legend()

    plt.tight_layout()
    plt.show()

def analyse_matches_pdb(filtered_pdb, filtered_af, pdb_to_af, jaccard_threshold = 0.5):
    pdb_status = (
        filtered_pdb["Status"]
        .value_counts()
        .reindex(["matched", "weak_match", "pdb_only"], fill_value=0)
    )
    pdb_status = pdb_status[pdb_status > 0]

    unique_stuctures = (
        filtered_pdb
        .groupby("Status")["ID"]
        .nunique()
    )

    print(f"Among {len(filtered_pdb)} PDB pockets after filtering there are:")
    for status, count in pdb_status.items():
        print(f"    {count} pockets with status '{status} in {unique_stuctures[status]} protein structures'")

    filtered_pdb["Jaccard"] = pd.to_numeric(
        filtered_pdb["Jaccard"],
        errors="coerce"
    )

    af_keys = set(
        zip(filtered_af["ID"], filtered_af["Pocket ID"])
    )

    pdb_matched = filtered_pdb[
        filtered_pdb["Status"].isin(["matched", "weak_match"])
        & (filtered_pdb["Jaccard"] >= jaccard_threshold)
    ].copy()

    pdb_matched["Partner retained"] = pdb_matched.apply(
        lambda r: (
            (pdb_to_af.get(r["ID"]), r["Match"])
            in af_keys
        ),
        axis=1
    )

    print(f"Number of PDB pockets with Jaccard >= {jaccard_threshold} and retained partner: {pdb_matched['Partner retained'].sum()}")

    return pdb_status, pdb_matched

def analyse_matches_af(filtered_pdb, filtered_af, af_to_pdb, jaccard_threshold = 0.5):
    af_status = (
        filtered_af["Status"]
        .value_counts()
        .reindex(["matched", "weak_match", "af_only"], fill_value=0)
    )
    af_status = af_status[af_status > 0]

    unique_models = (
        filtered_af
        .groupby("Status")["ID"]
        .nunique()
    )

    print(f"Among {len(filtered_af)} AlphaFold pockets after filtering there are:")
    for status, count in af_status.items():
        print(f"    {count} pockets with status '{status} in {unique_models[status]} protein models'")

    filtered_af["Jaccard"] = pd.to_numeric(
        filtered_af["Jaccard"],
        errors="coerce"
    )

    pdb_keys = set(
        zip(filtered_pdb["ID"], filtered_pdb["Pocket ID"])
    )

    af_matched = filtered_af[
        filtered_af["Status"].isin(["matched", "weak_match"])
        & (filtered_af["Jaccard"] >= jaccard_threshold)
    ].copy()

    af_matched["Partner retained"] = af_matched.apply(
        lambda r: (
            (af_to_pdb.get(r["ID"]), r["Match"])
            in pdb_keys
        ),
        axis=1
    )

    print(f"Number of AlphaFold pockets with Jaccard >= {jaccard_threshold} and retained partner: {af_matched['Partner retained'].sum()}")
    
    return af_status, af_matched

def identify_retained_pairs(retained_df, filtered_pdb, af_to_pdb):
    pdb_lookup = filtered_pdb.set_index(["ID", "Pocket ID"])

    rows = []

    for _, row in retained_df.iterrows():

        pdb_id = af_to_pdb[row["ID"]]
        pdb_pocket_id = int(row["Match"])

        try:
            pdb_row = pdb_lookup.loc[(pdb_id, pdb_pocket_id)]
        except KeyError:
            continue

        rows.append({
            "PDB ID": pdb_id,
            "AF ID": row["ID"],
            "PDB pocket ID": pdb_pocket_id,
            "AF pocket ID": int(row["Pocket ID"]),
            "Jaccard": row["Jaccard"],
            "Res": int(row["Res"]),
            "Mean pLDDT": row["Mean pLDDT"],
            "Score PDB": pdb_row["Score"],
            "Score AF": row["Score"],
            "Drug PDB": pdb_row["Drug"],
            "Drug AF": row["Drug"],
        })

    return pd.DataFrame(rows)

def analyze_origin(pockets_df, targets_df, mode="PDB"):
    if mode == "PDB":
        group_map = dict(
            zip(targets_df["PDB_ID"], targets_df["GROUP"])
        )
    else:
        group_map = dict(
            zip(targets_df["AF_ID"], targets_df["GROUP"])
        )        
    pockets_df["Group"] = pockets_df["ID"].map(group_map)

    drug_perc = pockets_df["Group"].sum() / len(pockets_df) * 100
    return drug_perc

def stats_corr(df, targets_df, stats="pocket"):
    matched_df = df[
        df["pair_status"].isin(["matched", "weak_match"])
    ].copy()
    group_map = dict(
        zip(targets_df["PDB_ID"], targets_df["GROUP"])
    )
    matched_df["Group"] = matched_df["pdb_id"].map(group_map)

    r, _ = pearsonr(
        matched_df[f"pdb_fpocket_{stats}_score"],
        matched_df[f"af_fpocket_{stats}_score"]
    )

    group1 = matched_df[matched_df["Group"] == 1]
    group0 = matched_df[matched_df["Group"] == 0]

    fig, ax = plt.subplots(figsize=(6, 6))

    ax.scatter(
        group1[f"pdb_fpocket_{stats}_score"],
        group1[f"af_fpocket_{stats}_score"],
        color="green",
        label="Group 1",
        alpha=0.7
    )

    ax.scatter(
        group0[f"pdb_fpocket_{stats}_score"],
        group0[f"af_fpocket_{stats}_score"],
        color="orange",
        label="Group 0",
        alpha=0.7
    )

    ax.text(
        0.05,
        0.95,
        f"Pearson r = {r:.3f}",
        transform=ax.transAxes,
        va="top",
        bbox=dict(boxstyle="round", alpha=0.3)
    )

    ax.legend()
    ax.set_xlabel(f"PDB {stats} score")
    ax.set_ylabel(f"AlphaFold {stats} score")
    ax.set_title(f"PDB vs AlphaFold {stats} scores")

    plt.tight_layout()
    plt.show()

    return r

def plddt_impact(pockets_af):
    plddt_stats = (
        pockets_af["Mean pLDDT"]
        .agg(["mean", "median", "max", "min", "std", "count"])
        .round(2)
    )
    plddt_stats = pd.DataFrame(
        [plddt_stats],
        index=["all"]
    )

    plddt_by_status = (
        pockets_af
        .groupby("Status")["Mean pLDDT"]
        .agg(["mean", "median", "max", "min", "std", "count"])
        .round(2)
    )
    plddt_by_status = pd.concat(
        [plddt_stats, plddt_by_status]
    )
    plddt_by_status["count"] = plddt_by_status["count"].astype(int)
    plddt_by_status = plddt_by_status.reindex(
        ["all", "matched", "weak_match", "af_only"]
    )

    # correlations
    matched_af = pockets_af[
        pockets_af["Status"].isin(["matched", "weak_match"])
    ].copy()
    matched_af["Jaccard"] = pd.to_numeric(
        matched_af["Jaccard"],
        errors="coerce"
    )
    matched_af["Mean pLDDT"] = pd.to_numeric(
        matched_af["Mean pLDDT"],
        errors="coerce"
    )
    
    r_score, _ = pearsonr(
        pockets_af["Mean pLDDT"],
        pockets_af["Score"]
    )
    r_drug, _ = pearsonr(
        pockets_af["Mean pLDDT"],
        pockets_af["Drug"]
    )
    r_jaccard, _ = pearsonr(
        matched_af["Mean pLDDT"],
        matched_af["Jaccard"]
    )

    r_dict = {
        "score" : r_score,
        "drug" : r_drug,
        "jaccard" : r_jaccard
    }

    print(f"Pearson correlation (pLDDT vs Score): {r_score:.3f}")
    print(f"Pearson correlation (pLDDT vs Drug): {r_drug:.3f}")
    print(f"Pearson correlation (pLDDT vs Jaccard) for matched and weak matched pairs: {r_jaccard:.3f}")

    return plddt_by_status, r_dict

def characteristics_corr(pockets_df, score = "Score"):
    characters = ["Vol", "Hydro", "Polar", "Charge"]
    corrs = {}

    for character in characters:
        r, _ = pearsonr(
            pockets_df[score],
            pockets_df[character]
        )
        corrs[character] = r

    # flex corr only considering PDB pockets
    pockets_flex = (
        pockets_df[pockets_df["Origin"] == "PDB"]
        .drop(columns="Origin")
        .reset_index(drop=True)
    ).copy()

    r_flex, _ = pearsonr(
        pockets_flex[score],
        pockets_flex["Flex"]
    )
    corrs["Flex"] = r_flex

    return corrs

def characteristics_summary(df):
    matched_df = df[
        df["pair_status"].isin(["matched", "weak_match"])
    ].copy()

    stats = ["volume", 
            "hydrophobicity_score",
            "polarity_score",
            "charge_score"]

    summary = []
    n = len(matched_df)

    for stat in stats:
        pdb_col = f"pdb_fpocket_{stat}"
        af_col = f"af_fpocket_{stat}"
        diff_col = f"delta_fpocket_{stat}"

        summary.append({
            "Descriptor": stat,
            "Mean PDB": matched_df[pdb_col].mean(),
            "Mean AF": matched_df[af_col].mean(),
            "Mean diff" : matched_df[diff_col].mean(),
            "% PDB > AF": 100 * (matched_df[pdb_col] > matched_df[af_col]).mean(),
            "% PDB = AF": 100 * (matched_df[pdb_col] == matched_df[af_col]).mean(),
            "% PDB < AF": 100 * (matched_df[pdb_col] < matched_df[af_col]).mean(),
        })

    summary_df = pd.DataFrame(summary).round(2)
    return summary_df

# visualization utils

def print_target_summary(pdb_id, af_id):
    TARGET_LIST_PATH = "targets/targets_list.csv"
    SEQ_RESULTS_PATH = "targets/sequence_compare_results.csv"
    RMSD_RESUTLS_PATH =  "targets/rmsd_results.csv"

    targets_df = pd.read_csv(TARGET_LIST_PATH)
    seq_results = pd.read_csv(SEQ_RESULTS_PATH)
    rmsd_results = pd.read_csv(RMSD_RESUTLS_PATH)

    # validate mapping
    target_row = targets_df[
        (targets_df["PDB_ID"] == pdb_id)
        & (targets_df["AF_ID"] == af_id)
    ]
    if target_row.empty:
        print(
            f"WARNING: {pdb_id} and {af_id} "
            "do not form a valid PDB–AF pair."
        )
        return
    target_row = target_row.iloc[0]

    # metadata
    protein_class = target_row["CLASS"]
    group = (
        "likely druggable"
        if target_row["GROUP"] == 1
        else "likely NON-druggable"
    )
    notes = target_row["NOTES"]

    # sequence alignment
    seq_row = seq_results[
        (seq_results["PDB_ID"] == pdb_id)
        & (seq_results["AF_ID"] == af_id)
    ].iloc[0]

    # structural alignment
    rmsd_row = rmsd_results[
        (rmsd_results["PDB_ID"] == pdb_id)
        & (rmsd_results["AF_ID"] == af_id)
    ].iloc[0]

    # output
    print("PROTEIN INFORMATION")
    print(f"  Class : {protein_class}")
    print(f"  Group : {group}")
    print(f"  Notes : {notes}")

    print("\nSequence alignment summary")
    print(f"  PDB sequence length : {seq_row['pdb_length']}")
    print(f"  AF sequence length  : {seq_row['af_length']}")
    print(f"  Identity            : {seq_row['identity']:.2f}%")
    print(f"  Coverage            : {seq_row['coverage']:.2f}%")

    print("\n3D alignment summary")
    print(f"  Matched Cα atoms : {rmsd_row['matched_ca_count']}")
    print(f"  RMSD             : {rmsd_row['RMSD']:.3f} Å")

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

def show_alignment(structure_pdb, structure_af, width=500, height=400, opacity=0.5):
    view = py3Dmol.view(
        width=width,
        height=height
    )

    with open(structure_pdb) as f:
        view.addModel(f.read(), "pdb")
    with open(structure_af) as f:
        view.addModel(f.read(), "pdb")

    view.setStyle(
        {"model": 0},
        {"cartoon": {"color": "violet"}}
    )
    view.addSurface(
        py3Dmol.MS,
        {
            "color": "violet",
            "opacity": opacity
        },
        {"model": 0}
    )
    view.setStyle(
        {"model": 1},
        {"cartoon": {"color": "orange"}}
    )
    view.addSurface(
        py3Dmol.MS,
        {
            "color": "orange",
            "opacity": opacity
        },
        {"model": 1}
    )

    view.setBackgroundColor("black")
    view.zoomTo()
    view.show()

def show_pocket_alignment(structure_pdb, structure_af,
                          pocket_files_pdb, pocket_files_af,
                          pdb_pocket, af_pocket,
                          width=500, height=400):

    view = py3Dmol.view(
        width=width,
        height=height
    )

    # show proteins
    with open(structure_pdb) as f:
        view.addModel(f.read(), "pdb")  # model 0
    with open(structure_af) as f:
        view.addModel(f.read(), "pdb")  # model 1

    view.setStyle(
        {"model": 0},
        {"cartoon": 
            {
            "color": "violet",
            "opacity": 0.5
            }
        }
    )
    view.setStyle(
        {"model": 1},
        {"cartoon": 
            {
            "color": "orange",
            "opacity": 0.5
            }
        }
    )

    # show PDB pocket
    pocket_vert_pdb = pocket_files_pdb[pdb_pocket]["vert"]
    pocket_atm_pdb = pocket_files_pdb[pdb_pocket]["atm"]

    with open(pocket_vert_pdb) as f:    
        view.addModel(f.read(), "pqr")  # model 2

    view.setStyle(
        {"model": 2},
        {"sphere": {"radius": 1.0, "color": "violet"}},
    )
    view.addSurface(
        py3Dmol.VDW,
        {
            "color": "violet",
            "opacity": 0.8
        },
        {"model": 2}
    )

    # show AF pocket
    pocket_vert_af = pocket_files_af[af_pocket]["vert"]
    pocket_atm_af = pocket_files_af[af_pocket]["atm"]

    with open(pocket_vert_af) as f:
        view.addModel(f.read(), "pqr")  # model 3

    view.setStyle(
        {"model": 3},
        {"sphere": {"radius": 1.0, "color": "orange"}},
    )
    view.addSurface(
        py3Dmol.VDW,
        {
            "color": "orange",
            "opacity": 0.8
        },
        {"model": 3}
    )

    # sticks
    with open(pocket_atm_pdb) as f:
        view.addModel(f.read(), "pdb")  # model 4
    with open(pocket_atm_af) as f:
        view.addModel(f.read(), "pdb")  # model 5

    view.setStyle(
        {"model": 4},
        {"stick": {"color": "violet"}}
    )
    view.setStyle(
        {"model": 5},
        {"stick": {"color": "orange"}}
    )

    view.setBackgroundColor("black")
    view.zoomTo({"model": 3})
    return view

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
    ANALYSIS_DIR = Path("analysis/pocket_comparison/outputs_fpocket/local")/f"{PDB_id}_vs_{AF_id}"
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
    ANALYSIS_DIR = Path("analysis/pocket_comparison/outputs_fpocket/local")/f"{PDB_id}_vs_{AF_id}"
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
