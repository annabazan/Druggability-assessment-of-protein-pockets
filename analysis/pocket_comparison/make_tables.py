#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd


INPUT_DIR = Path("analysis/pocket_comparison/outputs/global")
OUTPUT_DIR = Path("analysis/pocket_comparison/outputs/readable")


def existing_cols(df: pd.DataFrame, cols: list[str]) -> list[str]:
    return [c for c in cols if c in df.columns]


def round_numeric(df: pd.DataFrame, digits: int = 3) -> pd.DataFrame:
    out = df.copy()
    for col in out.select_dtypes(include="number").columns:
        out[col] = out[col].round(digits)
    return out


def add_interpretation_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    numeric_cols = [
        "jaccard",
        "shared_residues",
        "pdb_residue_count",
        "af_residue_count",
        "af_pocket_mean_plddt",
        "af_pocket_min_plddt",
        "pdb_fpocket_drug_score",
        "af_fpocket_drug_score",
        "pdb_fpocket_pocket_score",
        "af_fpocket_pocket_score",
    ]

    for col in numeric_cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")

    if {"shared_residues", "pdb_residue_count"}.issubset(out.columns):
        out["pdb_coverage"] = out["shared_residues"] / out["pdb_residue_count"]
        out["pdb_coverage"] = out["pdb_coverage"].replace([np.inf, -np.inf], np.nan)

    if {"shared_residues", "af_residue_count"}.issubset(out.columns):
        out["af_coverage"] = out["shared_residues"] / out["af_residue_count"]
        out["af_coverage"] = out["af_coverage"].replace([np.inf, -np.inf], np.nan)

    def match_quality(row):
        status = row.get("pair_status")
        j = row.get("jaccard", np.nan)

        if status == "matched":
            if pd.notna(j) and j >= 0.70:
                return "strong_match"
            if pd.notna(j) and j >= 0.50:
                return "good_match"
            return "accepted_match"

        if status == "weak_match":
            return "weak_overlap"

        return status

    out["match_quality"] = out.apply(match_quality, axis=1)

    return out


def make_global_pairs_readable(pairs: pd.DataFrame) -> pd.DataFrame:
    pairs = add_interpretation_columns(pairs)

    cols = [
        "target_id",
        "pdb_id",
        "af_id",
        "pair_status",
        "match_quality",
        "residue_mode_used",

        "pdb_pocket_id",
        "af_pocket_id",
        "pdb_pocket_num",
        "af_pocket_num",

        "jaccard",
        "pdb_coverage",
        "af_coverage",
        "shared_residues",
        "pdb_residue_count",
        "af_residue_count",

        "rmsd_common_ca_global",
        "rmsd_common_ca_local_aligned",

        "af_pocket_mean_plddt",
        "af_pocket_min_plddt",
        "af_pocket_fraction_plddt_lt_50",

        "pdb_fpocket_pocket_score",
        "af_fpocket_pocket_score",
        "delta_fpocket_pocket_score",
        "abs_delta_fpocket_pocket_score",

        "pdb_fpocket_drug_score",
        "af_fpocket_drug_score",
        "delta_fpocket_drug_score",
        "abs_delta_fpocket_drug_score",

        "pdb_fpocket_volume",
        "af_fpocket_volume",
        "delta_fpocket_volume",
    ]

    out = pairs[existing_cols(pairs, cols)].copy()

    sort_cols = existing_cols(
        out,
        ["target_id", "pair_status", "pdb_pocket_num", "af_pocket_num"],
    )

    out = out.sort_values(sort_cols, na_position="last")
    return round_numeric(out)


def make_matched_pockets_readable(pairs: pd.DataFrame) -> pd.DataFrame:
    pairs = add_interpretation_columns(pairs)

    out = pairs[pairs["pair_status"].isin(["matched", "weak_match"])].copy()

    cols = [
        "target_id",
        "pdb_id",
        "af_id",
        "pair_status",
        "match_quality",
        "residue_mode_used",

        "pdb_pocket_id",
        "af_pocket_id",
        "pdb_pocket_num",
        "af_pocket_num",

        "jaccard",
        "pdb_coverage",
        "af_coverage",
        "shared_residues",
        "pdb_residue_count",
        "af_residue_count",

        "rmsd_common_ca_local_aligned",

        "af_pocket_mean_plddt",
        "af_pocket_min_plddt",

        "pdb_fpocket_drug_score",
        "af_fpocket_drug_score",
        "delta_fpocket_drug_score",
        "abs_delta_fpocket_drug_score",

        "pdb_fpocket_pocket_score",
        "af_fpocket_pocket_score",
        "delta_fpocket_pocket_score",
        "abs_delta_fpocket_pocket_score",
    ]

    out = out[existing_cols(out, cols)]

    out = out.sort_values(
        ["target_id", "jaccard"],
        ascending=[True, False],
        na_position="last",
    )

    return round_numeric(out)


def make_target_summary_from_pairs(pairs: pd.DataFrame) -> pd.DataFrame:
    pairs = add_interpretation_columns(pairs)

    matched_like = pairs[pairs["pair_status"].isin(["matched", "weak_match"])].copy()
    matched_only = pairs[pairs["pair_status"] == "matched"].copy()

    counts = (
        pairs.groupby("target_id")["pair_status"]
        .value_counts()
        .unstack(fill_value=0)
    )

    for col in ["matched", "weak_match", "pdb_only", "af_only"]:
        if col not in counts.columns:
            counts[col] = 0

    base = pairs.groupby("target_id").agg(
        pdb_id=("pdb_id", "first"),
        af_id=("af_id", "first"),
        residue_mode_used=(
            "residue_mode_used",
            lambda x: ",".join(sorted(set(x.dropna().astype(str)))),
        ),
        n_pdb_pockets=(
            "pdb_pocket_id",
            lambda x: x.replace("", np.nan).dropna().nunique(),
        ),
        n_af_pockets=(
            "af_pocket_id",
            lambda x: x.replace("", np.nan).dropna().nunique(),
        ),
        max_jaccard=("jaccard", "max"),
        mean_af_plddt=("af_pocket_mean_plddt", "mean"),
        min_af_plddt=("af_pocket_min_plddt", "min"),
    )

    matched_like_stats = matched_like.groupby("target_id").agg(
        mean_jaccard_matched_like=("jaccard", "mean"),
        median_jaccard_matched_like=("jaccard", "median"),
        mean_rmsd_local=("rmsd_common_ca_local_aligned", "mean"),
        mean_abs_delta_drug_score=("abs_delta_fpocket_drug_score", "mean"),
        mean_abs_delta_pocket_score=("abs_delta_fpocket_pocket_score", "mean"),
    )

    matched_only_stats = matched_only.groupby("target_id").agg(
        mean_jaccard_matched_only=("jaccard", "mean"),
        median_jaccard_matched_only=("jaccard", "median"),
    )

    out = (
        base
        .join(counts[["matched", "weak_match", "pdb_only", "af_only"]], how="left")
        .join(matched_like_stats, how="left")
        .join(matched_only_stats, how="left")
    )

    out["matched_like"] = out["matched"] + out["weak_match"]

    out["strict_matched_fraction_pdb"] = (
        out["matched"] / out["n_pdb_pockets"].replace(0, np.nan)
    )

    out["matched_like_fraction_pdb"] = (
        out["matched_like"] / out["n_pdb_pockets"].replace(0, np.nan)
    )

    out["matched_like_fraction_af"] = (
        out["matched_like"] / out["n_af_pockets"].replace(0, np.nan)
    )

    ordered_cols = [
        "pdb_id",
        "af_id",
        "residue_mode_used",

        "n_pdb_pockets",
        "n_af_pockets",

        "matched",
        "weak_match",
        "matched_like",
        "pdb_only",
        "af_only",

        "strict_matched_fraction_pdb",
        "matched_like_fraction_pdb",
        "matched_like_fraction_af",

        "max_jaccard",
        "mean_jaccard_matched_only",
        "median_jaccard_matched_only",
        "mean_jaccard_matched_like",
        "median_jaccard_matched_like",

        "mean_rmsd_local",

        "mean_af_plddt",
        "min_af_plddt",

        "mean_abs_delta_drug_score",
        "mean_abs_delta_pocket_score",
    ]

    out = out[ordered_cols].reset_index()

    out = out.sort_values(
        ["strict_matched_fraction_pdb", "max_jaccard"],
        ascending=[True, True],
        na_position="first",
    )

    return round_numeric(out)


def make_problematic_targets(summary: pd.DataFrame) -> pd.DataFrame:
    out = summary.copy()

    # Prosta heurystyka: mało strict matched albo bardzo niski max Jaccard.
    out["problem_flag"] = ""

    out.loc[out["matched"] == 0, "problem_flag"] = "zero_strict_matches"
    out.loc[
        (out["matched"] > 0) & (out["strict_matched_fraction_pdb"] < 0.20),
        "problem_flag",
    ] = "low_strict_match_fraction"

    out.loc[
        (out["max_jaccard"] < 0.30) & (out["problem_flag"] == ""),
        "problem_flag",
    ] = "max_jaccard_below_match_threshold"

    out = out[out["problem_flag"] != ""].copy()

    return out.sort_values(
        ["matched", "max_jaccard", "strict_matched_fraction_pdb"],
        ascending=[True, True, True],
        na_position="first",
    )


def make_score_disagreements(pairs: pd.DataFrame) -> pd.DataFrame:
    pairs = add_interpretation_columns(pairs)

    out = pairs[pairs["pair_status"].isin(["matched", "weak_match"])].copy()

    cols = [
        "target_id",
        "pdb_id",
        "af_id",
        "pair_status",
        "match_quality",

        "pdb_pocket_id",
        "af_pocket_id",
        "jaccard",
        "pdb_coverage",
        "af_coverage",

        "af_pocket_mean_plddt",

        "pdb_fpocket_drug_score",
        "af_fpocket_drug_score",
        "delta_fpocket_drug_score",
        "abs_delta_fpocket_drug_score",

        "pdb_fpocket_pocket_score",
        "af_fpocket_pocket_score",
        "delta_fpocket_pocket_score",
        "abs_delta_fpocket_pocket_score",
    ]

    out = out[existing_cols(out, cols)]

    sort_by = []
    if "abs_delta_fpocket_drug_score" in out.columns:
        sort_by.append("abs_delta_fpocket_drug_score")
    if "abs_delta_fpocket_pocket_score" in out.columns:
        sort_by.append("abs_delta_fpocket_pocket_score")

    if sort_by:
        out = out.sort_values(sort_by, ascending=False, na_position="last")

    return round_numeric(out)


def make_low_plddt_high_score_cases(pairs: pd.DataFrame) -> pd.DataFrame:
    pairs = add_interpretation_columns(pairs)

    score_col = None
    for candidate in ["af_fpocket_drug_score", "af_fpocket_pocket_score"]:
        if candidate in pairs.columns:
            score_col = candidate
            break

    if score_col is None or "af_pocket_mean_plddt" not in pairs.columns:
        return pd.DataFrame()

    out = pairs[
        (pairs["af_pocket_mean_plddt"] < 70)
        & pairs[score_col].notna()
    ].copy()

    cols = [
        "target_id",
        "pdb_id",
        "af_id",
        "pair_status",
        "match_quality",
        "af_pocket_id",
        "af_pocket_num",
        "jaccard",
        "af_pocket_mean_plddt",
        "af_pocket_min_plddt",
        "af_pocket_fraction_plddt_lt_50",
        "af_fpocket_drug_score",
        "af_fpocket_pocket_score",
        "af_fpocket_volume",
    ]

    out = out[existing_cols(out, cols)]
    out = out.sort_values(score_col, ascending=False, na_position="last")

    return round_numeric(out)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    pairs_path = INPUT_DIR / "global_pocket_pairs.csv"

    if not pairs_path.exists():
        raise FileNotFoundError(f"Missing file: {pairs_path}")

    pairs = pd.read_csv(pairs_path)

    global_pairs_readable = make_global_pairs_readable(pairs)
    matched_pockets_readable = make_matched_pockets_readable(pairs)
    target_summary_readable = make_target_summary_from_pairs(pairs)
    problematic_targets = make_problematic_targets(target_summary_readable)
    score_disagreements = make_score_disagreements(pairs)
    low_plddt_high_score = make_low_plddt_high_score_cases(pairs)

    global_pairs_readable.to_csv(
        OUTPUT_DIR / "global_pocket_pairs_readable.csv",
        index=False,
    )

    matched_pockets_readable.to_csv(
        OUTPUT_DIR / "matched_pockets_readable.csv",
        index=False,
    )

    target_summary_readable.to_csv(
        OUTPUT_DIR / "target_summary_readable.csv",
        index=False,
    )

    problematic_targets.to_csv(
        OUTPUT_DIR / "problematic_targets_readable.csv",
        index=False,
    )

    score_disagreements.to_csv(
        OUTPUT_DIR / "score_disagreements_readable.csv",
        index=False,
    )

    low_plddt_high_score.to_csv(
        OUTPUT_DIR / "low_plddt_high_score_readable.csv",
        index=False,
    )

    print("Readable tables written to:")
    print(f"  {OUTPUT_DIR}")
    print()
    print("Files:")
    for path in sorted(OUTPUT_DIR.glob("*.csv")):
        print(f"  {path}")


if __name__ == "__main__":
    main()