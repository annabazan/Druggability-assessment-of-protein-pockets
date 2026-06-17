#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional, Any

import numpy as np
import pandas as pd
import scipy

from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.Align import PairwiseAligner

ResidueKey = Tuple[str, int, str]

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare fpocket pockets between PDB and AlphaFold structures."
    )

    parser.add_argument(
        "--targets",
        default="targets/targets_list.csv",
        help="CSV with at least PDB_ID and AF_ID columns.",
    )
    parser.add_argument(
        "--pdb-fpocket",
        default="pocket_detection/fpocket/pdb_out",
        help="Directory with fpocket outputs for PDB structures.",
    )
    parser.add_argument(
        "--af-fpocket",
        default="pocket_detection/fpocket/alpha_fold_out",
        help="Directory with fpocket outputs for AlphaFold structures.",
    )
    parser.add_argument(
        "--pdb-structures",
        default="targets/filtered_pdb",
        help="Directory with full PDB structures.",
    )
    parser.add_argument(
        "--af-structures",
        default="targets/3D_aligned_alpha_fold",
        help="Directory with full AlphaFold structures. pLDDT is read from B-factor field.",
    )
    parser.add_argument(
        "--out",
        default="analysis/pocket_comparison/outputs",
        help="Output directory.",
    )
    parser.add_argument(
        "--match-threshold",
        type=float,
        default=0.30,
        help="Jaccard threshold for accepted matched pockets.",
    )
    parser.add_argument(
        "--weak-threshold",
        type=float,
        default=0.10,
        help="Jaccard threshold for weak matches.",
    )
    parser.add_argument(
        "--residue-mode",
        choices=["auto", "chain", "number"],
        default="auto",
        help=(
            "Residue identity mode. "
            "'chain' uses chain+resseq+icode. "
            "'number' uses only residue number+icode. "
            "'auto' tries chain-aware first, then falls back to number-only if all overlaps are zero."
        ),
    )
    parser.add_argument(
        "--score-source",
        choices=["fpocket", "rescored"],
        default="fpocket",
        help="Source of score values used in statistics.",
    )
    parser.add_argument(
        "--pdb-rescored",
        default="pocket_detection/rescoring/fpocket_pdb_rescored_out",
        help="Directory with P2Rank/PRANK rescored CSV files for PDB fpocket pockets.",
    )
    parser.add_argument(
        "--af-rescored",
        default="pocket_detection/rescoring/fpocket_alpha_fold_rescored_out",
        help="Directory with P2Rank/PRANK rescored CSV files for AlphaFold fpocket pockets.",
    )

    return parser.parse_args()


def normalize_key(key: str) -> str:
    key = key.strip().lower()
    key = re.sub(r"^header\s*\d*\s*[-:]?\s*", "", key)
    key = key.replace("druggability", "drug")
    key = re.sub(r"[^a-z0-9]+", "_", key).strip("_")

    if key == "score":
        return "pocket_score"
    if "drug" in key and "score" in key:
        return "drug_score"
    if "pocket" in key and "score" in key:
        return "pocket_score"

    return key


def parse_float(text: str) -> Optional[float]:
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(text))
    if not match:
        return None
    try:
        return float(match.group(0))
    except ValueError:
        return None


def pocket_number_from_name(path: Path) -> Optional[int]:
    match = re.search(r"pocket(\d+)", path.name)
    if not match:
        return None
    return int(match.group(1))


def pocket_id_from_number(num: int) -> str:
    return f"pocket{num}"


def residue_key_from_pdb_line(line: str, mode: str) -> Optional[ResidueKey]:
    try:
        chain = line[21].strip() or "_"
        resseq = int(line[22:26].strip())
        icode = line[26].strip() or ""
    except Exception:
        return None

    if mode == "number":
        chain = "*"

    return chain, resseq, icode


def extract_sequence_and_residue_keys(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    sequence = []
    residue_keys = []

    seen = set()

    for model in structure:
        for chain in model:
            for residue in chain:

                if residue.id[0] != " ":
                    continue

                resname = residue.resname.upper()

                try:
                    aa = protein_letters_3to1[resname.capitalize()]
                except KeyError:
                    continue

                key = (
                    chain.id,
                    residue.id[1],
                    residue.id[2].strip()
                )

                if key in seen:
                    continue

                seen.add(key)

                sequence.append(aa)
                residue_keys.append(key)

    return "".join(sequence), residue_keys


def build_residue_mapping(
    pdb_keys,
    af_keys,
    alignment
):

    mapping = {}

    for pdb_block, af_block in zip(
        alignment.aligned[0],
        alignment.aligned[1]
    ):

        pdb_start, pdb_end = pdb_block
        af_start, af_end = af_block

        block_len = min(
            pdb_end - pdb_start,
            af_end - af_start
        )

        for i in range(block_len):

            mapping[
                pdb_keys[pdb_start + i]
            ] = af_keys[af_start + i]

    return mapping


def create_residue_mapping(
    pdb_structure_file,
    af_structure_file
):

    pdb_seq, pdb_keys = \
        extract_sequence_and_residue_keys(
            pdb_structure_file
        )

    af_seq, af_keys = \
        extract_sequence_and_residue_keys(
            af_structure_file
        )

    aligner = PairwiseAligner()

    alignment = aligner.align(
        pdb_seq,
        af_seq
    )[0]

    mapping = build_residue_mapping(
        pdb_keys,
        af_keys,
        alignment
    )

    return mapping


def read_residues_from_pdb(path: Path, mode: str) -> Set[ResidueKey]:
    residues: Set[ResidueKey] = set()

    if not path.exists():
        return residues

    with path.open("r", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                key = residue_key_from_pdb_line(line, mode)
                if key is not None:
                    residues.add(key)

    return residues


def read_ca_coordinates_and_bfactors(path: Path, mode: str):
    coords: Dict[ResidueKey, np.ndarray] = {}
    bfactors: Dict[ResidueKey, float] = {}

    if not path.exists():
        return coords, bfactors

    with path.open("r", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue

            key = residue_key_from_pdb_line(line, mode)
            if key is None:
                continue

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                b = float(line[60:66])
            except ValueError:
                continue

            coords[key] = np.array([x, y, z], dtype=float)
            bfactors[key] = b

    return coords, bfactors


def calculate_jaccard(a: Set[ResidueKey], b: Set[ResidueKey]) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def collect_common_ca_arrays(
    pdb_coords: Dict[ResidueKey, np.ndarray],
    af_coords: Dict[ResidueKey, np.ndarray],
    common_residues: Set[ResidueKey],
):
    xs = []
    ys = []

    for residue in sorted(common_residues):
        if residue in pdb_coords and residue in af_coords:
            xs.append(pdb_coords[residue])
            ys.append(af_coords[residue])

    if len(xs) < 3:
        return None, None

    return np.vstack(xs), np.vstack(ys)


def rmsd_between_arrays(x: np.ndarray, y: np.ndarray) -> float:
    diff = x - y
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def calculate_rmsd_global(
    pdb_coords: Dict[ResidueKey, np.ndarray],
    af_coords: Dict[ResidueKey, np.ndarray],
    common_residues: Set[ResidueKey],
) -> Optional[float]:
    x, y = collect_common_ca_arrays(pdb_coords, af_coords, common_residues)

    if x is None or y is None:
        return None

    return rmsd_between_arrays(x, y)


def calculate_rmsd_local_aligned(
    pdb_coords: Dict[ResidueKey, np.ndarray],
    af_coords: Dict[ResidueKey, np.ndarray],
    common_residues: Set[ResidueKey],
) -> Optional[float]:
    x, y = collect_common_ca_arrays(pdb_coords, af_coords, common_residues)

    if x is None or y is None:
        return None

    x_centroid = x.mean(axis=0)
    y_centroid = y.mean(axis=0)

    x_centered = x - x_centroid
    y_centered = y - y_centroid

    # Kabsch algorithm: align AF pocket coordinates y onto PDB pocket coordinates x.
    covariance = y_centered.T @ x_centered
    u, _, vt = np.linalg.svd(covariance)

    rotation = u @ vt

    # Prevent improper rotation/reflection.
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt

    y_aligned = y_centered @ rotation + x_centroid

    return rmsd_between_arrays(x, y_aligned)


def plddt_summary(
    af_bfactors: Dict[ResidueKey, float],
    residues: Set[ResidueKey],
) -> Dict[str, Any]:
    values = [af_bfactors[r] for r in residues if r in af_bfactors]

    if not values:
        return {
            "af_pocket_mean_plddt": np.nan,
            "af_pocket_median_plddt": np.nan,
            "af_pocket_min_plddt": np.nan,
            "af_pocket_std_plddt": np.nan,
            "af_pocket_fraction_plddt_ge_90": np.nan,
            "af_pocket_fraction_plddt_70_90": np.nan,
            "af_pocket_fraction_plddt_50_70": np.nan,
            "af_pocket_fraction_plddt_lt_50": np.nan,
            "af_pocket_plddt_count": 0,
            "af_pocket_plddt_warning": "missing",
        }

    arr = np.array(values, dtype=float)

    warning = ""
    if np.nanmin(arr) < 0 or np.nanmax(arr) > 100:
        warning = "values_outside_0_100_check_bfactor_field"

    return {
        "af_pocket_mean_plddt": float(np.nanmean(arr)),
        "af_pocket_median_plddt": float(np.nanmedian(arr)),
        "af_pocket_min_plddt": float(np.nanmin(arr)),
        "af_pocket_std_plddt": float(np.nanstd(arr)),
        "af_pocket_fraction_plddt_ge_90": float(np.mean(arr >= 90)),
        "af_pocket_fraction_plddt_70_90": float(np.mean((arr >= 70) & (arr < 90))),
        "af_pocket_fraction_plddt_50_70": float(np.mean((arr >= 50) & (arr < 70))),
        "af_pocket_fraction_plddt_lt_50": float(np.mean(arr < 50)),
        "af_pocket_plddt_count": int(len(arr)),
        "af_pocket_plddt_warning": warning,
    }


def parse_fpocket_descriptors(fpocket_out_dir: Path) -> Dict[int, Dict[str, float]]:
    descriptors: Dict[int, Dict[str, float]] = {}

    # 1. Parse *_info.txt if present.
    txt_files = sorted(fpocket_out_dir.glob("*_info.txt"))
    if not txt_files:
        txt_files = sorted(fpocket_out_dir.glob("*.txt"))

    for txt in txt_files:
        current_pocket: Optional[int] = None

        with txt.open("r", errors="ignore") as handle:
            for raw_line in handle:
                line = raw_line.strip()

                pocket_match = re.search(r"\bPocket\s+(\d+)\b", line, flags=re.IGNORECASE)
                if pocket_match:
                    current_pocket = int(pocket_match.group(1))
                    descriptors.setdefault(current_pocket, {})
                    continue

                if current_pocket is None:
                    continue

                if ":" not in line:
                    continue

                key, value = line.split(":", 1)
                value_num = parse_float(value)
                if value_num is None:
                    continue

                key_norm = normalize_key(key)
                descriptors.setdefault(current_pocket, {})[key_norm] = value_num

    # 2. Also parse pocket headers if descriptors are embedded there.
    for pocket_file in sorted((fpocket_out_dir / "pockets").glob("pocket*_atm.pdb")):
        pocket_num = pocket_number_from_name(pocket_file)
        if pocket_num is None:
            continue

        descriptors.setdefault(pocket_num, {})

        with pocket_file.open("r", errors="ignore") as handle:
            for line in handle:
                if line.startswith(("ATOM", "HETATM")):
                    break

                matches = re.findall(
                    r"([A-Za-z][A-Za-z0-9 /\._\-\(\)]*?)\s*:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)",
                    line,
                )

                for key, value in matches:
                    key_norm = normalize_key(key)
                    value_num = parse_float(value)
                    if value_num is not None:
                        descriptors[pocket_num][key_norm] = value_num

    return descriptors


def find_rescored_csv(rescored_root: Path, structure_id: str) -> Optional[Path]:
    if not rescored_root.exists():
        return None

    patterns = [
        f"{structure_id}.pdb_rescored.csv",
        f"{structure_id}_rescored.csv",
        f"*{structure_id}*rescored*.csv",
    ]

    for pattern in patterns:
        matches = sorted(rescored_root.rglob(pattern))
        if matches:
            return matches[0]

    return None


def normalize_csv_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [
        re.sub(r"[^a-z0-9]+", "_", str(col).strip().lower()).strip("_")
        for col in df.columns
    ]
    return df


def parse_rescored_descriptors(
    rescored_root: Path,
    structure_id: str,
) -> Dict[int, Dict[str, float]]:
    csv_path = find_rescored_csv(rescored_root, structure_id)
    if csv_path is None:
        return {}

    try:
        df = pd.read_csv(csv_path, comment="#", skipinitialspace=True)
    except Exception:
        return {}

    if df.empty:
        return {}

    df = normalize_csv_columns(df)

    if "old_rank" not in df.columns:
        return {}

    descriptors: Dict[int, Dict[str, float]] = {}

    for _, row in df.iterrows():
        if pd.isna(row.get("old_rank")):
            continue

        try:
            pocket_num = int(float(row["old_rank"]))
        except Exception:
            continue

        values: Dict[str, float] = {}

        for col, value in row.items():
            if pd.isna(value):
                continue

            try:
                numeric_value = float(value)
            except Exception:
                continue

            values[f"rescored_{col}"] = numeric_value

            if col == "score":
                values["pocket_score"] = numeric_value


            if col == "probability":
                values["drug_score"] = numeric_value

            if col in {"rank", "old_rank"}:
                values[col] = numeric_value

        if values:
            descriptors[pocket_num] = values

    return descriptors


def merge_rescored_descriptors(
    pockets: List[Dict[str, Any]],
    rescored_root: Path,
    structure_id: str,
    label: str,
) -> None:
    rescored = parse_rescored_descriptors(rescored_root, structure_id)

    if not rescored:
        raise FileNotFoundError(
            f"No rescored descriptors found for {label} {structure_id} in {rescored_root}"
        )

    for pocket in pockets:
        pocket_num = pocket["pocket_num"]
        if pocket_num in rescored:
            pocket["descriptors"].update(rescored[pocket_num])


def load_pockets(fpocket_out_dir: Path, mode: str) -> List[Dict[str, Any]]:
    pocket_dir = fpocket_out_dir / "pockets"
    descriptors = parse_fpocket_descriptors(fpocket_out_dir)

    pockets = []

    for pocket_file in sorted(pocket_dir.glob("pocket*_atm.pdb")):
        pocket_num = pocket_number_from_name(pocket_file)
        if pocket_num is None:
            continue

        residues = read_residues_from_pdb(pocket_file, mode)

        pockets.append(
            {
                "pocket_num": pocket_num,
                "pocket_id": pocket_id_from_number(pocket_num),
                "file": pocket_file,
                "residues": residues,
                "residue_count": len(residues),
                "descriptors": descriptors.get(pocket_num, {}),
            }
        )

    pockets.sort(key=lambda x: x["pocket_num"])
    return pockets


def build_similarity_matrix(
    pdb_pockets: List[Dict[str, Any]],
    af_pockets: List[Dict[str, Any]],
    residue_mapping
) -> np.ndarray:
    sim = np.zeros((len(pdb_pockets), len(af_pockets)), dtype=float)

    for i, p in enumerate(pdb_pockets):
        for j, a in enumerate(af_pockets):
            mapped_pdb = {
                residue_mapping[r]
                for r in p["residues"]
                if r in residue_mapping
            }

            sim[i, j] = calculate_jaccard(
                mapped_pdb,
                a["residues"]
            )
    return sim


def assign_unique_matches(similarity: np.ndarray) -> List[Tuple[int, int, float]]:
    if similarity.size == 0:
        return []

    try:
        from scipy.optimize import linear_sum_assignment

        row_ind, col_ind = linear_sum_assignment(-similarity)
        matches = [(int(i), int(j), float(similarity[i, j])) for i, j in zip(row_ind, col_ind)]
        matches.sort(key=lambda x: x[2], reverse=True)
        return matches

    except Exception:
        # Fallback: greedy unique matching.
        candidates = []
        n_rows, n_cols = similarity.shape
        for i in range(n_rows):
            for j in range(n_cols):
                candidates.append((float(similarity[i, j]), i, j))

        candidates.sort(reverse=True)

        used_rows = set()
        used_cols = set()
        matches = []

        for score, i, j in candidates:
            if i in used_rows or j in used_cols:
                continue
            used_rows.add(i)
            used_cols.add(j)
            matches.append((i, j, score))

        return matches


def add_prefixed_descriptors(row: Dict[str, Any], prefix: str, descriptors: Dict[str, float]) -> None:
    for key, value in descriptors.items():
        row[f"{prefix}_fpocket_{key}"] = value


def add_descriptor_deltas(
    row: Dict[str, Any],
    pdb_desc: Dict[str, float],
    af_desc: Dict[str, float],
) -> None:
    common_keys = set(pdb_desc) & set(af_desc)

    for key in common_keys:
        try:
            pdb_value = float(pdb_desc[key])
            af_value = float(af_desc[key])
        except Exception:
            continue

        row[f"delta_fpocket_{key}"] = af_value - pdb_value
        row[f"abs_delta_fpocket_{key}"] = abs(af_value - pdb_value)


def make_pair_rows(
    pdb_id: str,
    af_id: str,
    target_id: str,
    pdb_pockets: List[Dict[str, Any]],
    af_pockets: List[Dict[str, Any]],
    similarity: np.ndarray,
    pdb_coords: Dict[ResidueKey, np.ndarray],
    af_coords: Dict[ResidueKey, np.ndarray],
    af_bfactors: Dict[ResidueKey, float],
    match_threshold: float,
    weak_threshold: float,
    residue_mode_used: str,
    residue_mapping
) -> pd.DataFrame:
    raw_matches = assign_unique_matches(similarity)

    accepted_matches = []
    used_pdb = set()
    used_af = set()

    for i, j, score in raw_matches:
        if score < weak_threshold:
            continue

        status = "matched" if score >= match_threshold else "weak_match"
        accepted_matches.append((i, j, score, status))
        used_pdb.add(i)
        used_af.add(j)

    rows = []

    for i, j, score, status in accepted_matches:
        p = pdb_pockets[i]
        a = af_pockets[j]
        mapped_pdb = {
            residue_mapping[r]
            for r in p["residues"]
            if r in residue_mapping
        }

        common = mapped_pdb & a["residues"]

        row: Dict[str, Any] = {
            "target_id": target_id,
            "pdb_id": pdb_id,
            "af_id": af_id,
            "pair_status": status,
            "residue_mode_used": residue_mode_used,
            "pdb_pocket_id": p["pocket_id"],
            "af_pocket_id": a["pocket_id"],
            "pdb_pocket_num": p["pocket_num"],
            "af_pocket_num": a["pocket_num"],
            "jaccard": score,
            "shared_residues": len(common),
            "pdb_residue_count": p["residue_count"],
            "af_residue_count": a["residue_count"],
            "rmsd_common_ca_global": calculate_rmsd_global(pdb_coords, af_coords, common),
            "rmsd_common_ca_local_aligned": calculate_rmsd_local_aligned(pdb_coords, af_coords, common),
        }

        row.update(plddt_summary(af_bfactors, a["residues"]))

        add_prefixed_descriptors(row, "pdb", p["descriptors"])
        add_prefixed_descriptors(row, "af", a["descriptors"])
        add_descriptor_deltas(row, p["descriptors"], a["descriptors"])

        rows.append(row)

    for i, p in enumerate(pdb_pockets):
        if i in used_pdb:
            continue

        row = {
            "target_id": target_id,
            "pdb_id": pdb_id,
            "af_id": af_id,
            "pair_status": "pdb_only",
            "residue_mode_used": residue_mode_used,
            "pdb_pocket_id": p["pocket_id"],
            "af_pocket_id": "",
            "pdb_pocket_num": p["pocket_num"],
            "af_pocket_num": np.nan,
            "jaccard": 0.0,
            "shared_residues": 0,
            "pdb_residue_count": p["residue_count"],
            "af_residue_count": np.nan,
            "rmsd_common_ca_global": np.nan,
            "rmsd_common_ca_local_aligned": np.nan,
        }

        add_prefixed_descriptors(row, "pdb", p["descriptors"])
        rows.append(row)

    for j, a in enumerate(af_pockets):
        if j in used_af:
            continue

        row = {
            "target_id": target_id,
            "pdb_id": pdb_id,
            "af_id": af_id,
            "pair_status": "af_only",
            "residue_mode_used": residue_mode_used,
            "pdb_pocket_id": "",
            "af_pocket_id": a["pocket_id"],
            "pdb_pocket_num": np.nan,
            "af_pocket_num": a["pocket_num"],
            "jaccard": 0.0,
            "shared_residues": 0,
            "pdb_residue_count": np.nan,
            "af_residue_count": a["residue_count"],
            "rmsd_common_ca": np.nan,
        }

        row.update(plddt_summary(af_bfactors, a["residues"]))
        add_prefixed_descriptors(row, "af", a["descriptors"])
        rows.append(row)

    return pd.DataFrame(rows)


def safe_corr(df: pd.DataFrame, col_a: str, col_b: str) -> float:
    if col_a not in df.columns or col_b not in df.columns:
        return np.nan

    tmp = df[[col_a, col_b]].apply(pd.to_numeric, errors="coerce").dropna()

    if len(tmp) < 2:
        return np.nan

    if tmp[col_a].nunique() < 2 or tmp[col_b].nunique() < 2:
        return np.nan

    return float(tmp[col_a].corr(tmp[col_b]))


def safe_mean(df: pd.DataFrame, col: str) -> float:
    if col not in df.columns:
        return np.nan
    return float(pd.to_numeric(df[col], errors="coerce").mean())


def safe_median(df: pd.DataFrame, col: str) -> float:
    if col not in df.columns:
        return np.nan
    return float(pd.to_numeric(df[col], errors="coerce").median())


def make_protein_summary(pair_df: pd.DataFrame, pdb_id: str, af_id: str, target_id: str) -> pd.DataFrame:
    matched_like = pair_df[pair_df["pair_status"].isin(["matched", "weak_match"])].copy()
    matched = pair_df[pair_df["pair_status"] == "matched"].copy()

    n_pdb_pockets = int(pair_df["pdb_pocket_id"].replace("", np.nan).dropna().nunique())
    n_af_pockets = int(pair_df["af_pocket_id"].replace("", np.nan).dropna().nunique())

    top1_match = False
    top3_overlap_count = 0

    for _, row in matched_like.iterrows():
        if row.get("pdb_pocket_num") == 1 and row.get("af_pocket_num") == 1:
            top1_match = True

        try:
            if float(row.get("pdb_pocket_num")) <= 3 and float(row.get("af_pocket_num")) <= 3:
                top3_overlap_count += 1
        except Exception:
            pass

    summary = {
        "target_id": target_id,
        "pdb_id": pdb_id,
        "af_id": af_id,
        "n_pdb_pockets": n_pdb_pockets,
        "n_af_pockets": n_af_pockets,
        "n_matched": int((pair_df["pair_status"] == "matched").sum()),
        "n_weak_matches": int((pair_df["pair_status"] == "weak_match").sum()),
        "n_pdb_only": int((pair_df["pair_status"] == "pdb_only").sum()),
        "n_af_only": int((pair_df["pair_status"] == "af_only").sum()),
        "matched_fraction_pdb": float(len(matched_like) / n_pdb_pockets) if n_pdb_pockets else np.nan,
        "matched_fraction_af": float(len(matched_like) / n_af_pockets) if n_af_pockets else np.nan,
        "mean_jaccard_matched_like": safe_mean(matched_like, "jaccard"),
        "median_jaccard_matched_like": safe_median(matched_like, "jaccard"),
        "max_jaccard": safe_mean(pd.DataFrame({"jaccard": [pair_df["jaccard"].max()]}), "jaccard") if len(pair_df) else np.nan,
        "top1_pdb_matches_top1_af": bool(top1_match),
        "top3_overlap_count": int(top3_overlap_count),
        "mean_af_pocket_plddt": safe_mean(pair_df, "af_pocket_mean_plddt"),
        "median_af_pocket_plddt": safe_median(pair_df, "af_pocket_mean_plddt"),
        "fraction_af_pockets_mean_plddt_lt_70": float(
            np.mean(pd.to_numeric(pair_df.get("af_pocket_mean_plddt", pd.Series(dtype=float)), errors="coerce") < 70)
        )
        if "af_pocket_mean_plddt" in pair_df.columns
        else np.nan,
        "corr_pdb_af_pocket_score": safe_corr(matched_like, "pdb_fpocket_pocket_score", "af_fpocket_pocket_score"),
        "corr_pdb_af_drug_score": safe_corr(matched_like, "pdb_fpocket_drug_score", "af_fpocket_drug_score"),
        "mean_abs_delta_pocket_score": safe_mean(matched_like, "abs_delta_fpocket_pocket_score"),
        "mean_abs_delta_drug_score": safe_mean(matched_like, "abs_delta_fpocket_drug_score"),
    }

    return pd.DataFrame([summary])


def process_target_pair(
    row: pd.Series,
    args,
    out_root: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    pdb_id = str(row["PDB_ID"])
    af_id = str(row["AF_ID"])
    target_id = str(row["target_id"]) if "target_id" in row.index else f"{pdb_id}_vs_{af_id}"

    pdb_fpocket_dir = Path(args.pdb_fpocket) / f"{pdb_id}_out"
    af_fpocket_dir = Path(args.af_fpocket) / f"{af_id}_out"

    pdb_structure = Path(args.pdb_structures) / f"{pdb_id}.pdb"
    af_structure = Path(args.af_structures) / f"{af_id}.pdb"

    local_out = out_root / "local" / f"{pdb_id}_vs_{af_id}"
    local_out.mkdir(parents=True, exist_ok=True)

    if not pdb_fpocket_dir.exists():
        print(f"[WARN] Missing PDB fpocket output: {pdb_fpocket_dir}")
    if not af_fpocket_dir.exists():
        print(f"[WARN] Missing AF fpocket output: {af_fpocket_dir}")
    if not pdb_structure.exists():
        print(f"[WARN] Missing PDB structure: {pdb_structure}")
    if not af_structure.exists():
        print(f"[WARN] Missing AF structure: {af_structure}")

    modes_to_try = ["chain"]
    if args.residue_mode == "number":
        modes_to_try = ["number"]
    elif args.residue_mode == "auto":
        modes_to_try = ["chain", "number"]

    selected = None

    for mode in modes_to_try:
        pdb_pockets = load_pockets(pdb_fpocket_dir, mode)
        af_pockets = load_pockets(af_fpocket_dir, mode)
        residue_mapping = create_residue_mapping(pdb_structure, af_structure)
        similarity = build_similarity_matrix(pdb_pockets, af_pockets, residue_mapping)

        max_similarity = float(np.max(similarity)) if similarity.size else 0.0

        selected = (mode, pdb_pockets, af_pockets, similarity)

        if args.residue_mode != "auto":
            break

        if max_similarity > 0:
            break

    residue_mode_used, pdb_pockets, af_pockets, similarity = selected

    if args.score_source == "rescored":
        merge_rescored_descriptors(
            pdb_pockets,
            Path(args.pdb_rescored),
            pdb_id,
            "PDB",
        )
        merge_rescored_descriptors(
            af_pockets,
            Path(args.af_rescored),
            af_id,
            "AF",
        )

    pdb_coords, _ = read_ca_coordinates_and_bfactors(pdb_structure, residue_mode_used)
    af_coords, af_bfactors = read_ca_coordinates_and_bfactors(af_structure, residue_mode_used)

    pair_df = make_pair_rows(
        pdb_id=pdb_id,
        af_id=af_id,
        target_id=target_id,
        pdb_pockets=pdb_pockets,
        af_pockets=af_pockets,
        similarity=similarity,
        pdb_coords=pdb_coords,
        af_coords=af_coords,
        af_bfactors=af_bfactors,
        match_threshold=args.match_threshold,
        weak_threshold=args.weak_threshold,
        residue_mode_used=residue_mode_used,
        residue_mapping=residue_mapping
    )

    pair_df.to_csv(local_out / "pocket_pairs_detailed.csv", index=False)

    matrix_df = pd.DataFrame(
        similarity,
        index=[p["pocket_id"] for p in pdb_pockets],
        columns=[a["pocket_id"] for a in af_pockets],
    )
    matrix_df.to_csv(local_out / "jaccard_matrix.csv")

    summary_df = make_protein_summary(pair_df, pdb_id, af_id, target_id)
    summary_df.to_csv(local_out / "protein_summary.csv", index=False)

    print(
        f"[OK] {pdb_id} vs {af_id}: "
        f"{len(pdb_pockets)} PDB pockets, {len(af_pockets)} AF pockets, "
        f"{int((pair_df['pair_status'] == 'matched').sum())} matched, "
        f"mode={residue_mode_used}"
    )

    return pair_df, summary_df


def make_global_outputs(all_pairs: pd.DataFrame, all_summaries: pd.DataFrame, out_root: Path) -> None:
    global_dir = out_root / "global"
    global_dir.mkdir(parents=True, exist_ok=True)

    all_pairs.to_csv(global_dir / "global_pocket_pairs.csv", index=False)
    all_summaries.to_csv(global_dir / "global_protein_summary.csv", index=False)

    stats = []

    def add(metric: str, value: Any):
        stats.append({"metric": metric, "value": value})

    add("n_pair_rows", len(all_pairs))
    add("n_protein_pairs", len(all_summaries))
    add("n_matched", int((all_pairs["pair_status"] == "matched").sum()))
    add("n_weak_matches", int((all_pairs["pair_status"] == "weak_match").sum()))
    add("n_pdb_only", int((all_pairs["pair_status"] == "pdb_only").sum()))
    add("n_af_only", int((all_pairs["pair_status"] == "af_only").sum()))

    matched_like = all_pairs[all_pairs["pair_status"].isin(["matched", "weak_match"])].copy()

    add("mean_jaccard_matched_like", safe_mean(matched_like, "jaccard"))
    add("median_jaccard_matched_like", safe_median(matched_like, "jaccard"))
    add("mean_af_pocket_plddt", safe_mean(all_pairs, "af_pocket_mean_plddt"))
    add("median_af_pocket_plddt", safe_median(all_pairs, "af_pocket_mean_plddt"))
    add("corr_pdb_af_pocket_score", safe_corr(matched_like, "pdb_fpocket_pocket_score", "af_fpocket_pocket_score"))
    add("corr_pdb_af_drug_score", safe_corr(matched_like, "pdb_fpocket_drug_score", "af_fpocket_drug_score"))
    add("mean_abs_delta_pocket_score", safe_mean(matched_like, "abs_delta_fpocket_pocket_score"))
    add("mean_abs_delta_drug_score", safe_mean(matched_like, "abs_delta_fpocket_drug_score"))

    pd.DataFrame(stats).to_csv(global_dir / "global_statistics.csv", index=False)

    if "af_pocket_mean_plddt" in all_pairs.columns:
        score_col = None
        for candidate in ["af_fpocket_drug_score", "af_fpocket_pocket_score"]:
            if candidate in all_pairs.columns:
                score_col = candidate
                break

        if score_col is not None:
            cases = all_pairs.copy()
            cases["af_pocket_mean_plddt"] = pd.to_numeric(cases["af_pocket_mean_plddt"], errors="coerce")
            cases[score_col] = pd.to_numeric(cases[score_col], errors="coerce")

            low_plddt_high_score = cases[
                (cases["af_pocket_mean_plddt"] < 70) & cases[score_col].notna()
            ].sort_values(score_col, ascending=False)

            low_plddt_high_score.to_csv(global_dir / "low_plddt_high_score_cases.csv", index=False)

    if "abs_delta_fpocket_drug_score" in all_pairs.columns:
        all_pairs.sort_values("abs_delta_fpocket_drug_score", ascending=False).to_csv(
            global_dir / "top_drug_score_disagreements.csv",
            index=False,
        )

    if "abs_delta_fpocket_pocket_score" in all_pairs.columns:
        all_pairs.sort_values("abs_delta_fpocket_pocket_score", ascending=False).to_csv(
            global_dir / "top_pocket_score_disagreements.csv",
            index=False,
        )


def main():
    args = parse_args()

    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    targets_path = Path(args.targets)
    if not targets_path.exists():
        raise FileNotFoundError(f"Targets CSV not found: {targets_path}")

    targets = pd.read_csv(targets_path)

    required_cols = {"PDB_ID", "AF_ID"}
    missing = required_cols - set(targets.columns)
    if missing:
        raise ValueError(f"Targets CSV must contain columns: {required_cols}. Missing: {missing}")

    all_pair_tables = []
    all_summary_tables = []

    for _, row in targets.iterrows():
        pair_df, summary_df = process_target_pair(row, args, out_root)
        all_pair_tables.append(pair_df)
        all_summary_tables.append(summary_df)

    if all_pair_tables:
        all_pairs = pd.concat(all_pair_tables, ignore_index=True)
    else:
        all_pairs = pd.DataFrame()

    if all_summary_tables:
        all_summaries = pd.concat(all_summary_tables, ignore_index=True)
    else:
        all_summaries = pd.DataFrame()

    make_global_outputs(all_pairs, all_summaries, out_root)

    print()
    print("Analysis completed.")
    print(f"Local outputs:  {out_root / 'local'}")
    print(f"Global outputs: {out_root / 'global'}")


if __name__ == "__main__":
    main()