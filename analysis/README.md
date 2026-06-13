## PDB–AlphaFold pocket comparison

This step compares and analyses pockets detected by Fpocket in experimental PDB structures
and corresponding AlphaFold models.


The main script is:

```text
analysis/pocket_comparison/run_analysis.py
```

### Input requirements

Before running the analysis, make sure that Fpocket has already been run for both PDB
and AlphaFold structures.

Expected Fpocket output folders:

```text
pocket_detection/fpocket/pdb_out
pocket_detection/fpocket/alpha_fold_out
```

Expected target list:

```text
targets/targets_list.csv
```

The target list should contain at least:

```text
PDB_ID
AF_ID
```

The script also uses full PDB and AlphaFold structures to compute pocket RMSD and
AlphaFold pLDDT statistics.

Example structure folders:

```text
targets/PBD_raw
targets/3D_alignment/3D_aligned_alpha_fold
```

### Running the analysis

From the repository root, run:

```bash
python analysis/pocket_comparison/run_analysis.py \
  --pdb-structures targets/PBD_raw \
  --af-structures targets/3D_alignment/3D_aligned_alpha_fold
```

If your PDB structures are stored in a different folder, replace `targets/PBD_raw`
with the correct path.

Optional arguments:

```text
--targets             Path to targets_list.csv.
--pdb-fpocket         Path to Fpocket outputs for PDB structures.
--af-fpocket          Path to Fpocket outputs for AlphaFold structures.
--pdb-structures      Path to full PDB structures.
--af-structures       Path to full AlphaFold structures.
--out                 Output directory.
--match-threshold     Jaccard threshold for accepted matches. Default: 0.30.
--weak-threshold      Jaccard threshold for weak matches. Default: 0.10.
--residue-mode        Residue matching mode: auto, chain, or number.
```

Example with explicit thresholds:

```bash
python3 analysis/pocket_comparison/run_analysis.py \
  --pdb-structures targets/PBD_raw \
  --af-structures targets/3D_alignment/3D_aligned_alpha_fold \
  --match-threshold 0.30 \
  --weak-threshold 0.10
```

### Specification

For every `PDB_ID`–`AF_ID` pair, the script:

1. reads Fpocket pockets from PDB and AlphaFold outputs;
2. extracts pocket residues from `pocket*_atm.pdb` files;
3. compares every PDB pocket with every AlphaFold pocket using the Jaccard index;
4. assigns unique PDB–AlphaFold pocket matches;
5. labels pocket pairs as `matched`, `weak_match`, `pdb_only`, or `af_only`;
6. computes RMSD for common C-alpha atoms of matched pockets;
7. computes AlphaFold pocket pLDDT statistics from the B-factor field;
8. parses selected Fpocket descriptors, including `Pocket Score` and `Drug Score`;
9. writes local and global CSV output files.

Pocket similarity is computed as:

```text
Jaccard = shared pocket residues / union of pocket residues
```

Pocket pair status:

```text
matched      Jaccard >= match-threshold
weak_match   weak-threshold <= Jaccard < match-threshold
pdb_only     PDB pocket without an AlphaFold match
af_only      AlphaFold pocket without a PDB match
```

### Output files

The script writes results to:

```text
analysis/pocket_comparison/outputs
```

Local outputs are created separately for every PDB–AlphaFold pair (one for every protein):

```text
analysis/pocket_comparison/outputs/local/<PDB_ID>_vs_<AF_ID>/
```

Each local folder contains:

```text
pocket_pairs_detailed.csv
jaccard_matrix.csv
protein_summary.csv
```

Global outputs are written to:

```text
analysis/pocket_comparison/outputs/global/
```

Main global files:

```text
global_pocket_pairs.csv
global_protein_summary.csv
global_statistics.csv
low_plddt_high_score_cases.csv
top_drug_score_disagreements.csv
top_pocket_score_disagreements.csv
```

### Main output columns

`global_pocket_pairs.csv` contains pocket-level results.

Important columns:

```text
target_id
pdb_id
af_id
pair_status
residue_mode_used
pdb_pocket_id
af_pocket_id
jaccard
shared_residues
pdb_residue_count
af_residue_count
rmsd_common_ca_global
rmsd_common_ca_local_aligned
af_pocket_mean_plddt
af_pocket_median_plddt
af_pocket_min_plddt
af_pocket_fraction_plddt_lt_50
pdb_fpocket_pocket_score
af_fpocket_pocket_score
delta_fpocket_pocket_score
abs_delta_fpocket_pocket_score
pdb_fpocket_drug_score
af_fpocket_drug_score
delta_fpocket_drug_score
abs_delta_fpocket_drug_score
```

Column interpretation:

```text
jaccard
```

Residue-level overlap between a PDB pocket and an AlphaFold pocket.

```text
rmsd_common_ca_global
```

RMSD of common pocket C-alpha atoms in the current/global coordinate frame.

```text
rmsd_common_ca_local_aligned
```

Pocket-level RMSD after local Kabsch alignment of common C-alpha atoms.

```text
af_pocket_mean_plddt
```

Mean AlphaFold confidence score for residues in the AlphaFold pocket.

```text
delta_fpocket_drug_score
```

Difference between AlphaFold and PDB Fpocket Drug Score:

```text
AF Drug Score - PDB Drug Score
```

```text
abs_delta_fpocket_drug_score
```

Absolute difference in Fpocket Drug Score between matched PDB and AlphaFold pockets.

### Notes

The script compares Fpocket pockets only. P2Rank rescoring output is not required for
this analysis step.

If the script prints warnings such as:

````text
[WARN] Missing AF structure
````

the pocket overlap can still be computed from Fpocket pocket files, but RMSD and pLDDT
values for that structure may be missing.

If `residue_mode_used` is `number`, the script ignored chain identifiers when matching
residues. This usually means that chain IDs differ between the PDB and AlphaFold files,
so these cases should be interpreted more carefully.


After running `run_analysis.py`, more clear and readable summary tables can be generated with:

```bash
python analysis/pocket_comparison/make_tables.py

