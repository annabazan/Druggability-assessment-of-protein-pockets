## Introduction

This module implements a **preprocessing and validation pipeline** for protein structures, combining experimental data from the **Protein Data Bank (PDB)** with predicted models from the **AlphaFold Protein Structure Database**.

The goal is to ensure that both structure sources are **consistent, comparable, and restricted to equivalent sequence regions** prior to downstream analysis.

**The pipeline includes:**
- **automated** structure retrieval,
- structure **cleaning and standardization**,
- sequence-based **trimming of AlphaFold models**,
- **validation of sequence consistency** between PDB and AlphaFold.

## How to run

Ensure required Python dependencies are installed (`pandas`, `biopython`, `requests`).

The full **target preparation pipeline** should be executed from the `targets/` directory:

1. **Download structures**

   Download experimental structures from PDB: `python pdb_download.py`
   
   Download AlphaFold models: `python af_download.py`

2. **Preprocess structures**

    Clean PDB files (chain selection, remove non-protein residues): `python clear_pdb.py`

    Trim AlphaFold models to selected sequence ranges: `python cut_af.py`

3. **Validate sequence consistency**

    Run sequence comparison between PDB and AlphaFold:

    `python sequence_compare.py --pdb_dir filtered_pdb --af_dir cut_alpha_fold`


After running the pipeline, **the following directories will be created**:

- `pdb/` – downloaded experimental structures
- `alpha_fold/` – downloaded AlphaFold models
- `filtered_pdb/` – cleaned PDB structures
- `cut_alpha_fold/` – trimmed AlphaFold structures
- `alignment_results/` – sequence alignment reports


## File descriptions

1) `targets_list.csv`

    A curated list of protein targets used in the project, linking experimental structures from the Protein Data Bank (PDB) with corresponding AlphaFold models.

    The file contains **additional metadata** required for consistent preprocessing and comparison between structures:
    - **PDB structure ID** (`PDB_ID`) – identifier of the experimentally determined structure
    - **AlphaFold ID** (`AF_ID`) – corresponding UniProt/AlphaFold DB entry
    - **sequence range** (`START`-`END`) – defines the fragment of the AlphaFold model used in the analysis
    - **chain** (`CHAIN`) – specifies which chain from the experimental structure is considered
    - **protein class** (`CLASS`) – short description of the protein family or functional class

    This table serves as the **central reference for all downstream processing steps**, ensuring that both structure sources (PDB and AlphaFold) are aligned and comparable.

2) `pdb_download.py`

    A utility script for **automated downloading** of protein structures from the **Protein Data Bank (PDB)** based on entries listed in `targets_list.csv`.

    All downloaded structures are saved in the `pdb/` directory.

3) `af_download.py`

    A utility script for **automated downloading** of protein structure models from the **AlphaFold Protein Structure Database** based on entries listed in `targets_list.csv`.

    All downloaded structures are saved in the `alpha_fold/` directory.

4) `clear_pdb.py`

    A preprocessing script for **cleaning and standardizing PDB structures** based on definitions provided in `targets_list.csv`.

    **For each structure, the script:**
    - selects a **specific chain** (or all chains if `CHAIN` = '-'),
    - removes non-protein residues (e.g. water, ligands, heteroatoms),
    - preserves **only standard amino acids**,
    - filters `SEQRES` records to match the selected chain,
    - rewrites the structure into a clean, consistent PDB format.

    The resulting files contain **only relevant protein atoms and aligned sequence information**, making them suitable for downstream structural analysis. 
    
    Cleaned PDB files are saved in `filtered_pdb/` directory.

5) `cut_af.py`

    A preprocessing script for **trimming AlphaFold structures** to specific sequence ranges defined in `targets_list.csv`.

    **For each AlphaFold model, the script:**

    - selects only **residues within the specified sequence range** (`START`–`END`),
    - removes non-protein residues,
    - preserves **only standard amino acids** (`ATOM` records)
    - rebuilds the `SEQRES` section to match the trimmed structure,
    - outputs a clean and consistent PDB file.

    This ensures that AlphaFold models correspond to **the same sequence fragments** as their experimental counterparts.

    Trimmed PDB files are saved in `cut_alpha_fold/` directory.

6) `sequence_compare.py`

    A validation script for **comparing protein sequences** derived from experimental PDB structures and AlphaFold models.

    **For each target, the script:**

    - **extracts amino acid sequences** from PDB files,
    - performs a **substring check** (whether the PDB sequence is contained within the AlphaFold sequence),
    - runs a **local pairwise alignment** using the `BLOSUM62` substitution matrix,
    - calculates:
        - **sequence identity** (fraction of matching residues),
        - **coverage** (fraction of the PDB sequence aligned),
    - **classifies results** into three groups based on identity:
        - `low`: identity < **0.8**
        - `medium`: identity between **0.8** and **0.95**
        - `high`: identity > **0.95**,
    - **saves detailed alignment** outputs to corresponding directories.

    The script provides **a quantitative assessment of consistency** between experimental and predicted structures.

    **Usage:**

    `python sequence_compare.py --pdb_dir <pdb_dir> --af_dir <alpha_fold_dir>`

    **Output:**
    - alignment reports saved in:
        - `alignment_results/low_identity/`
        - `alignment_results/med_identity/`
        - `alignment_results/high_identity/`
    - summary statistics printed to stdout

    7) `3D_alignment/3D_alignment_superimposer.py`

    A 3d aligning script for **cut AlphaFold structures**.

    **For each AlphaFold model, the script:**

    - generates rotated pdb file in `/3D_aligned_alpha_fold/` using common residues and Superimposer
    - generates visualisations for the rotated version of the alpha fold structure (orange color) and experimental structure (violet color) in `/3D_alignment_visualisations/`

    This ensures that AlphaFold models position corresponds to their experimental counterpart.
