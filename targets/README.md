## Target Prepariation Instructions

This module implements a **preprocessing and validation pipeline** for protein structures, combining experimental data from the **Protein Data Bank (PDB)** with predicted models from the **AlphaFold Protein Structure Database**.

The goal is to ensure that both structure sources are **consistent, comparable, restricted to equivalent sequence regions, and where appropriate superimposed in 3D** prior to downstream analysis.

**The pipeline includes:**
- **automated** structure retrieval,
- structure **cleaning and standardization**,
- sequence-based **trimming** and/or **pLDDT filtering** of AlphaFold models,
- **validation of sequence consistency** between PDB and AlphaFold,
- **3D superposition of AlphaFold models onto their experimental counterparts** using common C-alpha residues.

## How to run

Ensure required Python dependencies are installed (`pandas`, `biopython`, `requests`).

The full **target preparation pipeline** should be executed from the `targets/` directory:

1. **Download structures**

   Download experimental structures from PDB: `python pdb_download.py --targets-file targets_list.csv --output-dir pdb`
   
   Download AlphaFold models: `python af_download.py --targets-file targets_list.csv --output-dir alpha_fold`

   Both download scripts are quiet by default. Add `--loud` for detailed per-target progress.

2. **Preprocess structures**

    Clean PDB files (chain selection, remove non-protein residues): `python clear_pdb.py`

    Trim AlphaFold models to selected sequence ranges: `python filter_af.py --mode range`

3. *Optional:* **Validate sequence consistency**

    Run sequence comparison between PDB and AlphaFold:

    `python sequence_compare.py --pdb_dir filtered_pdb --af_dir cut_alpha_fold`

4. **3D alignment**

    Align trimmed AlphaFold models to their corresponding experimental PDB structures:

    `python 3D_align_superimposer.py [--visualize] [--loud] [--direct-numbering]`

    The core alignment is always performed; optional PNG visualizations are generated only when `--visualize` is provided and PyMOL is available. The script uses sequence-based C-alpha matching by default and can overcome residue numbering differences. Use `--direct-numbering` only when the experimental and AlphaFold residue numbers are already aligned.

After running the pipeline, **the following directories will be created**:

- `pdb/` – downloaded experimental structures
- `alpha_fold/` – downloaded AlphaFold models
- `filtered_pdb/` – cleaned PDB structures
- `cut_alpha_fold/` (*default*) – trimmed AlphaFold structures (*in other modes:* `filtered_alpha_fold/` or `filtered_cut_alpha_fold/`)
- `alignment_results/` – sequence alignment reports
- `/3D_aligned_alpha_fold/` - rotated AlphaFold models
- `/3D_alignment_visualisations/` - PNG visualizations of aligned pairs (*optional:* when `--visualize` is enabled)


## File descriptions

1) `targets_list.csv`

    A curated list of protein targets used in the project, linking experimental structures from the Protein Data Bank (PDB) with corresponding AlphaFold models.

    The file contains **metadata required for consistent preprocessing, filtering, and comparison of protein structures**:

    * **record identifier** (`NR`) – unique target index used throughout the project
    * **PDB structure ID** (`PDB_ID`) – identifier of the experimentally determined structure
    * **AlphaFold ID** (`AF_ID`) – corresponding UniProt/AlphaFold DB entry
    * **sequence range** (`START`–`END`) – defines the fragment of the AlphaFold model used in the analysis
    * **chain** (`CHAIN`) – specifies which chain from the experimental structure is considered (`-` indicates that all chains are used)
    * **protein class** (`CLASS`) – short description of the protein family, function, or structural category
    * **group label** (`GROUP`) – target classification used in the study:
        - `1` – likely to contain druggable small-molecule binding pockets
        - `0` – control group; likely to lack druggable pockets or contain only small and/or poorly accessible pockets
    * **notes** (`NOTES`) – optional comments describing notable structural features, sequence discrepancies, insertions, or other observations relevant to downstream analysis

    This table serves as the **central reference for all preprocessing, filtering, validation, and comparative analysis steps**, ensuring that experimental PDB structures and AlphaFold models are processed consistently and remain directly comparable throughout the project.

    ---

2) `pdb_download.py`

    A utility script for **automated downloading** of protein structures from the **Protein Data Bank (PDB)** based on entries listed in `targets_list.csv`.

    All downloaded structures are saved in the `pdb/` directory.

    ---

3) `af_download.py`

    A utility script for **automated downloading** of protein structure models from the **AlphaFold Protein Structure Database** based on entries listed in `targets_list.csv`.

    All downloaded structures are saved in the `alpha_fold/` directory.

    ---

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

    ---

5) `filter_af.py`

    A preprocessing script for filtering AlphaFold structures based on:

    - sequence range (`START`–`END`)
    - per-residue confidence score (pLDDT)
    - or a combination of both

    **For each AlphaFold model, the script:**

    - selects residues within a specified sequence range (`range` mode),
    - computes average pLDDT per residue (from B-factor field),
    - removes residues below a specified confidence threshold (`plddt` mode),
    - preserves only standard amino acids (`ATOM` records),
    - rebuilds the `SEQRES` section to match the filtered structure,
    - optionally generates PyMOL visualizations comparing full and filtered structures.

    **Filtering modes:**
    - `range` – sequence-based trimming
    - `plddt` – confidence-based filtering
    - `complete` – combined filtering

    **Output directories:**
    - `cut_alpha_fold/`
    - `filtered_alpha_fold/`
    - `filtered_cut_alpha_fold/`

    **Usage:**

    `python filter_af.py --mode <mode> [--plddt <threshold>] [--visualize]`

    ---

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
    - detailed alignment reports and sequence summary saved in `alignment_results/`
    - summary statistics printed to stdout

    ---

7. `3D_alignment_superimposer.py`

    A structural alignment script for **superimposing filtered AlphaFold models onto their corresponding experimental PDB structures**.

    **Usage:**

    `python 3D_align_superimposer.py [--visualize] [--loud] [--direct-numbering]`

    - `--visualize` creates optional PyMOL PNG visualizations when PyMOL is installed.
    - `--loud` enables detailed alignment diagnostics and per-target status messages.
    - `--direct-numbering` attempts direct residue-number-based Cα matching; otherwise the script uses sequence-based matching by default.

    The script saves a summary report to `rmsd_results.csv` with the following columns:

    - `PDB_ID`
    - `AF_ID`
    - `ref_ca_count`
    - `sample_ca_count`
    - `matched_ca_count`
    - `method`
    - `RMSD`

    `ref_ca_count` and `sample_ca_count` are the numbers of Cα atoms extracted from the reference and AlphaFold structures, `matched_ca_count` is the number of Cα atoms used for superposition, and `method` indicates whether direct numbering or sequence alignment was used.

    **For each structure pair, the script:**

    * identifies matching residues shared between the experimental and predicted structures,
    * extracts corresponding Cα atoms,
    * computes an optimal rigid-body superposition using Biopython's `Superimposer`,
    * applies the calculated transformation to the AlphaFold model,
    * optionally generates PyMOL visualizations of the aligned structures when `--visualize` is used.
    * if PyMOL is not installed or visualization fails, the script still performs the core superposition and saves rotated AlphaFold models.

    The resulting aligned AlphaFold structures are placed in the same coordinate frame as their experimental counterparts, enabling direct structural comparison and downstream analyses.

    **Output directories:**

    * `3D_aligned_alpha_fold/` – transformed AlphaFold structures
    * `3D_alignment_visualisations/` – optional alignment visualizations (created only when `--visualize` is enabled)

    **Visualization convention:**

    * experimental structure – *violet*
    * AlphaFold model – *orange*

