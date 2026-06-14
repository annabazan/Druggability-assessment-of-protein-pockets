# Druggability assessment of protein pockets
Group project for the Jagiellonian University course **Machine Learning in Drug Design** in 2026, aimed at evaluating **protein pocket druggability** using **AlphaFold** structures.

---

#### Participants:
* Anna Bazan
* Julia Borowska
* Kacper Drozd

--- 

### Reproducing the Project Results

To reproduce the results presented in this project:

1. Ensure that both `fpocket` and `P2Rank` are installed and available in your environment.
2. Run the project pipeline:

   * Using the original `fpocket` scores:

     ```bash
     python run_pipeline.py [--visualize]
     ```
   * Using pocket scores rescored with `P2Rank`:

     ```bash
     python run_pipeline.py [--visualize] --score-source rescored
     ```

**Note:** The optional `--visualize` flag requires access to **PyMOL** for generating pocket visualizations.

---

### Result Analysis

The generated results can be further explored using the provided Jupyter notebooks:

* **`analysis.ipynb`** – comprehensive analysis of the obtained results, including various aspects of pocket detection and comparison between experimental and AlphaFold structures.

* **`notebook.ipynb`** – detailed inspection of individual protein pairs, with a particular focus on pocket-level analysis and visualization.
