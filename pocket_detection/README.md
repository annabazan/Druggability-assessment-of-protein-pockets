# Pocket Detection Instructions

This folder contains scripts and instructions for detecting protein pockets using **Fpocket** and **P2Rank**. We provide scripts to process both **AlphaFold models** and **PDB structures**.  

---

## Fpocket

1. Clone Fpocket (outside of this repo):
 
```bash
git clone https://github.com/Discngine/fpocket.git .
cd fpocket
sudo make install
```
2. In `pocket_detection/fpocket`:
```bash
chmod +x run_fpocket_pdb.sh 
chmod +x run_fpocket_alpha_fold.sh 
./run_fpocket_pdb.sh
./run_fpocket_alpha_fold.sh
```
3. You should have folders `pdb_out` and `alpha_fold_out` with results.
---

## P2Rank

1. Download P2Rank (outside of this repo), you already should have java:
 
```bash
wget https://github.com/rdk/p2rank/releases/latest/download/p2rank_2.5.1.tar.gz
tar -xzf p2rank_2.5.1.tar.gz
cd p2rank_2.5.1
chmod +x prank
```
2. To use P2Rank globally:

```bash
echo 'export PATH=$PATH:/full/path/to/p2rank_2.5.1' >> ~/.bashrc
source ~/.bashrc
```
3. In `pocket_detection/P2Rank`:
```bash
chmod +x run_P2Rank_pdb.sh 
chmod +x run_P2Rank_alpha_fold.sh 
./run_P2Rank_pdb.sh
./run_P2Rank_alpha_fold.sh
```
4. You should have folders `pdb_out` and `alpha_fold_out` with results.

# Comparison of pockets

Script `comparison.py` generates comparison results for each structure pair and visualisation for every matched pocket pair.

---

## P2Rank rescoring of Fpocket pockets

This step rescoring pockets already detected by **Fpocket** using **P2Rank/PRANK**.  
It does not detect new pockets. It only assigns a new P2Rank/PRANK score and ranking to existing Fpocket pockets.

1. Make sure that Fpocket has already been run and that the following folders exist:

```text
pocket_detection/fpocket/pdb_out
pocket_detection/fpocket/alpha_fold_out
```

2. In `pocket_detection/rescoring`:

```bash
bash run_rescoring.sh
```

3. You should have folders:

```text
fpocket_pdb_rescored_out
fpocket_alpha_fold_rescored_out
```

and dataset files:

```text
fpocket_pdb_rescore.ds
fpocket_alpha_fold_rescore.ds
```

The main output files are:

```text
<protein>.pdb_rescored.csv
<protein>.pdb_predictions.csv
```

The `*_rescored.csv` files contain the original Fpocket pockets with updated P2Rank/PRANK scores and ranks. Compare `old_rank` with `rank` to see how rescoring changed the original Fpocket ranking.