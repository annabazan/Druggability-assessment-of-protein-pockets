# Pocket Detection Instructions

This folder contains scripts and instructions for detecting protein pockets using **Fpocket** and **P2Rank**. We provide scripts to process both **AlphaFold models** and **raw PDB structures**.  

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
chmod +x run_fpocket_raw.sh 
chmod +x run_fpocket_alpha_fold.sh 
./run_fpocket_raw.sh
./run_fpocket_alpha_fold.sh
```
3. You should have folders `pdb_raw_out` and `alpha_fold_out` with results.
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
chmod +x run_P2Rank_raw.sh 
chmod +x run_P2Rank_alpha_fold.sh 
./run_P2Rank_raw.sh
./run_P2Rank_alpha_fold.sh
```
4. You should have folders `pdb_raw_out` and `alpha_fold_out` with results.
