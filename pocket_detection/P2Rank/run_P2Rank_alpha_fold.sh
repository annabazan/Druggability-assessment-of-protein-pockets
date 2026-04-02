#!/bin/bash

INPUT_DIR="../../targets/alpha_fold"
OUTPUT_DIR="./alpha_fold_out"

mkdir -p "$OUTPUT_DIR"

for pdb in "$INPUT_DIR"/*.pdb; do
    [ -e "$pdb" ] || { echo "No pdb files found"; break; }

    filename=$(basename "$pdb" .pdb)

    echo "Processing $filename..."

    prank predict -c alphafold -f "$pdb" -o "$OUTPUT_DIR/$filename"

done

echo "Done."