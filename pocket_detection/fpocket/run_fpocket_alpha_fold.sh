#!/bin/bash

INPUT_DIR="../../targets/3D_alignment/3D_aligned_alpha_fold"
OUTPUT_DIR="./alpha_fold_out"

mkdir -p "$OUTPUT_DIR"

for pdb in "$INPUT_DIR"/*.pdb; do
[ -e "$pdb" ] || { echo "No pdb files found"; break; }

filename=$(basename "$pdb" .pdb)

echo "Processing $filename..."

fpocket -f "$pdb"

if [ -d "$INPUT_DIR/${filename}_out" ]; then
    mv "$INPUT_DIR/${filename}_out" "$OUTPUT_DIR/"
fi


done

echo "Done."
