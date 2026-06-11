#!/bin/bash

INPUT_DIR="../../targets/3D_aligned_alpha_fold"
OUTPUT_DIR="./alpha_fold_out"

mkdir -p "$OUTPUT_DIR"

count=0
for pdb in "$INPUT_DIR"/*.pdb; do
    [ -e "$pdb" ] || { echo "No pdb files found in $INPUT_DIR"; break; }

    filename=$(basename "$pdb" .pdb)

    echo "Processing $filename..."

    fpocket -f "$pdb"

    if [ -d "$INPUT_DIR/${filename}_out" ]; then
        mv "$INPUT_DIR/${filename}_out" "$OUTPUT_DIR/"
    fi
    count=$((count+1))
done

echo "Processed $count AlphaFold files. Results in: $OUTPUT_DIR"