#!/bin/bash

INPUT_DIR="../../targets/filtered_pdb"
OUTPUT_DIR="./pdb_out"

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

echo "Processed $count PDB files. Results in: $OUTPUT_DIR"
