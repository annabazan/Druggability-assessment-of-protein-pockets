#!/bin/bash

INPUT_DIR="../../targets/filtered_pdb"
OUTPUT_DIR="./pdb_out"

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
