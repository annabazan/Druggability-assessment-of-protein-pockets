#!/bin/bash

INPUT_DIR="../../targets/pdb_raw"
OUTPUT_DIR="./pdb_raw_out"

mkdir -p "$OUTPUT_DIR"

for pdb in "$INPUT_DIR"/*.pdb; do
    [ -e "$pdb" ] || { echo "No pdb files found"; break; }

    filename=$(basename "$pdb" .pdb)

    echo "Processing $filename..."

    prank predict -f "$pdb" -o "$OUTPUT_DIR/$filename"

done

echo "Done."