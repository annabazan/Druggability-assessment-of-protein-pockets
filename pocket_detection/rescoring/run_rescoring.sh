#!/usr/bin/env bash
set -euo pipefail

# Rescore Fpocket-detected pockets with P2Rank/PRANK.
#
# This script assumes that Fpocket has already been run and that the following
# folders exist:
#
#   pocket_detection/fpocket/pdb_out
#   pocket_detection/fpocket/alpha_fold_out
#
# It builds PRANK/P2Rank rescoring dataset files and runs:
#
#   prank rescore ...
#
# for both PDB and AlphaFold structures.

THREADS=8
VISUALIZATIONS=0
AF_CONFIG="rescore_2024"
PRANK_BIN="${PRANK_BIN:-$HOME/tools/p2rank_2.5.1/prank}"

usage() {
    cat <<EOF_USAGE
Usage:
  bash pocket_detection/rescoring/run_rescoring.sh [options]

Options:
  --threads N              Number of threads for P2Rank/PRANK. Default: 8.
  --visualizations 0|1     Generate P2Rank visualization files. Default: 0.
  --prank-bin PATH         Path to P2Rank's prank binary.
                           Default: \$HOME/tools/p2rank_2.5.1/prank.
  --af-config NAME         Config for AlphaFold rescoring. Default: rescore_2024.
  --no-af-config           Use default rescoring config also for AlphaFold.
  -h, --help               Show this help.

Examples:
  bash pocket_detection/rescoring/run_rescoring.sh

  bash pocket_detection/rescoring/run_rescoring.sh --threads 8 --visualizations 1

  bash pocket_detection/rescoring/run_rescoring.sh --no-af-config
EOF_USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --visualizations)
            VISUALIZATIONS="$2"
            shift 2
            ;;
        --prank-bin)
            PRANK_BIN="$2"
            shift 2
            ;;
        --af-config)
            AF_CONFIG="$2"
            shift 2
            ;;
        --no-af-config)
            AF_CONFIG=""
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

PDB_INPUT="$REPO_ROOT/targets/pdb"
AF_INPUT="$REPO_ROOT/targets/alpha_fold"

FPOCKET_DIR="$REPO_ROOT/pocket_detection/fpocket"
RESCORING_DIR="$REPO_ROOT/pocket_detection/rescoring"

PDB_FPOCKET_OUT="$FPOCKET_DIR/pdb_out"
AF_FPOCKET_OUT="$FPOCKET_DIR/alpha_fold_out"

PDB_DS="$RESCORING_DIR/fpocket_pdb_rescore.ds"
AF_DS="$RESCORING_DIR/fpocket_alpha_fold_rescore.ds"

PDB_RESCORE_OUT="$RESCORING_DIR/fpocket_pdb_rescored_out"
AF_RESCORE_OUT="$RESCORING_DIR/fpocket_alpha_fold_rescored_out"

check_tools() {
    echo "=== Tool check ==="

    command -v java >/dev/null 2>&1 || {
        echo "ERROR: java not found."
        echo "Install Java, e.g.: sudo apt install -y openjdk-17-jre"
        exit 1
    }

    if [[ ! -x "$PRANK_BIN" ]]; then
        echo "ERROR: P2Rank prank binary not found or not executable:"
        echo "  $PRANK_BIN"
        echo
        echo "Use --prank-bin PATH or set PRANK_BIN."
        echo "Example:"
        echo "  bash pocket_detection/rescoring/run_rescoring.sh --prank-bin ~/tools/p2rank_2.5.1/prank"
        exit 1
    fi

    echo "java:  $(command -v java)"
    echo "prank: $PRANK_BIN"
    echo
}

build_rescore_ds() {
    local fpocket_out="$1"
    local input_dir="$2"
    local ds_file="$3"
    local label="$4"

    echo "=== Building rescoring dataset for $label ==="
    echo "Fpocket output: $fpocket_out"
    echo "Original PDBs:  $input_dir"
    echo "Dataset file:   $ds_file"

    if [[ ! -d "$fpocket_out" ]]; then
        echo "ERROR: fpocket output directory does not exist:"
        echo "  $fpocket_out"
        echo
        echo "Run Fpocket first."
        exit 1
    fi

    if [[ ! -d "$input_dir" ]]; then
        echo "ERROR: input structure directory does not exist:"
        echo "  $input_dir"
        exit 1
    fi

    cat > "$ds_file" <<EOF_DS
PARAM.PREDICTION_METHOD=fpocket
HEADER: prediction protein
EOF_DS

    local count=0

    while IFS= read -r -d '' pred; do
        local base
        base="$(basename "$pred" _out.pdb)"

        local protein="$input_dir/${base}.pdb"

        if [[ -f "$protein" ]]; then
            printf "%s %s\n" "$pred" "$protein" >> "$ds_file"
            count=$((count + 1))
        else
            echo "WARNING: missing original protein for:"
            echo "  prediction: $pred"
            echo "  expected:   $protein"
        fi
    done < <(find "$fpocket_out" -mindepth 2 -maxdepth 2 -name "*_out.pdb" -print0 | sort -z)

    if [[ "$count" -eq 0 ]]; then
        echo "ERROR: no valid fpocket prediction/protein pairs found for $label."
        echo
        echo "Expected pattern:"
        echo "  $fpocket_out/*_out/*_out.pdb"
        echo
        echo "Example:"
        echo "  $fpocket_out/1AAX_out/1AAX_out.pdb"
        exit 1
    fi

    echo "Pairs written: $count"
    echo "Preview:"
    head "$ds_file"
    echo
}

run_rescore() {
    local ds_file="$1"
    local output_dir="$2"
    local label="$3"
    local config="${4:-}"

    echo "=== Running P2Rank/PRANK rescoring for $label ==="
    echo "Dataset:        $ds_file"
    echo "Output:         $output_dir"
    echo "Threads:        $THREADS"
    echo "Visualizations: $VISUALIZATIONS"
    if [[ -n "$config" ]]; then
        echo "Config:         $config"
    else
        echo "Config:         default"
    fi

    rm -rf "$output_dir"

    local cmd=("$PRANK_BIN" rescore "$ds_file")

    if [[ -n "$config" ]]; then
        cmd+=(-c "$config")
    fi

    cmd+=(
        -threads "$THREADS"
        -visualizations "$VISUALIZATIONS"
        -o "$output_dir"
    )

    echo "Command:"
    printf ' %q' "${cmd[@]}"
    echo
    echo

    if ! "${cmd[@]}"; then
        if [[ "$label" == "AlphaFold" && -n "$config" ]]; then
            echo
            echo "WARNING: AlphaFold rescoring with config '$config' failed."
            echo "Retrying AlphaFold rescoring with default config..."
            echo

            rm -rf "$output_dir"

            "$PRANK_BIN" rescore "$ds_file" \
                -threads "$THREADS" \
                -visualizations "$VISUALIZATIONS" \
                -o "$output_dir"
        else
            echo "ERROR: rescoring failed for $label"
            exit 1
        fi
    fi

    echo "Done rescoring for $label"
    echo
}

main() {
    echo "=== Fpocket output rescoring with P2Rank/PRANK ==="
    echo "Repo root:       $REPO_ROOT"
    echo "Threads:         $THREADS"
    echo "Visualizations:  $VISUALIZATIONS"
    echo "PRANK binary:    $PRANK_BIN"
    echo

    check_tools

    mkdir -p "$RESCORING_DIR"

    build_rescore_ds "$PDB_FPOCKET_OUT" "$PDB_INPUT" "$PDB_DS" "PDB"
    build_rescore_ds "$AF_FPOCKET_OUT" "$AF_INPUT" "$AF_DS" "AlphaFold"

    run_rescore "$PDB_DS" "$PDB_RESCORE_OUT" "PDB" ""
    run_rescore "$AF_DS" "$AF_RESCORE_OUT" "AlphaFold" "$AF_CONFIG"

    echo "=== Final outputs ==="
    echo "PDB rescored output:"
    echo "  $PDB_RESCORE_OUT"
    echo
    echo "AlphaFold rescored output:"
    echo "  $AF_RESCORE_OUT"
    echo

    echo "=== CSV files ==="
    find "$RESCORING_DIR" -name "*.csv" | head -40
    echo

    echo "Rescoring pipeline completed."
}

main "$@"
