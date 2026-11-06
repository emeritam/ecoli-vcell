#!/usr/bin/env bash
# Usage: ./run_all.sh path/to/abundance_file [extra step-2 options]
# Example: ./run_all.sh li2014_tableS1.xlsx --value-col "MOPS complete" --skiprows 1
set -euo pipefail
cd "$(dirname "$0")"
[ $# -ge 1 ] || { echo "Give the abundance file as the first argument."; exit 1; }
python step1_fetch_proteome.py
python step2_add_abundance.py "$@"
python step3_pick_structures.py --top "${TOP:-0}"
python step4_download_structures.py
python step5_coarse_grain.py
