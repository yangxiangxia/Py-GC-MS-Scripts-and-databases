#!/usr/bin/env bash
# =============================================================================
# Purpose: Retrieve SMILES for reviewed identities using local libraries and PubChem.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${DATASET_ID+x}" ]]; then
  DATASET_ID_WAS_EXPLICIT="true"
else
  DATASET_ID_WAS_EXPLICIT="false"
fi
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/lib/shared/dataset_context.sh"
source "${SCRIPT_DIR}/config.sh"

echo "Dataset: ${DATASET_ID}"
echo "Input CDF dir: ${CDF_DIR}"
echo "Run dir: ${RUN_DIR}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/lib/shared/assert_context_matches.py" \
  --dataset-id "${DATASET_ID}" \
  --cdf-dir "${CDF_DIR}" \
  --work-dir "${WORK_DIR}"

MERGE_DIR="${WORK_DIR}/07_post_review_merge"
OUT_DIR="${WORK_DIR}/08_smiles"
mkdir -p "${OUT_DIR}"

# Keep local SMILES lookup order explicit: MoNA, then the configured NIST files.
SMILES_MSP_ARGS=(--msp "${MONA_FILE}")
IFS=';' read -r -a NIST_MSP_PATHS <<< "${NIST_FILES}"
for nist_msp in "${NIST_MSP_PATHS[@]}"; do
  SMILES_MSP_ARGS+=(--msp "${nist_msp}")
done

"${PYTHON_BIN}" "${SCRIPT_DIR}/lib/03_smiles/fill_smiles_by_compound_name.py" \
  --input "${MERGE_DIR}/final_compounds.csv" \
  --output "${OUT_DIR}/final_compounds_with_smiles.csv" \
  --name-column "final_identification" \
  --smiles-column "final_SMILES" \
  "${SMILES_MSP_ARGS[@]}" \
  --use-pubchem \
  --pubchem-cache "${OUT_DIR}/pubchem_smiles_cache.csv" \
  --pubchem-delay-sec "${PUBCHEM_DELAY_SEC}" \
  --pubchem-timeout-sec "${PUBCHEM_TIMEOUT_SEC}"

if [[ -f "${OUT_DIR}/final_compounds_with_smiles_source_summary.csv" ]]; then
  cp "${OUT_DIR}/final_compounds_with_smiles_source_summary.csv" "${OUT_DIR}/smiles_source_summary.csv"
fi
if [[ -f "${OUT_DIR}/final_compounds_with_smiles_summary.csv" ]]; then
  cp "${OUT_DIR}/final_compounds_with_smiles_summary.csv" "${OUT_DIR}/smiles_summary.csv"
fi

"${PYTHON_BIN}" "${SCRIPT_DIR}/lib/shared/processing_report.py" --run-dir "${WORK_DIR}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/lib/shared/export_readable_outputs.py" --run-dir "${WORK_DIR}" --output-dir "${MAIN_OUTPUT_DIR}" --foundin-min "${FOUNDIN_MIN}"
