#!/usr/bin/env bash
# =============================================================================
# Purpose: Run structural classification and export classified peak-area tables.
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

PROJECT_DIR="${PROJECT_DIR}" \
DATASET_ID="${DATASET_ID}" \
RUN_ROOT="${RUN_ROOT}" \
WORK_DIR="${WORK_DIR}" \
Rscript "${SCRIPT_DIR}/lib/04_classification/classification.R"

"${PYTHON_BIN}" "${SCRIPT_DIR}/lib/shared/processing_report.py" --run-dir "${WORK_DIR}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/lib/shared/export_readable_outputs.py" --run-dir "${WORK_DIR}" --output-dir "${MAIN_OUTPUT_DIR}" --foundin-min "${FOUNDIN_MIN}"
