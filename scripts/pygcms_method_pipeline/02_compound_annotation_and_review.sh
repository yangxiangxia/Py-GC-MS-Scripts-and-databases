#!/usr/bin/env bash
# =============================================================================
# Purpose: Generate putative compound annotations and finalize reviewed identities.
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

export DATASET_ID CDF_DIR RUN_DIR WORK_DIR MAIN_OUTPUT_DIR
export FOUNDIN_MIN COMPOUND_ID_PREFIX COMPOUND_ID_DIGITS POST_REVIEW_MERGE_RT_SEC
export POST_REVIEW_UNIDENTIFIED_COSINE POST_REVIEW_UNIDENTIFIED_MAX_JACCARD
export MONA_FILE NIST_FILE NIST_FILES TOP_HITS_PER_LIBRARY
export PRIMARY_TOP_HITS_PER_LIBRARY RESCUE_TOP_HITS_PER_LIBRARY
export ANNOTATION_DECISION_POLICY ANNOTATION_RECOMPUTE_FROM_CACHE
export ALKANE_RI_FILE SOM_REFERENCE_FILE RI_CACHE_DIR RI_SUPPORT_WINDOW RI_WEAK_WINDOW
export SPECTRAL_AUTO_THRESHOLD RI_QUERY_TIMEOUT_SEC RI_QUERY_DELAY_SEC RI_OFFLINE
export PUBCHEM_IDENTITY_CACHE_DIR PUBCHEM_QUERY_TIMEOUT_SEC PUBCHEM_QUERY_DELAY_SEC PUBCHEM_OFFLINE

echo "Dataset: ${DATASET_ID}"
echo "Input CDF dir: ${CDF_DIR}"
echo "Run dir: ${RUN_DIR}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/lib/02_annotation/compound_annotation_and_review.py"
