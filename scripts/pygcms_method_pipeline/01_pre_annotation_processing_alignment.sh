#!/usr/bin/env bash
# =============================================================================
# Purpose: Process CDF chromatograms through sample grouping, peak alignment,
# and reference-based retention index calculation.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${SCRIPT_DIR}/config.sh"

echo "Dataset: ${DATASET_ID}"
echo "Input CDF dir: ${CDF_DIR}"
echo "Run dir: ${RUN_DIR}"

GROUP_DIR="${WORK_DIR}/04_global_grouping"

mkdir -p "${WORK_DIR}/00_parameters" "${ALIGNMENT_DIR}" "${ALIGNMENT_PREBLOCK_DIR}" "${GROUP_DIR}" "${MAIN_OUTPUT_DIR}"

PROJECT_DIR="${PROJECT_DIR}" \
PIPELINE_DIR="${SCRIPT_DIR}/lib/01_preprocessing" \
CDF_DIR="${CDF_DIR}" \
OUT_DIR="${ALIGNMENT_DIR}" \
PREBLOCK_DIR="${ALIGNMENT_PREBLOCK_DIR}" \
ALIGNMENT_MODE="${ALIGNMENT_MODE}" \
ALIGNMENT_SINGLE_BLOCK_THRESHOLD="${ALIGNMENT_SINGLE_BLOCK_THRESHOLD}" \
ALIGNMENT_SINGLE_BLOCK_TIME_DIST_SEC="${ALIGNMENT_SINGLE_BLOCK_TIME_DIST_SEC}" \
N_SUBSET="all" \
MIN_BLOCK_SIZE="${ALIGNMENT_MIN_BLOCK_SIZE}" \
MAX_BLOCK_SIZE="${ALIGNMENT_MAX_BLOCK_SIZE}" \
TARGET_MAX_TIME_DIST="${ALIGNMENT_TARGET_MAX_TIME_DIST}" \
USE_EXISTING_PREBLOCK="${ALIGNMENT_USE_EXISTING_PREBLOCK}" \
USE_EXISTING_BLOCK_FOLDERS="${ALIGNMENT_USE_EXISTING_BLOCK_FOLDERS}" \
ALIGN_TIME_DIST="${ALIGNMENT_ALIGN_TIME_DIST}" \
MIN_PEAK_WIDTH="${MIN_PEAK_WIDTH}" \
MIN_PEAK_HEIGHT="${MIN_PEAK_HEIGHT}" \
AREA_FRACTION_FILTER="${AREA_FRACTION_FILTER}" \
NOISE_THRESHOLD="${NOISE_THRESHOLD}" \
ANALYSIS_START_MIN="${ANALYSIS_START_MIN}" \
ANALYSIS_END_MIN="${ANALYSIS_END_MIN}" \
MIN_SPECTRA_COR="${ALIGNMENT_MIN_SPECTRA_COR}" \
MZ_MIN="${ALIGNMENT_MZ_MIN}" \
MZ_MAX="${ALIGNMENT_MZ_MAX}" \
AVOID_PROCESSING_MZ="${ALIGNMENT_AVOID_PROCESSING_MZ}" \
POST_ALIGNMENT_FEATURE_CONSOLIDATION="${ALIGNMENT_POST_CONSOLIDATION}" \
FINAL_FEATURE_MERGE_RT_SEC="${ALIGNMENT_FINAL_FEATURE_MERGE_RT_SEC}" \
FINAL_FEATURE_MERGE_SPECTRAL_COS="${ALIGNMENT_FINAL_FEATURE_MERGE_SPECTRAL_COS}" \
FINAL_FEATURE_MERGE_BASE_MZ_TOL="${ALIGNMENT_FINAL_FEATURE_MERGE_BASE_MZ_TOL}" \
ALIGNID0_RESCUE_ENABLED="${ALIGNID0_RESCUE_ENABLED}" \
ALIGNID0_WITHIN_BLOCK_RT_SEC="${ALIGNID0_WITHIN_BLOCK_RT_SEC}" \
ALIGNID0_MIN_SPECTRAL_COSINE="${ALIGNID0_MIN_SPECTRAL_COSINE}" \
ALIGNID0_BASE_MZ_TOL="${ALIGNID0_BASE_MZ_TOL}" \
ALIGNID0_AMBIGUITY_MARGIN="${ALIGNID0_AMBIGUITY_MARGIN}" \
GLOBAL_MERGE_RT_SEC="${ALIGNMENT_GLOBAL_MERGE_RT_SEC}" \
GLOBAL_MERGE_SPECTRAL_COS="${ALIGNMENT_GLOBAL_MERGE_SPECTRAL_COS}" \
GLOBAL_MERGE_BASE_MZ_TOL="${ALIGNMENT_GLOBAL_MERGE_BASE_MZ_TOL}" \
RT_CORRECTION_METHOD="${ALIGNMENT_RT_CORRECTION_METHOD}" \
LANDMARK_IDENTITY_MEDIAN_SEC="${ALIGNMENT_LANDMARK_IDENTITY_MEDIAN_SEC}" \
LANDMARK_IDENTITY_MAX_SEC="${ALIGNMENT_LANDMARK_IDENTITY_MAX_SEC}" \
ALIGNMENT_LATE_ANCHOR_START_MIN="${ALIGNMENT_LATE_ANCHOR_START_MIN}" \
ALIGNMENT_LATE_ANCHOR_FULL_MIN="${ALIGNMENT_LATE_ANCHOR_FULL_MIN}" \
ALIGNMENT_LATE_ANCHOR_WARN_SEC="${ALIGNMENT_LATE_ANCHOR_WARN_SEC}" \
ALIGNMENT_LATE_ANCHOR_BLOCK_SEC="${ALIGNMENT_LATE_ANCHOR_BLOCK_SEC}" \
BLOCKWISE_OUT_DIR="${ALIGNMENT_DIR}" \
Rscript "${SCRIPT_DIR}/lib/01_preprocessing/00_run_alignment_pipeline.R"

selected_mode="$(<"${ALIGNMENT_PREBLOCK_DIR}/alignment_mode.txt")"

"${PYTHON_BIN}" "${SCRIPT_DIR}/lib/01_preprocessing/04_export_alignment_for_annotation.py" \
  --alignment-dir "${ALIGNMENT_DIR}" \
  --out-dir "${GROUP_DIR}" \
  --foundin-min "${FOUNDIN_MIN}"

ri_mapping_args=(
  --alignment-dir "${ALIGNMENT_DIR}"
  --output-dir "${REFERENCE_RI_OUTPUT_DIR}"
  --reference-is "${REFERENCE_RI_INTERNAL_STANDARD_RT_FILE}"
  --ri-ladder "${ALKANE_RI_FILE}"
  --foundin-min "${FOUNDIN_MIN}"
)
if [[ -n "${REFERENCE_RI_EXPECTED_BLOCK}" ]]; then
  ri_mapping_args+=(--expected-reference-block "${REFERENCE_RI_EXPECTED_BLOCK}")
fi
"${PYTHON_BIN}" "${SCRIPT_DIR}/lib/01_preprocessing/05_reference_ri_mapping.py" "${ri_mapping_args[@]}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/lib/shared/processing_report.py" --run-dir "${WORK_DIR}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/lib/shared/export_readable_outputs.py" --run-dir "${WORK_DIR}" --output-dir "${MAIN_OUTPUT_DIR}" --foundin-min "${FOUNDIN_MIN}"

mkdir -p "${RUN_ROOT}"
printf "%s\n" "${DATASET_ID}" > "${RUN_ROOT}/.last_dataset_id"

echo
echo "Alignment and reference-RI preprocessing complete."
echo "Alignment mode: ${selected_mode}"
echo "RI-bearing annotation input:"
echo "${REFERENCE_RI_OUTPUT_DIR}/features_for_annotation.csv"
echo "Readable annotation input:"
echo "${MAIN_OUTPUT_DIR}/features_for_annotation.csv"
echo
echo "Next run:"
echo "bash ${SCRIPT_DIR}/02_compound_annotation_and_review.sh"
