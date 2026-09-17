#!/usr/bin/env bash
# =============================================================================
# Purpose: Continue the matching dataset context across pipeline stages.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# =============================================================================


# Source this after config.sh in downstream pipeline steps.
# If DATASET_ID was not explicitly set, continue from the most recent 01 route
# only when it belongs to the same source dataset as the current config.

PROJECT_DIR_DEFAULT="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)}"
RUN_ROOT_DEFAULT="${RUN_ROOT:-${PROJECT_DIR_DEFAULT}/pygcms_method_outputs}"
LAST_DATASET_FILE="${RUN_ROOT_DEFAULT}/.last_dataset_id"

strip_route_suffix() {
  local value="$1"
  value="${value%_alignment}"
  printf "%s\n" "${value}"
}

if [[ "${DATASET_ID_WAS_EXPLICIT:-false}" != "true" && -f "${LAST_DATASET_FILE}" ]]; then
  LAST_DATASET_ID="$(tr -d '\r\n' < "${LAST_DATASET_FILE}")"
  CURRENT_SOURCE_DATASET_ID="$(strip_route_suffix "${DATASET_ID:-}")"
  LAST_SOURCE_DATASET_ID="$(strip_route_suffix "${LAST_DATASET_ID}")"

  if [[ -n "${LAST_DATASET_ID}" && "${LAST_SOURCE_DATASET_ID}" == "${CURRENT_SOURCE_DATASET_ID}" ]]; then
    DATASET_ID="${LAST_DATASET_ID}"
  fi
  export DATASET_ID
fi
