#!/usr/bin/env python3
# =============================================================================
# Purpose: Check that downstream inputs belong to the configured dataset.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# =============================================================================

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def norm(path: str | Path) -> str:
    return str(Path(path).expanduser().resolve())


def read_parameter_csv(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    return {
        str(row.get("parameter", "")).strip(): str(row.get("value", "")).strip()
        for row in rows
        if str(row.get("parameter", "")).strip()
    }


def read_preblock_cdf_dirs(preblock_dir: Path) -> set[str]:
    dirs: set[str] = set()
    for sample_list in sorted(preblock_dir.glob("auto_block_sample_list_*.csv")):
        with sample_list.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                cdf_file = row.get("cdf_file", "").strip()
                if cdf_file:
                    dirs.add(norm(Path(cdf_file).parent))
    return dirs


def fail(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stop downstream steps when existing work files came from another CDF directory."
    )
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--cdf-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    args = parser.parse_args()

    expected_cdf_dir = norm(args.cdf_dir)
    work_dir = args.work_dir

    observed_dirs: list[tuple[str, str]] = []

    alignment_params = work_dir / "04_alignment" / "00_run_parameters.csv"
    if alignment_params.exists():
        params = read_parameter_csv(alignment_params)
        if params.get("cdf_dir"):
            observed_dirs.append((str(alignment_params), norm(params["cdf_dir"])))

    preblock_dirs = read_preblock_cdf_dirs(work_dir / "03_alignment_preblock")
    for cdf_dir in sorted(preblock_dirs):
        observed_dirs.append((str(work_dir / "03_alignment_preblock"), cdf_dir))

    for deconv_params in [
        work_dir / "01_deconvolution" / "deconvolution_parameters.csv",
        work_dir / "audit_files" / "deconvolution_parameters.csv",
    ]:
        if deconv_params.exists():
            params = read_parameter_csv(deconv_params)
            if params.get("CDF_DIR"):
                observed_dirs.append((str(deconv_params), norm(params["CDF_DIR"])))

    mismatches = [
        (source, observed)
        for source, observed in observed_dirs
        if observed != expected_cdf_dir
    ]

    if mismatches:
        lines = [
            "Dataset context mismatch: existing work files were generated from another CDF directory.",
            f"Current DATASET_ID: {args.dataset_id}",
            f"Current CDF_DIR: {expected_cdf_dir}",
            "Observed previous CDF_DIR values:",
        ]
        lines.extend(f"  {source}: {observed}" for source, observed in mismatches)
        lines.extend([
            "",
            "Run the matching 01 preprocessing step again after archiving/removing the mismatched work directory,",
            "then rerun 02_compound_annotation_and_review.sh to create a fresh annotation_review_CHECK.csv.",
        ])
        fail("\n".join(lines))


if __name__ == "__main__":
    main()
