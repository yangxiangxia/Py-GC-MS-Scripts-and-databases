#!/usr/bin/env python3
# =============================================================================
# Purpose: Correct retention times using internal standards and calculate RI.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# =============================================================================

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import math
from pathlib import Path
import statistics
from typing import Iterable


@dataclass(frozen=True)
class InternalStandardShift:
    shifts_min: dict[str, float]
    robust_shift_min: float
    disagreement_sec: float
    status: str


@dataclass(frozen=True)
class LadderPoint:
    carbon_number: int
    ri: float
    rt_min: float


@dataclass(frozen=True)
class RIResult:
    measured_ri: float | None
    status: str
    confidence: str
    lower_carbon: int | None = None
    upper_carbon: int | None = None
    lower_rt_min: float | None = None
    upper_rt_min: float | None = None
    lower_ri: float | None = None
    upper_ri: float | None = None


@dataclass(frozen=True)
class MappingConfig:
    alignment_dir: Path
    output_dir: Path
    reference_is: Path
    ri_ladder: Path
    foundin_min: int = 3
    expected_reference_block: str = ""


def alignment_rt_confidence(metadata_row: dict[str, object], found_in: int) -> str:
    origins = {
        value.strip()
        for value in str(metadata_row.get("feature_origins", "") or "").split(";")
        if value.strip()
    }
    legacy_origins = {"erah_aligned", "alignid0_attached"}
    singleton_origins = {"singleton_cluster", "singleton_unmatched"}
    has_legacy = bool(origins & legacy_origins) or not origins
    has_attached_singleton = "alignid0_attached" in origins
    has_singleton = bool(origins & singleton_origins) or str(
        metadata_row.get("contains_singleton_derived", "")
    ).strip().lower() in {"1", "true", "yes"}
    if has_legacy and (has_singleton or has_attached_singleton):
        return "legacy_anchor_preserved_high"
    if has_legacy:
        return "legacy_aligned_high"
    if found_in < 3:
        return "singleton_derived_foundin_lt3_excluded"
    return "singleton_derived_reproducible_moderate"


def _strictly_increasing(values: Iterable[float]) -> bool:
    sequence = list(values)
    return all(right > left for left, right in zip(sequence, sequence[1:]))


def summarize_internal_standard_shifts(shifts_min: dict[str, float]) -> InternalStandardShift:
    if len(shifts_min) != 2:
        raise ValueError("Exactly two internal-standard shifts are required")
    values = list(shifts_min.values())
    if not all(math.isfinite(value) for value in values):
        raise ValueError("Internal-standard shifts must be finite")
    robust_shift = statistics.median(values)
    disagreement_sec = abs(values[0] - values[1]) * 60.0
    epsilon = 1e-9
    if disagreement_sec <= 10.0 + epsilon:
        status = "pass"
    elif disagreement_sec <= 20.0 + epsilon:
        status = "warning"
    else:
        status = "fail"
    return InternalStandardShift(
        shifts_min=dict(shifts_min),
        robust_shift_min=robust_shift,
        disagreement_sec=disagreement_sec,
        status=status,
    )


def interpolate_measured_ri(
    rt_min: float,
    ladder: list[LadderPoint],
) -> RIResult:
    ordered = sorted(ladder, key=lambda point: point.rt_min)
    if len(ordered) < 2:
        raise ValueError("At least two RI ladder points are required")
    if not _strictly_increasing(point.rt_min for point in ordered):
        raise ValueError("RI ladder RT values must be strictly increasing")
    if rt_min < ordered[0].rt_min or rt_min > ordered[-1].rt_min:
        return RIResult(None, "outside_calibration_range", "outside_calibration_range")

    confidence = "within_calibration_range"

    for point in ordered:
        if math.isclose(rt_min, point.rt_min, rel_tol=0.0, abs_tol=1e-12):
            return RIResult(
                point.ri,
                "ladder_anchor",
                confidence,
                point.carbon_number,
                point.carbon_number,
                point.rt_min,
                point.rt_min,
                point.ri,
                point.ri,
            )

    for lower, upper in zip(ordered, ordered[1:]):
        if lower.rt_min < rt_min < upper.rt_min:
            fraction = (rt_min - lower.rt_min) / (upper.rt_min - lower.rt_min)
            measured_ri = lower.ri + fraction * (upper.ri - lower.ri)
            return RIResult(
                measured_ri,
                "interpolated",
                confidence,
                lower.carbon_number,
                upper.carbon_number,
                lower.rt_min,
                upper.rt_min,
                lower.ri,
                upper.ri,
            )
    raise RuntimeError("Failed to bracket an RT that lies inside the RI ladder")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: list[str] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = []
        for row in rows:
            for key in row:
                if key not in fieldnames:
                    fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _float(value: object) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Expected a finite number, got {value!r}")
    return number


def _median_numeric(rows: list[dict[str, str]], column: str) -> float:
    values: list[float] = []
    for row in rows:
        try:
            values.append(_float(row[column]))
        except (KeyError, TypeError, ValueError):
            continue
    if not values:
        raise ValueError(f"No finite values found for {column}")
    return statistics.median(values)


def compute_internal_standard_shift(
    reference_rows: list[dict[str, str]],
    sample_rows: list[dict[str, str]],
    reference_block: str,
) -> tuple[InternalStandardShift, list[dict[str, object]]]:
    block_rows = [row for row in sample_rows if row.get("rt_block") == reference_block]
    if not block_rows:
        raise ValueError(f"No internal-standard rows found for reference block {reference_block}")
    reference_by_name = {
        row.get("standard", "").strip().lower(): _float(row["reference_rt_min"])
        for row in reference_rows
        if row.get("standard", "").strip()
    }
    definitions = [
        ("tetracosane-d50", "c24_rt"),
        ("chrysene-d12", "chrysene_rt"),
    ]
    shifts: dict[str, float] = {}
    audit: list[dict[str, object]] = []
    for standard, sample_column in definitions:
        if standard not in reference_by_name:
            raise ValueError(f"Missing reference RT for {standard}")
        block_median = _median_numeric(block_rows, sample_column)
        reference_rt = reference_by_name[standard]
        shift = reference_rt - block_median
        shifts[standard] = shift
        audit.append(
            {
                "standard": standard,
                "reference_block": reference_block,
                "reference_block_median_rt_min": block_median,
                "pygcms_reference_run_rt_min": reference_rt,
                "individual_shift_min": shift,
                "individual_shift_sec": shift * 60.0,
                "n_reference_block_samples": len(block_rows),
            }
        )
    summary = summarize_internal_standard_shifts(shifts)
    for row in audit:
        row["robust_shift_min"] = summary.robust_shift_min
        row["shift_disagreement_sec"] = summary.disagreement_sec
        row["shift_qc_status"] = summary.status
    return summary, audit


def _load_ladder(rows: list[dict[str, str]]) -> list[LadderPoint]:
    ladder = [
        LadderPoint(
            carbon_number=int(float(row["carbon_number"])),
            ri=_float(row["RI"]),
            rt_min=_float(row["RT_min"]),
        )
        for row in rows
    ]
    ladder.sort(key=lambda point: point.rt_min)
    if len({point.carbon_number for point in ladder}) != len(ladder):
        raise ValueError("RI ladder contains duplicate carbon numbers")
    if not _strictly_increasing(point.carbon_number for point in ladder):
        raise ValueError("RI ladder carbon numbers are not strictly increasing")
    if not _strictly_increasing(point.rt_min for point in ladder):
        raise ValueError("RI ladder RT values are not strictly increasing")
    return ladder


def _summary_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    return [dict(summary)]


def _write_failed_summary(
    output_dir: Path,
    summary: dict[str, object],
    failure_reason: str,
) -> dict[str, object]:
    summary.update({"status": "fail", "failure_reason": failure_reason})
    write_csv(output_dir / "shadow_run_summary.csv", _summary_rows(summary))
    return summary


def run_shadow_mapping(config: MappingConfig) -> dict[str, object]:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    filtered_output = (
        config.output_dir
        / f"foundin_ge{config.foundin_min}_features_measured_ri.csv"
    )
    ready_outputs = [
        config.output_dir / "all_global_features_measured_ri.csv",
        filtered_output,
        config.output_dir / "features_for_annotation.csv",
    ]
    for path in ready_outputs:
        if path.exists():
            path.unlink()

    required_alignment_files = {
        "metadata": config.alignment_dir / "02_global_feature_metadata.csv",
        "area": config.alignment_dir / "02_global_feature_area_matrix.csv",
        "membership": config.alignment_dir / "02_block_feature_to_global_feature_map.csv",
        "global_summary": config.alignment_dir / "02_global_merge_summary.csv",
        "corrections": config.alignment_dir / "02_rt_correction_by_block.csv",
        "internal_standards": config.alignment_dir / "00_internal_standard_rt_by_sample.csv",
    }
    for name, path in required_alignment_files.items():
        if not path.exists():
            return _write_failed_summary(
                config.output_dir,
                {"status": "fail"},
                f"missing_{name}_input",
            )

    metadata = read_csv(required_alignment_files["metadata"])
    area_rows = read_csv(required_alignment_files["area"])
    membership = read_csv(required_alignment_files["membership"])
    global_summary_rows = read_csv(required_alignment_files["global_summary"])
    corrections = read_csv(required_alignment_files["corrections"])
    sample_internal_standards = read_csv(required_alignment_files["internal_standards"])
    summary: dict[str, object] = {
        "status": "running",
        "n_global_features_input": len(metadata),
        "n_memberships_input": len(membership),
        "foundin_min": config.foundin_min,
    }
    if len(global_summary_rows) != 1:
        return _write_failed_summary(config.output_dir, summary, "invalid_global_merge_summary")
    reference_block = global_summary_rows[0].get("reference_block", "")
    summary["reference_block"] = reference_block
    if config.expected_reference_block and reference_block != config.expected_reference_block:
        return _write_failed_summary(config.output_dir, summary, "unexpected_reference_block")

    correction_by_block: dict[str, dict[str, str]] = {}
    alignment_qc: list[dict[str, object]] = []
    for row in corrections:
        block = row.get("rt_block", "")
        status = "pass"
        if not block or block in correction_by_block:
            status = "duplicate_or_missing_block_mapping"
        elif row.get("reference_block") != reference_block:
            status = "wrong_reference_block"
        elif block == reference_block and row.get("rt_correction_method") != "reference_identity":
            status = "reference_block_not_identity"
        correction_by_block[block] = row
        alignment_qc.append(
            {
                "rt_block": block,
                "reference_block": row.get("reference_block", ""),
                "rt_correction_method": row.get("rt_correction_method", ""),
                "input_qc_status": status,
            }
        )
    sample_blocks = {row.get("rt_block", "") for row in sample_internal_standards if row.get("rt_block")}
    missing_block_models = sorted(sample_blocks - set(correction_by_block))
    if missing_block_models:
        alignment_qc.append(
            {
                "rt_block": ";".join(missing_block_models),
                "reference_block": reference_block,
                "rt_correction_method": "",
                "input_qc_status": "missing_block_mapping",
            }
        )
    write_csv(config.output_dir / "reference_group_alignment_input_qc.csv", alignment_qc)
    if any(row["input_qc_status"] != "pass" for row in alignment_qc):
        return _write_failed_summary(config.output_dir, summary, "reference_group_alignment_input_qc_failed")

    reference_is_rows = read_csv(config.reference_is)
    try:
        shift, shift_audit = compute_internal_standard_shift(
            reference_is_rows,
            sample_internal_standards,
            reference_block,
        )
    except ValueError as error:
        summary["error_detail"] = str(error)
        return _write_failed_summary(config.output_dir, summary, "internal_standard_input_invalid")
    write_csv(config.output_dir / "internal_standard_reference_shift_audit.csv", shift_audit)
    summary.update(
        {
            "internal_standard_robust_shift_min": shift.robust_shift_min,
            "internal_standard_shift_disagreement_sec": shift.disagreement_sec,
            "internal_standard_shift_qc_status": shift.status,
        }
    )
    if shift.status == "fail":
        return _write_failed_summary(config.output_dir, summary, "internal_standard_shift_disagreement")

    ladder = _load_ladder(read_csv(config.ri_ladder))
    metadata_ids = [row.get("GlobalFeatureID", "") for row in metadata]
    area_by_id = {row.get("GlobalFeatureID", ""): row for row in area_rows}
    if not all(metadata_ids) or len(set(metadata_ids)) != len(metadata_ids):
        return _write_failed_summary(config.output_dir, summary, "invalid_global_feature_ids")
    if set(metadata_ids) != set(area_by_id):
        return _write_failed_summary(config.output_dir, summary, "metadata_area_id_mismatch")
    if any(row.get("GlobalFeatureID", "") not in area_by_id for row in membership):
        return _write_failed_summary(config.output_dir, summary, "membership_global_id_mismatch")

    coordinate_rows: list[dict[str, object]] = []
    sample_columns = [column for column in (area_rows[0].keys() if area_rows else []) if column != "GlobalFeatureID"]
    for row in metadata:
        feature_id = row["GlobalFeatureID"]
        reference_group_rt = _float(row["corrected_tmean"])
        reference_corrected_rt = reference_group_rt + shift.robust_shift_min
        found_in = sum(
            1
            for sample in sample_columns
            if _float(area_by_id[feature_id].get(sample, 0) or 0) > 0
        )
        coordinate_rows.append(
            {
                "GlobalFeatureID": feature_id,
                "reference_group_rt_min": reference_group_rt,
                "internal_standard_shift_min": shift.robust_shift_min,
                "reference_corrected_rt_min": reference_corrected_rt,
                "FoundIn": found_in,
            }
        )

    trace: list[dict[str, object]] = []
    results: list[dict[str, object]] = []
    for metadata_row, coordinate in zip(metadata, coordinate_rows):
        reference_corrected_rt = _float(coordinate["reference_corrected_rt_min"])
        ri_result = interpolate_measured_ri(reference_corrected_rt, ladder)
        trace_row = {
            "GlobalFeatureID": coordinate["GlobalFeatureID"],
            "reference_group_rt_min": coordinate["reference_group_rt_min"],
            "internal_standard_shift_min": coordinate["internal_standard_shift_min"],
            "reference_corrected_rt_min": reference_corrected_rt,
            "lower_alkane_carbon": ri_result.lower_carbon,
            "upper_alkane_carbon": ri_result.upper_carbon,
            "lower_alkane_rt_min": ri_result.lower_rt_min,
            "upper_alkane_rt_min": ri_result.upper_rt_min,
            "lower_alkane_ri": ri_result.lower_ri,
            "upper_alkane_ri": ri_result.upper_ri,
            "measured_RI": ri_result.measured_ri,
            "RI_status": ri_result.status,
            "RI_confidence": ri_result.confidence,
            "alignment_RT_confidence": alignment_rt_confidence(
                metadata_row, int(coordinate["FoundIn"])
            ),
            "FoundIn": coordinate["FoundIn"],
            "n_blocks": metadata_row.get("n_blocks", ""),
            "n_block_features": metadata_row.get("n_block_features", ""),
        }
        trace.append(trace_row)
        result_row = dict(metadata_row)
        result_row.update({key: value for key, value in trace_row.items() if key != "GlobalFeatureID"})
        results.append(result_row)

    write_csv(config.output_dir / "feature_rt_coordinate_trace.csv", trace)
    write_csv(config.output_dir / "all_global_features_measured_ri.csv", results)
    filtered_results = [row for row in results if int(row["FoundIn"]) >= config.foundin_min]
    write_csv(filtered_output, filtered_results)
    write_csv(config.output_dir / "features_for_annotation.csv", filtered_results)

    ri_calibration_qc = [
        {"metric": "reference_block", "value": reference_block, "status": "pass"},
        {"metric": "internal_standard_shift_disagreement_sec", "value": shift.disagreement_sec, "status": shift.status},
        {"metric": "non_null_measured_ri", "value": sum(row["measured_RI"] is not None for row in trace), "status": "reported"},
    ]
    write_csv(config.output_dir / "ri_calibration_qc.csv", ri_calibration_qc)
    summary.update(
        {
            "status": "pass",
            "failure_reason": "",
            "n_global_features_output": len(results),
            "n_foundin_filtered_output": len(filtered_results),
            "n_memberships_output": len(membership),
            "n_measured_ri": sum(row["measured_RI"] is not None for row in trace),
            "n_outside_calibration_range": sum(row["RI_status"] == "outside_calibration_range" for row in trace),
        }
    )
    write_csv(config.output_dir / "shadow_run_summary.csv", _summary_rows(summary))
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Apply internal-standard RT correction and calculate RI from the C7-C40 ladder."
    )
    parser.add_argument("--alignment-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--reference-is", required=True, type=Path)
    parser.add_argument("--ri-ladder", required=True, type=Path)
    parser.add_argument("--foundin-min", type=int, default=3)
    parser.add_argument("--expected-reference-block", default="")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = MappingConfig(
        alignment_dir=args.alignment_dir,
        output_dir=args.output_dir,
        reference_is=args.reference_is,
        ri_ladder=args.ri_ladder,
        foundin_min=args.foundin_min,
        expected_reference_block=args.expected_reference_block,
    )
    summary = run_shadow_mapping(config)
    print(f"Reference RI mapping status: {summary['status']}")
    print(f"Reference RI output: {config.output_dir}")
    if summary["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
