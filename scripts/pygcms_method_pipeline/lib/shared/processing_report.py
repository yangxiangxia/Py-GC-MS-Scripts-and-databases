#!/usr/bin/env python3
# =============================================================================
# Purpose: Summarize processing stages, retained areas, and annotation outputs.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# =============================================================================

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def read_first_existing(paths: list[Path]) -> list[dict[str, str]]:
    for path in paths:
        rows = read_csv(path)
        if rows:
            return rows
    return []


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
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


def fnum(value: object, default: float = 0.0) -> float:
    try:
        out = float(str(value).strip())
    except Exception:
        return default
    return out if math.isfinite(out) else default


def is_true(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def area_sum(rows: list[dict[str, str]], area_col: str = "Area") -> float:
    return sum(max(fnum(row.get(area_col)), 0.0) for row in rows)


def matrix_area(rows: list[dict[str, str]]) -> float:
    total = 0.0
    for row in rows:
        for key, value in row.items():
            if key in {"GlobalFeatureID", "CompoundID"}:
                continue
            total += max(fnum(value), 0.0)
    return total


STEP_LABELS = {
    "01_raw_deconvolution": "Raw deconvoluted peaks",
    "02_after_overdeconv_cleanup": "After over-deconvolution cleanup",
    "03_after_shoulder_merge": "After shoulder peak merging",
    "04_after_area_fraction_filter": "After within-sample area filtering",
    "05_after_global_grouping_all": "After cross-sample feature grouping",
    "06_after_foundin_filter": "After reproducibility filtering",
    "07_manual_review_table": "Manual annotation-review table",
    "08_after_manual_review_keep_identified": "After manual identification review",
    "09_after_post_review_same_name_rt_merge": "After same-name close-RT compound merging",
    "10_after_smiles": "After SMILES retrieval",
    "11_after_classification": "After structural classification",
}

WIDE_STEP_NAMES = {
    "01_raw_deconvolution": "raw_deconvolution",
    "02_after_overdeconv_cleanup": "overdeconv_cleanup",
    "03_after_shoulder_merge": "shoulder_merge",
    "04_after_area_fraction_filter": "area_filter",
    "05_after_global_grouping_all": "global_grouping",
    "06_after_foundin_filter": "foundin_filter",
    "09_after_post_review_same_name_rt_merge": "final_compounds",
}

SIMPLE_SAMPLE_SUMMARY_STEPS = [
    ("02_after_overdeconv_cleanup", "after_cleanup", "peak"),
    ("03_after_shoulder_merge", "after_shoulder_merge", "peak"),
    ("04_after_area_fraction_filter", "after_area_filter", "peak"),
    ("05_after_global_grouping_all", "after_global_grouping", "feature"),
    ("06_after_foundin_filter", "after_foundin_filter", "feature"),
    ("09_after_post_review_same_name_rt_merge", "after_final_compound_merge", "compound"),
]


def by_sample_peak(rows: list[dict[str, str]], step: str, raw_area_by_sample: dict[str, float]) -> list[dict[str, object]]:
    out = []
    samples = sorted({row.get("sample", "") for row in rows})
    for sample in samples:
        sample_rows = [row for row in rows if row.get("sample", "") == sample]
        area = area_sum(sample_rows)
        raw_area = raw_area_by_sample.get(sample, 0.0)
        out.append({
            "sample": sample,
            "step": step,
            "step_label": STEP_LABELS.get(step, step),
            "status": "complete",
            "unit": "peak",
            "n_retained": len(sample_rows),
            "area": round(area, 6),
            "area_pct_vs_raw_deconvolved": round(100 * area / raw_area, 3) if raw_area > 0 else "",
            "area_pct_vs_previous_step": "",
        })
    return out


def by_sample_matrix(rows: list[dict[str, str]], step: str, raw_area_by_sample: dict[str, float], unit: str) -> list[dict[str, object]]:
    if not rows:
        return []
    id_cols = {"GlobalFeatureID", "CompoundID"}
    sample_cols = [col for col in rows[0].keys() if col not in id_cols]
    out = []
    for sample in sample_cols:
        area = sum(max(fnum(row.get(sample)), 0.0) for row in rows)
        n_retained = sum(max(fnum(row.get(sample)), 0.0) > 0 for row in rows)
        raw_area = raw_area_by_sample.get(sample, 0.0)
        out.append({
            "sample": sample,
            "step": step,
            "step_label": STEP_LABELS.get(step, step),
            "status": "complete",
            "unit": unit,
            "n_retained": n_retained,
            "area": round(area, 6),
            "area_pct_vs_raw_deconvolved": round(100 * area / raw_area, 3) if raw_area > 0 else "",
            "area_pct_vs_previous_step": "",
        })
    return out


def add_previous_step_pct(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    previous_area_by_sample: dict[str, float] = {}
    out = []
    for row in rows:
        sample = str(row.get("sample", ""))
        area = fnum(row.get("area"))
        previous_area = previous_area_by_sample.get(sample)
        row["area_pct_vs_previous_step"] = (
            round(100 * area / previous_area, 3)
            if previous_area and previous_area > 0
            else ""
        )
        if area > 0:
            previous_area_by_sample[sample] = area
        out.append(row)
    return out


def build_sample_area_wide(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    samples: dict[str, dict[str, object]] = {}
    for row in rows:
        sample = str(row.get("sample", ""))
        step = str(row.get("step", ""))
        short_step = WIDE_STEP_NAMES.get(step)
        if not sample or not short_step:
            continue
        sample_row = samples.setdefault(sample, {"sample": sample})
        unit = str(row.get("unit", ""))
        count_col = "n_peaks" if unit == "peak" else f"n_{unit}s"
        sample_row[f"{short_step}_{count_col}"] = row.get("n_retained", "")
        sample_row[f"{short_step}_area"] = row.get("area", "")
        sample_row[f"{short_step}_pct_of_raw"] = row.get("area_pct_vs_raw_deconvolved", "")
        sample_row[f"{short_step}_pct_of_previous_step"] = row.get("area_pct_vs_previous_step", "")
    return [samples[sample] for sample in sorted(samples)]


def build_sample_area_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_sample_step = {
        (str(row.get("sample", "")), str(row.get("step", ""))): row
        for row in rows
    }
    samples = sorted({sample for sample, _step in by_sample_step})
    out: list[dict[str, object]] = []
    for sample in samples:
        raw = by_sample_step.get((sample, "01_raw_deconvolution"), {})
        summary: dict[str, object] = {
            "sample": sample,
            "raw_peak_count": raw.get("n_retained", ""),
            "raw_area": raw.get("area", ""),
        }
        for step, prefix, unit_name in SIMPLE_SAMPLE_SUMMARY_STEPS:
            row = by_sample_step.get((sample, step))
            if not row:
                continue
            summary[f"{prefix}_{unit_name}_count"] = row.get("n_retained", "")
            summary[f"{prefix}_area_pct_of_raw"] = row.get("area_pct_vs_raw_deconvolved", "")
        out.append(summary)
    return out


def step_row(step: str, unit: str, rows: list[dict[str, str]], raw_area: float, previous_area: float | None, area: float | None = None) -> dict[str, object]:
    if not rows:
        return {
            "step": step,
            "step_label": STEP_LABELS.get(step, step),
            "status": "not_run",
            "unit": unit,
            "n_retained": "",
            "total_area": "",
            "area_pct_vs_raw_deconvolved": "",
            "area_pct_vs_previous_step": "",
        }
    total_area = area_sum(rows) if area is None else area
    return {
        "step": step,
        "step_label": STEP_LABELS.get(step, step),
        "status": "complete",
        "unit": unit,
        "n_retained": len(rows),
        "total_area": round(total_area, 6),
        "area_pct_vs_raw_deconvolved": round(100 * total_area / raw_area, 3) if raw_area > 0 else "",
        "area_pct_vs_previous_step": round(100 * total_area / previous_area, 3) if previous_area and previous_area > 0 else "",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize Py-GC-MS method-pipeline peak, feature, compound, and area retention.")
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--review-dir", type=Path)
    parser.add_argument("--report-dir", type=Path)
    args = parser.parse_args()
    review_dir = args.review_dir if args.review_dir else args.run_dir / "06_manual_review"

    raw = read_first_existing([
        args.run_dir / "01_deconvolution" / "raw_deconvoluted_peaks.csv",
        args.run_dir / "01_deconvolution" / "deconvolved_peaks_raw.csv",
    ])
    clean = read_first_existing([
        args.run_dir / "02_cleanup" / "cleaned_peaks.csv",
        args.run_dir / "02_cleanup" / "peaks_after_overdeconv_cleanup.csv",
    ])
    shoulder = read_first_existing([
        args.run_dir / "02_cleanup" / "shoulder_merged_peaks.csv",
        args.run_dir / "02_cleanup" / "peaks_after_shoulder_merge.csv",
    ])
    filtered_files = sorted((args.run_dir / "03_area_filter").glob("filtered_peaks_*.csv"))
    filtered = read_first_existing([
        args.run_dir / "03_area_filter" / "area_filtered_peaks.csv",
        filtered_files[0] if filtered_files else args.run_dir / "03_area_filter" / "__missing__.csv",
    ])
    area_filter_audit = read_csv(
        args.run_dir / "04_alignment" / "01_pre_alignment_area_filter_audit.csv"
    )
    if area_filter_audit:
        filtered = []
        for row in area_filter_audit:
            if not is_true(row.get("retained")):
                continue
            normalized = dict(row)
            normalized["ID"] = normalized.get("peak_id", "")
            filtered.append(normalized)
    global_all = read_first_existing([
        args.run_dir / "audit_files" / "all_feature_area_matrix_before_reproducibility_filter.csv",
        args.run_dir / "04_global_grouping" / "global_feature_area_matrix_all.csv",
    ])
    global_foundin = sorted((args.run_dir / "04_global_grouping").glob("global_feature_area_matrix_foundin_ge*.csv"))
    foundin_matrix = read_first_existing([
        args.run_dir / "04_global_grouping" / "feature_area_matrix.csv",
        global_foundin[0] if global_foundin else args.run_dir / "04_global_grouping" / "__missing__.csv",
    ])
    review = read_csv(review_dir / "compound_identification_review_CHECK.csv")
    kept = read_csv(args.run_dir / "07_post_review_merge" / "checked_features_kept.csv")
    final_matrix = read_csv(args.run_dir / "07_post_review_merge" / "final_compound_area_matrix.csv")
    smiles = read_csv(args.run_dir / "08_smiles" / "final_compounds_with_smiles.csv")
    classified = read_csv(args.run_dir / "09_classification" / "classified_compounds.csv")

    raw_area = area_sum(raw)
    raw_area_by_sample: dict[str, float] = {}
    for row in raw:
        sample = row.get("sample", "")
        raw_area_by_sample[sample] = raw_area_by_sample.get(sample, 0.0) + max(fnum(row.get("Area")), 0.0)

    steps = []
    previous_area: float | None = None
    for name, unit, rows, area in [
        ("01_raw_deconvolution", "peak", raw, None),
        ("02_after_overdeconv_cleanup", "peak", clean, None),
        ("03_after_shoulder_merge", "peak", shoulder, None),
        ("04_after_area_fraction_filter", "peak", filtered, None),
        ("05_after_global_grouping_all", "feature", global_all, matrix_area(global_all)),
        ("06_after_foundin_filter", "feature", foundin_matrix, matrix_area(foundin_matrix)),
        ("07_manual_review_table", "feature", review, matrix_area(foundin_matrix) if review else 0.0),
        ("08_after_manual_review_keep_identified", "feature", kept, sum(fnum(row.get("feature_total_area")) for row in kept)),
        ("09_after_post_review_same_name_rt_merge", "compound", final_matrix, matrix_area(final_matrix)),
        ("10_after_smiles", "compound", smiles, matrix_area(final_matrix) if smiles else 0.0),
        ("11_after_classification", "compound", classified, matrix_area(final_matrix) if classified else 0.0),
    ]:
        row = step_row(name, unit, rows, raw_area, previous_area, area=area)
        steps.append(row)
        row_area = fnum(row["total_area"])
        if row_area > 0:
            previous_area = row_area

    by_sample = []
    by_sample.extend(by_sample_peak(raw, "01_raw_deconvolution", raw_area_by_sample))
    by_sample.extend(by_sample_peak(clean, "02_after_overdeconv_cleanup", raw_area_by_sample))
    by_sample.extend(by_sample_peak(shoulder, "03_after_shoulder_merge", raw_area_by_sample))
    by_sample.extend(by_sample_peak(filtered, "04_after_area_fraction_filter", raw_area_by_sample))
    by_sample.extend(by_sample_matrix(global_all, "05_after_global_grouping_all", raw_area_by_sample, "feature"))
    by_sample.extend(by_sample_matrix(foundin_matrix, "06_after_foundin_filter", raw_area_by_sample, "feature"))
    by_sample.extend(by_sample_matrix(final_matrix, "09_after_post_review_same_name_rt_merge", raw_area_by_sample, "compound"))
    by_sample = add_previous_step_pct(by_sample)

    report_dir = args.report_dir if args.report_dir else args.run_dir / "reports"
    write_csv(report_dir / "pipeline_summary.csv", steps)
    write_csv(report_dir / "sample_area_summary.csv", build_sample_area_summary(by_sample))
    write_csv(report_dir / "sample_step_area_summary.csv", by_sample)
    write_csv(report_dir / "sample_area_retention_wide.csv", build_sample_area_wide(by_sample))
    write_csv(report_dir / "pipeline_summary_by_sample.csv", by_sample)

    removed_rows = []
    for row in read_first_existing([
        args.run_dir / "audit_files" / "overdeconvolution_removed_peaks.csv",
        args.run_dir / "02_cleanup" / "overdeconvolution_removed_peaks.csv",
    ]):
        removed_rows.append({"removal_step": "overdeconvolution_cleanup", **row})
    for row in read_first_existing([
        args.run_dir / "audit_files" / "shoulder_merge_audit.csv",
        args.run_dir / "02_cleanup" / "shoulder_merge_audit.csv",
    ]):
        removed_rows.append({"removal_step": "shoulder_merge", **row})
    for row in area_filter_audit:
        if not is_true(row.get("retained")):
            removed_rows.append({
                "removal_step": "within_sample_area_fraction_filter",
                **row,
            })
    write_csv(args.run_dir / "audit_files" / "removed_peak_audit_combined.csv", removed_rows)

    print(f"Wrote processing report to {report_dir}")


if __name__ == "__main__":
    main()
