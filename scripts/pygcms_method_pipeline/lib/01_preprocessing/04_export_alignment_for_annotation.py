#!/usr/bin/env python3
# =============================================================================
# Purpose: Export global aligned features and their sample peak-area matrix.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# =============================================================================

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


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


def to_float(value: object) -> float:
    try:
        return float(value)
    except Exception:
        return 0.0


def normalize_metadata(
    metadata: list[dict[str, str]],
    matrix_rows: list[dict[str, str]],
    foundin_min: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    metadata_ids = [str(row.get("GlobalFeatureID", "")).strip() for row in metadata]
    matrix_ids = [str(row.get("GlobalFeatureID", "")).strip() for row in matrix_rows]
    if any(not fid for fid in metadata_ids) or any(not fid for fid in matrix_ids):
        raise ValueError("Global feature tables contain blank IDs")
    if len(metadata_ids) != len(set(metadata_ids)):
        raise ValueError("Global feature metadata contains duplicate IDs")
    if len(matrix_ids) != len(set(matrix_ids)):
        raise ValueError("Global feature area matrix contains duplicate IDs")
    if set(metadata_ids) != set(matrix_ids):
        only_metadata = sorted(set(metadata_ids) - set(matrix_ids))
        only_matrix = sorted(set(matrix_ids) - set(metadata_ids))
        raise ValueError(
            "Global feature metadata and area matrix ID sets differ: "
            f"metadata_only={only_metadata[:10]} matrix_only={only_matrix[:10]}"
        )
    matrix_by_id = {str(row.get("GlobalFeatureID", "")): row for row in matrix_rows}
    sample_cols = [col for col in (matrix_rows[0].keys() if matrix_rows else []) if col != "GlobalFeatureID"]

    normalized_meta: list[dict[str, object]] = []
    normalized_matrix: list[dict[str, object]] = []
    membership: list[dict[str, object]] = []

    for row in metadata:
        fid = str(row.get("GlobalFeatureID", ""))
        if not fid:
            continue
        area_row = matrix_by_id.get(fid, {"GlobalFeatureID": fid})
        areas = {sample: to_float(area_row.get(sample, 0)) for sample in sample_cols}
        found_in = sum(1 for value in areas.values() if value > 0)
        total_area = sum(areas.values())
        if found_in < foundin_min:
            continue

        out = dict(row)
        out["FoundIn"] = found_in
        out["n_samples"] = len(sample_cols)
        out["total_area"] = round(total_area, 6)
        normalized_meta.append(out)

        matrix_out: dict[str, object] = {"GlobalFeatureID": fid}
        matrix_out.update(areas)
        normalized_matrix.append(matrix_out)

        block_feature_ids = str(row.get("block_feature_ids", "")).strip()
        for block_feature_id in [part for part in block_feature_ids.split(";") if part]:
            membership.append({
                "GlobalFeatureID": fid,
                "block_feature_id": block_feature_id,
                "rt_blocks": row.get("rt_blocks", ""),
                "representative_block_feature_id": row.get("representative_block_feature_id", ""),
            })

    summary: list[dict[str, object]] = []
    all_foundins = []
    for row in metadata:
        fid = str(row.get("GlobalFeatureID", ""))
        area_row = matrix_by_id.get(fid, {"GlobalFeatureID": fid})
        all_foundins.append(sum(1 for sample in sample_cols if to_float(area_row.get(sample, 0)) > 0))
    for threshold in sorted(set([1, 2, 3, foundin_min])):
        ids = {
            str(row.get("GlobalFeatureID", ""))
            for row in metadata
            if sum(1 for sample in sample_cols if to_float(matrix_by_id.get(str(row.get("GlobalFeatureID", "")), {}).get(sample, 0)) > 0) >= threshold
        }
        total_area = sum(
            sum(to_float(matrix_by_id.get(fid, {}).get(sample, 0)) for sample in sample_cols)
            for fid in ids
        )
        summary.append({
            "foundin_min": threshold,
            "n_global_features": len(ids),
            "n_member_peaks": sum(1 for row in membership if str(row["GlobalFeatureID"]) in ids),
            "total_area": round(total_area, 6),
        })
    summary.append({
        "foundin_min": "all",
        "n_global_features": len(metadata),
        "n_member_peaks": sum(1 for row in metadata for part in str(row.get("block_feature_ids", "")).split(";") if part),
        "total_area": round(sum(sum(to_float(row.get(sample, 0)) for sample in sample_cols) for row in matrix_rows), 6),
        "median_foundin": sorted(all_foundins)[len(all_foundins) // 2] if all_foundins else 0,
    })
    return normalized_meta, normalized_matrix, membership, summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Export alignment global features to the standard annotation input layout.")
    parser.add_argument("--alignment-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--foundin-min", required=True, type=int)
    args = parser.parse_args()

    metadata_file = args.alignment_dir / "02_global_feature_metadata.csv"
    matrix_file = args.alignment_dir / "02_global_feature_area_matrix.csv"
    if not metadata_file.exists():
        raise SystemExit(f"Missing alignment metadata: {metadata_file}")
    if not matrix_file.exists():
        raise SystemExit(f"Missing alignment area matrix: {matrix_file}")

    metadata, matrix_rows, membership, summary = normalize_metadata(
        read_csv(metadata_file),
        read_csv(matrix_file),
        args.foundin_min,
    )

    write_csv(args.out_dir / "features_for_annotation.csv", metadata)
    write_csv(args.out_dir / "feature_area_matrix.csv", matrix_rows)
    write_csv(args.out_dir / "feature_membership.csv", membership)
    write_csv(args.out_dir / "reproducibility_filter_summary.csv", summary)
    print(f"Wrote alignment annotation inputs to {args.out_dir}")


if __name__ == "__main__":
    main()
