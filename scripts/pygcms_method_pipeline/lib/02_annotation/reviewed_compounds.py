#!/usr/bin/env python3
# =============================================================================
# Purpose: Finalize reviewed identifications and merge eligible compound rows.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# Responsibility: Validate review inputs and produce compound tables, areas, and membership.
# =============================================================================

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


BLANK_FINAL_VALUES = {"", "na", "n/a"}
UNKNOWN_FINAL_VALUES = {
    "unknown", "unknow", "unidentified", "not identified", "reject", "rejected"
}
CSV_ENCODINGS = ("utf-8-sig", "gb18030", "gbk")


def read_csv(path: Path) -> list[dict[str, str]]:
    last_error: UnicodeDecodeError | None = None
    for encoding in CSV_ENCODINGS:
        try:
            with path.open(newline="", encoding=encoding) as handle:
                return list(csv.DictReader(handle))
        except UnicodeDecodeError as exc:
            last_error = exc
    if last_error is not None:
        raise last_error
    return []


def validate_unique_nonblank_ids(
    rows: list[dict[str, str]],
    column: str,
    blank_description: str,
    duplicate_description: str,
) -> None:
    cleaned_ids: list[str] = []
    blank_positions: list[int] = []
    duplicate_ids: list[str] = []
    seen: set[str] = set()
    for csv_row_position, row in enumerate(rows, start=2):
        feature_id = str(row.get(column, "")).strip()
        cleaned_ids.append(feature_id)
        if not feature_id:
            blank_positions.append(csv_row_position)
        elif feature_id in seen and feature_id not in duplicate_ids:
            duplicate_ids.append(feature_id)
        seen.add(feature_id)
    if blank_positions:
        raise ValueError(
            blank_description
            + " at CSV row positions: "
            + ", ".join(str(position) for position in blank_positions)
        )
    if duplicate_ids:
        raise ValueError(
            duplicate_description + ": " + ", ".join(duplicate_ids)
        )
    for row, feature_id in zip(rows, cleaned_ids):
        row[column] = feature_id


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


def normalize_name(value: str) -> str:
    value = str(value or "").lower()
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return " ".join(value.split())


def final_inchikey(row: dict[str, object]) -> str:
    confirmed = str(row.get("final_InChIKey", "")).strip()
    if confirmed:
        return confirmed
    final_name = normalize_name(str(row.get("parsed_final_identification", "")))
    for name_col, key_col in (
        ("best_identification", "best_InChIKey"),
        ("library1_identification", "library1_InChIKey"),
        ("library2_identification", "library2_InChIKey"),
    ):
        if final_name and final_name == normalize_name(str(row.get(name_col, ""))):
            return str(row.get(key_col, "")).strip()
    return ""


def is_blank_final(value: str) -> bool:
    return normalize_name(value) in BLANK_FINAL_VALUES


def is_unknown_final(value: str) -> bool:
    return normalize_name(value) in UNKNOWN_FINAL_VALUES


def parse_final_formula_metadata(raw_name: str, raw_formula: str) -> tuple[str, str]:
    name = str(raw_name or "").strip()
    formula = str(raw_formula or "").strip()
    embedded_score = ""
    match = re.search(r"\(([^()]*)\)\s*$", name)
    if match:
        inside = match.group(1).strip()
        parts = [p.strip() for p in inside.split(";")]
        if len(parts) >= 2:
            embedded_score = parts[0]
            if not formula:
                formula = parts[1]
    return formula.strip(" ;,"), embedded_score


def feature_spectrum_summary(row: dict[str, object]) -> str:
    feature_id = str(row.get("GlobalFeatureID", "")).strip()
    base_mz = str(row.get("base_mz", "")).strip()
    ions = str(row.get("characteristic_ions", "")).strip()
    parts = []
    if base_mz:
        parts.append(f"base_mz={base_mz}")
    if ions:
        parts.append(f"ions={ions}")
    if not parts:
        return feature_id
    return f"{feature_id}({'; '.join(parts)})"


def first_nonblank(rows: list[dict[str, object]], column: str) -> str:
    for row in rows:
        value = str(row.get(column, "")).strip()
        if value:
            return value
    return ""


def joined_nonblank(rows: list[dict[str, object]], column: str) -> str:
    values: list[str] = []
    seen: set[str] = set()
    for row in rows:
        value = str(row.get(column, "")).strip()
        if value and value not in seen:
            values.append(value)
            seen.add(value)
    return ";".join(values)


def strongest_literature_support(rows: list[dict[str, object]]) -> tuple[str, str]:
    ranks = {"compound": 3, "homologous series": 2, "family": 1}
    levels = [str(row.get("literature_support_level", "")).strip() for row in rows]
    supported = [level for level in levels if level.lower() in ranks]
    if not supported:
        return "No", "Unavailable"
    level = max(supported, key=lambda value: ranks[value.lower()])
    return "Yes", level


def choose_primary_member(members: list[dict[str, object]]) -> dict[str, object]:
    """Use the most broadly observed reviewed feature as the compound representative."""
    return max(
        members,
        key=lambda m: (
            fnum(m.get("feature_found_in"), 0.0),
            fnum(m.get("feature_total_area"), 0.0),
            fnum(m.get("best_score"), 0.0),
        ),
    )


def parse_feature_spectrum(value: object) -> dict[int, float]:
    spectrum: dict[int, float] = {}
    for token in str(value or "").split():
        parts = token.split(",", 1)
        if len(parts) != 2:
            continue
        try:
            mz = int(round(float(parts[0])))
            intensity = float(parts[1])
        except (TypeError, ValueError):
            continue
        if intensity > 0 and math.isfinite(intensity):
            spectrum[mz] = spectrum.get(mz, 0.0) + intensity
    return spectrum


def spectral_cosine(
    left: dict[int, float], right: dict[int, float]
) -> float | None:
    if not left or not right:
        return None
    left_norm = math.sqrt(sum(value * value for value in left.values()))
    right_norm = math.sqrt(sum(value * value for value in right.values()))
    if left_norm <= 0 or right_norm <= 0:
        return None
    return sum(
        value * right.get(mz, 0.0) for mz, value in left.items()
    ) / (left_norm * right_norm)


def sample_presence(
    row: dict[str, object], sample_cols: list[str]
) -> set[str]:
    return {
        sample for sample in sample_cols
        if fnum(row.get(sample), 0.0) > 0
    }


def sample_jaccard(left: set[str], right: set[str]) -> float | None:
    union = left | right
    return len(left & right) / len(union) if union else None


def unidentified_pair_key(
    left: dict[str, object], right: dict[str, object]
) -> tuple[str, str]:
    return tuple(sorted((
        str(left["GlobalFeatureID"]),
        str(right["GlobalFeatureID"]),
    )))


def build_unidentified_groups(
    rows: list[dict[str, object]],
    spectra_by_id: dict[str, dict[int, float]],
    presence_by_id: dict[str, set[str]],
    rt_window_sec: float,
    min_cosine: float,
    max_jaccard: float,
) -> tuple[list[list[dict[str, object]]], list[dict[str, object]]]:
    ordered = sorted(
        rows,
        key=lambda row: (
            not math.isfinite(float(row["rt"])),
            float(row["rt"]),
            str(row["GlobalFeatureID"]),
        ),
    )
    pair_metrics: dict[tuple[str, str], dict[str, object]] = {}
    for index, left in enumerate(ordered):
        left_rt = float(left["rt"])
        left_id = str(left["GlobalFeatureID"])
        for right in ordered[index + 1:]:
            right_rt = float(right["rt"])
            right_id = str(right["GlobalFeatureID"])
            rt_diff_sec = abs(right_rt - left_rt) * 60.0
            if not math.isfinite(rt_diff_sec) or rt_diff_sec > rt_window_sec:
                continue
            cosine = spectral_cosine(
                spectra_by_id.get(left_id, {}),
                spectra_by_id.get(right_id, {}),
            )
            left_presence = presence_by_id.get(left_id, set())
            right_presence = presence_by_id.get(right_id, set())
            jaccard = sample_jaccard(left_presence, right_presence)
            key = unidentified_pair_key(left, right)
            pair_metrics[key] = {
                "GlobalFeatureID_1": left_id,
                "GlobalFeatureID_2": right_id,
                "final_identification": "Unidentified",
                "final_formula_1": "",
                "final_formula_2": "",
                "rt_1": left.get("corrected_RT_min", ""),
                "rt_2": right.get("corrected_RT_min", ""),
                "rt_diff_sec": round(rt_diff_sec, 6),
                "spectral_cosine": "" if cosine is None else round(cosine, 6),
                "n_sample_overlap": len(left_presence & right_presence),
                "n_sample_union": len(left_presence | right_presence),
                "sample_jaccard": "" if jaccard is None else round(jaccard, 6),
                "qualifies": bool(
                    cosine is not None
                    and cosine >= min_cosine
                    and jaccard is not None
                    and jaccard <= max_jaccard
                ),
            }

    groups: list[list[dict[str, object]]] = [[row] for row in ordered]
    qualifying_pairs = sorted(
        (metric for metric in pair_metrics.values() if metric["qualifies"]),
        key=lambda metric: (
            -float(metric["spectral_cosine"]),
            float(metric["rt_diff_sec"]),
            str(metric["GlobalFeatureID_1"]),
            str(metric["GlobalFeatureID_2"]),
        ),
    )
    for metric in qualifying_pairs:
        left_id = str(metric["GlobalFeatureID_1"])
        right_id = str(metric["GlobalFeatureID_2"])
        left_index = next(
            index for index, group in enumerate(groups)
            if any(str(row["GlobalFeatureID"]) == left_id for row in group)
        )
        right_index = next(
            index for index, group in enumerate(groups)
            if any(str(row["GlobalFeatureID"]) == right_id for row in group)
        )
        if left_index == right_index:
            continue
        proposed = sorted(
            groups[left_index] + groups[right_index],
            key=lambda row: (float(row["rt"]), str(row["GlobalFeatureID"])),
        )
        total_span_sec = (
            float(proposed[-1]["rt"]) - float(proposed[0]["rt"])
        ) * 60.0
        all_pairs_pass = all(
            bool(pair_metrics.get(unidentified_pair_key(left, right), {}).get("qualifies"))
            for member_index, left in enumerate(proposed)
            for right in proposed[member_index + 1:]
        )
        if total_span_sec <= rt_window_sec and all_pairs_pass:
            keep_index = min(left_index, right_index)
            drop_index = max(left_index, right_index)
            groups[keep_index] = proposed
            groups.pop(drop_index)

    group_by_id: dict[str, int] = {}
    for group_index, group in enumerate(groups):
        for row in group:
            group_by_id[str(row["GlobalFeatureID"])] = group_index
    merge_rule = (
        f"unidentified_total_span_within_{int(rt_window_sec)}sec_"
        f"full_cosine_ge_{min_cosine:.2f}_"
        f"sample_jaccard_le_{max_jaccard:.2f}_complete_link"
    )
    audit_rows: list[dict[str, object]] = []
    for metric in pair_metrics.values():
        left_id = str(metric["GlobalFeatureID_1"])
        right_id = str(metric["GlobalFeatureID_2"])
        group_index = group_by_id.get(left_id)
        if (
            metric["qualifies"]
            and group_index is not None
            and group_index == group_by_id.get(right_id)
            and len(groups[group_index]) > 1
        ):
            audit_rows.append({**metric, "merge_rule": merge_rule})
    return groups, audit_rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Apply manual review decisions and conservatively merge same-name close-RT compounds.")
    parser.add_argument("--review", required=True, type=Path)
    parser.add_argument("--area-matrix", required=True, type=Path)
    parser.add_argument("--feature-metadata", type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--rt-window-sec", type=float, default=60.0)
    parser.add_argument("--unidentified-cosine", type=float, default=0.90)
    parser.add_argument("--unidentified-max-jaccard", type=float, default=0.10)
    parser.add_argument("--compound-id-prefix", default="CMPD")
    parser.add_argument("--compound-id-digits", type=int, default=5)
    args = parser.parse_args()

    review_rows = read_csv(args.review)
    validate_unique_nonblank_ids(
        review_rows,
        "ID",
        "Reviewed rows have blank ID values",
        "Duplicate reviewed ID values",
    )
    review_columns = list(review_rows[0].keys()) if review_rows else []
    area_rows = read_csv(args.area_matrix)
    validate_unique_nonblank_ids(
        area_rows,
        "GlobalFeatureID",
        "Area-matrix rows have blank GlobalFeatureID values",
        "Duplicate area-matrix GlobalFeatureID values",
    )
    metadata_rows = read_csv(args.feature_metadata) if args.feature_metadata else []
    area_by_id = {row["GlobalFeatureID"]: row for row in area_rows}
    metadata_by_id = {
        str(row.get("GlobalFeatureID", "")).strip(): row
        for row in metadata_rows
        if str(row.get("GlobalFeatureID", "")).strip()
    }
    sample_cols = [c for c in (area_rows[0].keys() if area_rows else []) if c != "GlobalFeatureID"]

    missing_area_ids = list(dict.fromkeys(
        feature_id
        for row in review_rows
        if (feature_id := str(row.get("ID", "")).strip())
        and feature_id not in area_by_id
    ))
    if missing_area_ids:
        raise ValueError(
            "Reviewed feature IDs missing from area matrix: "
            + ", ".join(missing_area_ids)
        )

    kept: list[dict[str, object]] = []
    excluded: list[dict[str, object]] = []
    for row in review_rows:
        fid = row.get("ID", "").strip()
        if fid not in area_by_id:
            excluded.append({**row, "exclude_reason": "missing_area_matrix_row"})
            continue

        manual_name = (row.get("Final_Identification") or row.get("final_identification") or "").strip()
        manual_formula = (row.get("Formula") or row.get("final_formula") or "").strip()
        if manual_name:
            final_name = manual_name
            final_formula, embedded_score = parse_final_formula_metadata(
                manual_name,
                manual_formula,
            )
            source = "manual_final_identification"
        else:
            final_name, final_formula, embedded_score = "Unidentified", "", None
            source = "unidentified_after_review"

        if is_blank_final(final_name) or is_unknown_final(final_name):
            final_name, final_formula = "Unidentified", ""
            source = "unidentified_after_review"

        feature_found_in = sum(
            fnum(area_by_id[fid].get(sample), 0.0) > 0
            for sample in sample_cols
        )
        feature_total_area = sum(fnum(area_by_id[fid].get(sample, 0)) for sample in sample_cols)
        name_key = normalize_name(final_name)
        reported_before = row.get("reported_before") or row.get("Reported_in_SOM_literature", "")
        literature_support_level = row.get("literature_support_level") or row.get("Literature_Support", "")
        if final_name == "Unidentified":
            reported_before, literature_support_level = "No", "Unavailable"
        kept.append({
            **row,
            "GlobalFeatureID": fid,
            "corrected_RT_min": row.get("corrected_RT_min") or row.get("RT_min", ""),
            "measured_RI": row.get("measured_RI") or row.get("Measured_RI", ""),
            "reference_RI": row.get("reference_RI") or row.get("Reference_RI", ""),
            "base_mz": row.get("base_mz") or row.get("Base_mz", ""),
            "characteristic_ions": (
                row.get("characteristic_ions")
                or row.get("Top_8_Ions")
                or row.get("Top_5_Ions", "")
            ),
            "diagnostic_ions": (
                row.get("Observed_Diagnostic_Ions")
                or row.get("diagnostic_ions")
                or row.get("Matched_Diagnostic_Ions", "")
            ),
            "best_identification": row.get("best_identification") or row.get("Proposed_Identification", ""),
            "best_formula": row.get("best_formula") or row.get("Formula", ""),
            "best_library": row.get("best_library") or row.get("Source_Library", ""),
            "best_score": row.get("best_score") or row.get("Match_Score", ""),
            "reported_before": reported_before,
            "literature_support_level": literature_support_level,
            "parsed_final_identification": final_name,
            "parsed_final_formula": final_formula,
            "manual_embedded_score": embedded_score,
            "final_source": source,
            "name_key": name_key,
            "rt": fnum(row.get("corrected_RT_min") or row.get("RT_min"), float("nan")),
            "feature_found_in": feature_found_in,
            "feature_total_area": feature_total_area,
        })

    merge_pairs: list[dict[str, object]] = []
    grouped: dict[str, list[dict[str, object]]] = {}
    by_name: dict[str, list[dict[str, object]]] = defaultdict(list)
    identified_rows = [
        row for row in kept
        if row["parsed_final_identification"] != "Unidentified"
    ]
    unidentified_rows = [
        row for row in kept
        if row["parsed_final_identification"] == "Unidentified"
    ]
    for row in identified_rows:
        by_name[str(row["name_key"])].append(row)
    group_number = 0
    for rows in by_name.values():
        rows.sort(key=lambda row: (
            not math.isfinite(float(row["rt"])),
            float(row["rt"]),
            str(row["GlobalFeatureID"]),
        ))
        current: list[dict[str, object]] = []
        for row in rows:
            row_rt = float(row["rt"])
            span_sec = (row_rt - float(current[0]["rt"])) * 60.0 if current and math.isfinite(row_rt) else math.inf
            if current and span_sec > args.rt_window_sec:
                grouped[f"group_{group_number}"] = current
                group_number += 1
                current = []
            for left in current:
                rt_diff = abs(row_rt - float(left["rt"])) * 60.0
                merge_pairs.append({
                    "GlobalFeatureID_1": left["GlobalFeatureID"],
                    "GlobalFeatureID_2": row["GlobalFeatureID"],
                    "final_identification": left["parsed_final_identification"],
                    "final_formula_1": left["parsed_final_formula"],
                    "final_formula_2": row["parsed_final_formula"],
                    "rt_1": left.get("corrected_RT_min", ""),
                    "rt_2": row.get("corrected_RT_min", ""),
                    "rt_diff_sec": round(rt_diff, 6),
                    "merge_rule": f"same_manual_final_name_total_span_within_{int(args.rt_window_sec)}sec",
                })
            current.append(row)
        if current:
            grouped[f"group_{group_number}"] = current
            group_number += 1

    spectra_by_id = {
        feature_id: parse_feature_spectrum(row.get("Spectra", ""))
        for feature_id, row in metadata_by_id.items()
    }
    presence_by_id = {
        feature_id: sample_presence(area_by_id[feature_id], sample_cols)
        for feature_id in area_by_id
    }
    unidentified_groups, unidentified_pairs = build_unidentified_groups(
        unidentified_rows,
        spectra_by_id,
        presence_by_id,
        args.rt_window_sec,
        args.unidentified_cosine,
        args.unidentified_max_jaccard,
    )
    for members in unidentified_groups:
        grouped[f"group_{group_number}"] = members
        group_number += 1
    merge_pairs.extend(unidentified_pairs)

    compound_rows: list[dict[str, object]] = []
    compound_area_rows: list[dict[str, object]] = []
    membership_rows: list[dict[str, object]] = []
    for idx, members in enumerate(sorted(grouped.values(), key=lambda xs: min(float(x["rt"]) for x in xs)), start=1):
        compound_id = f"{args.compound_id_prefix}{idx:0{args.compound_id_digits}d}"
        area_row: dict[str, object] = {"CompoundID": compound_id}
        for sample in sample_cols:
            area_row[sample] = sum(fnum(area_by_id[str(m["GlobalFeatureID"])].get(sample, 0)) for m in members)
        total_area = sum(fnum(area_row[sample]) for sample in sample_cols)
        found_in = sum(fnum(area_row[sample]) > 0 for sample in sample_cols)
        primary_member = choose_primary_member(members)
        reported_before, literature_support_level = strongest_literature_support(members)
        source_ids = [str(m["GlobalFeatureID"]) for m in members]
        compound_rows.append({
            "CompoundID": compound_id,
            "source_feature_ids": ";".join(source_ids),
            "n_source_features": len(members),
            "post_review_same_name_rt_merge": "yes" if len(members) > 1 else "no",
            "representative_feature_id": primary_member.get("GlobalFeatureID", ""),
            "final_identification": str(primary_member["parsed_final_identification"]),
            "final_formula": str(primary_member["parsed_final_formula"]),
            "final_InChIKey": final_inchikey(primary_member),
            "original_RT_min": primary_member.get("Original_RT_min", ""),
            "corrected_RT_min": primary_member.get("corrected_RT_min", ""),
            "measured_RI": primary_member.get("measured_RI", ""),
            "reference_RI": primary_member.get("reference_RI", ""),
            "Delta_RI": primary_member.get("Delta_RI", ""),
            "representative_base_mz": primary_member.get("base_mz", ""),
            "representative_characteristic_ions": primary_member.get("characteristic_ions", ""),
            "diagnostic_ions": primary_member.get("diagnostic_ions", ""),
            "Proposed_Identification": primary_member.get("Proposed_Identification", ""),
            "Source_Library": primary_member.get("Source_Library", ""),
            "Match_Score": primary_member.get("Match_Score", ""),
            "Library1_MoNA_Top5": primary_member.get("Library1_MoNA_Top5", ""),
            "Library2_NIST_Top5": primary_member.get("Library2_NIST_Top5", ""),
            "RI_Matched_Candidate": primary_member.get("RI_Matched_Candidate", ""),
            "Identification_Evidence": primary_member.get("Identification_Evidence", ""),
            "Matched_Diagnostic_Ions": primary_member.get("Matched_Diagnostic_Ions", ""),
            "Reference_Diagnostic_Ions": primary_member.get("Reference_Diagnostic_Ions", ""),
            "Observed_Diagnostic_Ions": primary_member.get("Observed_Diagnostic_Ions", ""),
            "Diagnostic_Ion_Match": primary_member.get("Diagnostic_Ion_Match", ""),
            "Identification_Status": primary_member.get("Identification_Status", ""),
            "Alternative_Candidates": primary_member.get("Alternative_Candidates", ""),
            "Review_Notes": primary_member.get("Review_Notes", ""),
            "Reported_in_SOM_literature": reported_before,
            "Literature_Support": literature_support_level,
            "FoundIn": found_in,
            "total_area": round(total_area, 3),
        })
        compound_area_rows.append(area_row)
        for m in members:
            membership_row = {
                "CompoundID": compound_id,
                "GlobalFeatureID": m["GlobalFeatureID"],
                **{column: m.get(column, "") for column in review_columns},
                "FoundIn": m["feature_found_in"],
                "feature_total_area": round(float(m["feature_total_area"]), 3),
                "final_source": m["final_source"],
            }
            membership_rows.append(membership_row)

    n_unidentified = sum(
        row.get("parsed_final_identification") == "Unidentified" for row in kept
    )
    summary = [{
        "review_input_rows": len(review_rows),
        "kept_reviewed_features": len(kept),
        "kept_identified_features": len(kept) - n_unidentified,
        "kept_unidentified_features": n_unidentified,
        "excluded_missing_area_matrix_features": len(excluded),
        "merge_candidate_pairs": len(merge_pairs),
        "final_compounds": len(compound_rows),
        "merged_compound_groups": sum(int(row["n_source_features"]) > 1 for row in compound_rows),
        "rt_window_sec": args.rt_window_sec,
    }]

    write_csv(args.out_dir / "checked_features_kept.csv", kept)
    write_csv(args.out_dir / "checked_features_excluded.csv", excluded)
    write_csv(args.out_dir / "same_name_close_rt_merge_pairs.csv", merge_pairs)
    write_csv(args.out_dir / "final_compounds.csv", compound_rows)
    write_csv(args.out_dir / "final_compound_area_matrix.csv", compound_area_rows)
    write_csv(args.out_dir / "final_compound_membership.csv", membership_rows)
    write_csv(args.out_dir / "post_review_merge_summary.csv", summary)
    print(summary[0])
