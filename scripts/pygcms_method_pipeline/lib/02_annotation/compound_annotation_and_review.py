#!/usr/bin/env python3
# =============================================================================
# Purpose: Orchestrate compound annotation, review checkpoints, and final outputs.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# Responsibility: Coordinate focused annotation modules and shared report/export helpers.
# =============================================================================

from __future__ import annotations

import csv
import importlib.util
import os
import shutil
import sys
from pathlib import Path


def _load_module(path):
    """Load an adjacent stage/shared module without modifying the Python search path."""
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_MODULE_DIR = Path(__file__).resolve().parent
_EVIDENCE = _load_module(_MODULE_DIR / "annotation_evidence.py")
_SPECTRAL_SEARCH = _load_module(_MODULE_DIR / "spectral_search.py")
_REVIEWED_COMPOUNDS = _load_module(_MODULE_DIR / "reviewed_compounds.py")
_REPORT = _load_module(_MODULE_DIR.parent / "shared" / "processing_report.py")
_EXPORT = _load_module(_MODULE_DIR.parent / "shared" / "export_readable_outputs.py")

_feature_rt_for_annotation = _EVIDENCE._feature_rt_for_annotation
_RI_RUNNER = _EVIDENCE.main
_POST_RUNNER = _REVIEWED_COMPOUNDS.main
_REPORT_RUNNER = _REPORT.main
_EXPORT_RUNNER = _EXPORT.main


def _build_spectral_runner():
    return _SPECTRAL_SEARCH.build_runner(_feature_rt_for_annotation)


REVIEW_AUDIT_COLUMNS = [
    "ID", "corrected_RT_min", "original_RT_min", "measured_RI", "RI_calculation_status", "RI_confidence", "reference_RI",
    "reference_RT", "base_mz", "characteristic_ions", "FoundIn", "library1",
    "library1_identification", "library1_InChIKey", "library1_formula", "library1_score",
    "library1_top5", "library2", "library2_identification", "library2_InChIKey",
    "library2_formula", "library2_score", "library2_top5", "best_identification",
    "best_InChIKey", "best_formula", "best_library", "best_score", "library_agreement",
    "automated_status", "automated_spectral_eligible", "automated_RI_ambiguity_gap",
    "automated_identification", "automated_formula", "automated_match_score",
    "automated_reference_RI", "automated_delta_RI", "automated_base_peak_status",
    "automated_compound_ion_status", "automated_matched_compound_ions", "automated_family",
    "automated_family_ion_status", "automated_matched_family_ions", "automated_M_plus_status",
    "automated_literature_support", "automated_matched_diagnostic_ions",
    "automated_reference_diagnostic_ions", "automated_observed_diagnostic_ions",
    "automated_diagnostic_ion_match",
    "automated_identification_evidence", "homologous_series_class",
    "homologous_carbon_number", "homologous_RI_residual", "homologous_series_assignment",
    "homologous_series_status", "n_RI_supported_candidates", "RI_candidate_group",
    "final_identification", "final_formula", "final_match_factor", "confidence_after_review",
    "reported_before",
]

COMPACT_REVIEW_COLUMNS = [
    "ID", "RT_min", "Original_RT_min", "Measured_RI", "Top_8_Ions", "Base_mz",
    "Proposed_Identification", "Formula", "Source_Library", "Match_Score",
    "Library1_MoNA_Top5", "Library2_NIST_Top5",
    "RI_Matched_Candidate", "Reference_RI", "Delta_RI", "Identification_Evidence",
    "Matched_Diagnostic_Ions", "Reference_Diagnostic_Ions", "Observed_Diagnostic_Ions",
    "Diagnostic_Ion_Match", "Identification_Status", "Reported_in_SOM_literature", "Literature_Support",
    "Alternative_Candidates", "Final_Identification", "Review_Notes",
]

FORMAL_CANDIDATE_EVIDENCE_FIELDS = (
    "Proposed_Identification", "Formula", "Match_Score",
    "Reference_RI", "Delta_RI", "Identification_Evidence",
    "Matched_Diagnostic_Ions", "Reference_Diagnostic_Ions",
    "Observed_Diagnostic_Ions", "Diagnostic_Ion_Match",
    "Reported_in_SOM_literature", "Literature_Support",
    "Final_Identification",
)


def _clear_formal_candidate_evidence(row):
    for field in FORMAL_CANDIDATE_EVIDENCE_FIELDS:
        row[field] = ""
    return row


def _read_csv_rows(path):
    path = Path(path)
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _reported_in_som_literature(level):
    return "Yes" if _clean(level).lower() in {"compound", "homologous series", "family"} else "No"


def _write_csv_rows(path, rows, fieldnames):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _clean(value):
    return "" if value is None else str(value).strip()


_IDENTIFICATION_STATUS_REPAIRS = {
    "Putative Compound鈥擬anual review": "Putative Compound—Manual review",
    "Putative Compound__anual review": "Putative Compound—Manual review",
    "Putative Compound鈥擧igh confidence": "Putative Compound—High confidence",
    "Putative Compound__igh confidence": "Putative Compound—High confidence",
    "Manual Review鈥擲pectral below threshold": "Manual Review—Spectral below threshold",
    "Manual Review__pectral below threshold": "Manual Review—Spectral below threshold",
    "Ambiguous Candidates鈥擬anual Review": "Ambiguous Candidates—Manual Review",
    "Ambiguous Candidates__anual Review": "Ambiguous Candidates—Manual Review",
}


def normalize_identification_status(value: object) -> str:
    status = "" if value is None else str(value)
    return _IDENTIFICATION_STATUS_REPAIRS.get(status, status)


def _number_or_none(value):
    try:
        number = float(_clean(value))
    except (TypeError, ValueError):
        return None
    return number


def _r_numeric_text(value):
    text = _clean(value)
    number = _number_or_none(text)
    return text if number is None else format(number, ".15g")


def _normalized_review_name(value):
    import re

    return re.sub(r"[^a-z0-9]+", "", _clean(value).lower().replace("–", "-").replace("—", "-"))


def _is_noncompound_final_row(row):
    level = _clean(row.get("Final_Decision_Level") or row.get("final_decision_level")).lower()
    if level:
        return level not in {"compound", "specific compound", "具体化合物"}
    status = _clean(row.get("Identification_Status")).lower()
    return status in {
        "isomer group", "putative family", "homologous series",
        "family", "unidentified", "insufficient evidence",
    }


def _candidate_reference_ions(candidate):
    return _RI_RUNNER.parse_reference_ions(
        candidate.get("reference_diagnostic_ions_text")
        or candidate.get("expected_compound_diagnostic_ions")
        or ""
    )


def _final_candidate_matches(row, candidate):
    final_key = _clean(row.get("final_InChIKey")).upper()
    candidate_keys = {
        _clean(candidate.get("InChIKey")).upper(),
        _clean(candidate.get("resolved_InChIKey")).upper(),
    } - {""}
    if final_key:
        return final_key in candidate_keys
    final_name = _normalized_review_name(
        row.get("Final_Identification")
        or row.get("final_identification")
        or row.get("parsed_final_identification")
    )
    if not final_name:
        return False
    candidate_names = {
        _normalized_review_name(candidate.get("Name")),
        _normalized_review_name(candidate.get("reference_compound")),
        _normalized_review_name(candidate.get("resolved_name")),
    } - {""}
    return final_name in candidate_names


def apply_final_top5_diagnostic_evidence(
    reviewed_rows,
    candidate_rows,
    spectral_auto_threshold=850.0,
):
    candidates_by_feature = {}
    for candidate in candidate_rows:
        candidates_by_feature.setdefault(
            _clean(candidate.get("GlobalFeatureID")), []
        ).append(candidate)

    output = []
    for source in reviewed_rows:
        row = dict(source)
        row["Reference_Diagnostic_Ions"] = ""
        row["Observed_Diagnostic_Ions"] = ""
        final_name = _normalized_review_name(
            row.get("Final_Identification")
            or row.get("final_identification")
            or row.get("parsed_final_identification")
        )
        if _is_noncompound_final_row(row) or final_name in {
            "", "unknown", "unknow", "unidentified", "notidentified", "reject", "rejected"
        }:
            row["Diagnostic_Ion_Match"] = "Not applicable — non-compound level"
            output.append(row)
            continue

        matches = [
            candidate
            for candidate in candidates_by_feature.get(_clean(row.get("ID")), [])
            if _final_candidate_matches(row, candidate)
        ]
        if not matches:
            row["Diagnostic_Ion_Match"] = "Not evaluated — identity not linked"
            output.append(row)
            continue

        matching_scores = [
            score
            for candidate in matches
            if (score := _number_or_none(candidate.get("score_0_1000"))) is not None
        ]
        if matching_scores and max(matching_scores) < spectral_auto_threshold:
            best_score = max(matching_scores)
            _clear_formal_candidate_evidence(row)
            row["Final_Identification"] = "Unidentified"
            row["Identification_Status"] = "Manual Review—Spectral below threshold"
            row["Final_Decision_Level"] = "Unidentified"
            row["Final_Decision_Source"] = "Spectral threshold gate"
            row["Final_Confidence"] = ""
            row["final_InChIKey"] = ""
            row["Review_Notes"] = (
                "Candidate retained as suggestion only: maximum matching spectral "
                f"score {_r_numeric_text(best_score)} is below "
                f"{_r_numeric_text(spectral_auto_threshold)}."
            )
            output.append(row)
            continue

        matches.sort(key=lambda candidate: (
            not bool(_candidate_reference_ions(candidate)),
            _number_or_none(candidate.get("rank")) or float("inf"),
            -(_number_or_none(candidate.get("score_0_1000")) or 0),
        ))
        reference = _candidate_reference_ions(matches[0])
        row["Reference_Diagnostic_Ions"] = ";".join(str(mz) for mz in sorted(set(reference)))
        if not reference:
            row["Diagnostic_Ion_Match"] = "Not defined"
            output.append(row)
            continue

        # Top_8_Ions is retained for manual spectral review; formal compound-ion
        # support deliberately uses only the five most abundant measured ions.
        top5 = _RI_RUNNER.parse_top_ions(
            row.get("Top_8_Ions") or row.get("Top_5_Ions"),
            limit=5,
        )
        if not top5:
            row["Diagnostic_Ion_Match"] = "Not evaluated — measured Top5 unavailable"
            output.append(row)
            continue
        observed = _RI_RUNNER.top5_diagnostic_ions(reference, top5)
        row["Observed_Diagnostic_Ions"] = ";".join(str(mz) for mz in observed)
        row["Diagnostic_Ion_Match"] = _RI_RUNNER.diagnostic_ion_match_status(
            reference,
            observed,
        )
        output.append(row)
    return output


FINAL_FEATURE_COLUMNS = [
    *COMPACT_REVIEW_COLUMNS,
    "Final_Decision_Level", "Final_Decision_Source", "Final_Confidence",
    "final_InChIKey",
]


def write_final_feature_identifications_no_merge(source_path, output_path):
    rows = _read_csv_rows(source_path)
    feature_ids = [_clean(row.get("ID") or row.get("GlobalFeatureID")) for row in rows]
    if any(not feature_id for feature_id in feature_ids):
        raise ValueError("Final no-merge feature table contains a blank feature ID")
    if len(feature_ids) != len(set(feature_ids)):
        raise ValueError("Final no-merge feature table contains duplicate feature IDs")

    output_rows = []
    for source, feature_id in zip(rows, feature_ids):
        row = dict(source)
        row["ID"] = feature_id
        row["Top_8_Ions"] = _clean(
            row.get("Top_8_Ions") or row.get("Top_5_Ions")
        )
        row["Final_Identification"] = _clean(
            row.get("Final_Identification") or row.get("parsed_final_identification")
        )
        row["Formula"] = _clean(
            row.get("parsed_final_formula") or row.get("Formula")
        )
        output_rows.append(row)
    _write_csv_rows(output_path, output_rows, FINAL_FEATURE_COLUMNS)
    return output_rows


def _library_agreement(row):
    mona_score = _number_or_none(row.get("MoNA_GC_MS_score_0_1000"))
    nist_score = _number_or_none(row.get("NIST2020_EI_score_0_1000"))
    both = mona_score is not None and nist_score is not None and mona_score >= 850 and nist_score >= 850
    if not both:
        return "not_consistent_or_low_score"
    mona_key = _clean(row.get("MoNA_GC_MS_InChIKey"))
    nist_key = _clean(row.get("NIST2020_EI_InChIKey"))
    if mona_key and nist_key and mona_key == nist_key:
        return "same_inchikey_both_ge850"
    mona_name = _clean(row.get("MoNA_GC_MS_Name"))
    nist_name = _clean(row.get("NIST2020_EI_Name"))
    if mona_name and _normalized_review_name(mona_name) == _normalized_review_name(nist_name):
        return "same_name_both_ge850"
    mona_formula = _clean(row.get("MoNA_GC_MS_Formula"))
    if mona_formula and mona_formula == _clean(row.get("NIST2020_EI_Formula")):
        return "same_formula_name_diff_both_ge850"
    return "conflict_both_ge850"


def _format_top5(top_rows, library_name):
    grouped = {}
    rows = [row for row in top_rows if _clean(row.get("library")) == library_name]
    rows.sort(key=lambda row: (_clean(row.get("GlobalFeatureID")), _number_or_none(row.get("rank")) or float("inf")))
    for row in rows:
        grouped.setdefault(_clean(row.get("GlobalFeatureID")), []).append(row)
    result = {}
    for feature_id, hits in grouped.items():
        parts = []
        for hit in hits[:5]:
            parts.append(
                f"{_r_numeric_text(hit.get('rank'))}:{_clean(hit.get('Name'))} "
                f"({_r_numeric_text(hit.get('score_0_1000'))};{_clean(hit.get('Formula'))})"
            )
        result[feature_id] = " | ".join(parts)
    return result


def _reader_status(status):
    status = _clean(status)
    lower = status.lower()
    result = "Manual Review"
    if "insufficient" in lower or "unavailable" in lower:
        result = "Insufficient Evidence"
    if "ri_conflict" in lower or "ri_mismatch" in lower:
        result = "RI Conflict"
    if status in {"isomer_group", "isomer_group_or_manual_review"}:
        result = "Isomer Group"
    if status == "ambiguous_candidates_manual_review":
        result = "Ambiguous Candidates—Manual Review"
    if status == "manual_review_spectral_below_threshold":
        result = "Manual Review—Spectral below threshold"
    if status == "putative_best_candidate_medium_confidence":
        result = "Putative Compound—Medium confidence"
    if status == "putative_best_candidate_high_confidence":
        result = "Putative Compound—High confidence"
    if status in {
        "putative_level_2", "putative_compound_multiple_evidence",
        "putative_compound_ion_supported",
        "putative_compound_RI_selected",
    }:
        result = "Putative Compound"
    if status == "putative_family_ion_supported":
        result = "Putative Family"
    if status.startswith("homologous_series_"):
        result = "Homologous Series"
    return result


def _identification_evidence(row):
    status = _clean(row.get("automated_status"))
    delta = _number_or_none(row.get("automated_delta_RI"))
    ri_supported = bool(_clean(row.get("automated_reference_RI"))) and delta is not None and abs(delta) <= 20
    compound_supported = _clean(row.get("automated_compound_ion_status")) == "supported"
    family_supported = _clean(row.get("automated_family_ion_status")) == "supported"
    result = "Spectral Library"
    if ri_supported:
        result = "Spectral + RI"
    if compound_supported:
        result = "Spectral + Compound Ions"
    if ri_supported and compound_supported:
        result = "Spectral + RI + Compound Ions"
    if family_supported:
        result = "Spectral + Family Ions"
    if ri_supported and family_supported:
        result = "Spectral + RI + Family Ions"
    if status == "putative_compound_multiple_evidence":
        result = "Spectral + RI + Compound Ions"
    if status == "putative_compound_ion_supported":
        result = "Spectral + Compound Ions"
    if status == "putative_compound_RI_selected":
        result = "Spectral + RI"
    if status == "putative_family_ion_supported":
        result = "Spectral + RI + Family Ions" if ri_supported else "Spectral + Family Ions"
    if status.startswith("homologous_series_"):
        result = "Homologous RI + Diagnostic Ions"
    if status in {"isomer_group", "isomer_group_or_manual_review"}:
        result = "Spectral + RI (Isomer Unresolved)" if ri_supported else "Spectral Library (Isomer Unresolved)"
    return _clean(row.get("automated_identification_evidence")) or result


def _combine_candidates(existing, rejected):
    parts = []
    for value in _clean(existing).split("|"):
        value = value.strip()
        if value and value not in parts:
            parts.append(value)
    rejected = _clean(rejected)
    if rejected and rejected not in parts:
        parts.append(rejected)
    return " | ".join(parts)


def _format_ri_matched_candidates(rows, spectral_thresholds=None):
    spectral_thresholds = spectral_thresholds or {}
    grouped = {}
    for row in rows:
        feature_id = _clean(row.get("GlobalFeatureID"))
        name = _clean(row.get("Name"))
        delta = _number_or_none(row.get("delta_RI"))
        score = _number_or_none(row.get("score_0_1000"))
        ri_status = _clean(row.get("RI_status"))
        support = {
            "supported": ("RI supported", 0, 20),
            "weak_support": ("RI weak support", 1, 50),
        }.get(ri_status)
        if (
            not feature_id
            or not name
            or _clean(row.get("RI_calculation_status")) != "interpolated"
            or support is None
            or delta is None
            or abs(delta) > support[2]
        ):
            continue
        threshold = spectral_thresholds.get(feature_id, 850.0)
        if score is None:
            spectral_label = "spectral eligibility unavailable; suggestion only"
        elif score >= threshold:
            spectral_label = "spectral eligible"
        else:
            spectral_label = (
                f"spectral below {_r_numeric_text(threshold)}; suggestion only"
            )
        rendered = (
            f"[{support[0]}; {spectral_label}] {name} "
            f"(score={_r_numeric_text(row.get('score_0_1000'))}; "
            f"reference RI={_r_numeric_text(row.get('reference_RI'))}; "
            f"Delta RI={_r_numeric_text(row.get('delta_RI'))})"
        )
        grouped.setdefault(feature_id, []).append(
            (support[1], abs(delta), -(score or 0), name, rendered)
        )
    return {
        feature_id: " | ".join(item[4] for item in sorted(items))
        for feature_id, items in grouped.items()
    }


AUTOACCEPTED_STATUSES = {
    "Putative Compound",
    "Putative Compound—Medium confidence",
    "Putative Compound—High confidence",
    "Putative Family",
    "Homologous Series",
    "Isomer Group",
}


def partition_annotation_rows(rows):
    autoaccepted = []
    manual_review = []
    for source in rows:
        row = dict(source)
        proposed = _clean(row.get("Proposed_Identification"))
        if _clean(row.get("Identification_Status")) in AUTOACCEPTED_STATUSES and proposed:
            row["Final_Identification"] = proposed
            row["Review_Notes"] = ""
            autoaccepted.append(row)
        else:
            row["Final_Identification"] = ""
            manual_review.append(row)
    return autoaccepted, manual_review


def combine_autoaccepted_and_reviewed(autoaccepted_path, reviewed_path, output_path):
    accepted = _read_csv_rows(autoaccepted_path) if Path(autoaccepted_path).is_file() else []
    reviewed = _read_csv_rows(reviewed_path)
    reviewed_ids = {_clean(row.get("ID")) for row in reviewed}
    rows = [row for row in accepted if _clean(row.get("ID")) not in reviewed_ids] + reviewed
    columns = list(COMPACT_REVIEW_COLUMNS)
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    rows.sort(key=lambda row: (
        _number_or_none(row.get("RT_min")) is None,
        _number_or_none(row.get("RT_min")) or 0,
        _clean(row.get("ID")),
    ))
    _write_csv_rows(output_path, rows, columns)
    return rows


def build_complete_review_checkpoint(autoaccepted_path, manual_review_path, output_path):
    rows = combine_autoaccepted_and_reviewed(autoaccepted_path, manual_review_path, output_path)
    feature_ids = [_clean(row.get("ID")) for row in rows]
    if any(not feature_id for feature_id in feature_ids):
        raise ValueError("Complete review checkpoint contains a blank feature ID")
    if len(feature_ids) != len(set(feature_ids)):
        raise ValueError("Complete review checkpoint contains duplicate feature IDs")
    for row in rows:
        row["Review_Notes"] = ""
    _write_csv_rows(output_path, rows, COMPACT_REVIEW_COLUMNS)
    return rows


def build_annotation_review_tables(
    best_hits, top_hits, auto_identification, feature_metadata, out_dir,
    candidate_evidence=None,
):
    best_rows = _read_csv_rows(best_hits)
    auto_rows = _read_csv_rows(auto_identification)
    auto_by_id = {row.get("GlobalFeatureID", ""): row for row in auto_rows}
    spectral_thresholds = {
        _clean(row.get("GlobalFeatureID")): (
            _number_or_none(row.get("spectral_auto_threshold")) or 850.0
        )
        for row in auto_rows
        if _clean(row.get("GlobalFeatureID"))
    }
    metadata_by_id = {row.get("GlobalFeatureID", ""): row for row in _read_csv_rows(feature_metadata)}
    top_rows = _read_csv_rows(top_hits)
    ri_matched_by_id = _format_ri_matched_candidates(
        _read_csv_rows(candidate_evidence) if candidate_evidence else [],
        spectral_thresholds,
    )
    mona_top5 = _format_top5(top_rows, "MoNA_GC_MS")
    nist_top5 = _format_top5(top_rows, "NIST2020_EI")
    audit = []
    for best in best_rows:
        feature_id = _clean(best.get("GlobalFeatureID"))
        merged = dict(best)
        merged.update(auto_by_id.get(feature_id, {}))
        merged.update({"TopIons": metadata_by_id.get(feature_id, {}).get("TopIons", "")})
        row = {
            "ID": feature_id,
            "corrected_RT_min": _feature_rt_for_annotation({
                **merged, **metadata_by_id.get(feature_id, {}),
            }),
            "original_RT_min": merged.get("mean_original_tmean", ""),
            "measured_RI": merged.get("measured_RI", ""),
            "RI_calculation_status": merged.get("RI_calculation_status", ""),
            "RI_confidence": merged.get("RI_confidence", ""),
            "reference_RI": merged.get("reference_RI", ""),
            "reference_RT": merged.get("reference_RT", ""),
            "base_mz": merged.get("BaseMz", ""),
            "characteristic_ions": merged.get("TopIons", ""),
            "FoundIn": merged.get("FoundIn", ""),
            "library1": "MoNA",
            "library1_identification": merged.get("MoNA_GC_MS_Name", ""),
            "library1_InChIKey": merged.get("MoNA_GC_MS_InChIKey", ""),
            "library1_formula": merged.get("MoNA_GC_MS_Formula", ""),
            "library1_score": merged.get("MoNA_GC_MS_score_0_1000", ""),
            "library1_top5": mona_top5.get(feature_id, ""),
            "library2": "NIST",
            "library2_identification": merged.get("NIST2020_EI_Name", ""),
            "library2_InChIKey": merged.get("NIST2020_EI_InChIKey", ""),
            "library2_formula": merged.get("NIST2020_EI_Formula", ""),
            "library2_score": merged.get("NIST2020_EI_score_0_1000", ""),
            "library2_top5": nist_top5.get(feature_id, ""),
            "best_identification": merged.get("best_Name", ""),
            "best_InChIKey": merged.get("best_InChIKey", ""),
            "best_formula": merged.get("best_Formula", ""),
            "best_library": merged.get("best_source", ""),
            "best_score": merged.get("best_score_0_1000", ""),
            "library_agreement": _library_agreement(merged),
            "automated_status": merged.get("automated_status", ""),
            "automated_spectral_eligible": merged.get("automatic_identification_eligible", ""),
            "automated_RI_ambiguity_gap": merged.get("RI_ambiguity_gap", ""),
            "automated_identification": merged.get("selected_candidate", ""),
            "automated_formula": merged.get("selected_formula", ""),
            "automated_match_score": merged.get("selected_match_score", ""),
            "automated_reference_RI": merged.get("selected_reference_RI", ""),
            "automated_delta_RI": merged.get("selected_delta_RI", ""),
            "automated_base_peak_status": merged.get("selected_base_peak_status", ""),
            "automated_compound_ion_status": merged.get("selected_compound_ion_status", ""),
            "automated_matched_compound_ions": merged.get("selected_matched_compound_ions", ""),
            "automated_family": merged.get("selected_family", ""),
            "automated_family_ion_status": merged.get("selected_family_ion_status", ""),
            "automated_matched_family_ions": merged.get("selected_matched_family_ions", ""),
            "automated_M_plus_status": merged.get("selected_M_plus_status", ""),
            "automated_literature_support": merged.get("selected_literature_support_level", ""),
            "automated_matched_diagnostic_ions": merged.get("selected_matched_diagnostic_ions", ""),
            "automated_reference_diagnostic_ions": merged.get(
                "selected_reference_diagnostic_ions", ""
            ),
            "automated_observed_diagnostic_ions": merged.get(
                "selected_observed_diagnostic_ions", ""
            ),
            "automated_diagnostic_ion_match": merged.get(
                "selected_diagnostic_ion_match", ""
            ),
            "automated_identification_evidence": merged.get("selected_identification_evidence", ""),
            "homologous_series_class": merged.get("homologous_series_class", ""),
            "homologous_carbon_number": merged.get("homologous_carbon_number", ""),
            "homologous_RI_residual": merged.get("homologous_RI_residual", ""),
            "homologous_series_assignment": merged.get("homologous_series_assignment", ""),
            "homologous_series_status": merged.get("homologous_series_status", ""),
            "n_RI_supported_candidates": merged.get("n_RI_supported_candidates", ""),
            "RI_candidate_group": merged.get("candidate_group", ""),
            "final_identification": "", "final_formula": "", "final_match_factor": "",
            "confidence_after_review": "",
            "reported_before": _reported_in_som_literature(merged.get("selected_literature_support_level", "")),
        }
        numeric_columns = {
            "corrected_RT_min", "original_RT_min", "measured_RI", "reference_RI", "reference_RT", "base_mz",
            "FoundIn", "library1_score", "library2_score", "best_score",
            "automated_RI_ambiguity_gap", "automated_match_score", "automated_reference_RI",
            "automated_delta_RI", "homologous_carbon_number", "homologous_RI_residual",
            "n_RI_supported_candidates",
        }
        audit.append({
            key: _r_numeric_text(row.get(key, "")) if key in numeric_columns else _clean(row.get(key, ""))
            for key in REVIEW_AUDIT_COLUMNS
        })
    audit.sort(key=lambda row: (
        _number_or_none(row.get("corrected_RT_min")) is None,
        _number_or_none(row.get("corrected_RT_min")) or 0,
        row.get("ID", ""),
    ))
    compact = []
    for row in audit:
        status = _reader_status(row.get("automated_status"))
        homolog = status == "Homologous Series"
        isomer = status == "Isomer Group"
        rejected = status in {"RI Conflict", "Insufficient Evidence"}
        proposed = row.get("homologous_series_assignment") or row.get("automated_identification")
        if isomer:
            proposed = row.get("RI_candidate_group", "")
        formula = row.get("automated_formula", "")
        rejected_candidate = row.get("automated_identification", "") if rejected else ""
        is_autoaccepted = status in AUTOACCEPTED_STATUSES
        match_score = row.get("automated_match_score", "") or row.get("best_score", "")
        if not is_autoaccepted:
            proposed = ""
            formula = ""
            match_score = ""
        measured_ri = row.get("measured_RI", "")
        reference_ri = row.get("automated_reference_RI", "")
        delta_ri = row.get("automated_delta_RI", "")
        carbon = _number_or_none(row.get("homologous_carbon_number"))
        if homolog and carbon is not None:
            reference_ri = str(int(100 * carbon)) if (100 * carbon).is_integer() else str(100 * carbon)
        if homolog:
            delta_ri = row.get("homologous_RI_residual", "")
        alternatives = "" if isomer else _combine_candidates(row.get("RI_candidate_group", ""), rejected_candidate)
        if status not in AUTOACCEPTED_STATUSES and not alternatives:
            alternatives = _combine_candidates(row.get("library1_top5", ""), row.get("library2_top5", ""))
        compact_row = {
            "ID": row.get("ID", ""), "RT_min": row.get("corrected_RT_min", ""),
            "Original_RT_min": row.get("original_RT_min", ""),
            "Measured_RI": measured_ri, "Top_8_Ions": row.get("characteristic_ions", ""),
            "Base_mz": row.get("base_mz", ""), "Proposed_Identification": proposed,
            "Formula": formula, "Source_Library": row.get("best_library", ""),
            "Match_Score": match_score,
            "Library1_MoNA_Top5": row.get("library1_top5", ""),
            "Library2_NIST_Top5": row.get("library2_top5", ""),
            "RI_Matched_Candidate": ri_matched_by_id.get(row.get("ID", ""), ""),
            "Reference_RI": reference_ri,
            "Delta_RI": delta_ri, "Identification_Evidence": _identification_evidence(row),
            "Matched_Diagnostic_Ions": row.get("automated_matched_diagnostic_ions", ""),
            "Reference_Diagnostic_Ions": row.get(
                "automated_reference_diagnostic_ions", ""
            ),
            "Observed_Diagnostic_Ions": row.get(
                "automated_observed_diagnostic_ions", ""
            ),
            "Diagnostic_Ion_Match": row.get(
                "automated_diagnostic_ion_match", ""
            ),
            "Identification_Status": status,
            "Reported_in_SOM_literature": _reported_in_som_literature(row.get("automated_literature_support", "")),
            "Literature_Support": row.get("automated_literature_support", ""),
            "Alternative_Candidates": alternatives,
            "Final_Identification": "", "Review_Notes": "",
        }
        selected_score = _number_or_none(row.get("automated_match_score"))
        threshold = spectral_thresholds.get(row.get("ID", ""), 850.0)
        if (
            _clean(row.get("automated_status")) == "manual_review_spectral_below_threshold"
            or (selected_score is not None and selected_score < threshold)
        ):
            _clear_formal_candidate_evidence(compact_row)
        compact.append(compact_row)
    autoaccepted, manual_review = partition_annotation_rows(compact)
    out_dir = Path(out_dir)
    _write_csv_rows(out_dir / "annotation_evidence_full.csv", audit, REVIEW_AUDIT_COLUMNS)
    _write_csv_rows(out_dir / "annotation_autoaccepted.csv", autoaccepted, COMPACT_REVIEW_COLUMNS)
    _write_csv_rows(out_dir / "compound_identification_review.csv", manual_review, COMPACT_REVIEW_COLUMNS)
    agreement_counts = {}
    for row in audit:
        agreement_counts[row["library_agreement"]] = agreement_counts.get(row["library_agreement"], 0) + 1
    summary = [{"summary_type": "total", "category": "all_features", "n": len(audit)}]
    summary.extend(
        {"summary_type": "library_agreement", "category": key, "n": agreement_counts[key]}
        for key in sorted(agreement_counts)
    )
    summary.extend([
        {"summary_type": "workflow", "category": "automatically_accepted", "n": len(autoaccepted)},
        {"summary_type": "workflow", "category": "manual_review_required", "n": len(manual_review)},
        {"summary_type": "manual_columns", "category": "final_identification_blank", "n": len(manual_review)},
        {"summary_type": "manual_columns", "category": "review_notes_blank", "n": len(manual_review)},
    ])
    _write_csv_rows(out_dir / "compound_identification_review_summary.csv", summary, ["summary_type", "category", "n"])
    return manual_review, audit


def _invoke_runner(runner, argv):
    previous = sys.argv
    try:
        sys.argv = [previous[0], *[str(value) for value in argv]]
        return runner()
    finally:
        sys.argv = previous


def _copy_if_exists(source, target):
    source = Path(source)
    if source.is_file():
        Path(target).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


def _pipeline_paths():
    work = Path(os.environ["WORK_DIR"])
    main = Path(os.environ["MAIN_OUTPUT_DIR"])
    group = work / "04_global_grouping"
    production_ri_metadata = Path(os.environ.get("REFERENCE_RI_OUTPUT_DIR") or work / "05_reference_ri_mapping") / "features_for_annotation.csv"
    metadata_override = os.environ.get("ANNOTATION_FEATURE_METADATA", "").strip()
    matrix_override = os.environ.get("ANNOTATION_FEATURE_AREA_MATRIX", "").strip()
    metadata = (
        Path(metadata_override)
        if metadata_override
        else (
            production_ri_metadata
            if production_ri_metadata.is_file()
            else group / "features_for_annotation.csv"
        )
    )
    matrix = (
        Path(matrix_override)
        if matrix_override
        else group / "feature_area_matrix.csv"
    )
    nearby = work / "audit_files" / "similar_feature_pairs_after_grouping.csv"
    return work, main, metadata, matrix, nearby


def enrich_review_with_confirmed_identity(reviewed_path, audit_path, output_path):
    reviewed_rows = _read_csv_rows(reviewed_path)
    audit_by_id = {row.get("ID", ""): row for row in _read_csv_rows(audit_path)}
    fieldnames = list(reviewed_rows[0].keys()) if reviewed_rows else []
    if "final_InChIKey" not in fieldnames:
        fieldnames.append("final_InChIKey")
    for row in reviewed_rows:
        final_name = _clean(row.get("Final_Identification") or row.get("final_identification"))
        normalized_final = _normalized_review_name(final_name)
        confirmed_key = ""
        if normalized_final and normalized_final not in {
            "unknown", "unknow", "unidentified", "notidentified", "reject", "rejected"
        }:
            audit = audit_by_id.get(_clean(row.get("ID")), {})
            for name_column, key_column in (
                ("best_identification", "best_InChIKey"),
                ("library1_identification", "library1_InChIKey"),
                ("library2_identification", "library2_InChIKey"),
            ):
                if normalized_final == _normalized_review_name(audit.get(name_column, "")):
                    confirmed_key = _clean(audit.get(key_column, ""))
                    if confirmed_key:
                        break
        row["final_InChIKey"] = confirmed_key
    _write_csv_rows(output_path, reviewed_rows, fieldnames)
    return reviewed_rows


def _validate_authoritative_review_csv_headers(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        header = next(csv.reader(handle), [])
    seen = set()
    duplicates = []
    blank_positions = []
    for column_position, name in enumerate(header, start=1):
        if not name.strip():
            blank_positions.append(column_position)
            continue
        if name in seen and name not in duplicates:
            duplicates.append(name)
        seen.add(name)
    problems = []
    if blank_positions:
        problems.append(
            "blank headers at 1-based column positions: "
            + ", ".join(str(position) for position in blank_positions)
        )
    if duplicates:
        problems.append("duplicate headers: " + ", ".join(duplicates))
    if problems:
        raise ValueError(
            "Authoritative review CSV has " + "; ".join(problems)
        )


def prepare_review_for_finalization(
    reviewed_path, audit_path, output_path, feature_metadata_path=None,
):
    _validate_authoritative_review_csv_headers(reviewed_path)
    rows = enrich_review_with_confirmed_identity(
        reviewed_path,
        audit_path,
        output_path,
    )
    original_rt_by_id = {}
    if feature_metadata_path and Path(feature_metadata_path).is_file():
        original_rt_by_id = {
            _clean(row.get("GlobalFeatureID")): _clean(row.get("mean_original_tmean"))
            for row in _read_csv_rows(feature_metadata_path)
            if _clean(row.get("GlobalFeatureID"))
        }
    columns = []
    for row in rows:
        row["Top_8_Ions"] = _clean(
            row.get("Top_8_Ions") or row.get("Top_5_Ions")
        )
        row.pop("Top_5_Ions", None)
        if not _clean(row.get("Original_RT_min")):
            row["Original_RT_min"] = original_rt_by_id.get(_clean(row.get("ID")), "")
        row["Identification_Status"] = normalize_identification_status(
            row.get("Identification_Status", "")
        )
        for key in row:
            if key not in columns:
                columns.append(key)
    if "Original_RT_min" in columns:
        columns.remove("Original_RT_min")
        rt_index = columns.index("RT_min") + 1 if "RT_min" in columns else 1
        columns.insert(rt_index, "Original_RT_min")
    if "Top_8_Ions" in columns:
        columns.remove("Top_8_Ions")
        ion_index = columns.index("Measured_RI") + 1 if "Measured_RI" in columns else 1
        columns.insert(ion_index, "Top_8_Ions")
    _write_csv_rows(output_path, rows, columns)
    return rows


def select_phase(main_output_dir):
    from pathlib import Path

    return (
        "review_finalization"
        if Path(main_output_dir, "annotation_review.csv").is_file()
        else "automatic_annotation"
    )


def run_automatic_annotation():
    work, main_output, metadata, matrix, nearby = _pipeline_paths()
    existing_check = main_output / "annotation_review_CHECK.csv"
    existing_work_check = work / "06_manual_review" / "compound_identification_review_CHECK.csv"
    if existing_check.is_file() or existing_work_check.is_file():
        print("\nManual-review checkpoint already exists; automatic annotation was not rerun or overwritten.")
        print(existing_check if existing_check.is_file() else existing_work_check)
        print("Save the completed table as:")
        print(main_output / "annotation_review.csv")
        print("Then run Stage 02 again.")
        return 0
    if not metadata.is_file() or not matrix.is_file():
        raise FileNotFoundError("Annotation inputs are missing; run Stage 01 first")
    annotation = work / "05_annotation"
    review = work / "06_manual_review"
    annotation.mkdir(parents=True, exist_ok=True)
    review.mkdir(parents=True, exist_ok=True)
    os.environ.update({
        "FEATURE_METADATA": str(metadata), "FEATURE_QC": str(metadata),
        "FEATURE_AREA_MATRIX": str(matrix), "RESIDUAL_FLAGS": str(nearby),
        "SPECTRAL_LIBRARY_MODE": "both", "ANNOTATION_OUT_DIR": str(annotation),
        "ANNOTATION_FOUNDIN_GT": str(int(os.environ.get("FOUNDIN_MIN", "3")) - 1),
    })
    _build_spectral_runner()()
    generated = {
        "04_foundin_gt3_mona_nist_best_hits.csv": "annotation_best_hits.csv",
        "04_foundin_gt3_mona_nist_top_hits.csv": "annotation_top_hits.csv",
        "04_foundin_gt3_mona_nist_annotation_summary.csv": "annotation_summary.csv",
        "04_foundin_gt3_mona_nist_library_search_summary.csv": "annotation_library_search_summary.csv",
    }
    for old, new in generated.items():
        source = annotation / old
        target = annotation / new
        if target.exists():
            target.unlink()
        source.replace(target)
    ri_args = [
        "--top-hits", annotation / "annotation_top_hits.csv",
        "--feature-metadata", metadata, "--alkanes", os.environ["ALKANE_RI_FILE"],
        "--reference-workbook", os.environ["SOM_REFERENCE_FILE"],
        "--candidate-output", annotation / "annotation_candidate_evidence.csv",
        "--decision-output", annotation / "annotation_automated_identification.csv",
        "--cache-dir", os.environ["RI_CACHE_DIR"],
        "--pubchem-cache-dir", os.environ["PUBCHEM_IDENTITY_CACHE_DIR"],
        "--support-window", os.environ["RI_SUPPORT_WINDOW"],
        "--weak-window", os.environ["RI_WEAK_WINDOW"],
        "--spectral-auto-threshold", os.environ["SPECTRAL_AUTO_THRESHOLD"],
        "--decision-policy", os.environ["ANNOTATION_DECISION_POLICY"],
        "--primary-top-n", os.environ["PRIMARY_TOP_HITS_PER_LIBRARY"],
        "--rescue-top-n", os.environ["RESCUE_TOP_HITS_PER_LIBRARY"],
        "--timeout", os.environ["RI_QUERY_TIMEOUT_SEC"], "--delay", os.environ["RI_QUERY_DELAY_SEC"],
        "--pubchem-timeout", os.environ["PUBCHEM_QUERY_TIMEOUT_SEC"],
        "--pubchem-delay", os.environ["PUBCHEM_QUERY_DELAY_SEC"],
    ]
    if os.environ.get("RI_OFFLINE", "false").lower() == "true":
        ri_args.append("--offline")
    if os.environ.get("PUBCHEM_OFFLINE", "false").lower() == "true":
        ri_args.append("--pubchem-offline")
    _invoke_runner(_RI_RUNNER, ri_args)
    build_annotation_review_tables(
        annotation / "annotation_best_hits.csv", annotation / "annotation_top_hits.csv",
        annotation / "annotation_automated_identification.csv", metadata, review,
        annotation / "annotation_candidate_evidence.csv",
    )
    check = review / "compound_identification_review_CHECK.csv"
    build_complete_review_checkpoint(
        review / "annotation_autoaccepted.csv",
        review / "compound_identification_review.csv",
        check,
    )
    _copy_if_exists(check, main_output / "annotation_review_CHECK.csv")
    _copy_if_exists(review / "annotation_autoaccepted.csv", main_output / "annotation_autoaccepted.csv")
    _copy_if_exists(review / "annotation_evidence_full.csv", main_output / "annotation_evidence_full.csv")
    print("\nAutomatic annotation completed. Only unresolved features require manual review:")
    print(main_output / "annotation_review_CHECK.csv")
    print("Save the reviewed table as:")
    print(main_output / "annotation_review.csv")
    print("Then run Stage 02 again.")
    return 0


def recompute_automatic_annotation_from_cached_hits():
    work, main_output, metadata, _matrix, _nearby = _pipeline_paths()
    annotation = work / "05_annotation"
    review = work / "06_manual_review"
    best_hits = annotation / "annotation_best_hits.csv"
    top_hits = annotation / "annotation_top_hits.csv"
    missing = [path for path in (metadata, best_hits, top_hits) if not path.is_file()]
    if missing:
        rendered = "\n".join(str(path) for path in missing)
        raise FileNotFoundError(
            "Cached annotation inputs are missing; run automatic Stage 02 first:\n"
            f"{rendered}"
        )

    candidate_evidence = annotation / "annotation_candidate_evidence.csv"
    automated_identification = annotation / "annotation_automated_identification.csv"
    ri_args = [
        "--top-hits", top_hits,
        "--feature-metadata", metadata, "--alkanes", os.environ["ALKANE_RI_FILE"],
        "--reference-workbook", os.environ["SOM_REFERENCE_FILE"],
        "--candidate-output", candidate_evidence,
        "--decision-output", automated_identification,
        "--cache-dir", os.environ["RI_CACHE_DIR"],
        "--pubchem-cache-dir", os.environ["PUBCHEM_IDENTITY_CACHE_DIR"],
        "--support-window", os.environ["RI_SUPPORT_WINDOW"],
        "--weak-window", os.environ["RI_WEAK_WINDOW"],
        "--spectral-auto-threshold", os.environ["SPECTRAL_AUTO_THRESHOLD"],
        "--decision-policy", os.environ["ANNOTATION_DECISION_POLICY"],
        "--primary-top-n", os.environ["PRIMARY_TOP_HITS_PER_LIBRARY"],
        "--rescue-top-n", os.environ["RESCUE_TOP_HITS_PER_LIBRARY"],
        "--timeout", os.environ["RI_QUERY_TIMEOUT_SEC"],
        "--delay", os.environ["RI_QUERY_DELAY_SEC"],
        "--pubchem-timeout", os.environ["PUBCHEM_QUERY_TIMEOUT_SEC"],
        "--pubchem-delay", os.environ["PUBCHEM_QUERY_DELAY_SEC"],
    ]
    if os.environ.get("RI_OFFLINE", "false").lower() == "true":
        ri_args.append("--offline")
    if os.environ.get("PUBCHEM_OFFLINE", "false").lower() == "true":
        ri_args.append("--pubchem-offline")

    _invoke_runner(_RI_RUNNER, ri_args)
    build_annotation_review_tables(
        best_hits,
        top_hits,
        automated_identification,
        metadata,
        review,
        candidate_evidence,
    )
    check = review / "compound_identification_review_CHECK.csv"
    build_complete_review_checkpoint(
        review / "annotation_autoaccepted.csv",
        review / "compound_identification_review.csv",
        check,
    )
    _copy_if_exists(check, main_output / "annotation_review_CHECK.csv")
    _copy_if_exists(
        review / "annotation_autoaccepted.csv",
        main_output / "annotation_autoaccepted.csv",
    )
    _copy_if_exists(
        review / "annotation_evidence_full.csv",
        main_output / "annotation_evidence_full.csv",
    )
    print("\nAutomatic identification recomputed from cached Top20 spectral hits.")
    print(f"Review checkpoint: {main_output / 'annotation_review_CHECK.csv'}")
    print(f"Candidate evidence: {candidate_evidence}")
    print(f"Automated decisions: {automated_identification}")
    return 0


def finalize_reviewed_annotation():
    work, main_output, metadata, matrix, _nearby = _pipeline_paths()
    reviewed = main_output / "annotation_review.csv"
    if not reviewed.is_file():
        raise FileNotFoundError(f"Missing reviewed table: {reviewed}")
    if not matrix.is_file():
        raise FileNotFoundError(f"Missing feature area matrix: {matrix}")
    merge_dir = work / "07_post_review_merge"
    audit = work / "06_manual_review" / "annotation_evidence_full.csv"
    if not audit.is_file():
        audit = main_output / "annotation_evidence_full.csv"
    enriched_review = merge_dir / "review_for_finalization_enriched.csv"
    prepare_review_for_finalization(
        reviewed,
        audit,
        enriched_review,
        metadata,
    )
    _invoke_runner(_POST_RUNNER, [
        "--review", enriched_review, "--area-matrix", matrix, "--out-dir", merge_dir,
        "--feature-metadata", metadata,
        "--rt-window-sec", os.environ["POST_REVIEW_MERGE_RT_SEC"],
        "--unidentified-cosine", os.environ["POST_REVIEW_UNIDENTIFIED_COSINE"],
        "--unidentified-max-jaccard", os.environ["POST_REVIEW_UNIDENTIFIED_MAX_JACCARD"],
        "--compound-id-prefix", os.environ["COMPOUND_ID_PREFIX"],
        "--compound-id-digits", os.environ["COMPOUND_ID_DIGITS"],
    ])
    final_feature_table = main_output / "final_feature_identifications_no_merge.csv"
    write_final_feature_identifications_no_merge(
        merge_dir / "checked_features_kept.csv",
        final_feature_table,
    )
    _invoke_runner(_REPORT_RUNNER, ["--run-dir", work])
    _invoke_runner(_EXPORT_RUNNER, [
        "--run-dir", work, "--output-dir", main_output,
        "--foundin-min", os.environ.get("FOUNDIN_MIN", "3"),
    ])
    print("\nManual review finalized.")
    print(f"Final compound table: {main_output / 'final_compounds.csv'}")
    print(f"Final no-merge feature table: {final_feature_table}")
    return 0


def main():
    phase = select_phase(Path(os.environ["MAIN_OUTPUT_DIR"]))
    if os.environ.get("ANNOTATION_RECOMPUTE_FROM_CACHE", "false").lower() == "true":
        recompute_automatic_annotation_from_cached_hits()
        return 0
    if phase == "review_finalization":
        return finalize_reviewed_annotation()
    return run_automatic_annotation()


if __name__ == "__main__":
    raise SystemExit(main())
