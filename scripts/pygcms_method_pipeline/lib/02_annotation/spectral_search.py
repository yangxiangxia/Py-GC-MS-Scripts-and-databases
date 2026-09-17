#!/usr/bin/env python3
# =============================================================================
# Purpose: Search MoNA and NIST EI libraries against aligned feature spectra.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# Responsibility: Bind run settings, rank library candidates, and write spectral search tables.
# =============================================================================

from __future__ import annotations


def build_runner(_feature_rt_for_annotation):
    """Bind current input paths/settings and the shared RT selector for one search."""

    import csv
    import gc
    import heapq
    import math
    import os
    import re
    from collections import defaultdict
    from pathlib import Path
    from typing import Dict, Iterable, List, Tuple


    # Stage 02 supplies current run inputs; library defaults live in config.sh.
    PROJECT = Path(os.environ.get("PROJECT_DIR", Path(__file__).resolve().parents[4]))
    FEATURE_METADATA = Path(os.environ["FEATURE_METADATA"])
    FEATURE_QC = Path(os.environ["FEATURE_QC"])
    FEATURE_AREA_MATRIX = Path(os.environ["FEATURE_AREA_MATRIX"])
    RESIDUAL_FLAGS = Path(os.environ["RESIDUAL_FLAGS"])

    MONA_FILE = Path(os.environ["MONA_FILE"])
    NIST_FILES = [Path(value) for value in os.environ["NIST_FILES"].split(";") if value.strip()]
    SPECTRAL_LIBRARY_MODE = os.environ.get("SPECTRAL_LIBRARY_MODE", "both").strip().lower()

    OUT_DIR = Path(os.environ["ANNOTATION_OUT_DIR"])
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    MZ_MIN = 46
    MZ_MAX = 650
    FOUNDIN_MIN_EXCLUSIVE = int(os.environ.get("ANNOTATION_FOUNDIN_GT", "3"))
    QUERY_TOP_N_IONS = int(os.environ.get("QUERY_TOP_N_IONS", "20"))
    TOP_HITS_PER_LIBRARY = int(os.environ.get("TOP_HITS_PER_LIBRARY", "5"))
    MAX_POSTINGS_PER_MZ = int(os.environ.get("MAX_POSTINGS_PER_MZ", "120000"))

    Spectrum = Dict[int, float]

    KEEP_FIELDS = {
        "Name", "Synon", "DB#", "CAS", "CASNO", "Formula", "MW", "ExactMass",
        "InChIKey", "InChI", "SMILES", "Spectrum_type", "Instrument_type",
        "Instrument", "Ion_mode", "Comments", "Comment", "Splash", "Retention_index",
    }


    def library_jobs_for_mode(mode: str, mona_file: Path, nist_files: list[Path]) -> list[tuple[list[Path], str, bool]]:
        jobs = {
            "mona": [([mona_file], "MoNA_GC_MS", True)],
            "nist": [(nist_files, "NIST2020_EI", False)],
            "both": [
                ([mona_file], "MoNA_GC_MS", True),
                (nist_files, "NIST2020_EI", False),
            ],
        }
        if mode not in jobs:
            raise ValueError("SPECTRAL_LIBRARY_MODE must be one of: mona, nist, both")
        return jobs[mode]


    def add_peak(spec: Spectrum, mz: float, intensity: float) -> None:
        if not math.isfinite(mz) or not math.isfinite(intensity) or intensity <= 0:
            return
        mz_i = int(round(mz))
        if MZ_MIN <= mz_i <= MZ_MAX:
            spec[mz_i] = spec.get(mz_i, 0.0) + float(intensity)


    def normalize(spec: Spectrum) -> Spectrum:
        norm = math.sqrt(sum(v * v for v in spec.values()))
        if norm <= 0:
            return {}
        return {mz: inten / norm for mz, inten in spec.items()}


    def keep_top_ions(spec: Spectrum, n: int = QUERY_TOP_N_IONS) -> Spectrum:
        if not spec or len(spec) <= n:
            return spec
        top = dict(sorted(spec.items(), key=lambda kv: (-kv[1], kv[0]))[:n])
        return normalize(top)


    def spectral_audit_metrics(
        query: Spectrum,
        library: Spectrum,
        measured_base_mz: int | None,
        spectral_cosine: float | None = None,
    ) -> dict[str, object]:
        """Describe query coverage by a library spectrum without affecting ranking."""
        query_mz = set(query)
        library_mz = set(library)
        matched_mz = query_mz & library_mz
        query_ion_count = len(query_mz)
        query_intensity = sum(query.values())
        matched_intensity = sum(query[mz] for mz in matched_mz)
        matched_ion_fraction = len(matched_mz) / query_ion_count if query_ion_count else 0.0
        matched_intensity_fraction = matched_intensity / query_intensity if query_intensity > 0 else 0.0
        base_peak_matched = measured_base_mz is not None and measured_base_mz in library_mz
        possible_coelution = (
            spectral_cosine is not None
            and spectral_cosine >= 0.80
            and matched_intensity_fraction < 0.60
        )
        return {
            "matched_ion_count": len(matched_mz),
            "query_ion_count": query_ion_count,
            "matched_ion_fraction": matched_ion_fraction,
            "matched_intensity_fraction": matched_intensity_fraction,
            "base_peak_matched": base_peak_matched,
            "possible_coelution": possible_coelution,
        }


    def parse_peak_pairs(line: str) -> List[Tuple[float, float]]:
        out: List[Tuple[float, float]] = []
        for mz, inten in re.findall(r"(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s*;?", line):
            try:
                out.append((float(mz), float(inten)))
            except ValueError:
                pass
        return out


    def is_ei_gc_like(fields: dict) -> bool:
        text = " ".join(
            str(fields.get(k, ""))
            for k in ("Spectrum_type", "Instrument_type", "Instrument", "Ion_mode", "Comments", "Comment")
        ).lower()
        spectrum_type = str(fields.get("Spectrum_type", "")).lower()
        instrument_type = str(fields.get("Instrument_type", "")).lower()
        if "ms2" in spectrum_type or "lc-" in instrument_type or "esi" in instrument_type:
            return False
        return (
            "ei" in text
            or "gc-" in text
            or "gc/ms" in text
            or "gc-ms" in text
            or "electron impact" in text
        )


    def emit_record(fields: dict, spec: Spectrum, library: str, source_file: str, record_id: int, require_ei_gc: bool) -> dict | None:
        if not fields.get("Name") or not spec:
            return None
        if require_ei_gc and not is_ei_gc_like(fields):
            return None
        rec = dict(fields)
        match = re.search(r"(?:^|\s)SemiStdNP=([-+]?\d+(?:\.\d+)?)", str(fields.get("Retention_index", "")), re.I)
        rec["SemiStdNP_RI"] = float(match.group(1)) if match else None
        rec["spectrum"] = normalize(spec)
        rec["n_peaks"] = len(spec)
        rec["library"] = library
        rec["source_file"] = source_file
        rec["library_record_id"] = record_id
        return rec


    def parse_msp_records(path: Path, library: str, require_ei_gc: bool = False) -> Iterable[dict]:
        fields: dict = {}
        spec: Spectrum = {}
        in_peaks = False
        record_id = 0

        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for raw in handle:
                line = raw.strip()
                if not line:
                    rec = emit_record(fields, spec, library, path.name, record_id, require_ei_gc)
                    if rec is not None:
                        yield rec
                    fields = {}
                    spec = {}
                    in_peaks = False
                    continue

                if ":" in line and not re.match(r"^\d", line):
                    key, value = line.split(":", 1)
                    key = key.strip()
                    value = value.strip()
                    if key.lower() == "num peaks":
                        in_peaks = True
                        record_id += 1
                        continue
                    if key == "Synon" and key in fields:
                        fields[key] = fields[key] + " | " + value
                    elif key in KEEP_FIELDS:
                        fields[key] = value
                    continue

                if in_peaks or re.match(r"^\d", line):
                    for mz, inten in parse_peak_pairs(line):
                        add_peak(spec, mz, inten)

        rec = emit_record(fields, spec, library, path.name, record_id, require_ei_gc)
        if rec is not None:
            yield rec


    def parse_feature_spectrum(value: str) -> Spectrum:
        spec: Spectrum = {}
        for token in str(value or "").split():
            if "," not in token:
                continue
            mz_s, int_s = token.split(",", 1)
            try:
                add_peak(spec, float(mz_s), float(int_s))
            except ValueError:
                continue
        return normalize(spec)


    def compact_hit(rec: dict, score: float, audit_metrics: dict[str, object] | None = None) -> dict:
        hit = {
            "score": score,
            "Name": rec.get("Name", ""),
            "InChIKey": rec.get("InChIKey", ""),
            "CAS": clean_cas(rec),
            "CASNO": rec.get("CASNO", ""),
            "Formula": rec.get("Formula", ""),
            "MW": rec.get("MW", ""),
            "source_file": rec.get("source_file", ""),
            "library_record_id": rec.get("library_record_id", ""),
            "SemiStdNP_RI": rec.get("SemiStdNP_RI"),
        }
        hit.update(audit_metrics or {})
        return hit


    def compound_identity(rec: dict) -> tuple[str, str]:
        inchikey = str(rec.get("InChIKey", "")).strip().upper()
        if inchikey:
            return ("inchikey", inchikey)
        cas = clean_cas(rec).strip()
        if cas:
            return ("cas", cas)
        name_formula = f"{norm_text(rec.get('Name', ''))}|{str(rec.get('Formula', '')).strip()}"
        if name_formula != "|":
            return ("name_formula", name_formula)
        return ("record", f"{rec.get('source_file', '')}:{rec.get('library_record_id', '')}")


    def stream_search_library(features: list[dict], paths: list[Path], library: str, require_ei_gc: bool, top_n: int) -> tuple[dict[tuple[str, str], dict], list[dict], dict]:
        query_index: dict[int, list[tuple[int, float]]] = defaultdict(list)
        feature_ids: list[str] = []
        for idx, feature in enumerate(features):
            feature_ids.append(feature["GlobalFeatureID"])
            for mz, weight in feature["query_spectrum"].items():
                query_index[mz].append((idx, weight))

        candidates: dict[int, dict[tuple[str, str], dict]] = {idx: {} for idx in range(len(features))}
        compound_meta: dict[tuple[str, str], dict] = {}
        total = 0
        searched = 0
        searched_files = []
        for path in paths:
            searched_files.append(str(path))
            for rec in parse_msp_records(path, library, require_ei_gc=require_ei_gc):
                total += 1
                identity = compound_identity(rec)
                meta = compound_meta.setdefault(identity, {"ri": [], "files": set(), "records": 0})
                meta["records"] += 1
                meta["files"].add(rec.get("source_file", ""))
                if rec.get("SemiStdNP_RI") is not None:
                    meta["ri"].append(float(rec["SemiStdNP_RI"]))
                if len(rec["spectrum"]) < 3:
                    continue
                searched += 1
                scores: dict[int, float] = defaultdict(float)
                for mz, lib_weight in rec["spectrum"].items():
                    for q_idx, q_weight in query_index.get(mz, ()):
                        scores[q_idx] += lib_weight * q_weight
                for q_idx, score in scores.items():
                    if score <= 0:
                        continue

                    def audited_hit() -> dict:
                        base_value = features[q_idx].get("BaseMz")
                        try:
                            measured_base_mz = int(round(float(base_value)))
                        except (TypeError, ValueError):
                            measured_base_mz = None
                        audit_metrics = spectral_audit_metrics(
                            features[q_idx]["query_spectrum"],
                            rec["spectrum"],
                            measured_base_mz,
                            spectral_cosine=score,
                        )
                        return compact_hit(rec, score, audit_metrics)

                    current = candidates[q_idx].get(identity)
                    if current is not None:
                        if score > current["score"]:
                            candidates[q_idx][identity] = audited_hit()
                        continue
                    if len(candidates[q_idx]) < top_n:
                        candidates[q_idx][identity] = audited_hit()
                    else:
                        worst_key = min(candidates[q_idx], key=lambda key: candidates[q_idx][key]["score"])
                        if score > candidates[q_idx][worst_key]["score"]:
                            del candidates[q_idx][worst_key]
                            candidates[q_idx][identity] = audited_hit()
                if searched % 50000 == 0:
                    print(f"  {library}: searched {searched} records", flush=True)

        best: dict[tuple[str, str], dict] = {}
        top_rows: list[dict] = []
        for idx, feature in enumerate(features):
            fid = feature_ids[idx]
            hits = sorted(candidates[idx].items(), key=lambda item: item[1]["score"], reverse=True)
            for rank, (identity, hit) in enumerate(hits, start=1):
                score = hit["score"]
                meta = compound_meta[identity]
                ri_values = sorted(meta["ri"])
                hit["embedded_SemiStdNP_values"] = ";".join(f"{v:g}" for v in ri_values)
                hit["embedded_SemiStdNP_count"] = len(ri_values)
                hit["embedded_SemiStdNP_RI"] = "" if not ri_values else (ri_values[len(ri_values)//2] if len(ri_values)%2 else (ri_values[len(ri_values)//2-1]+ri_values[len(ri_values)//2])/2)
                hit["embedded_SemiStdNP_min"] = "" if not ri_values else min(ri_values)
                hit["embedded_SemiStdNP_max"] = "" if not ri_values else max(ri_values)
                hit["nist_hit_record_count"] = meta["records"]
                hit["nist_source_files"] = ";".join(sorted(meta["files"]))
                if rank == 1:
                    best[(fid, library)] = hit
                top_rows.append({
                    "GlobalFeatureID": fid,
                    "FoundIn": feature.get("FoundIn", ""),
                    "corrected_tmean": _feature_rt_for_annotation(feature),
                    "BaseMz": feature.get("BaseMz", ""),
                    "library": library,
                    "rank": rank,
                    "score_0_1000": round(score * 1000, 1),
                    "spectral_cosine": round(score, 6),
                    "matched_ion_count": hit.get("matched_ion_count", ""),
                    "query_ion_count": hit.get("query_ion_count", ""),
                    "matched_ion_fraction": round(float(hit.get("matched_ion_fraction", 0)), 6),
                    "matched_intensity_fraction": round(float(hit.get("matched_intensity_fraction", 0)), 6),
                    "base_peak_matched": hit.get("base_peak_matched", ""),
                    "possible_coelution": hit.get("possible_coelution", ""),
                    "Name": hit.get("Name", ""),
                    "InChIKey": hit.get("InChIKey", ""),
                    "CAS": hit.get("CAS", ""),
                    "Formula": hit.get("Formula", ""),
                    "MW": hit.get("MW", ""),
                    "source_file": hit.get("source_file", ""),
                    "record_id": hit.get("library_record_id", ""),
                    "embedded_SemiStdNP_RI": hit.get("embedded_SemiStdNP_RI", ""),
                    "embedded_SemiStdNP_values": hit.get("embedded_SemiStdNP_values", ""),
                    "embedded_SemiStdNP_count": hit.get("embedded_SemiStdNP_count", ""),
                    "embedded_SemiStdNP_min": hit.get("embedded_SemiStdNP_min", ""),
                    "embedded_SemiStdNP_max": hit.get("embedded_SemiStdNP_max", ""),
                    "nist_hit_record_count": hit.get("nist_hit_record_count", ""),
                    "nist_source_files": hit.get("nist_source_files", ""),
                })

        try:
            display_path = ";".join(str(path.resolve().relative_to(PROJECT)) for path in paths)
        except ValueError:
            display_path = ";".join(str(path) for path in paths)
        summary = {
            "library": library,
            "file": display_path,
            "records_loaded_or_passed_filter": total,
            "records_searched": searched,
            "ei_gc_filter": require_ei_gc,
        }
        return best, top_rows, summary


    def clean_cas(rec: dict) -> str:
        return rec.get("CAS") or rec.get("CASNO") or ""


    def norm_text(value: str) -> str:
        return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


    def hit_columns(hit: dict | None, prefix: str) -> dict:
        if hit is None:
            return {
                f"{prefix}_score_0_1000": "",
                f"{prefix}_spectral_cosine": "",
                f"{prefix}_matched_ion_count": "",
                f"{prefix}_query_ion_count": "",
                f"{prefix}_matched_ion_fraction": "",
                f"{prefix}_matched_intensity_fraction": "",
                f"{prefix}_base_peak_matched": "",
                f"{prefix}_possible_coelution": "",
                f"{prefix}_Name": "",
                f"{prefix}_InChIKey": "",
                f"{prefix}_CAS": "",
                f"{prefix}_Formula": "",
                f"{prefix}_MW": "",
                f"{prefix}_source_file": "",
                f"{prefix}_record_id": "",
            }
        return {
            f"{prefix}_score_0_1000": round(float(hit["score"]) * 1000, 1),
            f"{prefix}_spectral_cosine": round(float(hit["score"]), 6),
            f"{prefix}_matched_ion_count": hit.get("matched_ion_count", ""),
            f"{prefix}_query_ion_count": hit.get("query_ion_count", ""),
            f"{prefix}_matched_ion_fraction": round(float(hit.get("matched_ion_fraction", 0)), 6),
            f"{prefix}_matched_intensity_fraction": round(float(hit.get("matched_intensity_fraction", 0)), 6),
            f"{prefix}_base_peak_matched": hit.get("base_peak_matched", ""),
            f"{prefix}_possible_coelution": hit.get("possible_coelution", ""),
            f"{prefix}_Name": hit.get("Name", ""),
            f"{prefix}_InChIKey": hit.get("InChIKey", ""),
            f"{prefix}_CAS": clean_cas(hit),
            f"{prefix}_Formula": hit.get("Formula", ""),
            f"{prefix}_MW": hit.get("MW", ""),
            f"{prefix}_source_file": hit.get("source_file", ""),
            f"{prefix}_record_id": hit.get("library_record_id", ""),
        }


    def best_overall(mona: dict | None, nist: dict | None) -> tuple[str, dict | None]:
        candidates = [("MoNA_GC_MS", mona), ("NIST2020_EI", nist)]
        candidates = [(src, hit) for src, hit in candidates if hit is not None]
        if not candidates:
            return "", None
        return max(candidates, key=lambda item: float(item[1]["score"]))


    def confidence(score: float | None, agreement: str) -> str:
        if score is None:
            return "no_hit"
        if score >= 0.90 and agreement in {"same_inchikey", "same_name"}:
            return "high_library_agreement"
        if score >= 0.90:
            return "high_single_library"
        if score >= 0.80:
            return "medium"
        if score >= 0.70:
            return "low"
        return "very_low"


    def load_residual_flags() -> dict[str, str]:
        flags: dict[str, set[str]] = defaultdict(set)
        if not RESIDUAL_FLAGS.exists():
            return {}
        with RESIDUAL_FLAGS.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                risk = row.get("risk_class") or row.get("relation") or "residual_pair"
                for col in ("GlobalFeatureID_1", "GlobalFeatureID_2"):
                    fid = row.get(col)
                    if fid:
                        flags[fid].add(risk)
        return {fid: ";".join(sorted(vals)) for fid, vals in flags.items()}


    def load_feature_qc_rows() -> list[dict]:
        if FEATURE_QC.exists():
            with FEATURE_QC.open(newline="", encoding="utf-8-sig") as handle:
                return list(csv.DictReader(handle))
        if not FEATURE_AREA_MATRIX.exists():
            raise FileNotFoundError(f"Missing FEATURE_QC and FEATURE_AREA_MATRIX: {FEATURE_QC}, {FEATURE_AREA_MATRIX}")
        rows: list[dict] = []
        with FEATURE_AREA_MATRIX.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            sample_cols = [c for c in (reader.fieldnames or []) if c != "GlobalFeatureID"]
            for row in reader:
                vals = []
                for col in sample_cols:
                    try:
                        vals.append(float(row.get(col, 0) or 0))
                    except ValueError:
                        vals.append(0.0)
                nonzero = [v for v in vals if v > 0]
                rows.append({
                    "GlobalFeatureID": row["GlobalFeatureID"],
                    "FoundIn": len(nonzero),
                    "total_area": sum(vals),
                    "median_nonzero_area": "" if not nonzero else sorted(nonzero)[len(nonzero) // 2],
                })
        return rows


    def main() -> None:
        with FEATURE_METADATA.open(newline="", encoding="utf-8-sig") as handle:
            metadata = {row["GlobalFeatureID"]: row for row in csv.DictReader(handle)}
        qc_rows = load_feature_qc_rows()

        features = []
        for row in qc_rows:
            found_in = int(float(row.get("FoundIn", 0) or 0))
            if found_in <= FOUNDIN_MIN_EXCLUSIVE:
                continue
            fid = row["GlobalFeatureID"]
            merged = dict(metadata[fid])
            merged.update(row)
            merged["query_spectrum"] = keep_top_ions(parse_feature_spectrum(merged.get("Spectra", "")))
            merged["query_n_peaks"] = len(merged["query_spectrum"])
            features.append(merged)

        residual_flags = load_residual_flags()
        best: dict[tuple[str, str], dict] = {}
        top_rows: list[dict] = []
        lib_summary: list[dict] = []

        library_jobs = library_jobs_for_mode(SPECTRAL_LIBRARY_MODE, MONA_FILE, NIST_FILES)
        print(f"Spectral library mode: {SPECTRAL_LIBRARY_MODE}", flush=True)
        for paths, lib_name, require_ei_gc in library_jobs:
            missing = [path for path in paths if not path.exists()]
            if missing:
                print(f"Skipping missing {lib_name}: {missing}", flush=True)
                lib_summary.append({
                    "library": lib_name,
                    "file": ";".join(map(str, missing)),
                    "records_loaded_or_passed_filter": 0,
                    "records_searched": 0,
                    "ei_gc_filter": require_ei_gc,
                    "status": "missing_file",
                })
                continue
            print(f"Streaming {lib_name}: {len(paths)} file(s)", flush=True)
            lib_best, lib_top_rows, summary = stream_search_library(
                features,
                paths,
                lib_name,
                require_ei_gc=require_ei_gc,
                top_n=TOP_HITS_PER_LIBRARY,
            )
            best.update(lib_best)
            top_rows.extend(lib_top_rows)
            summary["status"] = "searched"
            lib_summary.append(summary)
            gc.collect()

        include_n_blocks = any("n_blocks" in feature for feature in features)
        include_n_block_features = any("n_block_features" in feature for feature in features)
        out_rows: list[dict] = []
        for feature in features:
            fid = feature["GlobalFeatureID"]
            mona = best.get((fid, "MoNA_GC_MS"))
            nist = best.get((fid, "NIST2020_EI"))
            mona_key = (mona or {}).get("InChIKey", "")
            nist_key = (nist or {}).get("InChIKey", "")
            same_inchikey = bool(mona_key and nist_key and mona_key == nist_key)
            same_name = bool(
                norm_text((mona or {}).get("Name", "")) and
                norm_text((mona or {}).get("Name", "")) == norm_text((nist or {}).get("Name", ""))
            )
            same_formula = bool((mona or {}).get("Formula", "") and (mona or {}).get("Formula", "") == (nist or {}).get("Formula", ""))
            agreement = (
                "same_inchikey" if same_inchikey else
                "same_name" if same_name else
                "same_formula_only" if same_formula else
                "conflict_or_single_library"
            )
            best_source, best_hit = best_overall(mona, nist)
            best_score = float(best_hit["score"]) if best_hit is not None else None
            row = {
                "GlobalFeatureID": fid,
                "FoundIn": feature.get("FoundIn", ""),
                "corrected_tmean": _feature_rt_for_annotation(feature),
                "mean_original_tmean": feature.get("mean_original_tmean", ""),
                "BaseMz": feature.get("BaseMz", ""),
            }
            for optional_column in ("measured_RI", "reference_RI", "reference_RT"):
                if optional_column in feature:
                    row[optional_column] = feature.get(optional_column, "")
            if include_n_blocks:
                row["n_blocks"] = feature.get("n_blocks", "")
            if include_n_block_features:
                row["n_block_features"] = feature.get("n_block_features", "")
            row.update({
                "query_n_peaks": feature.get("query_n_peaks", ""),
                "best_source": best_source,
                "best_score_0_1000": "" if best_score is None else round(best_score * 1000, 1),
                "best_spectral_cosine": "" if best_score is None else round(best_score, 6),
                "best_Name": "" if best_hit is None else best_hit.get("Name", ""),
                "best_InChIKey": "" if best_hit is None else best_hit.get("InChIKey", ""),
                "best_CAS": "" if best_hit is None else clean_cas(best_hit),
                "best_Formula": "" if best_hit is None else best_hit.get("Formula", ""),
                "best_MW": "" if best_hit is None else best_hit.get("MW", ""),
                "agreement_flag": agreement,
                "confidence": confidence(best_score, agreement),
                "residual_duplicate_flag": residual_flags.get(fid, ""),
            })
            row.update(hit_columns(mona, "MoNA_GC_MS"))
            row.update(hit_columns(nist, "NIST2020_EI"))
            out_rows.append(row)

        best_file = OUT_DIR / "04_foundin_gt3_mona_nist_best_hits.csv"
        top_file = OUT_DIR / "04_foundin_gt3_mona_nist_top_hits.csv"
        summary_file = OUT_DIR / "04_foundin_gt3_mona_nist_annotation_summary.csv"
        lib_file = OUT_DIR / "04_foundin_gt3_mona_nist_library_search_summary.csv"

        with best_file.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(out_rows[0].keys()))
            writer.writeheader()
            writer.writerows(out_rows)
        with top_file.open("w", newline="", encoding="utf-8") as handle:
            top_fields = list(top_rows[0].keys()) if top_rows else [
                "GlobalFeatureID", "FoundIn", "corrected_tmean", "BaseMz", "library",
                "rank", "score_0_1000", "spectral_cosine", "Name", "InChIKey",
                "CAS", "Formula", "MW", "source_file", "record_id",
            ]
            writer = csv.DictWriter(handle, fieldnames=top_fields)
            writer.writeheader()
            writer.writerows(top_rows)

        summary = {
            "n_input_global_features": len(qc_rows),
            "foundin_filter": f"FoundIn>{FOUNDIN_MIN_EXCLUSIVE}",
            "n_features_searched": len(out_rows),
            "n_mona_best_hits": sum(bool(r["MoNA_GC_MS_Name"]) for r in out_rows),
            "n_nist_best_hits": sum(bool(r["NIST2020_EI_Name"]) for r in out_rows),
            "n_best_score_ge_900": sum(float(r["best_score_0_1000"] or 0) >= 900 for r in out_rows),
            "n_best_score_ge_800": sum(float(r["best_score_0_1000"] or 0) >= 800 for r in out_rows),
            "n_best_score_ge_700": sum(float(r["best_score_0_1000"] or 0) >= 700 for r in out_rows),
            "n_same_inchikey": sum(str(r["agreement_flag"]) == "same_inchikey" for r in out_rows),
            "n_same_name": sum(str(r["agreement_flag"]) == "same_name" for r in out_rows),
            "n_same_formula_only": sum(str(r["agreement_flag"]) == "same_formula_only" for r in out_rows),
            "n_with_residual_duplicate_flag": sum(bool(r["residual_duplicate_flag"]) for r in out_rows),
        }
        with summary_file.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summary.keys()))
            writer.writeheader()
            writer.writerow(summary)
        with lib_file.open("w", newline="", encoding="utf-8") as handle:
            lib_fields = []
            for row in lib_summary:
                for key in row:
                    if key not in lib_fields:
                        lib_fields.append(key)
            writer = csv.DictWriter(handle, fieldnames=lib_fields)
            writer.writeheader()
            writer.writerows(lib_summary)

        print("Wrote", best_file)
        print("Wrote", top_file)
        print("Wrote", summary_file)
        print(summary)
    return main
