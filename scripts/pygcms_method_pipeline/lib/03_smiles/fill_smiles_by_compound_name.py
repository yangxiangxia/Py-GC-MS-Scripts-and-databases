#!/usr/bin/env python3
# =============================================================================
# Purpose: Resolve reviewed identities to SMILES, prioritizing InChIKey matches.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# =============================================================================

"""Fill SMILES in a compound table by matching compound names to local sources."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import ssl
import time
from pathlib import Path
from typing import Dict, Iterable, List, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


PROJECT = Path(os.environ.get("PROJECT_DIR", Path(__file__).resolve().parents[4]))

NAME_COLUMN_CANDIDATES = [
    "Compound",
    "compound",
    "Name",
    "name",
    "Best_compound",
    "Best_name",
    "compound_name",
    "CompoundName",
]

INCHIKEY_COLUMN_CANDIDATES = [
    "final_InChIKey", "final_inchikey", "InChIKey", "inchikey",
    "best_InChIKey", "best_inchikey",
]

CSV_NAME_SMILES_PAIRS = [
    ("Compound", "SMILES"),
    ("compound", "SMILES"),
    ("Name", "SMILES"),
    ("name", "SMILES"),
    ("Best_compound", "Best_SMILES"),
    ("Best_name", "Best_SMILES"),
    ("MoNA_GC_MS_compound", "MoNA_GC_MS_SMILES"),
    ("MoNA_GC_MS_Name", "MoNA_GC_MS_SMILES"),
    ("MassBank_EI_GCMS_compound", "MassBank_EI_GCMS_SMILES"),
    ("MassBank_NIST_Name", "MassBank_NIST_SMILES"),
    ("NIST2023_EI_compound", "NIST2023_EI_SMILES"),
]


def normalize_name(value: str) -> str:
    value = str(value or "").strip().lower()
    value = value.replace("\u2010", "-").replace("\u2011", "-").replace("\u2013", "-").replace("\u2014", "-")
    value = re.sub(r"\s+", " ", value)
    value = value.strip(" \t\r\n\"'")
    return value


def compact_name_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", normalize_name(value))


def clean_smiles(value: str) -> str:
    value = str(value or "").strip().strip('"').strip("'")
    if not value or value.upper() in {"N/A", "NA", "NULL", "NONE"}:
        return ""
    return value


def extract_smiles_from_text(text: str) -> str:
    patterns = (
        r'(?:^|["\s])SMILES=([^"]+)',
        r'(?:^|["\s])computed SMILES=([^"]+)',
        r'(?:^|["\s])smiles=([^"]+)',
        r'(?:^|["\s])computed smiles=([^"]+)',
    )
    for pattern in patterns:
        match = re.search(pattern, text or "")
        if match:
            smiles = clean_smiles(match.group(1))
            if smiles:
                return smiles
    return ""


def add_mapping(mapping: Dict[str, Tuple[str, str]], name: str, smiles: str, source: str) -> None:
    smiles = clean_smiles(smiles)
    if not name or not smiles:
        return
    for key in (normalize_name(name), compact_name_key(name)):
        if key and key not in mapping:
            mapping[key] = (smiles, source)


def normalize_inchikey(value: str) -> str:
    return str(value or "").strip().upper()


def is_unidentified_name(value: str) -> bool:
    return normalize_name(value) in {
        "", "unknown", "unknow", "unidentified", "not identified", "reject", "rejected"
    }


def add_inchikey_mapping(mapping: Dict[str, Tuple[str, str]], inchikey: str, smiles: str, source: str) -> None:
    key = normalize_inchikey(inchikey)
    smiles = clean_smiles(smiles)
    if key and smiles and key not in mapping:
        mapping[key] = (smiles, source)


def read_csv_mapping(path: Path, mapping: Dict[str, Tuple[str, str]]) -> int:
    if not path.exists():
        return 0
    added = 0
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return 0
        for row in reader:
            before = len(mapping)
            for name_col, smiles_col in CSV_NAME_SMILES_PAIRS:
                if name_col in row and smiles_col in row:
                    add_mapping(mapping, row.get(name_col, ""), row.get(smiles_col, ""), str(path))
            added += len(mapping) - before
    return added


def split_msp_records(path: Path) -> Iterable[List[str]]:
    record: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n\r")
            if not line.strip():
                if record:
                    yield record
                    record = []
                continue
            record.append(line)
    if record:
        yield record


def parse_msp_fields(record: List[str]) -> Dict[str, str]:
    fields: Dict[str, str] = {}
    comments: List[str] = []
    synons: List[str] = []
    for line in record:
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()
        if key.lower() == "num peaks":
            break
        key_lower = key.lower()
        if key_lower == "smiles":
            fields["SMILES"] = value
        elif key_lower == "comments":
            comments.append(value)
        elif key_lower == "comment":
            comments.append(value)
        elif key_lower == "synon":
            synons.append(value)
        elif key not in fields:
            fields[key] = value
    if comments:
        fields["Comments"] = " ".join(comments)
    if synons:
        fields["Synon"] = " | ".join(synons)
    return fields


def read_msp_mapping(path: Path, mapping: Dict[str, Tuple[str, str]], inchikey_mapping: Dict[str, Tuple[str, str]] | None = None) -> int:
    if not path.exists():
        return 0
    added = 0
    for record in split_msp_records(path):
        fields = parse_msp_fields(record)
        smiles = clean_smiles(fields.get("SMILES", "")) or extract_smiles_from_text(fields.get("Comments", ""))
        if not smiles:
            continue
        before = len(mapping)
        before_keys = len(inchikey_mapping or {})
        add_mapping(mapping, fields.get("Name", ""), smiles, str(path))
        for synon in str(fields.get("Synon", "")).split("|"):
            add_mapping(mapping, synon.strip(), smiles, str(path))
        if inchikey_mapping is not None:
            add_inchikey_mapping(inchikey_mapping, fields.get("InChIKey", ""), smiles, str(path))
        added += (len(mapping) - before) + (len(inchikey_mapping or {}) - before_keys)
    return added


def build_mapping(annotation_csvs: List[Path], msp_files: List[Path]) -> Tuple[Dict[str, Tuple[str, str]], Dict[str, Tuple[str, str]], List[dict]]:
    mapping: Dict[str, Tuple[str, str]] = {}
    inchikey_mapping: Dict[str, Tuple[str, str]] = {}
    summary: List[dict] = []
    for path in annotation_csvs:
        added = read_csv_mapping(path, mapping)
        summary.append({"source_type": "annotation_csv", "source": str(path), "new_name_keys_added": added})
    for path in msp_files:
        added = read_msp_mapping(path, mapping, inchikey_mapping)
        summary.append({"source_type": "msp", "source": str(path), "new_name_keys_added": added})
    return mapping, inchikey_mapping, summary


def detect_name_column(fieldnames: List[str], requested: str | None) -> str:
    if requested:
        if requested not in fieldnames:
            raise ValueError(f"Name column '{requested}' not found. Available columns: {', '.join(fieldnames)}")
        return requested
    for col in NAME_COLUMN_CANDIDATES:
        if col in fieldnames:
            return col
    raise ValueError(
        "Could not auto-detect compound-name column. Use --name-column. "
        f"Available columns: {', '.join(fieldnames)}"
    )


def lookup_smiles(mapping: Dict[str, Tuple[str, str]], name: str) -> Tuple[str, str, str]:
    for key, match_type in ((normalize_name(name), "normalized_name"), (compact_name_key(name), "compact_name")):
        if key in mapping:
            smiles, source = mapping[key]
            return smiles, source, match_type
    return "", "", ""


def lookup_smiles_by_inchikey(mapping: Dict[str, Tuple[str, str]], inchikey: str) -> Tuple[str, str, str]:
    key = normalize_inchikey(inchikey)
    if key in mapping:
        smiles, source = mapping[key]
        return smiles, source, "inchikey"
    return "", "", ""


def lookup_local_smiles(name_mapping, inchikey_mapping, name: str, inchikey: str) -> Tuple[str, str, str]:
    if normalize_inchikey(inchikey):
        hit = lookup_smiles_by_inchikey(inchikey_mapping, inchikey)
        if hit[0]:
            return hit
    return lookup_smiles(name_mapping, name)


def pubchem_name_variants(name: str) -> List[str]:
    name = str(name or "").strip()
    normalized = normalize_name(name)
    variants = [name]
    alias_map = {
        "meta xylene": "m-Xylene",
        "meta-xylene": "m-Xylene",
        "meta xylol": "m-Xylene",
        "ortho xylene": "o-Xylene",
        "ortho-xylene": "o-Xylene",
        "para xylene": "p-Xylene",
        "para-xylene": "p-Xylene",
    }
    if normalized in alias_map:
        variants.append(alias_map[normalized])
    out: List[str] = []
    seen = set()
    for variant in variants:
        key = normalize_name(variant)
        if key and key not in seen:
            seen.add(key)
            out.append(variant)
    return out


def read_pubchem_cache(path: Path) -> Dict[str, dict]:
    cache: Dict[str, dict] = {}
    if not path.exists():
        return cache
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            query = normalize_name(row.get("query", ""))
            if query:
                if row.get("status") == "ok" and not (
                    clean_smiles(row.get("canonical_smiles", "")) or
                    clean_smiles(row.get("isomeric_smiles", ""))
                ):
                    row["status"] = "needs_refresh"
                cache[query] = row
    return cache


def write_pubchem_cache(path: Path, cache: Dict[str, dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "query", "status", "cid", "canonical_smiles", "isomeric_smiles",
        "inchikey", "error"
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for query in sorted(cache):
            row = {key: cache[query].get(key, "") for key in fieldnames}
            writer.writerow(row)


def _pubchem_lookup(identifier: str, namespace: str, cache: Dict[str, dict], delay_sec: float, timeout_sec: float) -> dict:
    query = normalize_name(identifier)
    cache_key = query if namespace == "name" else f"{namespace}:{query}"
    if not query:
        return {"query": query, "status": "empty_query", "error": "empty query"}
    if cache_key in cache:
        cached = cache[cache_key]
        # A previous timeout or server error must not prevent a later retry.
        if cached.get("status") in {"ok", "not_found"}:
            return cached

    encoded = quote(identifier.strip(), safe="")
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{namespace}/"
        f"{encoded}/property/CanonicalSMILES,IsomericSMILES,InChIKey/JSON"
    )
    row = {
        "query": cache_key,
        "status": "",
        "cid": "",
        "canonical_smiles": "",
        "isomeric_smiles": "",
        "inchikey": "",
        "error": "",
    }
    try:
        req = Request(url, headers={"User-Agent": "gcms-smiles-fill/1.0"})
        try:
            with urlopen(req, timeout=timeout_sec) as response:
                data = json.loads(response.read().decode("utf-8"))
        except URLError as exc:
            if "CERTIFICATE_VERIFY_FAILED" not in str(exc):
                raise
            context = ssl._create_unverified_context()
            with urlopen(req, timeout=timeout_sec, context=context) as response:
                data = json.loads(response.read().decode("utf-8"))
        props = data.get("PropertyTable", {}).get("Properties", [])
        if not props:
            row["status"] = "not_found"
        else:
            first = props[0]
            connectivity_smiles = clean_smiles(first.get("ConnectivitySMILES", ""))
            pubchem_smiles = clean_smiles(first.get("SMILES", ""))
            row.update({
                "status": "ok",
                "cid": str(first.get("CID", "")),
                "canonical_smiles": clean_smiles(first.get("CanonicalSMILES", "")) or connectivity_smiles or pubchem_smiles,
                "isomeric_smiles": clean_smiles(first.get("IsomericSMILES", "")) or pubchem_smiles or connectivity_smiles,
                "inchikey": first.get("InChIKey", ""),
            })
    except HTTPError as exc:
        row["status"] = "not_found" if exc.code == 404 else "http_error"
        row["error"] = f"HTTP {exc.code}: {exc.reason}"
    except (URLError, TimeoutError, json.JSONDecodeError) as exc:
        row["status"] = "request_error"
        row["error"] = str(exc)

    cache[cache_key] = row
    if delay_sec > 0:
        time.sleep(delay_sec)
    return row


def pubchem_lookup(name: str, cache: Dict[str, dict], delay_sec: float, timeout_sec: float) -> dict:
    return _pubchem_lookup(name, "name", cache, delay_sec, timeout_sec)


def pubchem_inchikey_lookup(inchikey: str, cache: Dict[str, dict], delay_sec: float, timeout_sec: float) -> dict:
    return _pubchem_lookup(normalize_inchikey(inchikey), "inchikey", cache, delay_sec, timeout_sec)


def smiles_from_pubchem_row(row: dict, smiles_type: str) -> str:
    if row.get("status") != "ok":
        return ""
    if smiles_type == "canonical":
        return clean_smiles(row.get("canonical_smiles", ""))
    return clean_smiles(row.get("isomeric_smiles", "")) or clean_smiles(row.get("canonical_smiles", ""))


def default_output_path(input_path: Path) -> Path:
    return input_path.with_name(input_path.stem + "_with_smiles" + input_path.suffix)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="CSV table containing final compound names.")
    parser.add_argument("--output", default="", help="Output CSV path. Defaults to *_with_smiles.csv.")
    parser.add_argument("--name-column", default="", help="Compound-name column. Auto-detected if omitted.")
    parser.add_argument("--inchikey-column", default="", help="InChIKey column. Auto-detected when available.")
    parser.add_argument("--smiles-column", default="SMILES", help="Column to fill or create.")
    parser.add_argument("--overwrite-existing", action="store_true", help="Replace existing non-empty SMILES values.")
    parser.add_argument("--annotation-csv", action="append", default=[], help="Annotation CSV source; may be repeated.")
    parser.add_argument("--msp", action="append", default=[], help="MSP source; may be repeated.")
    parser.add_argument("--use-pubchem", action="store_true", help="Query PubChem by compound name when local sources do not match.")
    parser.add_argument(
        "--pubchem-cache",
        default=str(PROJECT / "pubchem_smiles_cache.csv"),
        help="CSV cache for PubChem name lookups.",
    )
    parser.add_argument(
        "--pubchem-smiles-type",
        choices=("isomeric", "canonical"),
        default="isomeric",
        help="Which PubChem SMILES to fill when both are available.",
    )
    parser.add_argument("--pubchem-delay-sec", type=float, default=0.12, help="Delay after uncached PubChem requests.")
    parser.add_argument("--pubchem-timeout-sec", type=float, default=5, help="Timeout for each PubChem request.")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output) if args.output else default_output_path(input_path)

    annotation_csvs = [Path(p) for p in args.annotation_csv]
    msp_files = [Path(p) for p in args.msp]

    mapping, inchikey_mapping, source_summary = build_mapping(annotation_csvs, msp_files)
    pubchem_cache_path = Path(args.pubchem_cache)
    pubchem_cache = read_pubchem_cache(pubchem_cache_path) if args.use_pubchem else {}

    with input_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"No header found in {input_path}")
        fieldnames = list(reader.fieldnames)
        rows = list(reader)

    name_col = detect_name_column(fieldnames, args.name_column or None)
    if args.inchikey_column:
        if args.inchikey_column not in fieldnames:
            raise ValueError(f"InChIKey column '{args.inchikey_column}' not found")
        inchikey_col = args.inchikey_column
    else:
        inchikey_col = next((col for col in INCHIKEY_COLUMN_CANDIDATES if col in fieldnames), "")
    if args.smiles_column not in fieldnames:
        fieldnames.append(args.smiles_column)
    source_col = f"{args.smiles_column}_source"
    match_col = f"{args.smiles_column}_match_type"
    for col in (source_col, match_col):
        if col not in fieldnames:
            fieldnames.append(col)

    n_filled = 0
    n_filled_pubchem = 0
    n_existing_kept = 0
    n_unmatched = 0
    n_skipped_unidentified = 0
    for row in rows:
        if is_unidentified_name(row.get(name_col, "")):
            row[args.smiles_column] = ""
            row[source_col] = ""
            row[match_col] = "skipped_unidentified"
            n_skipped_unidentified += 1
            continue
        existing = clean_smiles(row.get(args.smiles_column, ""))
        if existing and not args.overwrite_existing:
            n_existing_kept += 1
            row[args.smiles_column] = existing
            row.setdefault(source_col, "existing")
            row.setdefault(match_col, "existing")
            continue
        inchikey = row.get(inchikey_col, "") if inchikey_col else ""
        smiles, source, match_type = lookup_local_smiles(mapping, inchikey_mapping, row.get(name_col, ""), inchikey)
        if not smiles and args.use_pubchem:
            if normalize_inchikey(inchikey):
                pc_row = pubchem_inchikey_lookup(
                    inchikey, pubchem_cache, args.pubchem_delay_sec, args.pubchem_timeout_sec,
                )
                write_pubchem_cache(pubchem_cache_path, pubchem_cache)
                smiles = smiles_from_pubchem_row(pc_row, args.pubchem_smiles_type)
                if smiles:
                    source = f"PubChem CID {pc_row.get('cid', '')}"
                    match_type = f"pubchem_inchikey_{args.pubchem_smiles_type}"
                    n_filled_pubchem += 1
            for query_name in pubchem_name_variants(row.get(name_col, "")):
                if smiles:
                    break
                pc_row = pubchem_lookup(
                    query_name,
                    pubchem_cache,
                    args.pubchem_delay_sec,
                    args.pubchem_timeout_sec,
                )
                write_pubchem_cache(pubchem_cache_path, pubchem_cache)
                smiles = smiles_from_pubchem_row(pc_row, args.pubchem_smiles_type)
                if smiles:
                    source = f"PubChem CID {pc_row.get('cid', '')}"
                    match_type = f"pubchem_name_{args.pubchem_smiles_type}"
                    if normalize_name(query_name) != normalize_name(row.get(name_col, "")):
                        match_type += "_alias"
                    n_filled_pubchem += 1
                    break
        if smiles:
            row[args.smiles_column] = smiles
            row[source_col] = source
            row[match_col] = match_type
            n_filled += 1
        else:
            row[args.smiles_column] = existing
            row[source_col] = row.get(source_col, "")
            row[match_col] = row.get(match_col, "")
            n_unmatched += 1

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    summary_path = output_path.with_name(output_path.stem + "_summary.csv")
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames_summary = ["metric", "value"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames_summary)
        writer.writeheader()
        writer.writerows([
            {"metric": "input", "value": str(input_path)},
            {"metric": "output", "value": str(output_path)},
            {"metric": "name_column", "value": name_col},
            {"metric": "smiles_column", "value": args.smiles_column},
            {"metric": "inchikey_column", "value": inchikey_col},
            {"metric": "n_rows", "value": len(rows)},
            {"metric": "n_filled", "value": n_filled},
            {"metric": "n_filled_pubchem", "value": n_filled_pubchem},
            {"metric": "n_existing_kept", "value": n_existing_kept},
            {"metric": "n_unmatched", "value": n_unmatched},
            {"metric": "n_skipped_unidentified", "value": n_skipped_unidentified},
            {"metric": "n_mapping_keys", "value": len(mapping)},
            {"metric": "n_inchikey_mapping_keys", "value": len(inchikey_mapping)},
            {"metric": "use_pubchem", "value": args.use_pubchem},
            {"metric": "pubchem_cache", "value": str(pubchem_cache_path)},
        ])

    source_summary_path = output_path.with_name(output_path.stem + "_source_summary.csv")
    with source_summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source_type", "source", "new_name_keys_added"])
        writer.writeheader()
        writer.writerows(source_summary)

    if args.use_pubchem:
        write_pubchem_cache(pubchem_cache_path, pubchem_cache)

    print(f"Wrote: {output_path}")
    print(f"Wrote: {summary_path}")
    print(f"Wrote: {source_summary_path}")
    print({
        "n_rows": len(rows),
        "n_filled": n_filled,
        "n_filled_pubchem": n_filled_pubchem,
        "n_existing_kept": n_existing_kept,
        "n_unmatched": n_unmatched,
        "n_skipped_unidentified": n_skipped_unidentified,
        "n_mapping_keys": len(mapping),
    })


if __name__ == "__main__":
    main()
