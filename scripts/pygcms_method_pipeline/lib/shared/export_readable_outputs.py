#!/usr/bin/env python3
# =============================================================================
# Purpose: Export selected workflow outputs to the user-facing output directory.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# =============================================================================

from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path


def copy_if_exists(src: Path, dst: Path) -> bool:
    if not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def first_match(folder: Path, pattern: str) -> Path | None:
    hits = sorted(folder.glob(pattern)) if folder.exists() else []
    return hits[0] if hits else None


def first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def write_readme(out_dir: Path, copied: list[str]) -> None:
    readme = out_dir / "README.md"
    readme.write_text(
        "\n".join([
            "# Main Py-GC-MS Outputs",
            "",
            "This folder contains the key tables needed for review, downstream analysis, and figures.",
            "Intermediate and audit files are kept separately in `work_files` when needed.",
            "",
            "- `features_for_annotation.csv`: reproducible features used for annotation.",
            "- `feature_area_matrix.csv`: feature-by-sample peak-area matrix after reproducibility filtering.",
            "- `annotation_best_hits.csv`: best library hits before manual review.",
            "- `annotation_candidate_evidence.csv`: every retained Top-5 candidate with measured/reference RI evidence.",
            "- `annotation_automated_identification.csv`: conservative feature-level RI-guided decisions.",
            "- `annotation_autoaccepted.csv`: unambiguous automated compound-, homologous-series-, and family-level assignments.",
            "- `annotation_review_CHECK.csv`: compact manual-review table containing only unresolved features.",
            "- `annotation_evidence_full.csv`: complete annotation evidence retained for audit and troubleshooting.",
            "- `annotation_review.csv`: protected reviewed table name; if present, downstream finalization uses this first.",
            "- `final_compounds.csv`: reviewed final compound table.",
            "- `final_compound_area_matrix.csv`: final compound-by-sample peak-area matrix.",
            "- `final_compound_membership.csv`: mapping from final compounds back to reviewed feature-level rows.",
            "- `same_name_close_rt_merge_pairs.csv`: feature pairs merged after review.",
            "- `post_review_merge_summary.csv`: counts after manual review and post-review merging.",
            "- `checked_features_kept.csv` / `checked_features_excluded.csv`: reviewed feature-level rows kept or excluded.",
            "- `final_compounds_with_smiles.csv`: final compounds after SMILES retrieval.",
            "- `classified_compounds.csv`: final compounds with structural categories.",
            "- `classified_area_table.csv`: classified compound table with sample areas.",
            "- `classification_summary.csv`: compound counts by structural category.",
            "- `pipeline_summary.csv`: total peak/feature/compound retention across steps.",
            "- `sample_area_summary.csv`: simplified per-sample area-retention table for checking.",
            "- `sample_area_retention_wide.csv`: per-sample area retention in one row per sample.",
            "- `sample_step_area_summary.csv`: per-sample retained area after each processing step.",
            "- `sample_area_filter_summary.csv`: within-sample area-filter summary.",
            "",
            "## Files copied in this run",
            "",
            *[f"- `{path}`" for path in copied],
            "",
        ]),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Create easy-to-read copies of key Py-GC-MS pipeline outputs.")
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--foundin-min", type=int, default=2)
    parser.add_argument("--annotation-dir", type=Path)
    parser.add_argument("--review-dir", type=Path)
    parser.add_argument("--report-dir", type=Path)
    args = parser.parse_args()

    run_dir = args.run_dir
    out_dir = args.output_dir if args.output_dir else run_dir / "main_outputs"
    annotation_dir = args.annotation_dir if args.annotation_dir else run_dir / "05_annotation"
    review_dir = args.review_dir if args.review_dir else run_dir / "06_manual_review"
    report_dir = args.report_dir if args.report_dir else run_dir / "reports"
    copied: list[str] = []

    mappings: list[tuple[Path | None, str]] = [
        (run_dir / "03_area_filter" / "sample_area_filter_summary.csv", "sample_area_filter_summary.csv"),
        (first_existing([
            Path(os.environ.get("REFERENCE_RI_OUTPUT_DIR") or run_dir / "05_reference_ri_mapping") / "features_for_annotation.csv",
            run_dir / "04_global_grouping" / "features_for_annotation.csv",
            run_dir / "04_global_grouping" / f"global_feature_metadata_foundin_ge{args.foundin_min}.csv",
        ]), "features_for_annotation.csv"),
        (first_existing([
            run_dir / "04_global_grouping" / "feature_area_matrix.csv",
            run_dir / "04_global_grouping" / f"global_feature_area_matrix_foundin_ge{args.foundin_min}.csv",
        ]), "feature_area_matrix.csv"),
        (annotation_dir / "annotation_best_hits.csv", "annotation_best_hits.csv"),
        (annotation_dir / "annotation_candidate_evidence.csv", "annotation_candidate_evidence.csv"),
        (annotation_dir / "annotation_automated_identification.csv", "annotation_automated_identification.csv"),
        (review_dir / "annotation_autoaccepted.csv", "annotation_autoaccepted.csv"),
        (review_dir / "compound_identification_review_CHECK.csv", "annotation_review_CHECK.csv"),
        (review_dir / "annotation_evidence_full.csv", "annotation_evidence_full.csv"),
        (run_dir / "07_post_review_merge" / "final_compounds.csv", "final_compounds.csv"),
        (run_dir / "07_post_review_merge" / "final_compound_area_matrix.csv", "final_compound_area_matrix.csv"),
        (run_dir / "07_post_review_merge" / "final_compound_membership.csv", "final_compound_membership.csv"),
        (run_dir / "07_post_review_merge" / "same_name_close_rt_merge_pairs.csv", "same_name_close_rt_merge_pairs.csv"),
        (run_dir / "07_post_review_merge" / "post_review_merge_summary.csv", "post_review_merge_summary.csv"),
        (run_dir / "07_post_review_merge" / "checked_features_kept.csv", "checked_features_kept.csv"),
        (run_dir / "07_post_review_merge" / "checked_features_excluded.csv", "checked_features_excluded.csv"),
        (run_dir / "08_smiles" / "final_compounds_with_smiles.csv", "final_compounds_with_smiles.csv"),
        (run_dir / "09_classification" / "classified_compounds.csv", "classified_compounds.csv"),
        (run_dir / "09_classification" / "classified_area_table.csv", "classified_area_table.csv"),
        (run_dir / "09_classification" / "classification_summary.csv", "classification_summary.csv"),
        (report_dir / "pipeline_summary.csv", "pipeline_summary.csv"),
        (report_dir / "sample_area_summary.csv", "sample_area_summary.csv"),
        (report_dir / "sample_area_retention_wide.csv", "sample_area_retention_wide.csv"),
        (report_dir / "sample_step_area_summary.csv", "sample_step_area_summary.csv"),
    ]

    for src, rel_dst in mappings:
        if src is None:
            continue
        if copy_if_exists(src, out_dir / rel_dst):
            copied.append(rel_dst)

    write_readme(out_dir, copied)
    print(f"Wrote main outputs to {out_dir}")


if __name__ == "__main__":
    main()
