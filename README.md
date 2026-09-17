# Soil Py-GC-MS workflow

Scripts for processing soil organic matter Py-GC-MS data, including peak deconvolution and alignment, compound annotation, SMILES retrieval, and chemical classification.

Code accompanying the manuscript **“A reproducible workflow for high-throughput pyrolysis gas chromatography–mass spectrometry analysis of soil organic matter”** by Xiangxia Yang, Zhujun Wang, and Kazuo Isobe.

## Requirements

Run the scripts in a Bash terminal with Python 3.10 or newer and R available as `python3` and `Rscript`.

```bash
python3 -m pip install "openpyxl>=3.1,<4"
```

Install the R packages from an R session:

```r
install.packages(c("erah", "ncdf4", "dplyr", "readr", "ggplot2", "rcdk", "tidyverse"))
```

Chemical classification uses `rcdk`, which requires Java. Check the installation with `Rscript -e 'library(rcdk)'`. Annotation and SMILES retrieval use internet access for NIST WebBook and PubChem queries.

## Input files and paths

Open a terminal in the repository root. Set the dataset name and replace the three input paths below with your own locations:

```bash
export DATASET_ID="case2_alignment"
export CDF_DIR_OVERRIDE="/path/to/case2 spectrum"
export MONA_FILE="/path/to/MoNA-GC-MS.msp"
export NIST_FILES="/path/to/nist-main.MSP;/path/to/nist-replicates.MSP"
```

- `DATASET_ID` names the output folder. Use a new name when processing another dataset or comparing parameter settings.
- `CDF_DIR_OVERRIDE` points to the folder containing the sample CDF chromatograms.
- `MONA_FILE` points to the MoNA EI spectral library in MSP format.
- `NIST_FILES` contains one or more NIST EI library files. Separate multiple paths with semicolons inside the quotation marks.

For library preparation, see [Exporting the NIST Mass Spectral Library to MSP Format](NIST_to_MSP.md).

Case 2 raw chromatograms are packaged separately as `Case2_raw_data.zip`. Extract its `case2 spectrum/` folder into the repository root, or set `CDF_DIR_OVERRIDE` to the extracted folder. Prepare the EI spectral libraries separately using [MoNA](https://mona.fiehnlab.ucdavis.edu/) and your [NIST library installation](https://www.nist.gov/srd/nist-standard-reference-database-1a).

Keep these environment settings in the same terminal for all four stages. To use another Python installation, also set `PYTHON_BIN` to its executable path.

The repository includes the calibration and compound-reference files used by the workflow. For a different analytical setup, replace calibration inputs with matching measurements while retaining their column names:

| Setting in `config.sh` | Default input | Purpose |
|---|---|---|
| `ALKANE_RI_FILE` | `alkane_RI.csv` | Alkane carbon numbers, RI values, and retention times |
| `REFERENCE_RI_INTERNAL_STANDARD_RT_FILE` | `reference_internal_standard_rts.csv` | Internal-standard reference retention times |
| `SOM_REFERENCE_FILE` | `SOM_PyGCMS_EI_pipeline_reference_library_v1.0.xlsx` | Compound, diagnostic-ion, and literature references |

### Preparing the two calibration tables

Both CSV files are stored in the repository root. They were prepared from separate standard measurements before sample processing; the four-stage workflow reads these tables without regenerating them.

**`alkane_RI.csv` — n-alkane RI calibration**

1. Deconvolute the direct GC-MS C7–C40 standard run (`C7-C40-3uL-2.cdf`) and identify the n-alkanes using their elution order and EI spectra, including molecular ions where available.
2. Retain one confirmed RT per carbon number. The supplied table contains 34 selected component RTs, rounded to three decimal places in minutes.
3. Export `carbon_number`, `RI`, and `RT_min`, assigning `RI = 100 × carbon_number`; `source` records the identification notes. The workflow interpolates between these RT–RI points to calculate sample RI.

**`reference_internal_standard_rts.csv` — internal-standard RT references**

1. Deconvolute the separate Py-GC-MS reference run containing the internal standards (`C7-C40-IS.cdf`) with eRah.
2. Select tetracosane-d50 and chrysene-d12 using their expected RT regions and characteristic ions at m/z 66 and 240. Cross-check the selected component RTs against the corresponding extracted-ion chromatogram (EIC) apex RTs.
3. Save the selected eRah component RTs in `reference_rt_min`: **50.1545 min** for tetracosane-d50 and **51.9803 min** for chrysene-d12. These are the values used for RT correction; `eic_rt_min` and `erah_eic_delta_sec` record the cross-check only.

For new reference measurements, replace the calibration values while keeping the column names and internal-standard names unchanged. The `source` and `source_cdf` fields are provenance records, not file paths opened during sample processing.

Sample grouping searches for tetracosane-d50 at m/z 66 around 50.7 min and chrysene-d12 at m/z 240 around 52.7 min. These expected RTs locate the peaks; they are separate from the reference RTs used for correction. If the peak locations change, update the `find_eic_apex()` calls in `scripts/pygcms_method_pipeline/lib/01_preprocessing/01_internal_standard_guided_sample_grouping.R`.

## Processing parameters

The main default settings in [`config.sh`](scripts/pygcms_method_pipeline/config.sh) are:

| Parameter | Default | What it controls |
|---|---|---|
| `MIN_PEAK_HEIGHT` | `8000` | Minimum peak height for detection |
| `MIN_PEAK_WIDTH` | `2.5` | Minimum peak width in seconds |
| `NOISE_THRESHOLD` | `1000` | Noise threshold for deconvolution |
| `ANALYSIS_START_MIN`, `ANALYSIS_END_MIN` | `3`, `90` | Retention-time interval to process, in minutes |
| `AREA_FRACTION_FILTER` | `0.005` | Retain peaks contributing at least 0.5% of the sample's total peak area after peak cleanup |
| `FOUNDIN_MIN` | `3` | Minimum number of samples in which a feature must occur |
| `ALIGNMENT_MODE` | `auto` | Use one alignment group for up to 25 samples; otherwise group samples by internal-standard RT drift |
| `ALIGNMENT_MIN_BLOCK_SIZE`, `ALIGNMENT_MAX_BLOCK_SIZE` | `10`, `25` | Target sample-group sizes for blockwise alignment |
| `ALIGNMENT_ALIGN_TIME_DIST` | `60` | Within-group alignment RT window in seconds |
| `ALIGNMENT_MIN_SPECTRA_COR` | `0.90` | Minimum spectral correlation for alignment |
| `ALIGNMENT_MZ_MIN`, `ALIGNMENT_MZ_MAX` | `46`, `650` | Mass range used for processing |
| `ALIGNMENT_AVOID_PROCESSING_MZ` | `74:75,147:149,207,281` | Ions excluded during processing; a colon denotes an inclusive range |
| `ALIGNMENT_GLOBAL_MERGE_RT_SEC` | `60` | RT window for merging features across sample groups, in seconds |
| `ALIGNMENT_RT_CORRECTION_METHOD` | `adaptive_landmark` | Adaptive RT correction using matched spectral landmarks across sample groups |
| `TOP_HITS_PER_LIBRARY` | `20` | Number of spectral candidates retained from each library |
| `PRIMARY_TOP_HITS_PER_LIBRARY`, `RESCUE_TOP_HITS_PER_LIBRARY` | `5`, `20` | Evaluate the top 5 candidates first, then use top-20 rescue for unresolved features |
| `SPECTRAL_AUTO_THRESHOLD` | `850` | Spectral score required for automatic annotation eligibility, on a 0–1000 scale |
| `RI_SUPPORT_WINDOW`, `RI_WEAK_WINDOW` | `20`, `50` | Absolute RI-difference limits for supported and weakly supported candidates |
| `POST_REVIEW_MERGE_RT_SEC` | `60` | RT window for merging reviewed features with the same identity, in seconds |

For another dataset, adjust the detection thresholds, RT/mass ranges, and minimum sample occurrence to match the acquisition and study design. Edit the default values in `config.sh`, or override individual settings before running:

```bash
export MIN_PEAK_HEIGHT="8000"
export AREA_FRACTION_FILTER="0.005"
export FOUNDIN_MIN="3"
```

Exported values take precedence over the defaults in `config.sh`.

## Run the workflow

Run the following commands from the repository root. The files in `scripts/pygcms_method_pipeline/lib/` are called by these entry scripts and should remain in place.

The implementation follows the same processing order:

| Folder under `lib/` | Responsibility |
|---|---|
| `01_preprocessing/` | Sample grouping, deconvolution, area filtering, alignment, RT correction, and RI calculation |
| `02_annotation/` | `spectral_search.py` searches libraries; `annotation_evidence.py` evaluates candidates; `reviewed_compounds.py` merges reviewed features; `compound_annotation_and_review.py` coordinates Stage 2 and creates review tables |
| `03_smiles/` | SMILES retrieval for reviewed identities |
| `04_classification/` | Structural classification and final area tables |
| `shared/` | Dataset checks, processing summaries, and output export used across stages |

### 1. Deconvolution, alignment, and RI calculation

```bash
bash scripts/pygcms_method_pipeline/01_pre_annotation_processing_alignment.sh
```

This stage groups samples, deconvolves and aligns peaks, and applies the area and sample-occurrence filters. Globally aligned retention times are corrected using the internal-standard reference RTs, then used directly for RI interpolation against `alkane_RI.csv`.

Inside `lib/01_preprocessing/`, steps `01`–`05` follow the processing order. `00_run_alignment_pipeline.R` calls steps `01`–`03`; the Stage 1 shell script then calls `04` and `05`. Module `06` supplies functions used within step `02`. These files are called automatically by the Stage 1 shell script.

| File prefix or folder | Role |
|---|---|
| `00_` | Coordinate sample grouping and alignment (`01`–`03`) |
| `01_` | Group samples using internal-standard RT drift |
| `02_` | Deconvolute peaks, clean and filter them by area, and align within sample groups |
| `03_` | Correct RT differences between groups and align features globally |
| `04_` | Apply the minimum sample-occurrence filter and export annotation inputs |
| `05_` | Correct RT against internal-standard references and calculate RI from the alkane table |
| `06_` | Match and merge unaligned peaks (`AlignID=0`); called within step `02` |

### 2. Compound annotation and manual review

```bash
bash scripts/pygcms_method_pipeline/02_compound_annotation_and_review.sh
```

This stage searches both spectral libraries, combines spectral and RI evidence, and prepares the annotation review table.

Open `pygcms_method_outputs/case2_alignment/main_outputs/annotation_review_CHECK.csv` (replace `case2_alignment` with your dataset ID). Review the candidate names, spectral scores, RI evidence, and diagnostic ions. Fill or correct `Final_Identification`; use `Unidentified` where no identity can be assigned, and record comments in `Review_Notes`. Keep all feature rows, IDs, and column names, including the automatically filled assignments.

When changing `Final_Identification`, also update `Formula` to the reviewed identity, or clear it if uncertain. The script does not infer a new formula from the edited name. Spectral scores, reference RI, and diagnostic evidence remain records of the original candidate assessment.

Save the completed table as **`annotation_review.csv`** in the same folder, then run Stage 2 again to finalize the identities and merge eligible features:

```bash
bash scripts/pygcms_method_pipeline/02_compound_annotation_and_review.sh
```

### 3. Retrieve SMILES

```bash
bash scripts/pygcms_method_pipeline/03_retrieve_smiles.sh
```

This stage retrieves SMILES for the reviewed compound names from the local libraries and PubChem.

Local lookup prioritizes InChIKey matches over name matches. Within each lookup type, it checks MoNA before the configured NIST files, then queries PubChem if no local match is found.

### 4. Chemical classification

```bash
bash scripts/pygcms_method_pipeline/04_chemical_classification.sh
```

This stage assigns structural categories from SMILES and exports the classified peak-area table.

## Outputs

Results are saved under `pygcms_method_outputs/<DATASET_ID>/main_outputs/`. Intermediate files are stored in the adjacent `work_files/` folder.

In the annotation and review tables, `corrected_RT_min` or `RT_min` reports the internal-standard-corrected RT used to calculate the measured RI.

| Output | Contents |
|---|---|
| `features_for_annotation.csv`, `feature_area_matrix.csv` | Aligned features and sample peak areas retained for annotation |
| `annotation_review_CHECK.csv` | Annotation table for manual review |
| `final_compounds.csv`, `final_compound_area_matrix.csv` | Final identities and compound peak areas after review and merging |
| `final_compounds_with_smiles.csv` | Final compound identities with retrieved SMILES |
| `classified_area_table.csv` | Compound identities, structural categories, and sample peak areas for downstream analysis |
