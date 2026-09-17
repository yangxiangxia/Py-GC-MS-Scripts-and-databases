#!/usr/bin/env Rscript
# =============================================================================
# Purpose: Classify reviewed compounds using explicit structural and name rules.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# =============================================================================


# Inputs prepared by the preceding workflow stages.
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = getwd())
DATASET_ID <- Sys.getenv("DATASET_ID", unset = "case2_alignment")
RUN_ROOT <- Sys.getenv("RUN_ROOT", unset = file.path(PROJECT_DIR, "pygcms_method_outputs"))
WORK_DIR <- Sys.getenv("WORK_DIR", unset = file.path(RUN_ROOT, DATASET_ID, "work_files"))

input_file <- file.path(WORK_DIR, "08_smiles", "final_compounds_with_smiles.csv")
area_file <- file.path(WORK_DIR, "07_post_review_merge", "final_compound_area_matrix.csv")
out_dir <- file.path(WORK_DIR, "09_classification")
intermediate_dir <- file.path(out_dir, "_intermediate")
classification_rules_label <- "integrated rules in lib/04_classification/classification.R"
pipeline_out_dir <- out_dir

if (!file.exists(input_file)) stop("Missing SMILES table: ", input_file)
if (!file.exists(area_file)) stop("Missing final area matrix: ", area_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

prepared_input <- file.path(intermediate_dir, "classification_input_prepared.csv")
prepared_area <- file.path(intermediate_dir, "classification_area_matrix_prepared.csv")

compounds <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!"CompoundID" %in% names(compounds)) stop("Input table must contain CompoundID")
if (!"GlobalFeatureID" %in% names(compounds)) compounds$GlobalFeatureID <- compounds$CompoundID
write.csv(compounds, prepared_input, row.names = FALSE, na = "")

area <- read.csv(area_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!"CompoundID" %in% names(area)) stop("Area matrix must contain CompoundID")
if (!"GlobalFeatureID" %in% names(area)) area$GlobalFeatureID <- area$CompoundID
area <- area[, c("GlobalFeatureID", "CompoundID", setdiff(names(area), c("GlobalFeatureID", "CompoundID"))), drop = FALSE]
write.csv(area, prepared_area, row.names = FALSE, na = "")

suppressPackageStartupMessages({
  library(rcdk)
  library(tidyverse)
})

# Structural classification and detailed audit outputs.
input_file <- prepared_input
area_file <- prepared_area
out_dir <- intermediate_dir
classified_file <- file.path(out_dir, "final_foundin_gt3_classified_compounds.csv")
complete_file <- file.path(out_dir, "final_foundin_gt3_classified_compound_area_table.csv")
complete_xlsx_file <- file.path(out_dir, "final_foundin_gt3_classified_compound_area_table.xlsx")
summary_file <- file.path(out_dir, "final_foundin_gt3_classification_summary.csv")

# Classification labels and structural thresholds.
CARB_LABEL <- Sys.getenv("CARB_LABEL", unset = "Carbohydrates")
# Minimum number of consecutive acyclic sp3 carbons that marks a compound as a
# long-chain (lipid/wax-like) aliphatic rather than a sugar pyrolysis fragment.
LONG_CHAIN_MIN_C <- as.integer(Sys.getenv("LONG_CHAIN_MIN_C", unset = "8"))

df <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!"GlobalFeatureID" %in% names(df) && "ID" %in% names(df)) {
  df$GlobalFeatureID <- df$ID
}
n <- nrow(df)
cols <- names(df)

pick_first <- function(cands) {
  hit <- cands[cands %in% cols]
  if (!length(hit)) NA_character_ else hit[1]
}

smiles_col <- pick_first(c("final_SMILES", "SMILES", "Smiles", "smiles"))
name_col <- pick_first(c("final_identification", "Name", "Name.1", "Compound", "Compound_name", "CompoundName", "PeakName"))
mf_col <- pick_first(c(
  "final_score_0_1000", "best_score_0_1000",
  "best_score",
  "MatchFactor", "Match.factor", "match_factor", "match.factor",
  "MF", "Similarity", "similarity", "SI", "sim", "Match"
))

if (is.na(smiles_col)) stop("No SMILES column found.")
if (is.na(name_col)) stop("No compound name column found.")

smiles <- trimws(df[[smiles_col]])
smiles[smiles %in% c("", "NA", "Na", "na", "N/A", "n/a", "Unidentified", "unidentified")] <- NA_character_
match_factor <- if (!is.na(mf_col)) suppressWarnings(as.numeric(df[[mf_col]])) else rep(NA_real_, n)
# Classification is intentionally structure/name based. Match scores are kept
# as audit metadata only and do not gate structural category assignment.
low_conf_id <- rep(FALSE, n)

isT <- function(x) !is.na(x) & x
safe_grepl <- function(pat, x, ...) {
  out <- rep(FALSE, length(x))
  idx <- which(!is.na(x))
  if (length(idx)) out[idx] <- grepl(pat, x[idx], ...)
  out
}

mol_all <- vector("list", n)
parse_ok <- rep(FALSE, n)

prepare_mol <- function(m) {
  tryCatch(rcdk::set.atom.types(m), error = function(e) NULL)
  tryCatch(rcdk::do.aromaticity(m), error = function(e) NULL)
  m
}

idx_smiles <- which(!is.na(smiles))
for (i in idx_smiles) {
  m <- tryCatch(rcdk::parse.smiles(smiles[i])[[1]], error = function(e) NULL)
  if (!is.null(m)) {
    mol_all[[i]] <- prepare_mol(m)
    parse_ok[i] <- TRUE
  }
}

mol_ok <- parse_ok & !vapply(mol_all, is.null, logical(1))
mols_valid <- mol_all[mol_ok]

smart_match_vec <- function(smarts) {
  out <- rep(NA, n)
  if (!length(mols_valid)) return(out)
  out[mol_ok] <- tryCatch(
    rcdk::matches(smarts, mols_valid),
    error = function(e) rep(FALSE, sum(mol_ok))
  )
  out
}

smart_match_any <- function(sv) {
  out <- rep(FALSE, n)
  for (s in sv) out <- out | isT(smart_match_vec(s))
  out
}

atom_symbols <- vector("list", n)
elem_set <- vector("list", n)
carbon_count <- rep(NA_integer_, n)
oxygen_count <- rep(NA_integer_, n)
nitrogen_count <- rep(NA_integer_, n)

for (i in seq_len(n)) {
  if (mol_ok[i]) {
    syms <- vapply(rcdk::get.atoms(mol_all[[i]]), rcdk::get.symbol, character(1))
    atom_symbols[[i]] <- syms
    elem_set[[i]] <- unique(syms)
    carbon_count[i] <- sum(syms == "C")
    oxygen_count[i] <- sum(syms == "O")
    nitrogen_count[i] <- sum(syms == "N")
  } else {
    atom_symbols[[i]] <- character(0)
    elem_set[[i]] <- character(0)
  }
}

has_N <- map_lgl(elem_set, ~ "N" %in% .x)
has_O <- map_lgl(elem_set, ~ "O" %in% .x)
# S, P, halogens, metals, and any other elements outside the target C/H/O/N
# chemistry are explicitly retained as structurally valid but miscellaneous.
has_non_target_element <- map_lgl(
  elem_set,
  ~ length(.x) > 0 && any(!.x %in% c("C", "H", "O", "N"))
)
only_CH <- map_lgl(elem_set, ~ length(.x) > 0 && all(.x %in% c("C", "H")))
only_CHO <- map_lgl(elem_set, ~ length(.x) > 0 && all(.x %in% c("C", "H", "O")))

has_double_bond <- safe_grepl("=", smiles, fixed = TRUE)
has_triple_bond <- safe_grepl("#", smiles, fixed = TRUE)

is_aromatic <- smart_match_vec("a")
has_benzene <- smart_match_vec("c1ccccc1")

# --- Robust benzene / aromaticity detection -------------------------------
# Some library SMILES arrive in Kekule form (upper-case C with explicit "=")
# or in non-standard notation, and CDK aromaticity perception can then fail,
# leaving a genuine aromatic ring flagged as non-aromatic. That is what pushed
# xylene / styrene / indene / fluorene / methyl-indene into the "alkene" bucket.
# has_benzene_kekule matches a 6-membered all-carbon ring with alternating
# double bonds (i.e. a benzene ring that was NOT aromatized). Combined with the
# aromatic-symbol match above it covers a benzene ring in either representation.
has_benzene_kekule <- smart_match_any(c("C1=CC=CC=C1"))
has_benzene_robust <- isT(has_benzene) | isT(has_benzene_kekule)
is_aromatic_robust <- isT(is_aromatic) | has_benzene_robust

has_phenol_on_benzene <- smart_match_any(c(
  "[OX2H]c1ccccc1", "c1cc([OX2H])ccc1",
  "c1ccc([OX2H])cc1", "c1c([OX2H])cccc1"
))

has_lignin_same_ring <- smart_match_any(c(
  "[OX2H]c1c([OD2;!R][C,c])cccc1",
  "[OX2H]c1cc([OD2;!R][C,c])ccc1",
  "[OX2H]c1ccc([OD2;!R][C,c])cc1",
  "c1([OD2;!R][C,c])c([OX2H])cccc1",
  "c1([OD2;!R][C,c])cc([OX2H])ccc1",
  "c1([OD2;!R][C,c])ccc([OX2H])cc1"
))

has_aromatic_oxidized_sidechain <- smart_match_any(c(
  "c[CH2][OX2H]", "c[CX3H1](=O)", "c[CX3](=O)[OX2H]", "c[CX3](=O)[#6]"
))

has_phenylpropanoid_like <- smart_match_any(c(
  "c[CX3]=[CX3][CX3]", "c[CX3]=[CX3][CX3](=O)",
  "c[CX4][CX3]=[CX3]", "c[CX4][CX4][CX3](=O)"
))

is_lignin <- only_CHO & has_O & is_aromatic_robust & (
  isT(has_lignin_same_ring) |
    (isT(has_phenol_on_benzene) & has_aromatic_oxidized_sidechain) |
    has_phenylpropanoid_like
)

is_phenolic <- isT(has_phenol_on_benzene) & !is_lignin

is_n_compound <- has_N

has_furan <- smart_match_any(c("o1cccc1", "c1ccoc1", "[o;r5]", "[OX2;r5]", "C1=CC=CO1", "O1C=CC=C1"))
has_pyran <- smart_match_any(c("O1CCCCC1", "O1C=CCCC1", "O1CC=CCC1", "O1C=CC=CC1", "o1ccccc1",
                               "[OX2;r6]", "[o;r6]", "O=C1C=COC=C1"))
has_cyclopentenone <- smart_match_any(c("C1=CC(=O)CC1", "O=C1C=CCC1", "O=C1CC=CC1"))
has_lactone <- smart_match_vec("[C;R](=O)O[C;R]")
has_furfural_like <- smart_match_any(c("c1ccoc1C=O", "O=CC1=COC=C1"))
has_hmf_like <- smart_match_any(c("c1cc(CO)oc1C=O", "O=CC1=COC(CO)=C1"))

# Oxygen functional groups (defined here so the small-fragment carbohydrate rule
# below can use them; also reused for lipid subclassing further down).
has_carboxylic_acid <- smart_match_any(c("[CX3](=O)[OX2H]", "[CX3](=O)[OX1-]"))
has_ester <- smart_match_any(c("[CX3](=O)[OX2][#6]"))
has_alcohol <- smart_match_any(c("[CX4][OX2H]"))
has_aldehyde <- smart_match_any(c("[CX3H1](=O)"))
has_ketone <- smart_match_any(c("[#6][CX3](=O)[#6]"))

has_cyclopentanone <- smart_match_any(c("O=C1CCCC1", "C1CCC(=O)C1"))
has_cyclopentanedione <- smart_match_any(c("O=C1C(=O)CCC1", "O=C1CCC(=O)C1", "OC1=C(O)CCC1"))
has_polyol <- smart_match_any(c("[CX4;!R]([OX2H])[CX4;!R][OX2H]", "[CX4]([OX2H])[CX4]([OX2H])"))
has_anhydrosugar <- smart_match_any(c(
  "[OX2H][CX4]1[CX4]([OX2H])[CX4]([OX2H])[CX4]2[OX2][CX4]1[OX2]2",
  "[CX4;R]([OX2H])[CX4;R]([OX2H])[CX4;R][OX2;R]"
))

# --- Long-chain (lipid / wax-like) aliphatic detection ---------------------
# The single most useful discriminator between a sugar pyrolysis fragment and a
# lipid/wax-derived oxygenate is the presence of a long unbranched acyclic
# methylene run. Eight consecutive acyclic sp3 carbons cannot arise from a
# cellulose/hemicellulose monomer but is the defining feature of fatty acids,
# fatty alcohols, long-chain methyl ketones (2-tritriacontanone) and long-chain
# alkylfurans (2-tetradecylfuran).
long_chain_smarts <- paste0("[CX4;!R]", strrep("[CX4;!R]", LONG_CHAIN_MIN_C - 1))
has_long_alkyl_chain <- smart_match_any(long_chain_smarts)

is_small_lactone <- isT(has_lactone) & !is.na(carbon_count) & carbon_count <= 8
is_large_lactone <- isT(has_lactone) & !is.na(carbon_count) & carbon_count >= 9

# Any non-aromatic C/H/O species that carries a long methylene run, or is simply
# large and aliphatic, is lipid/wax-like and must NOT be called a carbohydrate.
# It is routed to the broad "Aliphatic compounds" class instead (detailed label
# "Other aliphatic compounds"), together with fatty acids and fatty alcohols.
is_long_chain_aliphatic_CHO <- only_CHO & has_O & !has_N & !has_benzene_robust & (
  isT(has_long_alkyl_chain) |
    (!is.na(carbon_count) & carbon_count >= 14)
)

# --- Carbohydrate decision tree -------------------------------------------
# A compound is called a carbohydrate pyrolysis product only if it POSITIVELY
# shows a sugar-derived motif. The former blanket rule ("everything that
# contains only C, H and O is Carb") is removed: it swept fatty acids,
# long-chain methyl ketones and long alkylfurans into the sugar pool.
#
#   1. furan / pyran ring, but NOT benzo-fused (-> oxygenated aromatic) and NOT
#      carrying a long alkyl chain (-> lipid-derived alkylfuran)
#   2. furfural / HMF / methylfurfural type aldehydes
#   3. cyclopentenone / cyclopentanone / cyclopentanedione (levoglucosenone
#      secondary products)
#   4. small lactones (<= C8), e.g. 2(5H)-furanone, angelicalactone
#   5. vicinal-diol / anhydrosugar motifs (levoglucosan, glycolaldehyde dimers)
#   6. small (<= C6) acyclic oxygenates - acetic acid, hydroxyacetone,
#      2,3-butanedione, etc. These are genuinely ambiguous (hemicellulose acetyl
#      groups vs. lipid/protein fragments); they are kept as carbohydrate to
#      stay comparable with published Py-GC/MS schemes, but are flagged via
#      flag_ambiguous_small_oxygenate so the assignment can be audited or
#      reversed downstream.
is_carb_ring_marker <- (
  (isT(has_furan) & !has_benzene_robust & !isT(has_long_alkyl_chain)) |
    (isT(has_pyran) & !has_benzene_robust) |
    isT(has_cyclopentenone) | isT(has_cyclopentanone) | isT(has_cyclopentanedione) |
    isT(has_furfural_like) | isT(has_hmf_like) |
    is_small_lactone
)
is_carb_polyol_marker <- only_CHO & (isT(has_polyol) | isT(has_anhydrosugar)) &
  !is.na(carbon_count) & carbon_count <= 12
is_small_acyclic_oxygenate <- only_CHO & has_O & !has_N & !has_benzene_robust &
  !isT(has_long_alkyl_chain) &
  !is.na(carbon_count) & carbon_count <= 6 &
  (isT(has_carboxylic_acid) | isT(has_aldehyde) | isT(has_ketone) |
     isT(has_alcohol) | isT(has_ester))

is_carbohydrate_derived <- (is_carb_ring_marker | is_carb_polyol_marker | is_small_acyclic_oxygenate) &
  !is_long_chain_aliphatic_CHO
# Audit flag: small oxygenates admitted on convention rather than on a
# sugar-specific motif (acetic acid, 1-penten-3-one, ...).
flag_ambiguous_small_oxygenate <- is_small_acyclic_oxygenate & !is_carb_ring_marker & !is_carb_polyol_marker

carbohydrate_subclass <- case_when(
  !is_carbohydrate_derived ~ NA_character_,
  isT(has_furfural_like) ~ "furfural-like",
  isT(has_hmf_like) ~ "HMF-like",
  isT(has_anhydrosugar) ~ "anhydrosugar-like",
  isT(has_polyol) ~ "polyol / sugar fragment",
  is_small_lactone ~ "small lactone",
  isT(has_furan) ~ "furan-like",
  isT(has_pyran) ~ "pyran-like",
  isT(has_cyclopentenone) | isT(has_cyclopentanone) | isT(has_cyclopentanedione) ~ "cyclopentanone/-enone-like",
  is_small_acyclic_oxygenate ~ "small acyclic oxygenate (ambiguous origin)",
  TRUE ~ NA_character_
)

# Aliphatic hydrocarbons are reported under the broad "Aliphatic compounds"
# class. Any alkane not covered by the short/long n-alkane bins is kept as
# "Other aliphatic compounds" rather than a separate "Other alkanes" bucket.

is_nonarom_hc <- only_CH & !is_aromatic_robust
is_alkane <- rep(FALSE, n)
is_alkene <- rep(FALSE, n)

for (i in idx_smiles) {
  if (!is_nonarom_hc[i] || is.na(smiles[i]) || smiles[i] == "") next
  if (!grepl("=", smiles[i]) && !grepl("#", smiles[i])) is_alkane[i] <- TRUE
  if (grepl("=", smiles[i])) is_alkene[i] <- TRUE
}

is_short_chain_alkane <- is_alkane & !is.na(carbon_count) & carbon_count >= 7 & carbon_count <= 22
is_long_chain_alkane <- is_alkane & !is.na(carbon_count) & carbon_count >= 23 & carbon_count <= 33
is_hydrocarbon <- is_alkane | is_alkene

# "Other aliphatic compounds" is the detailed bin that now also receives the
# OXYGENATED aliphatics which used to be mis-filed as carbohydrates: fatty
# acids, fatty alcohols, long-chain methyl ketones (2-tritriacontanone) and
# long-chain alkylfurans (2-tetradecylfuran), plus any short alkane outside the
# C7-C22 window. The broad class stays "Aliphatic compounds", so the requested
# 4-way detail (Short-chain alkanes / Long-chain alkanes / Alkenes / Other
# aliphatic compounds) is preserved.
is_aliphatic_oxygenate <- only_CHO & has_O & !has_N & !has_benzene_robust & !is_carbohydrate_derived
is_other_aliphatic <- is_aliphatic_oxygenate |
  (is_alkane & !is_short_chain_alkane & !is_long_chain_alkane)

hydrocarbon_subclass <- case_when(
  is_short_chain_alkane ~ "short-chain alkane C7-C22",
  is_long_chain_alkane ~ "long-chain alkane C23-C33",
  is_aliphatic_oxygenate & isT(has_carboxylic_acid) & isT(has_long_alkyl_chain) ~ "long-chain fatty acid (lipid-derived)",
  is_aliphatic_oxygenate & isT(has_ester) & isT(has_long_alkyl_chain) ~ "wax ester / fatty acid ester (lipid-derived)",
  is_aliphatic_oxygenate & isT(has_ketone) & isT(has_long_alkyl_chain) ~ "long-chain ketone (wax-derived)",
  is_aliphatic_oxygenate & isT(has_alcohol) & isT(has_long_alkyl_chain) ~ "long-chain alcohol (lipid-derived)",
  is_aliphatic_oxygenate & isT(has_furan) & isT(has_long_alkyl_chain) ~ "long-chain alkylfuran (lipid-derived)",
  is_aliphatic_oxygenate ~ "other aliphatic oxygenate",
  is_alkane ~ "other aliphatic compound",
  is_alkene ~ "alkene / unsaturated hydrocarbon",
  TRUE ~ NA_character_
)

is_arom_hc <- only_CH & is_aromatic_robust
ring_labels <- rep(NA_real_, n)
if (length(idx_smiles) > 0) {
  s <- smiles[idx_smiles]
  ring_labels[idx_smiles] <- (str_count(s, "(?<!%)\\d") + str_count(s, "%\\d\\d")) / 2
}
n_arom_atoms <- rep(NA_integer_, n)
if (length(idx_smiles) > 0) {
  n_arom_atoms[idx_smiles] <- str_count(smiles[idx_smiles], "[cnosp]")
}
n_rings <- ring_labels

# --- PAH definition: CONVENTIONAL (operational) scope ----------------------
# Strict IUPAC usage reserves "PAH" for ortho-fused, fully aromatic ring
# systems (naphthalene, phenanthrene, pyrene). This pipeline deliberately
# follows the WIDER convention used throughout the Py-GC/MS and environmental
# literature - the same convention behind the US-EPA 16 priority PAHs, which
# already include the partly hydrogenated species fluorene, acenaphthene and
# acenaphthylene. Under that convention a PAH is any hydrocarbon carrying two
# or more rings of which at least one is aromatic, i.e.:
#   (a) ortho-fused aromatics            naphthalene, methylnaphthalenes,
#                                        phenanthrene, anthracene, pyrene
#   (b) cyclopenta-fused / hydroaromatic indene, methylindene, indane,
#                                        fluorene, acenaphthene
#   (c) non-fused ring assemblies        biphenyl, terphenyl, binaphthyl
# This keeps the classification internally consistent: previously fluorene was
# a PAH while indene (the same kind of sp3-bridged part-aromatic bicyclic) was
# demoted to "monocyclic", which was indefensible. The three sub-types remain
# separable via aromatic_subclass and the flag_* columns, so a stricter,
# fused-only PAH sum can still be recomputed downstream at any time.
# Report the definition in the methods, e.g.: "PAHs were defined as compounds
# containing two or more rings, including partially hydrogenated and
# cyclopenta-fused aromatics (e.g. indene, fluorene) and non-fused biaryls
# (e.g. biphenyl)."
has_fused_polyaromatic <- smart_match_any(c(
  "c1ccc2ccccc2c1",            # naphthalene core (two ortho-fused aromatic rings)
  "c1ccc2cc3ccccc3cc2c1",      # linear/angular tri-aromatic (anthracene/phenanthrene)
  "C1=CC=C2C=CC=CC2=C1",       # Kekule naphthalene core
  "C1=CC2=CC=CC=C2C=C1"        # Kekule naphthalene core (alt.)
))
# Benzene ortho-fused to a 5-membered carbocycle. [#6] (any carbon, aromatic or
# not) lets one pattern cover indene, methylindene, indane, fluorene and the
# acenaphthene 5-ring in both aromatic and Kekule notation.
has_cyclopenta_fused_arom <- smart_match_any(c(
  "c1ccc2c(c1)[#6][#6][#6]2",
  "c1ccc2c(c1)[#6][#6]2",
  "C1=CC=C2C(=C1)[#6][#6][#6]2",
  "C1=CC=C2C(=C1)[#6][#6]2"
))
# Non-fused ring assemblies: two aromatic ring atoms joined by a bond that is
# itself not part of any ring ("!@"). Matches biphenyl / phenylnaphthalene /
# terphenyl but NOT toluene or styrene (their side chains are aliphatic).
has_biaryl <- smart_match_any(c(
  "[cR]!@[cR]",
  "c1ccccc1-c1ccccc1",
  "C1=CC=C(C=C1)C1=CC=CC=C1"
))
# Ring-count fallback, independent of SMARTS/aromaticity perception: a ring
# closure digit pair == one ring. An aromatic hydrocarbon with >= 2 rings is
# polycyclic by definition, whatever the fusion type.
has_multiple_rings <- !is.na(n_rings) & n_rings >= 2

is_fused_PAH <- is_arom_hc & isT(has_fused_polyaromatic)
is_hydroaromatic_PAH <- is_arom_hc & !is_fused_PAH & isT(has_cyclopenta_fused_arom)
is_biaryl_PAH <- is_arom_hc & !is_fused_PAH & !is_hydroaromatic_PAH & isT(has_biaryl)
is_other_polycyclic_arom_hc <- is_arom_hc & !is_fused_PAH & !is_hydroaromatic_PAH &
  !is_biaryl_PAH & has_multiple_rings

# PAH (conventional scope) = fused + hydroaromatic + biaryl + any other
# multi-ring aromatic hydrocarbon. MAH is then strictly single-ring.
is_PAH_candidate <- is_fused_PAH | is_hydroaromatic_PAH | is_biaryl_PAH | is_other_polycyclic_arom_hc
is_mono_arom_hc <- is_arom_hc & !is_PAH_candidate
is_other_oxygenated_aromatic <- is_aromatic_robust & has_O & only_CHO & !is_lignin & !is_phenolic & !is_PAH_candidate
is_other_aromatic <- is_aromatic_robust & !is_lignin & !is_phenolic & !is_n_compound & !is_mono_arom_hc & !is_PAH_candidate

aromatic_subclass <- case_when(
  is_fused_PAH ~ "PAH - ortho-fused polyaromatic",
  is_hydroaromatic_PAH ~ "PAH - cyclopenta-fused / hydroaromatic",
  is_biaryl_PAH ~ "PAH - non-fused biaryl / ring assembly",
  is_other_polycyclic_arom_hc ~ "PAH - other polycyclic aromatic hydrocarbon",
  is_mono_arom_hc ~ "monocyclic aromatic hydrocarbon",
  is_other_oxygenated_aromatic ~ "oxygenated aromatic",
  is_other_aromatic ~ "other aromatic",
  TRUE ~ NA_character_
)

# Name-based safety net: rescue clear manual identifications when SMILES are
# missing or when the structural rules miss a known Py-GC/MS marker. This is
# intentionally conservative; true Unknown/Unidentified rows stay unidentified.
name_grepl <- function(pat, x) {
  out <- rep(FALSE, length(x))
  idx <- which(!is.na(x))
  if (length(idx)) out[idx] <- grepl(pat, x[idx], perl = TRUE)
  out
}
compound_name_lc <- tolower(trimws(as.character(df[[name_col]])))
unknown_name <- is.na(compound_name_lc) |
  compound_name_lc %in% c("", "na", "n/a", "unknown", "unknow", "unidentified", "not identified")
has_usable_name <- !unknown_name

name_carbon <- rep(NA_integer_, n)
name_carbon_match <- stringr::str_match(compound_name_lc, "\\b(?:n[- ]?)?c(\\d{1,2})\\b")
idx_name_carbon <- which(!is.na(name_carbon_match[, 2]))
if (length(idx_name_carbon)) {
  name_carbon[idx_name_carbon] <- suppressWarnings(as.integer(name_carbon_match[idx_name_carbon, 2]))
}

lignin_name_patterns <- c(
  "eugen", "guaiac", "syring", "vanill", "conifer", "sinap",
  "ferulic", "ferulate", "feruloyl",
  "coumaryl", "coumaric", "coumarate", "coumaroyl",
  "methoxyphenol", "dimethoxyphenol", "trimethoxyphenol",
  "cinnamyl", "cinnamaldehyde", "cinnamic", "cinnamate"
)
name_lignin_related <- name_grepl(paste(lignin_name_patterns, collapse = "|"), compound_name_lc)
name_phenolic <- name_grepl("phenol|cresol", compound_name_lc)
# Some finalized isomer groups span a benzofuran and a vinylphenol candidate.
# Their common defensible structural level is oxygenated aromatic; the name of
# the phenolic candidate must not force the unresolved group into Phenolics.
name_mixed_oxygenated_aromatic_isomer <- has_usable_name &
  name_grepl("benzofuran", compound_name_lc) &
  name_grepl("phenol", compound_name_lc)
name_n_compound <- name_grepl("nitrile|pyridine|pyrrole|indole|imidazole|pyrazine|pyrimidine|amine|amide|aniline", compound_name_lc)
# --- Long-chain (lipid / wax) name-net ------------------------------------
# Evaluated BEFORE the carbohydrate net and given priority over it in the final
# case_when. This is what stops "C11-C12 fatty acid", "2-tetradecylfuran" and
# "2-tritriacontanone" from being filed as sugars when their SMILES is missing.
long_alkyl_stem <- paste0(
  "(oct|non|dec|undec|dodec|tridec|tetradec|pentadec|hexadec|heptadec|octadec|nonadec|",
  "eicos|heneicos|docos|tricos|tetracos|pentacos|hexacos|heptacos|octacos|nonacos|",
  "triacont|hentriacont|dotriacont|tritriacont|tetratriacont|pentatriacont)"
)
name_long_chain_aliphatic <- has_usable_name & (
  name_grepl("fatty acid|fatty alcohol|palmitic|stearic|oleic|linoleic|linolenic|myristic|lauric|arachidic|behenic|lignoceric|wax ester", compound_name_lc) |
    name_grepl(paste0("\\b\\S*", long_alkyl_stem, "an(oic|ol|al|one|amide)?\\b"), compound_name_lc) |
    name_grepl(paste0(long_alkyl_stem, "yl"), compound_name_lc)
) & !name_grepl("benzene|phenol|cresol|guaiac|syring|vanill|pyridine|pyrrole|indole|nitrile|amine|aniline|naphthalene|phenanthrene|biphenyl", compound_name_lc)

# --- Carbohydrate name-net ------------------------------------------------
# Positive sugar-pyrolysis motifs only. The old net additionally listed
# "fatty acid|palmitic|stearic|..." because every C/H/O compound was Carb by
# convention; those entries are now handled by name_long_chain_aliphatic.
# Benzofurans stay excluded (oxygenated aromatics).
name_carbohydrate <- has_usable_name & name_grepl(
  paste(
    "furfural|hydroxymethylfurfural|levoglucosan|levoglucosenone|anhydroglucose|anhydrosugar",
    "furanone|furfuryl|furan|pyranone|maltol",
    "cyclopentanedione|cyclopentenolone|cyclopentenone|cyclopenten-[0-9]+-one|cyclopentanone",
    "acetic acid|formic acid|propionic acid|propanoic acid|butyric acid|butanoic acid",
    "hydroxyacetone|acetol|glycolaldehyde|butanedione|glucose|xylose|mannose|galactose",
    sep = "|"
  ),
  compound_name_lc
) &
  !name_grepl("benzofuran|dibenzofuran", compound_name_lc) &
  !name_long_chain_aliphatic
name_short_chain_alkane <- has_usable_name & (
  name_grepl("\\b(heptane|octane|nonane|decane|undecane|dodecane|tridecane|tetradecane|pentadecane|hexadecane|heptadecane|octadecane|nonadecane|eicosane|heneicosane|docosane)\\b", compound_name_lc) |
    (!is.na(name_carbon) & name_carbon >= 7 & name_carbon <= 22 & name_grepl("alkane", compound_name_lc))
)
name_long_chain_alkane <- has_usable_name & (
  name_grepl("\\b(tricosane|tetracosane|pentacosane|hexacosane|heptacosane|octacosane|nonacosane|triacontane|hentriacontane|dotriacontane|tritriacontane)\\b", compound_name_lc) |
    (!is.na(name_carbon) & name_carbon >= 23 & name_carbon <= 33 & name_grepl("alkane", compound_name_lc))
)
# Conventional PAH name-net. Now ALSO carries the cyclopenta-fused /
# hydroaromatic species (indene, indane, acenaphthene) and the non-fused ring
# assemblies (biphenyl, terphenyl), which previously routed to MAH via the old
# name_other_arom_hc rule. "indene" also catches "1H-Indene, 3-methyl-" and
# other methyl-indene isomers reported by NIST.
name_pah <- has_usable_name & name_grepl(
  paste(
    "naphthalene|phenanthrene|anthracene|pyrene|fluorene|fluoranthene|chrysene",
    "acenaphthylene|acenaphthene|perylene|triphenylene|coronene|benzo\\[",
    "indene|indane|indan\\b|biphenyl|diphenyl|terphenyl|azulene|biphenylene|stilbene",
    sep = "|"
  ),
  compound_name_lc
) & !name_grepl(
  "indole|benzofuran|dibenzofuran|phenol|cresol|nitrile|amine|amide|aniline|ether|oxide|sulfide|sulfone|methanone|ketone|carboxaldehyde|carboxylic|quinone",
  compound_name_lc
)
name_mono_arom_hc <- has_usable_name & name_grepl("benzene|toluene|xylene|ethylbenzene|propylbenzene|methylbenzene|styrene|vinylbenzene|cumene|cymene|mesitylene", compound_name_lc) &
  !name_grepl("phenol|benzofuran|dibenzofuran|nitrile|amine|amide|aniline|naphthalene|phenanthrene|anthracene|pyrene|biphenyl|diphenyl|indene|indane", compound_name_lc)
# Oxygenated aromatics (benzofurans, aryl aldehydes/ketones) -> detailed
# "Oxygenated aromatic compounds". Note that indene / biphenyl are NO LONGER
# routed here or to MAH: under the conventional definition adopted above they
# are PAHs and are captured by name_pah.
name_oxy_aromatic <- has_usable_name & name_grepl("benzofuran|dibenzofuran|benzaldehyde|acetophenone|benzoic|anisole|veratrole|chromene|chromone|xanthene", compound_name_lc)
# name_alkene must NOT swallow aromatic-hydrocarbon names that simply end in
# "-ene" (xylene, styrene, indene, fluorene, naphthalene, ...). The greedy
# "\\b[a-z]*ene\\b" pattern used to do exactly that, which is why those aromatics
# were mislabelled as alkenes. Explicitly exclude aromatic / heteroaromatic
# names AND anything already claimed by an aromatic name-net.
name_alkene <- has_usable_name & name_grepl("alkene|\\b[a-z]*ene\\b", compound_name_lc) &
  !name_grepl("benzene|toluene|xylene|styrene|cumene|cymene|mesitylene|indene|indane|fluorene|fluoranthene|naphthalene|phenanthrene|anthracene|pyrene|chrysene|acenaphthylene|acenaphthene|azulene|stilbene|biphenyl|diphenyl|terphenyl|benzofuran|dibenzofuran|phenol|cresol|guaiac|syring|vanill|aniline|pyridine|pyrrole|furan|thiophene|imidazole|pyrazine|indole", compound_name_lc) &
  !name_pah & !name_mono_arom_hc & !name_oxy_aromatic

missing_or_invalid_structure <- is.na(smiles) | !mol_ok
is_unidentified <- missing_or_invalid_structure & !has_usable_name
cat_final <- rep(NA_character_, n)
cat_final[is_unidentified] <- "Unidentified"

remaining <- is.na(cat_final)
cat_final[remaining] <- case_when(
  (is_lignin[remaining] | name_lignin_related[remaining]) ~ "Lignin-derived compounds",
  name_mixed_oxygenated_aromatic_isomer[remaining] ~ "Oxygenated aromatic compounds",
  (is_phenolic[remaining] | name_phenolic[remaining]) ~ "Phenolic compounds",
  (is_n_compound[remaining] | name_n_compound[remaining]) ~ "N-containing compounds",
  # Lipid / wax-like long-chain aliphatic oxygenates are resolved BEFORE the
  # carbohydrate branch: fatty acids, fatty alcohols, long-chain methyl ketones
  # (2-tritriacontanone) and long alkylfurans (2-tetradecylfuran) are aliphatic,
  # not sugar-derived.
  (is_long_chain_aliphatic_CHO[remaining] | name_long_chain_aliphatic[remaining]) ~ "Other aliphatic compounds",
  (is_carbohydrate_derived[remaining] | name_carbohydrate[remaining]) ~ CARB_LABEL,
  # Aromatic hydrocarbons are resolved BEFORE the alkane/alkene branches so that
  # a benzene/styrene/xylene/indene/fluorene ring can never be captured by an
  # "-ene" name rule or a mis-perceived double bond.
  # PAH (conventional scope) = any aromatic hydrocarbon with >= 2 rings, i.e.
  # ortho-fused polyaromatics PLUS hydroaromatics (indene, methyl-indene,
  # fluorene) PLUS non-fused biaryls (biphenyl). MAH is strictly single-ring.
  (is_PAH_candidate[remaining] | name_pah[remaining]) ~ "Polycyclic aromatic hydrocarbons",
  (is_mono_arom_hc[remaining] | name_mono_arom_hc[remaining]) ~ "Monocyclic aromatic hydrocarbons",
  (is_other_oxygenated_aromatic[remaining] | name_oxy_aromatic[remaining]) ~ "Oxygenated aromatic compounds",
  (is_short_chain_alkane[remaining] | name_short_chain_alkane[remaining]) ~ "Short-chain alkanes",
  (is_long_chain_alkane[remaining] | name_long_chain_alkane[remaining]) ~ "Long-chain alkanes",
  (is_alkene[remaining] | name_alkene[remaining]) ~ "Alkenes",
  # Remaining non-aromatic C/H/O species and off-window alkanes. The former
  # blanket rule "any leftover C/H/O compound is a carbohydrate" is DELETED:
  # a compound now has to earn the carbohydrate label from a positive sugar
  # motif, otherwise it lands here (aliphatic) or in "Others".
  (is_other_aliphatic[remaining] | is_alkane[remaining]) ~ "Other aliphatic compounds",
  missing_or_invalid_structure[remaining] ~ "Unidentified",
  has_non_target_element[remaining] ~ "Others",
  # Any other valid structure not captured by a predefined class is also
  # structurally miscellaneous.
  TRUE ~ "Others"
)

# Two-level taxonomy: cat_final is the fine-grained (detailed) label; the broad
# label groups the detailed classes into the top-level categories used for
# reporting (restores broad_structural_category / detailed_structural_category).
detailed_structural_category <- cat_final
broad_structural_category <- case_when(
  cat_final == "Lignin-derived compounds" ~ "Lignin-derived compounds",
  cat_final == "Phenolic compounds" ~ "Phenolic compounds",
  cat_final == "N-containing compounds" ~ "N-containing compounds",
  cat_final == CARB_LABEL ~ CARB_LABEL,
  cat_final %in% c("Short-chain alkanes", "Long-chain alkanes", "Other aliphatic compounds", "Alkenes") ~ "Aliphatic compounds",
  # Single aromatic umbrella (excludes lignin/phenolic, which keep their own broad
  # classes): MAH + PAH + oxygenated aromatics all roll up to "Other aromatic
  # compounds".
  cat_final %in% c("Monocyclic aromatic hydrocarbons", "Polycyclic aromatic hydrocarbons", "Oxygenated aromatic compounds") ~ "Other aromatic compounds",
  cat_final == "Unidentified" ~ "Unidentified",
  TRUE ~ "Others"
)

classified <- df %>%
  mutate(
    !!smiles_col := smiles,
    match_factor_used = match_factor,
    low_confidence_annotation = low_conf_id,
    broad_structural_category = broad_structural_category,
    detailed_structural_category = detailed_structural_category,
    `structural category` = cat_final,
    carbohydrate_subclass = carbohydrate_subclass,
    hydrocarbon_subclass = hydrocarbon_subclass,
    aromatic_subclass = aromatic_subclass,
    flag_valid_smiles = mol_ok,
    flag_lignin_relaxed = is_lignin,
    flag_name_lignin_related = name_lignin_related,
    flag_phenolic = is_phenolic,
    flag_n_compound = is_n_compound,
    flag_carbohydrate_derived = is_carbohydrate_derived,
    flag_hydrocarbon = is_hydrocarbon,
    flag_alkane = is_alkane,
    flag_alkene = is_alkene,
    flag_short_chain_alkane = is_short_chain_alkane,
    flag_long_chain_alkane = is_long_chain_alkane,
    flag_other_aromatic = is_other_aromatic,
    # PAH sub-flags: keep the strict (fused-only) subset separable so a
    # narrower PAH sum can be recomputed from the same output file.
    flag_fused_PAH = is_fused_PAH,
    flag_hydroaromatic_PAH = is_hydroaromatic_PAH,
    flag_biaryl_PAH = is_biaryl_PAH,
    flag_other_polycyclic_aromatic_hydrocarbon = is_other_polycyclic_arom_hc,
    flag_PAH_candidate = is_PAH_candidate,
    flag_mono_aromatic_hydrocarbon = is_mono_arom_hc,
    flag_other_oxygenated_aromatic = is_other_oxygenated_aromatic,
    flag_ambiguous_small_oxygenate = flag_ambiguous_small_oxygenate,
    flag_long_chain_aliphatic_CHO = is_long_chain_aliphatic_CHO,
    flag_aliphatic_oxygenate = is_aliphatic_oxygenate,
    flag_long_alkyl_chain = isT(has_long_alkyl_chain),
    flag_other_aliphatic = is_other_aliphatic,
    n_rings_from_smiles = n_rings,
    carbon_count = carbon_count,
    oxygen_count = oxygen_count,
    nitrogen_count = nitrogen_count
  )

area <- read.csv(area_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!"GlobalFeatureID" %in% names(area)) names(area)[1] <- "GlobalFeatureID"
area_final <- area %>% semi_join(classified %>% select(GlobalFeatureID), by = "GlobalFeatureID")

complete <- classified %>%
  select(
    GlobalFeatureID,
    compound_name = all_of(name_col),
    final_SMILES = all_of(smiles_col),
    broad_structural_category,
    detailed_structural_category,
    `structural category`,
    carbohydrate_subclass,
    hydrocarbon_subclass,
    aromatic_subclass,
    final_score_0_1000 = any_of("final_score_0_1000"),
    confidence = any_of("confidence"),
    final_database = any_of("final_database"),
    final_Formula = any_of("final_Formula"),
    final_MW = any_of("final_MW"),
    carbon_count,
    oxygen_count,
    nitrogen_count
  ) %>%
  left_join(area_final, by = "GlobalFeatureID")

summary_df <- classified %>%
  count(broad_structural_category, detailed_structural_category, name = "n_compounds") %>%
  arrange(desc(n_compounds), broad_structural_category, detailed_structural_category)

# Audit table for the PAH class: lets you check at a glance which compounds
# entered PAH through the conventional (non-strict) part of the definition.
pah_audit <- classified %>%
  filter(detailed_structural_category == "Polycyclic aromatic hydrocarbons") %>%
  select(
    GlobalFeatureID,
    compound_name = all_of(name_col),
    aromatic_subclass,
    flag_fused_PAH, flag_hydroaromatic_PAH, flag_biaryl_PAH,
    flag_other_polycyclic_aromatic_hydrocarbon,
    n_rings_from_smiles
  )
pah_audit_file <- file.path(out_dir, "final_foundin_gt3_PAH_definition_audit.csv")

write.csv(classified, classified_file, row.names = FALSE, na = "")
write.csv(complete, complete_file, row.names = FALSE, na = "")
write.csv(summary_df, summary_file, row.names = FALSE, na = "")
write.csv(pah_audit, pah_audit_file, row.names = FALSE, na = "")
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list(
      complete_compound_area_table = complete,
      classified_compounds = classified,
      classification_summary = summary_df,
      PAH_definition_audit = pah_audit
    ),
    complete_xlsx_file
  )
}

cat("\n=== Reviewed Compound Classification Summary ===\n")
cat("Input rows:", n, "\n")
cat("SMILES valid:", sum(mol_ok), "/", n, "\n")
cat("Match factor column:", ifelse(is.na(mf_col), "none", mf_col), "\n")
cat("Score gate: disabled; classification uses SMILES and compound name only\n")
cat("PAH definition: CONVENTIONAL (>=2 rings; incl. hydroaromatics and biaryls)\n")
cat("  ortho-fused PAH        :", sum(is_fused_PAH), "\n")
cat("  hydroaromatic PAH      :", sum(is_hydroaromatic_PAH), "\n")
cat("  non-fused biaryl PAH   :", sum(is_biaryl_PAH), "\n")
cat("  other polycyclic AH    :", sum(is_other_polycyclic_arom_hc), "\n")
cat("Carbohydrate label in use:", CARB_LABEL, "\n")
cat("Carbohydrate rule: positive sugar motif required (blanket C/H/O rule removed)\n")
cat("  carbohydrate-assigned    :", sum(is_carbohydrate_derived), "\n")
cat("  of which ambiguous small oxygenates:", sum(flag_ambiguous_small_oxygenate), "\n")
cat("  long-chain aliphatic C/H/O moved out of Carb:", sum(is_long_chain_aliphatic_CHO), "\n")
cat("Area matrix rows retained:", nrow(area_final), "\n\n")
print(summary_df)
cat("\nClassified compounds:", classified_file, "\n")
cat("Complete compound-area table:", complete_file, "\n")
if (file.exists(complete_xlsx_file)) cat("Complete compound-area workbook:", complete_xlsx_file, "\n")
cat("Summary:", summary_file, "\n")
cat("PAH definition audit:", pah_audit_file, "\n")
out_dir <- pipeline_out_dir

# Final tables for downstream analysis.
classified_raw <- read.csv(
  file.path(intermediate_dir, "final_foundin_gt3_classified_compounds.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!all(c("broad_structural_category", "detailed_structural_category") %in% names(classified_raw))) {
  stop(
    "Classification rules script did not output broad_structural_category and detailed_structural_category: ",
    classification_rules_label
  )
}

if ("structural category" %in% names(classified_raw)) {
  mismatch <- !is.na(classified_raw[["structural category"]]) &
    classified_raw[["structural category"]] != "" &
    classified_raw[["structural category"]] != classified_raw$detailed_structural_category
  if (any(mismatch)) {
    warning(
      "The legacy structural category differs from detailed_structural_category for ",
      sum(mismatch),
      " rows. Keeping broad/detail from the rules script unchanged."
    )
  }
}

if (!file.exists(file.path(intermediate_dir, "final_foundin_gt3_classification_summary.csv"))) {
  warning("Classification summary was not written by: ", classification_rules_label)
}

drop_cols <- c(
  "structural category",
  "carbohydrate_subclass",
  "lipid_subclass",
  "hydrocarbon_subclass",
  "aromatic_subclass",
  "carbon_count",
  "oxygen_count",
  "nitrogen_count",
  "match_factor_used",
  "low_confidence_annotation"
)

internal_cols <- c(
  "member_original_best_identifications",
  "final_SMILES_source",
  "final_SMILES_match_type",
  "flag_valid_smiles",
  "flag_lignin_relaxed",
  "flag_name_lignin_related",
  "flag_phenolic",
  "flag_n_compound",
  "flag_carbohydrate_derived",
  "flag_lipid_like",
  "flag_hydrocarbon",
  "flag_alkane",
  "flag_alkene",
  "flag_short_chain_alkane",
  "flag_long_chain_alkane",
  "flag_other_aromatic",
  "flag_PAH_candidate",
  "flag_mono_aromatic_hydrocarbon",
  "flag_other_oxygenated_aromatic"
)

if (!"compound_name" %in% names(classified_raw)) {
  name_col <- if ("final_identification" %in% names(classified_raw)) "final_identification" else "GlobalFeatureID"
  classified_raw$compound_name <- classified_raw[[name_col]]
}
if (!"CompoundID" %in% names(classified_raw)) classified_raw$CompoundID <- classified_raw$GlobalFeatureID
if (!"final_SMILES" %in% names(classified_raw)) classified_raw$final_SMILES <- ""

remove_cols <- intersect(c(drop_cols, internal_cols), names(classified_raw))
classified_out <- classified_raw[, setdiff(names(classified_raw), remove_cols), drop = FALSE]

area_table_out <- merge(
  classified_out[, c("GlobalFeatureID", "CompoundID", "compound_name", "final_SMILES", "broad_structural_category", "detailed_structural_category"), drop = FALSE],
  area,
  by = c("GlobalFeatureID", "CompoundID"),
  all.y = TRUE,
  sort = FALSE
)

preferred_front <- c(
  "GlobalFeatureID", "CompoundID", "compound_name", "final_SMILES",
  "broad_structural_category", "detailed_structural_category"
)
reorder_cols <- function(df) {
  front <- preferred_front[preferred_front %in% names(df)]
  df[, c(front, setdiff(names(df), front)), drop = FALSE]
}
classified_out <- reorder_cols(classified_out)
area_table_out <- reorder_cols(area_table_out)

summary_broad <- as.data.frame(table(classified_out$broad_structural_category), stringsAsFactors = FALSE)
names(summary_broad) <- c("class", "n_compounds")
summary_broad$level <- "broad"
summary_detail <- as.data.frame(table(classified_out$detailed_structural_category), stringsAsFactors = FALSE)
names(summary_detail) <- c("class", "n_compounds")
summary_detail$level <- "detail"
summary_out <- rbind(summary_broad, summary_detail)
summary_out <- summary_out[, c("level", "class", "n_compounds")]
summary_out <- summary_out[order(summary_out$level, -summary_out$n_compounds, summary_out$class), ]

write.csv(classified_out, file.path(out_dir, "classified_compounds.csv"), row.names = FALSE, na = "")
write.csv(area_table_out, file.path(out_dir, "classified_area_table.csv"), row.names = FALSE, na = "")
write.csv(summary_out, file.path(out_dir, "classification_summary.csv"), row.names = FALSE, na = "")

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list(
      classified_area_table = area_table_out,
      classified_compounds = classified_out,
      classification_summary = summary_out
    ),
    file.path(out_dir, "classified_area_table.xlsx")
  )
}

cat("Classification complete:\n")
cat("Rules: ", classification_rules_label, "\n", sep = "")
cat("Main table: ", file.path(out_dir, "classified_area_table.csv"), "\n", sep = "")
cat("Summary: ", file.path(out_dir, "classification_summary.csv"), "\n", sep = "")
