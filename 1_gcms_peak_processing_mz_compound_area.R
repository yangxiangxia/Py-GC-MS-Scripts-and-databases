## =============================================================================
## PyGCMS Compound Classification - Conservative Evidence-Based Approach
## SMILES-based identification with strict, rule-based criteria
## =============================================================================
##
## =========================
## Classification rules (overview)
## =========================
## The script assigns ONE final label (PyroClass) per compound using a strict
## priority order (top = highest priority).
##
## YOUR REQUIRED PRIORITY (STRICT):
## 1) Lignin and derivatives
## 2) Phenolic compounds
## 3) Locked ClassyFire categories [HARD LOCK]:
##    - Alkaloids and derivatives (superclass)
##    - Nucleosides/nucleotides (superclass)
##    - Lipids and lipid-like molecules (superclass)
##    - Carbohydrates and carbohydrate conjugates (subclass under Organic oxygen compounds)
##    - Unidentified (superclass, this section should be set by yourself, the match factor less than 80)
## 4) Benzene derivatives
## 5) Polycyclic aromatic hydrocarbons
## 6) Monocyclic aromatic hydrocarbons
## 7) Short-chain alkanes (C7-C22)
## 8) Long-chain alkanes (C23-C33)
## 9) n-Alkenes
## 10) Additional Lipids (structure-based extension)
## 11) N-heterocyclic compounds
## 12) Non-heterocyclic N compounds
## 13) Additional Carbohydrates (structure-based extension)
## 14) Others (including small alkanes C<7, large alkanes C>33, etc.)
##
## IMPORTANT:
## - HARD LOCK applies ONLY AFTER lignin + phenol rules.
##   So: if a compound is "Alkaloids and derivatives" but also matches Lignin
##   rule, it will be assigned "Lignin and derivatives" (higher priority).
##   If it matches Phenol rule, it will be assigned "Phenolic compounds".
##   Otherwise, it keeps its locked category and cannot be overridden by
##   any lower-priority rules.
## - Locked categories include superclasses (Alkaloids, Nucleosides, Lipids, 
##   Unidentified) and specific subclass (Carbohydrates).
##
## =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rcdk)
  library(stringr)
})

## -----------------------------
## Input / Output
## -----------------------------
input_file  <- "peaks_area_with_classification.csv"
output_file <- "peaks_area_with_classification_revised.csv"

df <- read.csv(input_file, stringsAsFactors = FALSE)
n  <- nrow(df)
cols <- names(df)

pick_first <- function(cands) {
  hit <- cands[cands %in% cols]
  if (length(hit) == 0) NA_character_ else hit[1]
}

smiles_col  <- pick_first(c("SMILES","Smiles","smiles"))
name_col    <- pick_first(c("Name","Name.1","Compound","Compound_name","CompoundName","PeakName"))
formula_col <- pick_first(c("Formula","Formula.1","MolecularFormula","formula"))

superclass_col <- pick_first(c("superclass","Superclass"))
class_col      <- pick_first(c("class","Class"))
subclass_col   <- pick_first(c("subclass","Subclass"))

if (is.na(smiles_col)) stop("No SMILES column found.")

smiles <- df[[smiles_col]]
nm     <- if (!is.na(name_col)) df[[name_col]] else rep(NA_character_, n)
form   <- if (!is.na(formula_col)) df[[formula_col]] else rep(NA_character_, n)

sc_raw <- if (!is.na(superclass_col)) df[[superclass_col]] else rep(NA_character_, n)
cl <- if (!is.na(class_col))      df[[class_col]]      else rep(NA_character_, n)
sb <- if (!is.na(subclass_col))   df[[subclass_col]]   else rep(NA_character_, n)

## Normalize superclass strings to avoid whitespace issues
sc <- ifelse(is.na(sc_raw), NA_character_, str_squish(sc_raw))

isT <- function(x) !is.na(x) & x

## -----------------------------
## (0) Locked ClassyFire superclasses and specific subclasses (HARD LOCK, but AFTER Lignin + Phenol)
## -----------------------------
locked_superclass <- c(
  "Alkaloids and derivatives",
  "Nucleosides, nucleotides, and analogues",
  "Lipids and lipid-like molecules",
  "Unidentified"
)
is_locked_superclass <- !is.na(sc) & sc %in% locked_superclass

## Locked carbohydrates (specific subclass under Organic oxygen compounds)
is_locked_carbohydrate <- (!is.na(sc) & sc == "Organic oxygen compounds" &
                             !is.na(sb) & sb == "Carbohydrates and carbohydrate conjugates")

## -----------------------------
## (1) Formula parsing
## -----------------------------
parse_formula_counts <- function(formula_str) {
  if (is.na(formula_str) || formula_str == "") return(integer(0))
  m <- str_match_all(formula_str, "([A-Z][a-z]?)([0-9]*)")[[1]]
  if (nrow(m) == 0) return(integer(0))
  elems <- m[,2]
  nums  <- m[,3]; nums[nums == ""] <- "1"
  out <- as.integer(nums); names(out) <- elems
  out
}

## -----------------------------
## (2) Parse SMILES -> molecules (rcdk)
## -----------------------------
mol_all  <- vector("list", n)
parse_ok <- rep(FALSE, n)

prepare_mol <- function(m) {
  tryCatch(rcdk::set.atom.types(m), error = function(e) NULL)
  tryCatch(rcdk::do.aromaticity(m), error = function(e) NULL)
  m
}

idx_smiles <- which(!is.na(smiles) & smiles != "")
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
  if (length(mols_valid) == 0) return(out)
  m <- tryCatch(rcdk::matches(smarts, mols_valid),
                error = function(e) rep(FALSE, length(mols_valid)))
  out[mol_ok] <- m
  out
}

## -----------------------------
## (3) Element detection (has_N / has_O / only_CH / only_CHO)
## -----------------------------
elem_set <- vector("list", n)
for (i in 1:n) {
  if (!is.na(form[i]) && form[i] != "") {
    cnt <- parse_formula_counts(form[i])
    elem_set[[i]] <- names(cnt)
  } else if (mol_ok[i]) {
    elem_set[[i]] <- unique(vapply(rcdk::get.atoms(mol_all[[i]]),
                                   rcdk::get.symbol, character(1)))
  } else {
    elem_set[[i]] <- character(0)
  }
}

has_N   <- map_lgl(elem_set, ~ "N" %in% .x)
has_O   <- map_lgl(elem_set, ~ "O" %in% .x)  # require true oxygen presence
only_CH <- map_lgl(elem_set, ~ length(.x) > 0 && all(.x %in% c("C","H")))
only_CHO <- map_lgl(elem_set, ~ length(.x) > 0 && all(.x %in% c("C","H","O")))

## -----------------------------
## (4) Aromatic + lignin-related structural detection
## -----------------------------
is_aromatic     <- smart_match_vec("a")
has_benzene     <- smart_match_vec("c1ccccc1")
has_phenolic_OH <- smart_match_vec("[OX2H]c")

## Aryl-ether (wide rule): Ar–O–(C or aromatic C), O must be exocyclic (!R)
has_aryl_ether_wide <- smart_match_vec("c[OD2;!R][C,c]") | smart_match_vec("c[OD2;!R]c")
has_aryl_ether <- has_aryl_ether_wide

## Specific LgC structures (exact SMILES hits)
## Define the three specific SMILES (source:pubchem) that should be classified as LgC:
## 1) C=CC1=CC=C(C=C1)O      - 4-Vinylphenol
## 2) CC1=CC=C(C=C1)O        - 4-Methylphenol/p-Cresol/Para-Cresol
## 3) CCC1=CC=C(C=C1)O       - 4-Ethylphenol/Para-Ethylphenol/P-Ethylphenol
specific_lgc_smiles <- c(
  "C=CC1=CC=C(C=C1)O",   # 4-Vinylphenol
  "CC1=CC=C(C=C1)O",     # 4-Methylphenol
  "CCC1=CC=C(C=C1)O"     # 4-Ethylphenol
)

is_specific_lgc <- rep(FALSE, n)
if (length(idx_smiles) > 0) {
  normalized_specific <- tolower(trimws(specific_lgc_smiles))
  for (i in idx_smiles) {
    s <- tolower(trimws(smiles[i]))
    if (!is.na(s) && s != "" && s %in% normalized_specific) {
      is_specific_lgc[i] <- TRUE
    }
  }
}

## Lignin flag: (benzene + phenolic OH + aryl-ether) OR exact SMILES hit
LgC_flag <- (isT(has_benzene) & isT(is_aromatic) & isT(has_phenolic_OH) & isT(has_aryl_ether)) | is_specific_lgc

## Phenol rule (note: lignin takes precedence later)
is_phenol_rule <- isT(has_benzene) & isT(is_aromatic) & isT(has_phenolic_OH)

## -----------------------------
## (5) PAH vs monocyclic aromatic hydrocarbons (hydrocarbons only)
## -----------------------------
norm_name <- function(x) {
  tolower(x) |> str_replace_all("[^a-z0-9]", "")
}

EPA16 <- c(
  "Naphthalene","Acenaphthylene","Acenaphthene","Fluorene","Phenanthrene","Anthracene",
  "Fluoranthene","Pyrene","Benzo(a)anthracene","Chrysene","Benzo(b)fluoranthene",
  "Benzo(k)fluoranthene","Benzo(a)pyrene","Dibenzo(a,h)anthracene",
  "Benzo(ghi)perylene","Indeno(1,2,3-cd)pyrene"
)
EPA16_hit <- norm_name(ifelse(is.na(nm), "", nm)) %in% norm_name(EPA16)

## Ring-closure pairs from SMILES
ring_labels <- rep(NA_integer_, n)
if (length(idx_smiles) > 0) {
  s <- smiles[idx_smiles]
  n_single <- str_count(s, "(?<!%)\\d")
  n_multi  <- str_count(s, "%\\d\\d")
  ring_labels[idx_smiles] <- n_single + n_multi
}
ring_pairs <- ring_labels / 2

## Aromatic atom count from SMILES (fallback)
n_arom_atoms_smiles <- rep(NA_integer_, n)
if (length(idx_smiles) > 0) {
  n_arom_atoms_smiles[idx_smiles] <- str_count(smiles[idx_smiles], "[cnosp]")
}

## Aromatic hydrocarbons (ONLY C and H)
is_arom_hc <- only_CH & isT(is_aromatic)

## PAH rule within aromatic hydrocarbons (ONLY C and H)
is_PAH <- is_arom_hc & (
  EPA16_hit |
    (!is.na(ring_pairs) & ring_pairs >= 2) |
    (!is.na(n_arom_atoms_smiles) & n_arom_atoms_smiles > 6)
)

## Monocyclic aromatic hydrocarbons (subset of aromatic hydrocarbons)
is_mono_arom_hc <- is_arom_hc & !is_PAH

## -----------------------------
## (6) Benzene derivatives (oxygenated non-phenolic aromatics)
## -----------------------------
## NOTE: require has_O so aromatic hydrocarbons cannot be captured here
is_benzene_derivative <- isT(has_benzene) &
  isT(is_aromatic) &
  has_O &
  only_CHO &
  !isT(has_phenolic_OH)

## -----------------------------
## (7) Non-aromatic hydrocarbons classification (SMILES-based)
## -----------------------------
is_nonarom_hc <- only_CH & !isT(is_aromatic)

## Carbon atom count from SMILES (for chain-length bins)
carbon_count <- rep(NA_integer_, n)
if (length(idx_smiles) > 0) {
  for (i in idx_smiles) {
    s <- smiles[i]
    if (is.na(s) || s == "") next
    s_clean <- s
    s_clean <- str_remove_all(s_clean, "(?<!%)\\d+")
    s_clean <- str_remove_all(s_clean, "%\\d\\d")
    carbon_count[i] <- str_count(s_clean, "C")
  }
}

## Identify alkanes (no double/triple bonds)
is_alkane <- rep(FALSE, n)
if (length(idx_smiles) > 0) {
  for (i in idx_smiles) {
    if (!is_nonarom_hc[i]) next
    s <- smiles[i]
    if (is.na(s) || s == "") next
    if (!grepl("=", s) && !grepl("#", s)) is_alkane[i] <- TRUE
  }
}

## Identify alkenes (contain double bonds)
is_alkene <- rep(FALSE, n)
if (length(idx_smiles) > 0) {
  for (i in idx_smiles) {
    if (!is_nonarom_hc[i]) next
    s <- smiles[i]
    if (is.na(s) || s == "") next
    if (grepl("=", s)) is_alkene[i] <- TRUE
  }
}

## Classify alkanes by chain length
is_short_chain_alkane <- is_alkane & !is.na(carbon_count) & carbon_count >= 7  & carbon_count <= 22
is_long_chain_alkane  <- is_alkane & !is.na(carbon_count) & carbon_count >= 23 & carbon_count <= 33

## Small alkanes (C < 7) - should not be classified as carbohydrates
is_small_alkane <- is_alkane & !is.na(carbon_count) & carbon_count < 7

## -----------------------------
## (8) Carbohydrates + lactones + Lipids extension
## -----------------------------
has_furan <- smart_match_vec("o1cccc1")
has_pyran <- smart_match_vec("O1CCCCC1")
has_cyclopentenone <- smart_match_vec("C1=CC(=O)CC1")
has_lactone <- smart_match_vec("[C;R](=O)O[C;R]")

## Furfural-like and HMF-like structures
has_furfural_like <- smart_match_vec("c1ccoc1C=O")  # Furfural: furan ring with aldehyde
has_hmf_like <- smart_match_vec("c1cc(CO)oc1C=O")   # 5-HMF: furan ring with CH2OH and aldehyde

is_small_lactone <- isT(has_lactone) & (is.na(carbon_count) | carbon_count <= 8)
is_large_lactone <- isT(has_lactone) & !is.na(carbon_count) & carbon_count >= 9

is_carb_marker <- isT(has_furan) | isT(has_pyran) | isT(has_cyclopentenone) | is_small_lactone

## Exclude small alkanes from short_CHO to prevent misclassification
is_short_CHO <- only_CHO & !is_carb_marker & !isT(has_benzene) &
  (is.na(carbon_count) | carbon_count <= 6) & !is_small_alkane

is_medium_CHO_lipid <- only_CHO & !is_carb_marker & !isT(has_benzene) &
  !is.na(carbon_count) & carbon_count >= 7 & carbon_count <= 8

is_long_CHO_lipid <- only_CHO & !is_carb_marker &
  !is.na(carbon_count) & carbon_count >= 9

is_carbohydrate <- is_carb_marker | is_short_CHO

## Furfural-like or HMF-like structures (CHO only) - structure-based extension
is_furfural_hmf <- only_CHO & (isT(has_furfural_like) | isT(has_hmf_like))

## Additional structure-based carbohydrates (not locked by ClassyFire)
is_structural_carb <- is_carbohydrate | is_furfural_hmf

## Lipids and lipid-like molecules (structure-based extension only)
## ClassyFire Lipids are already locked in priority 3
## Here we only add additional lipids identified by structure
is_structural_lipid <- is_medium_CHO_lipid | is_long_CHO_lipid | is_large_lactone

## -----------------------------
## (9) Nitrogen compounds split
## -----------------------------
has_aromatic_N <- smart_match_vec("[n]")
has_ring_N     <- smart_match_vec("[N;R]")

is_N_heterocycles <- has_N & (isT(has_aromatic_N) | isT(has_ring_N))
is_other_Ntg      <- has_N & !is_N_heterocycles

## -----------------------------
## (10) FINAL classification (single label, STRICT priority)
## -----------------------------
PyroClass_final <- case_when(
  ## 1) Lignin and derivatives (highest priority)
  LgC_flag ~ "Lignin and derivatives",
  
  ## 2) Phenolic compounds
  is_phenol_rule ~ "Phenolic compounds",
  
  ## 3) Locked ClassyFire categories (HARD LOCK, after lignin + phenol)
  !is.na(sc) & sc == "Nucleosides, nucleotides, and analogues" ~ "Nucleosides/nucleotides",
  is_locked_superclass & sc != "Nucleosides, nucleotides, and analogues" ~ sc,
  is_locked_carbohydrate ~ "Carbohydrates",
  
  ## 4) Benzene derivatives
  is_benzene_derivative ~ "Benzene derivatives",
  
  ## 5-6) Aromatic hydrocarbons (hydrocarbons only)
  is_PAH ~ "Polycyclic aromatic hydrocarbons",
  is_mono_arom_hc ~ "Monocyclic aromatic hydrocarbons",
  
  ## 7-9) Non-aromatic hydrocarbons (strict carbon number ranges)
  is_short_chain_alkane ~ "Short-chain alkanes",
  is_long_chain_alkane  ~ "Long-chain alkanes",
  is_alkene ~ "n-Alkenes",
  
  ## 10) Additional Lipids (structure-based, ClassyFire Lipids already locked at priority 3)
  is_structural_lipid ~ "Lipids and lipid-like molecules",
  
  ## 11-12) Nitrogen groups
  is_N_heterocycles ~ "N-heterocyclic compounds",
  is_other_Ntg ~ "Non-heterocyclic N compounds",
  
  ## 13) Additional Carbohydrates (structure-based, ClassyFire Carbohydrates already locked at priority 3)
  is_structural_carb ~ "Carbohydrates",
  
  ## 14) Fallback (includes small alkanes C<7, large alkanes C>33, and other unclassified compounds)
  TRUE ~ "Others"
)

## -----------------------------
## (11) Diagnostics (priority integrity checks)
## -----------------------------
cat("\n========================================\n")
cat("PyGCMS Classification Diagnostics\n")
cat("========================================\n")
cat("Total rows:", n, "\n")
cat("SMILES non-empty:", length(idx_smiles), "/", n, "\n")
cat("SMILES parsed successfully:", sum(parse_ok), "/", n, "\n")
cat("Molecules valid:", sum(mol_ok), "/", n, "\n\n")

cat("--- ClassyFire Lipids (locked at priority 3) ---\n")
is_classyfire_lipid <- !is.na(sc) & sc == "Lipids and lipid-like molecules"
cat("ClassyFire Lipids total:", sum(is_classyfire_lipid, na.rm = TRUE), "\n")
cat("Structural lipids (medium CHO, C7-C8):", sum(is_medium_CHO_lipid, na.rm = TRUE), "\n")
cat("Structural lipids (long CHO, C≥9):", sum(is_long_CHO_lipid, na.rm = TRUE), "\n")
cat("Structural lipids (large lactones, C≥9):", sum(is_large_lactone, na.rm = TRUE), "\n")
cat("Structural lipids total:", sum(is_structural_lipid, na.rm = TRUE), "\n")
cat("All Lipids (ClassyFire + structural):", sum(is_classyfire_lipid | is_structural_lipid, na.rm = TRUE), "\n\n")

cat("--- Carbohydrates identification ---\n")
cat("Structure-based carbohydrates:", sum(is_carbohydrate, na.rm = TRUE), "\n")
cat("Furfural/HMF-like:", sum(is_furfural_hmf, na.rm = TRUE), "\n")
cat("Structural carbohydrates total:", sum(is_structural_carb, na.rm = TRUE), "\n")
cat("ClassyFire carbohydrates (locked at priority 3):", sum(is_locked_carbohydrate, na.rm = TRUE), "\n")
cat("All Carbohydrates (ClassyFire + structural):", sum(is_locked_carbohydrate | is_structural_carb, na.rm = TRUE), "\n\n")

cat("--- Alkanes classification ---\n")
cat("Short-chain alkanes (C7-C22):", sum(is_short_chain_alkane, na.rm = TRUE), "\n")
cat("Long-chain alkanes (C23-C33):", sum(is_long_chain_alkane, na.rm = TRUE), "\n")
cat("Small alkanes (C<7) -> Other:", sum(is_small_alkane, na.rm = TRUE), "\n")
cat("Large alkanes (C>33) -> Other:", sum(is_alkane & !is.na(carbon_count) & carbon_count > 33, na.rm = TRUE), "\n\n")

cat("N-heterocyclic compounds:", sum(PyroClass_final == "N-heterocyclic compounds", na.rm = TRUE), "\n")
cat("Non-heterocyclic N compounds:", sum(PyroClass_final == "Non-heterocyclic N compounds", na.rm = TRUE), "\n")

cat("--- Locked ClassyFire superclasses and specific subclasses ---\n")
cat("Locked superclasses total:", sum(is_locked_superclass, na.rm = TRUE), "\n")
cat("  Alkaloids and derivatives:", sum(!is.na(sc) & sc == "Alkaloids and derivatives", na.rm = TRUE), "\n")
cat("  Nucleosides/nucleotides:", sum(!is.na(sc) & sc == "Nucleosides, nucleotides, and analogues", na.rm = TRUE), "\n")
cat("  Lipids and lipid-like molecules:", sum(is_classyfire_lipid, na.rm = TRUE), "\n")
cat("  Unidentified:", sum(!is.na(sc) & sc == "Unidentified", na.rm = TRUE), "\n")
cat("Locked carbohydrates (subclass):", sum(is_locked_carbohydrate, na.rm = TRUE), "\n")
cat("Total locked categories:", sum(is_locked_superclass | is_locked_carbohydrate, na.rm = TRUE), "\n\n")

## Check how many locked compounds are overridden by higher-priority rules (expected behavior)
cat("--- Locked overridden by higher priorities (expected) ---\n")
all_locked <- is_locked_superclass | is_locked_carbohydrate
cat("Locked but final = Lignin and derivatives:",
    sum(all_locked & PyroClass_final == "Lignin and derivatives", na.rm = TRUE), "\n")
cat("Locked but final = Phenolic compounds:",
    sum(all_locked & PyroClass_final == "Phenolic compounds", na.rm = TRUE), "\n\n")

## Check how many locked compounds are incorrectly overridden by LOWER priorities (should be 0)
cat("--- Locked incorrectly overridden by LOWER priorities (should be 0) ---\n")
locked_kept_superclass <- (is_locked_superclass & 
                             ((sc == "Nucleosides, nucleotides, and analogues" & PyroClass_final == "Nucleosides/nucleotides") |
                                (sc != "Nucleosides, nucleotides, and analogues" & PyroClass_final == sc)))
locked_kept_carb <- is_locked_carbohydrate & (PyroClass_final == "Carbohydrates")
locked_kept <- locked_kept_superclass | locked_kept_carb
locked_overridden_lower <- all_locked &
  !(PyroClass_final %in% c("Lignin and derivatives", "Phenolic compounds")) &
  !locked_kept

cat("Locked overridden by LOWER priorities:", sum(locked_overridden_lower, na.rm = TRUE), "\n")
if (sum(locked_overridden_lower, na.rm = TRUE) > 0) {
  cat("WARNING: Some locked categories were incorrectly overridden:\n")
  print(table(PyroClass_final[locked_overridden_lower]))
}
cat("\n")

## Check ClassyFire Lipids retention (should be 100% now)
cat("--- ClassyFire Lipids retention check ---\n")
classyfire_lipids_kept <- is_classyfire_lipid & PyroClass_final == "Lipids and lipid-like molecules"
classyfire_lipids_lost <- is_classyfire_lipid & PyroClass_final != "Lipids and lipid-like molecules"
cat("ClassyFire Lipids kept:", sum(classyfire_lipids_kept, na.rm = TRUE), "/", sum(is_classyfire_lipid, na.rm = TRUE), "\n")
cat("ClassyFire Lipids lost to other categories:", sum(classyfire_lipids_lost, na.rm = TRUE), "(should be 0)\n")
if (sum(classyfire_lipids_lost, na.rm = TRUE) > 0) {
  cat("WARNING: Some ClassyFire Lipids were overridden. Lost to:\n")
  print(table(PyroClass_final[classyfire_lipids_lost]))
}
cat("\n")

## Check ClassyFire Carbohydrates retention (should be 100% now)
cat("--- ClassyFire Carbohydrates retention check ---\n")
classyfire_carb_kept <- is_locked_carbohydrate & PyroClass_final == "Carbohydrates"
classyfire_carb_lost <- is_locked_carbohydrate & PyroClass_final != "Carbohydrates"
cat("ClassyFire Carbohydrates kept:", sum(classyfire_carb_kept, na.rm = TRUE), "/", sum(is_locked_carbohydrate, na.rm = TRUE), "\n")
cat("ClassyFire Carbohydrates lost to other categories:", sum(classyfire_carb_lost, na.rm = TRUE), "(should be 0)\n")
if (sum(classyfire_carb_lost, na.rm = TRUE) > 0) {
  cat("WARNING: Some ClassyFire Carbohydrates were overridden. Lost to:\n")
  print(table(PyroClass_final[classyfire_carb_lost]))
}
cat("\n")

cat("========================================\n")
cat("Final Classification Counts:\n")
cat("========================================\n")
print(table(PyroClass_final, useNA = "ifany"))
cat("\n")

## -----------------------------
## (12) Output
## -----------------------------
out <- df %>% mutate(PyroClass = PyroClass_final)
write.csv(out, output_file, row.names = FALSE)
