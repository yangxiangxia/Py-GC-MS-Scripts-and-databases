# This workflow outlines the process for handling GC-MS (Gas Chromatography-Mass Spectrometry) data, 
# including peak deconvolution, alignment, compound annotation, and classification. 
# The eRah package is used for peak processing, while ClassyFire is utilized for compound classification. 
# MassBank is used as the primary database for compound annotation. 
# eRah is employed for peak deconvolution, alignment, and compound annotation.
# ClassyFire is used for compound classification.
#Please note that the current version of the MassBank database does not include mass spectral data for the internal standards N-tetracosane-d50 and chrysene-d12. 
#The base-peak m/z of n-tetracosane-d50 is 66 (sometimes 50), and the base-peak m/z of chrysene-d12 is 240.
# MassBank_NIST_20241126_rev.msp include mass spectral data for the internal standards N-tetracosane-d50 and chrysene-d12. 
# Currently, InChIKeys from the MassBank database are used to retrieve compound classifications. 
#This script provides results for both peak area and peak height.
######


# Loading Required Packages
library(erah) #https://cran.r-project.org/web/packages/erah/vignettes/erah-manual.html
library(tidyr) 
library(data.table)
library(dplyr)
library(ggplot2)
library(classyfireR) #https://aberhrml.github.io/classyfireR/articles/Introduction_to_classyfireR.html
library(purrr) #https://purrr.tidyverse.org
library(tibble)
library(stringr)
library(future)
# Define the base directory for data storage and retrieval.
base_dir <- "~/Desktop/larch" #file path

##############################
# 1. Loading GC-MS Data
##############################
# GC-MS data is stored in CDF format. The .cdf files are retrieved from the specified directory and 
# processed using the eRah package to create an Instrumental Table.

# Retrieving CDF format data
cdf_dir <- file.path(base_dir, "gcms_cdf_data/")
cdf_files <- list.files(cdf_dir, pattern="\\.cdf$", full.names=TRUE)

# Creating phenotype information
instrumental <- createInstrumentalTable(cdf_files)
sample_groups <- rep("Sample", length(cdf_files))  # tentatively, all samples are set to be "Sample" group
phenotype <- createPhenoTable(cdf_files, sample_groups)

# Creating an eRah experiment object
exp <- newExp(instrumental=instrumental, phenotype=phenotype, info="GCMS Untargeted Analysis")


#######################################
# 2. Loading Database Data
#######################################

# MassBank_NIST.msp: https://github.com/MassBank/MassBank-data/releases/tag/2024.11

MassBank_NIST_file_path <- file.path(base_dir, "MassBank_NIST_20241126_rev.msp") 

# Importing the MassBank MSP database
MassBank_NIST_db <- importMSP(file = MassBank_NIST_file_path, 
                              DB.name = "MassBank_NIST",              # Database name
                              DB.version = "2024.11", # Database version
                              DB.info = "https://github.com/MassBank/MassBank-data/releases/tag/2024.11"  # Database information
                              )

# Save and load
save(MassBank_NIST_db, file = file.path(base_dir, "MassBank_NIST_db.rda"))
load(file.path(base_dir, "MassBank_NIST_db.rda"))

# Parsing MSP data
# A custom function, read_msp_with_structures(), extracts Name, Formula, InChIKey, InChI, and SMILES
# from the MSP file and converts it into a data frame.
read_msp_with_structures <- function(msp_file) {
  lines <- readLines(msp_file)
  compounds <- list()
  
  current_name <- NULL
  current_inchikey <- NULL
  current_formula <- NULL
  current_smiles <- NULL
  current_inchi <- NULL 
  
  for (line in lines) {

    if (grepl("^Name:", line)) {
      current_name <- gsub("^Name: ", "", line)
    }
    if (grepl("^Formula:", line)) {
      current_formula <- gsub("^Formula: ", "", line)
    }
    if (grepl("^InChIKey:", line)) {
      current_inchikey <- gsub("^InChIKey: ", "", line)
    }
    if (grepl("^InChI:", line)) {
      current_inchi <- gsub("^InChI: ", "", line)
    }
    if (grepl("^SMILES:", line)) {
      current_smiles <- gsub("^SMILES: ", "", line)
    }
    
    if (line == "") {  
      if (!is.null(current_name) && !is.null(current_inchikey)) {
        compounds <- append(compounds, list(list(Name = current_name, 
                                                 InChIKey = current_inchikey,
                                                 Formula = current_formula,
                                                 SMILES = current_smiles,
                                                 InChI = ifelse(is.null(current_inchi), NA, current_inchi)))) 
      }
      
      current_name <- current_inchikey <- current_formula <- current_smiles <- current_inchi <- NULL
    }
  }
  
  # Convert to data frame
  if (length(compounds) > 0) {
    msp_df <- do.call(rbind, lapply(compounds, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  } else {
    stop("No valid compounds found in the MSP file.")
  }
  
  return(msp_df)
}

# Call a custom function to read MSP data,for the same name and formula, keep the one with the most complete structural information
MassBank_NIST_data <- read_msp_with_structures(MassBank_NIST_file_path)

MassBank_NIST_lookup <- MassBank_NIST_data %>%
  mutate(across(c(Name, Formula, InChIKey, InChI, SMILES),
                ~ na_if(str_trim(.), ""))) %>%
  mutate(struct_score = rowSums(across(c(InChIKey, InChI, SMILES), ~ !is.na(.)))) %>%
  group_by(Name, Formula) %>%
  slice_max(order_by = struct_score, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Name, Formula, InChIKey, InChI, SMILES)

#####################################
# 3. Peak Deconvolution
#####################################
plan(future::multisession,workers = 4)
# Setting deconvolution parameters
dec_params <- setDecPar(min.peak.width = 1.8, # Min peak width (sec)
                        min.peak.height = 1000, # Minimum peak height (intensity threshold), should be >5 times of noise
                        noise.threshold = 200,  # Noise threshold
                        avoid.processing.mz = c(0),
                        analysis.time = c(0, 92)) # Analysis time (min)

# Executing peak deconvolution
exp_dec <- deconvolveComp(exp, decParameters = dec_params)

# Save and load
save(exp_dec, file = file.path(base_dir, "exp_dec.rda"))
load(file.path(base_dir, "exp_dec.rda"))

# Extract the list of deconvoluted peaks
peak_list <- exp_dec@Data@FactorList

# Sample names (you can modify them according to your actual sample order)
sample_names <- names(peak_list)
if (is.null(sample_names)) {
  sample_names <- paste0("Sample_", seq_along(peak_list))
}

# Combine all sample peaks into one data frame and add a "sample" column
all_peaks <- bind_rows(
  lapply(seq_along(peak_list), function(i) {
    df <- as.data.frame(peak_list[[i]])
    df$sample <- sample_names[i]
    return(df)
  })
)

head(all_peaks)

# Save to a CSV file
write.csv(all_peaks, file = file.path(base_dir, "1_deconvoluted_peaks.csv"), row.names = FALSE)

## Detecting standard compound peaks
## Specific standard compounds (N-Tetracosane-d50, Chrysene-d12) are identified based on:
## base m/z values and expected retention time (RT)
## A verification step ensures that all samples contain these standard compounds. 
## If not detected, remeasuremant of samples is recommended.

std_info <- data.frame(
            Name = c("N-Tetracosane-d50", "Chrysene-d12"),
            Target_mz = c(66, 240),  # base m/z values
            Expected_RT = c(50.7, 52.7)  # expectedRT (decide based on raw data)
            )

peak_list <- exp_dec@Data@FactorList  # peak list
sample_names <- names(peak_list)  # sample name

std_dec <- data.frame(
      Sample = sample_names,
      N_Tetracosane_d50 = FALSE, Chrysene_d12 = FALSE, # presence or absence
      N_Tetracosane_d50_RT = NA, Chrysene_d12_RT = NA, # retention time
      N_Tetracosane_d50_Height = NA, Chrysene_d12_Height = NA # height
      )

## detecting standard compound peaks in all samples
for (i in seq_along(sample_names)) {
  sample_name <- sample_names[i]
  sample_peaks <- peak_list[[sample_name]] 
  
  for (j in seq_along(std_info$Name)) {
    std_name <- std_info$Name[j]
    target_mz <- std_info$Target_mz[j]
    expected_rt <- std_info$Expected_RT[j]
    
    # detecting peaks with m/z amd RT of standard compounds
    match <- sample_peaks %>%
      filter(abs(RT - expected_rt) < 0.5) %>%  # RT ± 0.5 min
      filter(grepl(paste0("\\b", target_mz, ","), Spectra))  # m/z
    
    if (nrow(match) > 0) {
      peak_height <- sum(as.numeric(match$`Peak Height`), na.rm = TRUE)   # sum of peak height
      peak_rt <- mean(match$RT, na.rm = TRUE)  # average RT
      
      if (std_name == "N-Tetracosane-d50") {
        std_dec$N_Tetracosane_d50[i] <- TRUE
        std_dec$N_Tetracosane_d50_RT[i] <- peak_rt
        std_dec$N_Tetracosane_d50_Height[i] <- peak_height
      } else if (std_name == "Chrysene-d12") {
        std_dec$Chrysene_d12[i] <- TRUE
        std_dec$Chrysene_d12_RT[i] <- peak_rt
        std_dec$Chrysene_d12_Height[i] <- peak_height
      }
    }
  }
}

print(std_dec)

## Evaluating retention time variability
## The maximum and minimum retention times (RT) for the standard compounds are calculated. 
## This RT variation (in seconds) is used to determine an appropriate max.time.dist parameter for alignment.
RT_diff_N_Tetracosane_d50 <- (max(std_dec$N_Tetracosane_d50_RT, na.rm = TRUE) - 
                                min(std_dec$N_Tetracosane_d50_RT, na.rm = TRUE)) * 60  # min -> sec
RT_diff_Chrysene_d12 <- (max(std_dec$Chrysene_d12_RT, na.rm = TRUE) - 
                           min(std_dec$Chrysene_d12_RT, na.rm = TRUE)) * 60  # min -> sec
print(paste("N_Tetracosane_d50:", RT_diff_N_Tetracosane_d50, "sec"))
print(paste("Chrysene_d12:", RT_diff_Chrysene_d12, "sec"))


#############################
# 4. Peak Alignment
#############################
# Setting alignment parameters
align_params <- setAlPar(min.spectra.cor = 0.80,  # Minimum spectral correlation (0 to 1)
                         max.time.dist = 120,　# Maximum retention time tolerance (sec)
                         mz.range = 46:650) # m/z range for analysis.

# Executing peak alignment
exp_align <- alignComp(exp_dec, alParameters = align_params) #Small number of samples
#exp_align <- alignComp(exp_dec, alParameters = align_params, blocks.size = 10 ) #large number of samples#7
# Save and load
save(exp_align, file = file.path(base_dir, "exp_align.rda"))
load(file.path(base_dir, "exp_align.rda"))

## Verifying standard compound alignment
## Standard compounds should have the same AlignID across all samples. 
## If discrepancies are found, the data quality should be re-evaluated and remeasurement of samples is recommended.

aligned_peaks <- exp_align@Results@Alignment
sample_names <- colnames(aligned_peaks)[-(1:5)] 

std_align <- data.frame(
        Sample = sample_names,
        N_Tetracosane_d50 = FALSE, Chrysene_d12 = FALSE, # presence/absence
        N_Tetracosane_d50_AlignID = NA, Chrysene_d12_AlignID = NA, # AlignID
        N_Tetracosane_d50_RT = NA, Chrysene_d12_RT = NA, # retention time
        N_Tetracosane_d50_Height = NA, Chrysene_d12_Height = NA # height
        )

## detecting standard compound peaks in all samples
for (i in seq_along(sample_names)) {
  sample_name <- sample_names[i]
  
  for (j in seq_along(std_info$Name)) {
    std_name <- std_info$Name[j]
    target_mz <- std_info$Target_mz[j]
    expected_rt <- std_info$Expected_RT[j]
    
    # detecting peaks with m/z amd RT of standard compounds
    match <- aligned_peaks %>%
      filter(abs(tmean - expected_rt) < 0.5) %>%  # RT ± 0.5 min
      filter(grepl(paste0("\\b", target_mz, ","), Spectra))  # m/z
    
    if (nrow(match) > 0) {
      peak_height <- sum(match[[sample_name]], na.rm = TRUE)  # Sum of peak height
      peak_rt <- mean(match$tmean, na.rm = TRUE)  # average RT
      align_id <- match$AlignID[1]  # first AlignID
      
      if (peak_height > 0) {  
        if (std_name == "N-Tetracosane-d50") {
          std_align$N_Tetracosane_d50[i] <- TRUE
          std_align$N_Tetracosane_d50_AlignID[i] <- align_id
          std_align$N_Tetracosane_d50_RT[i] <- peak_rt
          std_align$N_Tetracosane_d50_Height[i] <- peak_height
        } else if (std_name == "Chrysene-d12") {
          std_align$Chrysene_d12[i] <- TRUE
          std_align$Chrysene_d12_AlignID[i] <- align_id
          std_align$Chrysene_d12_RT[i] <- peak_rt
          std_align$Chrysene_d12_Height[i] <- peak_height
        }
      }
    }
  }
}

print(std_align)
write.csv(aligned_peaks, file = file.path(base_dir, "2.1_aligned_peaks_height.csv"), row.names = FALSE)

# Extracting the area using the alignList function
aligned_areas <- alignList(exp_align, by.area = TRUE)

# If aligned_areas is a list, convert it to a data frame
aligned_area_df <- as.data.frame(aligned_areas)

# View data structure
head(aligned_area_df)

# Save the data frame containing Spectra and Area as a CSV file
write.csv(aligned_area_df, file = file.path(base_dir, "2.2_aligned_peaks_area.csv"), row.names = FALSE)

##############################
# 5. Compound Annotation
##############################
## Identifying Compounds
## Annotate peaks by matching them with the MassBank_NIST database.
## Only the highest MatchFactor compound candidate is retained.
exp_annot <- identifyComp(exp_align, 
                          id.database = MassBank_NIST_db,  # MSP database
                          mz.range = 46:650,
                          n.putative = 1
)

# Save and load
save(exp_annot, file = file.path(base_dir, "exp_annot.rda"))
load(file.path(base_dir, "exp_annot.rda"))

# Load aligned peaks (if not in memory already)
aligned_peaks <- exp_align@Results@Alignment

# Load annotation result (if not in memory already)
annot_result <- idList(exp_annot, id.database = MassBank_NIST_db)

# Select the relevant columns in annot_result (keep only the columns you need)
annot_result_relevant <- annot_result %>%
  dplyr::select(AlignID, Name.1, MatchFactor.1, DB.Id.1, CAS.1, Formula.1)

# Merge annotation information
aligned_peaks_annotated <- aligned_peaks %>%
  left_join(annot_result_relevant, by = "AlignID")

# View the merged data frame
head(aligned_peaks_annotated)

# Export results
write.csv(aligned_peaks_annotated, 
          file = file.path(base_dir, "3_aligned_peaks_height_with_annotation.csv"), 
          row.names = FALSE)

# View the exported data frame
head(aligned_peaks_annotated)

####################################################
# 6.Duplicate peak consolidation
####################################################
##extract the base m/z from a Spectra string
get_base_mz <- function(spectra_str) {
  if (is.na(spectra_str) || spectra_str == "") return(NA)
  
  # Split the string into "mz,intensity" pairs
  pairs <- strsplit(spectra_str, " ")[[1]]
  mz_int_matrix <- do.call(rbind, strsplit(pairs, ","))  # Create a 2-column character matrix
  
  # Convert to numeric values
  mz <- as.numeric(mz_int_matrix[, 1])
  intensity <- as.numeric(mz_int_matrix[, 2])
  
  if (all(is.na(intensity))) return(NA)
  
  # Find the m/z value with the highest intensity (base peak)
  base_mz <- mz[which.max(intensity)]
  return(base_mz)
}

# Apply the function to the Spectra column to extract base m/z for each aligned peak
aligned_peaks_annotated$base_mz <- sapply(aligned_peaks_annotated$Spectra, get_base_mz)

# save the updated table with annotation and base m/z
write.csv(aligned_peaks_annotated,
          file = file.path(base_dir, "4.1_aligned_peaks_height_with_annotation_and_base_mz.csv"),
          row.names = FALSE)

head(aligned_peaks_annotated[, c("AlignID", "Spectra", "base_mz")])


# Merge Spectra, base_mz, Name.1, MatchFactor.1, DB.Id.1, CAS.1, Formula.1 into the area table
aligned_area_df <- read.csv(file.path(base_dir, "2.2_aligned_peaks_area.csv"))
aligned_peaks_annotated <- read.csv(file.path(base_dir, "4.1_aligned_peaks_height_with_annotation_and_base_mz.csv"))

# Make sure the annotation data frame
head(aligned_peaks_annotated)

# Extract annotation data
annot_result <- idList(exp_annot, id.database = MassBank_NIST_db)

# Extract relevant columns
aligned_peaks_annotated_relevant <- aligned_peaks_annotated %>%
  dplyr::select(AlignID, Spectra, base_mz, Name.1, MatchFactor.1, DB.Id.1, CAS.1, Formula.1) 

# Merge annotation data and m/z data into an area data frame
final_df <- aligned_area_df %>%
  left_join(aligned_peaks_annotated_relevant, by = "AlignID")

# Export the merged data frame as a CSV file
write.csv(final_df, file = file.path(base_dir, "4.2_aligned_peaks_area_with_annotation_base_mz.csv"), row.names = FALSE)

# Confirm data
head(final_df)

# Loading data
df <- read.csv(file.path(base_dir, "4.2_aligned_peaks_area_with_annotation_base_mz.csv"))

non_sample_cols <- c("AlignID", "Factor", "Spectra", "tmean", "base_mz", 
                     "MatchFactor.1", "Name.1", "Formula.1", "CAS.1", "DB.Id.1", "FoundIn")

sample_cols <- setdiff(names(df), non_sample_cols)

# Filter out sample columns that contain only numeric data
sample_cols <- sample_cols[sapply(df[, sample_cols], is.numeric)]

print(sample_cols)

## Grouping: by base_mz and tmean ± time window (the time setting is consistent with the alignment time distance)

time_window <- 120 / 60

df <- df %>% arrange(Name.1, base_mz, tmean)

grouped_rows <- list()
used_indices <- rep(FALSE, nrow(df))

# ===== Requires the name of a compound that triggers 50/66 fusions. =====
special_names <- c("Tetracosane-d50")
# ============================================

for (i in seq_len(nrow(df))) {
  if (used_indices[i]) next 
  
  ref_name <- df$Name.1[i]
  ref_mz   <- df$base_mz[i]
  ref_rt   <- df$tmean[i]
  
  same_name  <- df$Name.1 == ref_name
  exact_mz   <- df$base_mz == ref_mz
  time_close <- abs(df$tmean - ref_rt) <= time_window
  not_used   <- !used_indices
  
  # =====Only Tetracosane-d50 allows for 50/66 interoperability. =====
  special_cross_mz_ok_scalar <- (ref_name %in% special_names) && (ref_mz %in% c(50, 66))
  special_cross_mz_vec <- special_cross_mz_ok_scalar & (df$base_mz %in% c(50, 66))
  mz_ok <- exact_mz | special_cross_mz_vec
  # ============================================
  
  group_idx <- which(same_name & mz_ok & time_close & not_used)
  
  grouped_rows[[length(grouped_rows) + 1]] <- group_idx
  used_indices[group_idx] <- TRUE
}


nonzero_sample_count <- function(group_df, sample_cols) {
  sample_matrix <- group_df[, sample_cols, drop = FALSE]
  sample_matrix <- as.data.frame(sample_matrix)
  sample_matrix[] <- lapply(sample_matrix, function(x) as.numeric(as.character(x)))
  
  if (is.null(dim(sample_matrix)) || length(dim(sample_matrix)) < 2) {
    sample_matrix <- matrix(as.numeric(sample_matrix), ncol = length(sample_cols))
  }
  
  sample_has_peak <- colSums(sample_matrix > 0, na.rm = TRUE) > 0
  return(sum(sample_has_peak))
}

filtered_peaks <- lapply(grouped_rows, function(group) {
  group_df <- df[group, ]
  
  best <- group_df[which.max(group_df$MatchFactor.1), ]
  
  for (col in sample_cols) {
    best[[col]] <- max(as.numeric(as.character(group_df[[col]])), na.rm = TRUE)
  }
  
  best$FoundIn <- nonzero_sample_count(group_df, sample_cols)
  
  # ===== If the group is a 50/66 parallel group of Tetracosane-d50, 
  #then the base_mz value will be uniformly retained as 66 in the output. =====
  if (best$Name.1 %in% special_names && all(group_df$base_mz %in% c(50, 66))) {
    best$base_mz <- 66
  }
  # =====================================================================================
  
  return(best)
})

filtered_df <- bind_rows(filtered_peaks)

write.csv(
  filtered_df,
  file = file.path(base_dir, "5_filtered_annotated_peaks_area.csv"),
  row.names = FALSE
)

# Integrating annotation results
#The annotation results are merged with MassBank_NIST_lookup to append chemical information (Name, Formula, InChIKey, InChI, SMILES).

annot_result <- idList(exp_annot, id.database = MassBank_NIST_db) %>%
  mutate(Name.1   = stringr::str_trim(Name.1),
         Formula.1 = stringr::str_trim(Formula.1)) %>%
  left_join(MassBank_NIST_lookup %>% mutate(Name = str_trim(Name), Formula = str_trim(Formula)),
            by = c("Name.1" = "Name", "Formula.1" = "Formula"))

##Merge structure info into filtered_df
# Load the filtered peak table
# Load MassBank data (with structural information)
# It should contain columns: Name, Formula, InChIKey, InChI, SMILES
# Merge structural information: Based on compound name and molecular formula
# Append structure information to filtered_df
filtered_df <- read.csv(file.path(base_dir, "5_filtered_annotated_peaks_area.csv"), check.names = FALSE) %>%
  mutate(`Name.1`   = stringr::str_trim(`Name.1`),
         `Formula.1` = stringr::str_trim(`Formula.1`)) %>%
  left_join(MassBank_NIST_lookup %>% mutate(Name = str_trim(Name), Formula = str_trim(Formula)),
            by = c("Name.1" = "Name", "Formula.1" = "Formula"))

write.csv(filtered_df,
          file = file.path(base_dir, "6_filtered_peaks_area_with_InChIKey.csv"),
          row.names = FALSE)

## Confirming annotation of standard compounds
## The AlignID of the standard compounds is used to filter the annotation results, ensuring correct identification.
std_annot <- annot_result %>%
  filter(AlignID %in% unique(std_align$N_Tetracosane_d50_AlignID) | AlignID %in% unique(std_align$Chrysene_d12_AlignID))
print(std_annot)

for (id in unique(std_annot$AlignID)) {
  message(sprintf("Plotting profile for AlignID: %s", id))
  plotProfile(exp_annot, id)
}

# Visualizing initial raw peaks "Match factor" distribution
annot_result <- annot_result %>%
  mutate(MatchFactor.1 = as.numeric(MatchFactor.1))

ggplot(annot_result, aes(x = MatchFactor.1, y = FoundIn)) +
  geom_point() +
  labs(x = "Match Factor", y = "FoundIn") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
    axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)
  )

# Visualizing Match Factor distribution for filtered peaks
filtered_df <- read.csv(file.path(base_dir, "5_filtered_annotated_peaks_area.csv"))
filtered_df <- filtered_df %>%
  mutate(MatchFactor.1 = as.numeric(MatchFactor.1),
         FoundIn       = as.numeric(FoundIn))

p2 <- ggplot(filtered_df, aes(x = MatchFactor.1, y = FoundIn)) +
  geom_point() +                                  
  labs(x = "Match Factor",                       
       y = "Number of Samples (FoundIn)") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14)
  )

print(p2)

####################################################
# 7. Compound Classification (for filtered compounds)
#####################################################

# Using ClassyFire for Compound Classification
## Currently, InChIKeys from the annotation results are used to retrieve compound classifications via get_classification().
## Classification levels include Kingdom, Superclass, Class, and Subclass, Level 5-7.
## But, InChIKeys are missing for some compounds despite the presence of InChI and SMILES.
## So, Should consider using SMILES or needs to get InChIKeys for those samples.


#Load filtered compound table with structure
filtered_df_with_structure <- read.csv(file.path(base_dir, "6_filtered_peaks_area_with_InChIKey.csv"))

#  Extract InChIKeys from filtered table
filtered_inchikeys <- filtered_df_with_structure$InChIKey %>%
  unique() %>%
  na.omit()

cat("🔍 Number of unique InChIKeys to classify:", length(filtered_inchikeys), "\n")

# ClassyFire query classification
filtered_classification_list <- purrr::map(filtered_inchikeys, get_classification)
filtered_classification_list <- filtered_classification_list[!sapply(filtered_classification_list, is.null)]

# Extract categorical fields
extract_classification <- function(class_obj) {
  if (is.null(class_obj)) return(tibble())  # skip empty
  
  inchikey <- gsub("^InChIKey=", "", class_obj@meta$inchikey)
  smiles <- class_obj@meta$smiles
  classification_tbl <- class_obj@classification
  
  safe_extract <- function(tbl, level) {
    val <- tbl$Classification[tbl$Level == level]
    if (length(val) > 0) return(val) else return(NA)
  }
  
  tibble(
    InChIKey = inchikey,
    SMILES = smiles,
    Kingdom = safe_extract(classification_tbl, "kingdom"),
    Superclass = safe_extract(classification_tbl, "superclass"),
    Class = safe_extract(classification_tbl, "class"),
    Subclass = safe_extract(classification_tbl, "subclass"),
    Level_5 = safe_extract(classification_tbl, "level 5"),
    Level_6 = safe_extract(classification_tbl, "level 6"),
    Level_7 = safe_extract(classification_tbl, "level 7")
  )
}

#Summarize all classification results
filtered_classification <- map_dfr(filtered_classification_list, extract_classification)

cat("number of successfully classified compounds:", nrow(filtered_classification), "\n")

#Merge the classification information back into the filtered data
filtered_df_classified <- filtered_df_with_structure %>%
  left_join(filtered_classification, by = "InChIKey")

#Save the filter results with classification information
write.csv(filtered_df_classified,
          file = file.path(base_dir, "7_filtered_peaks_area_with_classification.csv"),
          row.names = FALSE)

## Create summary classification table
# Make sure two data are loaded/generated:
# 1. filtered_df_with_structure: filtered peak table + InChIKey
# 2. filtered_classification: classification information obtained from get_classification()

compound_summary <- filtered_df_with_structure %>%
  dplyr::select(AlignID, Name = Name.1, MatchFactor = MatchFactor.1, InChIKey) %>%
  left_join(
    filtered_classification,
    by = "InChIKey"
  ) %>%
  distinct()  #remove duplicate rows

# Export as a classification summary table
write.csv(compound_summary,
          file = file.path(base_dir, "8_classified_compounds_summary.csv"),
          row.names = FALSE)



#######################################################
# 8.Structural refinement
#######################################################
#In this step, you can use "7_filtered_peaks_area_with_classification.csv"
#Change the identification result of compounds with MF below 80 to "Unidentified", and set the corresponding SMILES to NA.

## Required packages ----
library(rcdk)
library(tidyverse)

## Input/output file paths ----
input_file  <- "7.1_filtered_peaks_area_with_classification.csv"
output_file <- "7.1_filtered_peaks_area_with_LgC_class_revised.csv"

## Load data ----
df <- read.csv(input_file, stringsAsFactors = FALSE)

## Extract SMILES ----
smiles <- df$SMILES

## Parse SMILES (rcdk::parse.smiles), skipping NA/empty ----
valid_idx    <- !is.na(smiles) & smiles != ""
smiles_valid <- smiles[valid_idx]

mols    <- parse.smiles(smiles_valid)  # list of IAtomContainer
n_total <- nrow(df)

## Helper function wrapping rcdk::matches ----
## Returns a logical vector with the same length as df
smart_match_vec <- function(smarts, mol_list, valid_index, n_total) {
  out <- rep(FALSE, n_total)
  if (length(mol_list) == 0) return(out)
  m <- rcdk::matches(smarts, mol_list)  # TRUE/FALSE for valid SMILES only
  out[valid_index] <- m
  out
}

## 1. Detect aromatic atoms (broad) ----
## "a" matches aromatic atoms in SMARTS
is_aromatic <- smart_match_vec("a", mols, valid_idx, n_total)

## 1'. Detect benzene ring (C6 aromatic ring) ----
## Lignin is defined as containing a benzene ring + phenolic OH + aryl-ether
is_aryl_benzene <- smart_match_vec("c1ccccc1", mols, valid_idx, n_total)

## 2. Detect phenolic OH ----
## [OX2H]c = hydroxyl attached to an aromatic carbon (Ar–OH)
has_phenolic_OH <- smart_match_vec("[OX2H]c", mols, valid_idx, n_total)

## 3. Detect aryl-ether / methoxy groups ----
## Lignin definition requires Ar–O–CH3 (e.g., guaiacol).
has_aryl_ether_core <- smart_match_vec("cO[CH3]", mols, valid_idx, n_total)
# has_aryl_ether_wide <- smart_match_vec("cO[C,c]", mols, valid_idx, n_total)  # optional broader rule

## Use methoxy-only version for strict lignin detection
has_aryl_ether <- has_aryl_ether_core
## To broaden detection, uncomment:
## has_aryl_ether <- has_aryl_ether_core | has_aryl_ether_wide

## 4. LgC / non-LgC aromatic / non-aromatic classification ----
LgC_flag <- is_aryl_benzene & has_phenolic_OH & has_aryl_ether
non_LgC_aromatics_flag <- is_aromatic & !LgC_flag

df$Aromatic_group <- NA_character_
df$Aromatic_group[!is_aromatic]             <- "non-aromatic"
df$Aromatic_group[non_LgC_aromatics_flag]   <- "non-LgC-aromatics"
df$Aromatic_group[LgC_flag]                 <- "LgC"

## (Optional) Additional phenol group breakdown ----
PhC_flag <- is_aromatic & has_phenolic_OH & !has_aryl_ether

df$Phenol_group <- NA_character_
df$Phenol_group[LgC_flag]                       <- "LgC"
df$Phenol_group[PhC_flag]                       <- "PhC"
df$Phenol_group[is_aromatic & !has_phenolic_OH] <- "other-aromatics"
df$Phenol_group[!is_aromatic]                   <- "non-aromatic"

## 5. Lignin vs non-lignin aromatic grouping ----
df$Lignin_nonLignin_group <- NA_character_
df$Lignin_nonLignin_group[df$Aromatic_group == "LgC"]               <- "Lignin and derivatives"
df$Lignin_nonLignin_group[df$Aromatic_group == "non-LgC-aromatics"] <- "Non-lignin aromatics"
df$Lignin_nonLignin_group[df$Aromatic_group == "non-aromatic"]      <- "non-aromatic"

## 6. Subclassification for non-aromatic compounds ----
df$NonArom_subgroup <- NA_character_

## 6-1) Non-aromatic Organic oxygen compounds → O-rich
df$NonArom_subgroup[
  df$Aromatic_group == "non-aromatic" &
    df$Superclass == "Organic oxygen compounds"
] <- "non-aromatic O-rich compounds"

## 6-2) Non-aromatic Organoheterocyclic compounds → heterocycles
df$NonArom_subgroup[
  df$Aromatic_group == "non-aromatic" &
    df$Superclass == "Organoheterocyclic compounds"
] <- "non-aromatic heterocycles"

## 6-3) All other non-aromatic
df$NonArom_subgroup[
  df$Aromatic_group == "non-aromatic" &
    is.na(df$NonArom_subgroup)
] <- "other non-aromatic"

## 6-4) Upgrade carbohydrates to polysaccharide-derived
## (from within non-aromatic O-rich compounds)
df$NonArom_subgroup[
  df$Aromatic_group == "non-aromatic" &
    df$Superclass == "Organic oxygen compounds" &
    df$Subclass == "Carbohydrates and carbohydrate conjugates"
] <- "polysaccharide-derived"

## 7. Create final Superclass_updated classification ----
df$Superclass_updated <- case_when(
  ## 7-1) N-containing classes
  df$Superclass == "Organic nitrogen compounds" ~ "Organic nitrogen compounds",
  df$Superclass == "Nucleosides, nucleotides, and analogues" ~ "Nucleosides/nucleotides",
  df$Superclass == "Alkaloids and derivatives" ~ "Alkaloids and derivatives",
  
  ## 7-2) Polysaccharide-derived
  df$Aromatic_group == "non-aromatic" &
    df$NonArom_subgroup == "polysaccharide-derived" ~ "polysaccharide-derived",
  
  ## 7-3) Other major biochemical groups
  df$Superclass == "Organic acids and derivatives" ~ "Organic acids and derivatives",
  df$Superclass == "Lipids and lipid-like molecules" ~ "Lipids and lipid-like molecules",
  df$Superclass %in% c("Hydrocarbons", "Hydrocarbon derivatives") ~ "Hydrocarbons",
  df$Superclass == "Organohalogen compounds" ~ "Organohalogen compounds",
  df$Superclass == "Organophosphorus compounds" ~ "Organophosphorus compounds",
  df$Superclass == "Organosulfur compounds" ~ "Organosulfur compounds",
  df$Superclass == "Organometallic compounds" ~ "Organometallic compounds",
  df$Superclass == "Unidentified" ~ "Unidentified",
  ## 7-4) Aromatic groups
  df$Aromatic_group == "LgC"               ~ "Lignin and derivatives",
  df$Aromatic_group == "non-LgC-aromatics" ~ "Non-lignin aromatics",
  
  ## 7-5) Non-aromatic subgroups
  df$Aromatic_group == "non-aromatic" &
    df$NonArom_subgroup == "non-aromatic O-rich compounds" ~ "non-aromatic O-rich compounds",
  df$Aromatic_group == "non-aromatic" &
    df$NonArom_subgroup == "non-aromatic heterocycles" ~ "non-aromatic heterocycles",
  df$Aromatic_group == "non-aromatic" &
    df$NonArom_subgroup == "other non-aromatic" ~ "non-aromatic phenylpropanoid and polyketide",
  
  ## 7-6) Fallback
  TRUE ~ "Unclassified"
)

## 8. Optional summary tables ----
cat("Superclass_updated counts:\n")
print(with(df, table(Superclass_updated)))

cat("\nSuperclass vs Superclass_updated:\n")
print(with(df, table(Superclass, Superclass_updated)))

cat("\nSuperclass × Aromatic_group:\n")
print(with(df, table(Superclass, Aromatic_group)))

cat("\nNonArom_subgroup counts:\n")
print(with(df, table(NonArom_subgroup)))

## 9. Write output ----
write.csv(df, output_file, row.names = FALSE)


#####################################
# 9. Saving the Complete Environment
#####################################
save.image(file = file.path(base_dir, "gsms_peak_processing.RData"))

