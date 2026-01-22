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
#exp_align <- alignComp(exp_dec, alParameters = align_params) #Small number of samples
exp_align <- alignComp(exp_dec, alParameters = align_params, blocks.size = 10 ) #large number of samples#7
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
# 7. Compound for ClassyFire Classification (for filtered compounds)
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



####################################
# 8.Structural refinement
####################################
#In this step, you can use "7_filtered_peaks_area_with_classification.csv"
#Change the identification result of compounds with MF below 80 to "Unidentified", and set the corresponding SMILES to NA.
 
## =============================================================================
## PyGCMS Compound Classification - Conservative Evidence-Based Approach
## SMILES-based identification with strict, rule-based criteria
## =============================================================================
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
structural_category_final <- case_when(
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

cat("N-heterocyclic compounds:", sum(structural_category_final == "N-heterocyclic compounds", na.rm = TRUE), "\n")
cat("Non-heterocyclic N compounds:", sum(structural_category_final == "Non-heterocyclic N compounds", na.rm = TRUE), "\n")

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
    sum(all_locked & structural_category_final == "Lignin and derivatives", na.rm = TRUE), "\n")
cat("Locked but final = Phenolic compounds:",
    sum(all_locked & structural_category_final == "Phenolic compounds", na.rm = TRUE), "\n\n")




cat("Final Classification Counts:\n")
cat("========================================\n")
print(table(structural_category_final, useNA = "ifany"))
cat("\n")

## -----------------------------
## (12) Output
## -----------------------------
out <- df %>% mutate(`structural category` = structural_category_final)
write.csv(out, output_file, row.names = FALSE)
