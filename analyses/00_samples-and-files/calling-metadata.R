# Load necessary tidyverse packages
library(tidyverse)

# List of samples
# This vector contains all sequenced samples in the APA project 
# Format: {individual}-{type} where type can be:
# - N: Normal
# - T: Tumor 
# - TS: Tumor Secondary
# - F: Functional
# - NF: Non-Functional
samples <- c(
  "APA10-N", "APA10-T", "APA12-N", "APA12-T", "APA13-N", "APA13-T",
  "APA15-N", "APA15-T", "APA17-N", "APA17-T", "APA18-N", "APA18-T",
  "APA19-N", "APA19-T", "APA1-N", "APA1-T", "APA23-N", "APA23-TS",
  "APA25-N", "APA25-T", "APA29-N", "APA29-T", "APA2-T", "APA31-N",
  "APA31-T", "APA32-N", "APA32-TS", "APA34-N", "APA34-T", "APA36-N",
  "APA36-T", "APA38-N", "APA38-T", "APA3-N", "APA3-T", "APA41-N",
  "APA41-T", "APA43-N", "APA43-T", "APA44-N", "APA44-T", "APA45-N",
  "APA45-T", "APA47-N", "APA47-T", "APA48-N", "APA48-T", "APA49-N",
  "APA49-T", "APA51-T", "APA52-N", "APA52-T", "APA53-N", "APA53-T",
  "APA55-T", "APA56-N", "APA56-T", "APA58-N", "APA58-T", "APA59-N",
  "APA59-T", "APA61-N", "APA61-T", "APA61-TS", "APA69-N", "APA69-T",
  "APA72-N", "APA72-T", "APA74-N", "APA74-T", "APA75-N", "APA75-T",
  "APA78-N", "APA78-T", "APA79-N", "APA79-T", "APA81-N", "APA81-T",
  "APA86-N", "APA86-T", "APA87-N", "APA87-T", "iAPA1-F", "iAPA1-NF",
  "iAPA1-N", "iAPA2-F", "iAPA2-NF", "iAPA2-N", "iAPA3-F", "iAPA3-NF",
  "iAPA3-N", "iAPA5-F", "iAPA5-NF", "iAPA5-N", "iAPA6-F", "iAPA6-NF",
  "iAPA6-N"
)

#' Extract individual and sample type from sample identifier.
#'
#' @param sample A character string with the sample identifier (e.g., "APA10-N")
#' @return A vector with [1] individual ID and [2] sample type
extract_info <- function(sample) {
  individual <- str_remove(sample, "-.*")
  type <- str_extract(sample, "(?<=-).+")
  return(c(individual, type))
}

#------------------------------------------------------------------------------
# Create a sample inclusion/exclusion status table
#------------------------------------------------------------------------------

# Define sample status table with detailed exclusion reasons
sample_status <- tribble(
  ~individual, ~exclude, ~reason,
  "APA2",     TRUE,     "MAPM - unselective dissection",
  "APA3",     TRUE,     "Other sample contamination in tumor",
  "APA32",    TRUE,     "Hyperplasia - No clear tumor, micronoduli in normal tissue",
  "APA38",    TRUE,     "Sample switch with unknown sample",
  "APA45",    TRUE,     "Duplicate of iAPA5",
  "APA56",    TRUE,     "High contamination with other sample",
  "APA69",    TRUE,     "MAPM - unselective dissection"
)

# Additional specific samples to exclude (not entire individuals)
specific_samples_to_exclude <- c("APA61-TS")

# Convert to tibble with all individuals
all_individuals <- tibble(
  individual = unique(str_remove(samples, "-.*"))
) %>%
  left_join(sample_status, by = "individual") %>%
  mutate(
    exclude = replace_na(exclude, FALSE),
    reason = replace_na(reason, "")
  )

# Extract samples and individuals to exclude from the status table
samples_to_exclude <- c(
  samples[str_remove(samples, "-.*") %in% filter(all_individuals, exclude)$individual],
  specific_samples_to_exclude
)
individuals_to_exclude <- filter(all_individuals, exclude)$individual

# Define analyses types to exclude (based on type combinations)
analyses_types_to_exclude <- c(
  "NvsN",   # Normal vs Normal (same individual)
  "NvsF",   # Normal vs Functional
  "NvsNF",  # Normal vs Non-Functional
  "FvsF",   # Functional vs Functional
  "NFvsNF", # Non-Functional vs Non-Functional
  "NvsT",   # Normal vs Tumor
  "NvsTS"   # Normal vs Tumor Secondary
)

# Normal sample types to identify in the filtering step
normal_types <- c("N", "NS")

#------------------------------------------------------------------------------
# Phase 1: Generate pairwise combinations
#------------------------------------------------------------------------------

# Create a tibble with all sample combinations
pairwise_analyses <- crossing(sample1 = samples, sample2 = samples) %>%
  # Remove self-comparisons
  filter(sample1 != sample2) %>%
  # Extract individual and type information
  mutate(
    individual1 = map_chr(sample1, ~ extract_info(.)[1]),
    individual2 = map_chr(sample2, ~ extract_info(.)[1]),
    type1 = map_chr(sample1, ~ extract_info(.)[2]),
    type2 = map_chr(sample2, ~ extract_info(.)[2])
  ) %>%
  # Filter to keep only pairs from the same individual but different sample types
  filter(individual1 == individual2 & type1 != type2) %>%
  # Create bam file names (assuming they match sample names, adjust if needed)
  mutate(
    bam1_file_basename = sample1,
    bam2_file_basename = sample2,
    analysis = str_c(type1, "vs", type2)
  ) %>%
  # Filter out unwanted analysis types
  filter(!analysis %in% analyses_types_to_exclude) %>%
  # Select and order columns
  select(
    sample1, sample2, bam1_file_basename, bam2_file_basename,
    individual1, individual2, analysis
  )

#------------------------------------------------------------------------------
# Generate tumor-only analyses
#------------------------------------------------------------------------------

# Create a tibble from samples, extract sample types
sample_info <- tibble(sample = samples) %>%
  mutate(
    individual = map_chr(sample, ~ extract_info(.)[1]),
    type = map_chr(sample, ~ extract_info(.)[2])
  )

# Extract tumor samples (anything that is not a normal sample)
tumor_samples <- sample_info %>% 
  filter(!type %in% normal_types) %>% 
  pull(sample)

# Create "Tumor only" analysis for tumor samples as a tibble
to_analysis <- tibble(
  sample1 = tumor_samples,
  sample2 = NA_character_,
  bam1_file_basename = tumor_samples,
  bam2_file_basename = NA_character_,
  # Using more specific individual identifiers for NF/F samples to avoid conflicts
  individual1 = map_chr(tumor_samples, function(sample) {
    info <- extract_info(sample)
    ind <- info[1]
    type <- info[2]
    # For samples with types F or NF, include the type in the individual identifier
    if (type %in% c("F", "NF")) {
      return(paste0(ind, "-", type))
    } else {
      return(ind)
    }
  }),
  individual2 = NA_character_,
  analysis = "To"
)

# Combine pairwise and tumor-only analyses
initial_metadata <- bind_rows(pairwise_analyses, to_analysis)

#------------------------------------------------------------------------------
# Phase 2: Apply corrections for known sample issues
#------------------------------------------------------------------------------

#' Apply corrections to sample metadata based on known issues.
#'
#' @param metadata_df A dataframe containing sample metadata
#' @return A dataframe with corrected values for all affected columns
apply_corrections <- function(metadata_df) {
  # Create a copy of the input dataframe to avoid direct mutation during processing
  corrected_df <- metadata_df
  
  #----------------------------------------------------------------------------
  # Correction 1: APA12/APA13 Tumor Switch
  # Issue: The tumor samples for patients APA12 and APA13 were swapped
  # Fix: We need to swap sample names and bam_file_basenames to use the correct
  #      genomic data, but keep individual identifiers matching the metadata
  #----------------------------------------------------------------------------
  
  # First, gather the original values before any changes
  # Extract the original basenames associated with each tumor sample
  apa12t_basename <- metadata_df %>% 
    filter(sample1 == "APA12-T" | sample2 == "APA12-T") %>%
    pull(bam1_file_basename) %>%
    unique() %>%
    first()
  
  apa13t_basename <- metadata_df %>% 
    filter(sample1 == "APA13-T" | sample2 == "APA13-T") %>%
    pull(bam1_file_basename) %>%
    unique() %>%
    first()
  
  # Identify all rows that need to be modified for the APA12/13 tumor swap
  # For sample1
  apa12_sample1_rows <- which(corrected_df$sample1 == "APA12-T")
  apa13_sample1_rows <- which(corrected_df$sample1 == "APA13-T")
  
  # For sample2
  apa12_sample2_rows <- which(corrected_df$sample2 == "APA12-T")
  apa13_sample2_rows <- which(corrected_df$sample2 == "APA13-T")
  
  # Apply corrections for sample1
  if (length(apa12_sample1_rows) > 0) {
    # For APA12-T, change to APA13-T and swap bam_file_basename
    corrected_df$sample1[apa12_sample1_rows] <- "APA13-T"
    corrected_df$bam1_file_basename[apa12_sample1_rows] <- apa13t_basename
    # Keep individual1 as APA12 (unchanged) to maintain the metadata organization
  }
  
  if (length(apa13_sample1_rows) > 0) {
    # For APA13-T, change to APA12-T and swap bam_file_basename 
    corrected_df$sample1[apa13_sample1_rows] <- "APA12-T"
    corrected_df$bam1_file_basename[apa13_sample1_rows] <- apa12t_basename
    # Keep individual1 as APA13 (unchanged) to maintain the metadata organization
  }
  
  # Apply corrections for sample2
  if (length(apa12_sample2_rows) > 0) {
    # For APA12-T, change to APA13-T and swap bam_file_basename
    corrected_df$sample2[apa12_sample2_rows] <- "APA13-T"
    corrected_df$bam2_file_basename[apa12_sample2_rows] <- apa13t_basename
    # Keep individual2 unchanged
  }
  
  if (length(apa13_sample2_rows) > 0) {
    # For APA13-T, change to APA12-T and swap bam_file_basename
    corrected_df$sample2[apa13_sample2_rows] <- "APA12-T"
    corrected_df$bam2_file_basename[apa13_sample2_rows] <- apa12t_basename
    # Keep individual2 unchanged
  }
  
  # Note: By changing both the sample names AND bam_file_basename but keeping
  # the individual columns unchanged, we ensure the correct BAM files are used
  # while maintaining the individual's identity for output naming
  
  #----------------------------------------------------------------------------
  # Correction 2: APA58 T/N Switch
  # Issue: The tumor and normal samples for patient APA58 were swapped
  # Fix: We swap sample designations and ensure analysis type is TvsN
  #----------------------------------------------------------------------------
  
  # Identify all rows involving APA58 samples
  apa58_rows <- which(
    str_detect(corrected_df$sample1, "^APA58") | 
    (!is.na(corrected_df$sample2) & str_detect(corrected_df$sample2, "^APA58"))
  )
  
  # Apply corrections for each affected row
  for (row_idx in apa58_rows) {
    sample1 <- corrected_df$sample1[row_idx]
    sample2 <- corrected_df$sample2[row_idx]
    
    # Handle sample1 swap
    if (sample1 == "APA58-N") {
      # Change N to T in sample name and BAM file basename
      corrected_df$sample1[row_idx] <- "APA58-T"
      corrected_df$bam1_file_basename[row_idx] <- "APA58-T"
    } else if (sample1 == "APA58-T") {
      # Change T to N in sample name and BAM file basename
      corrected_df$sample1[row_idx] <- "APA58-N"
      corrected_df$bam1_file_basename[row_idx] <- "APA58-N"
    }
    
    # Handle sample2 swap (if it exists)
    if (!is.na(sample2)) {
      if (sample2 == "APA58-N") {
        # Change N to T in sample name and BAM file basename
        corrected_df$sample2[row_idx] <- "APA58-T"
        corrected_df$bam2_file_basename[row_idx] <- "APA58-T"
      } else if (sample2 == "APA58-T") {
        # Change T to N in sample name and BAM file basename
        corrected_df$sample2[row_idx] <- "APA58-N"
        corrected_df$bam2_file_basename[row_idx] <- "APA58-N"
      }
      
      # For APA58, ALWAYS ensure analysis is TvsN (not NvsT)
      if (str_detect(corrected_df$sample1[row_idx], "^APA58") && 
          str_detect(corrected_df$sample2[row_idx], "^APA58") &&
          corrected_df$analysis[row_idx] == "NvsT") {
        
        # Force analysis to be TvsN
        corrected_df$analysis[row_idx] <- "TvsN"
      }
    }
  }
  
  #----------------------------------------------------------------------------
  # Correction 3: iAPA5/iAPA6 Complex Sample Swap
  # Issue: There's a complex swap between iAPA5 and iAPA6 samples:
  # - iAPA5 group has: iAPA6-F, iAPA6-N, iAPA5-NF
  # - iAPA6 group has: iAPA5-F, iAPA5-N, iAPA6-NF
  # Fix: Assign the correct individual identifiers for each sample
  #----------------------------------------------------------------------------
  
  # Update individual identifiers for all iAPA5/iAPA6 sample combinations
  # First pass - identify all relevant rows
  iapa_rows <- which(
    str_detect(corrected_df$sample1, "^iAPA[56]") | 
    (!is.na(corrected_df$sample2) & str_detect(corrected_df$sample2, "^iAPA[56]"))
  )
  
  # Apply specific corrections based on the correct sample groupings
  for (row_idx in iapa_rows) {
    sample1 <- corrected_df$sample1[row_idx]
    sample2 <- corrected_df$sample2[row_idx]
    
    # Check if this is a tumor-only analysis (sample2 is NA)
    is_tumor_only <- is.na(sample2)
    
    # Handle tumor-only analyses
    if (is_tumor_only) {
      # For tumor-only (To) analyses, we need specific individual identifiers
      # that include the type to ensure unique output basenames
      if (sample1 == "iAPA6-F") {
        # iAPA6-F belongs to iAPA5 group
        corrected_df$individual1[row_idx] <- "iAPA5-F"
      }
      else if (sample1 == "iAPA6-N") {
        # iAPA6-N belongs to iAPA5 group
        corrected_df$individual1[row_idx] <- "iAPA5-N"
      }
      else if (sample1 == "iAPA5-NF") {
        # iAPA5-NF belongs to iAPA5 group
        corrected_df$individual1[row_idx] <- "iAPA5-NF"
      }
      else if (sample1 == "iAPA5-F") {
        # iAPA5-F belongs to iAPA6 group
        corrected_df$individual1[row_idx] <- "iAPA6-F"
      }
      else if (sample1 == "iAPA5-N") {
        # iAPA5-N belongs to iAPA6 group
        corrected_df$individual1[row_idx] <- "iAPA6-N"
      }
      else if (sample1 == "iAPA6-NF") {
        # iAPA6-NF belongs to iAPA6 group
        corrected_df$individual1[row_idx] <- "iAPA6-NF"
      }
    } 
    else {
      # Handle pairwise analyses
      
      # Process sample1
      if (sample1 == "iAPA6-F" || sample1 == "iAPA6-N" || sample1 == "iAPA5-NF") {
        # These samples belong to the iAPA5 group
        corrected_df$individual1[row_idx] <- "iAPA5"
      }
      else if (sample1 == "iAPA5-F" || sample1 == "iAPA5-N" || sample1 == "iAPA6-NF") {
        # These samples belong to the iAPA6 group
        corrected_df$individual1[row_idx] <- "iAPA6"
      }
      
      # Process sample2
      if (sample2 == "iAPA6-F" || sample2 == "iAPA6-N" || sample2 == "iAPA5-NF") {
        # These samples belong to the iAPA5 group
        corrected_df$individual2[row_idx] <- "iAPA5"
      }
      else if (sample2 == "iAPA5-F" || sample2 == "iAPA5-N" || sample2 == "iAPA6-NF") {
        # These samples belong to the iAPA6 group
        corrected_df$individual2[row_idx] <- "iAPA6"
      }
    }
  }
  
  #----------------------------------------------------------------------------
  # Correction 4: Standardize analysis names
  # Issue: Some samples use "TS" (Tumor Secondary) instead of "T", resulting in
  # analysis names like "TSvsN" instead of the standard "TvsN"
  # Fix: Change "TSvsN" to "TvsN" for consistency across all tumor/normal analyses
  #----------------------------------------------------------------------------
  
  # Identify rows with "TSvsN" analysis
  tsvs_rows <- which(corrected_df$analysis == "TSvsN")
  
  # Change "TSvsN" to "TvsN" for consistency
  if (length(tsvs_rows) > 0) {
    corrected_df$analysis[tsvs_rows] <- "TvsN"
  }
  
  # Clean up any temporary columns
  corrected_df <- corrected_df %>%
    select(-matches("^type[12]$"), -matches("^orig_"))
  
  return(corrected_df)
}

# Apply corrections
corrected_metadata <- apply_corrections(initial_metadata)

#------------------------------------------------------------------------------
# Phase 3: Apply filters to exclude problematic samples/individuals
#------------------------------------------------------------------------------

#' Apply filters to exclude problematic samples and individuals.
#'
#' @param metadata_df A dataframe containing sample metadata
#' @param samples_to_exclude Vector of sample IDs to exclude
#' @param individuals_to_exclude Vector of individual IDs to exclude
#' @return A filtered dataframe
apply_filters <- function(metadata_df, samples_to_exclude, individuals_to_exclude) {
  filtered_df <- metadata_df %>%
    # Filter out rows with excluded samples
    filter(!sample1 %in% samples_to_exclude) %>%
    filter(is.na(sample2) | !sample2 %in% samples_to_exclude) %>%
    # Filter out rows with excluded individuals
    filter(!individual1 %in% individuals_to_exclude)
  
  return(filtered_df)
}

# Apply filters
final_metadata <- apply_filters(corrected_metadata, samples_to_exclude, 
                               individuals_to_exclude)

#------------------------------------------------------------------------------
# Filename Conflict Detection
#------------------------------------------------------------------------------

# Add a temporary column with potential output basenames used by the Snakemake pipeline
final_metadata_with_basenames <- final_metadata %>%
  mutate(output_basename = str_c(individual1, "_", analysis))

# Count occurrences of each output basename
output_basename_counts <- final_metadata_with_basenames %>%
  count(output_basename) %>%
  arrange(desc(n))

# Find conflicts (basenames with more than one occurrence)
conflicting_basenames <- output_basename_counts %>%
  filter(n > 1)

# Check if conflicts exist and take appropriate action
if (nrow(conflicting_basenames) > 0) {
  # Extract the list of conflicting basenames
  conflict_list <- conflicting_basenames %>% 
    pull(output_basename)
  
  # Get the full metadata rows for conflicting basenames for debugging
  conflict_details <- final_metadata_with_basenames %>%
    filter(output_basename %in% conflict_list) %>%
    arrange(output_basename)
  
  # Generate error message
  cat("\n============== ERROR: FILENAME CONFLICTS DETECTED ==============\n")
  cat("The following output basenames would be duplicated in the Snakemake pipeline:\n\n")
  
  # Print each conflicting basename and count
  conflicting_basenames %>%
    mutate(conflict_message = str_c(output_basename, " (", n, " occurrences)")) %>%
    pull(conflict_message) %>%
    cat(sep = "\n")
  
  # Print detailed information about the conflicting rows
  cat("\n\nConflicting metadata rows:\n\n")
  print(conflict_details)
  
  # Stop script execution
  stop("Potential filename conflicts detected! Please review metadata corrections and filtering logic.")
} else {
  cat("\n============== INFO: No output filename conflicts detected ==============\n")
  cat("All output basenames ({individual1}_{analysis}) are unique in the final metadata table.\n")
}

#------------------------------------------------------------------------------
# Output Generation
#------------------------------------------------------------------------------

# Create a detailed status table for all samples
sample_inclusion_status <- sample_info %>%
  left_join(select(all_individuals, individual, exclude, reason), 
            by = c("individual")) %>%
  mutate(
    status = if_else(exclude, "exclude", "include"),
    issue_details = case_when(
      individual == "APA12" & type == "T" ~ "Tumor sample swapped with APA13-T",
      individual == "APA13" & type == "T" ~ "Tumor sample swapped with APA12-T",
      individual == "APA58" ~ "Tumor and normal samples swapped",
      sample == "APA61-TS" ~ "Sample sequenced twice (APA61-T is from normal FFPE slides)",
      TRUE ~ reason
    )
  ) %>%
  select(sample, individual, type, status, issue_details)

# Write the sample status table
write_tsv(sample_inclusion_status, "sample_inclusion_status.tsv")

# Ensure row names are not included in the output (removes automatic numeric indexing)
# Also ensure consistent row ordering to prevent display inconsistencies
final_metadata <- final_metadata %>%
  arrange(individual1, analysis, sample1) %>%
  as_tibble()  # Convert to tibble to ensure no row names

# Write the final result to a TSV file
write_tsv(final_metadata, "calling_metadata.tsv", na = "")

# Print a summary
cat("=============== Calling Metadata Summary ===============\n")
cat("Generated calling metadata file with", nrow(final_metadata), "sample combinations\n")
cat("Included", n_distinct(final_metadata$individual1), "unique individuals\n")
cat("Excluded samples:", str_c(samples_to_exclude, collapse = ", "), "\n")
cat("Excluded individuals:", str_c(individuals_to_exclude, collapse = ", "), "\n")
cat("Excluded analyses:", str_c(analyses_types_to_exclude, collapse = ", "), "\n")

# Count analysis types for reference
analysis_counts <- final_metadata %>%
  count(analysis) %>%
  arrange(desc(n))

cat("\nAnalysis type distribution:\n")
print(analysis_counts)

# Sample inclusion/exclusion summary
exclusion_summary <- sample_inclusion_status %>%
  count(status) %>%
  arrange(desc(n))

cat("\nSample inclusion/exclusion summary:\n")
print(exclusion_summary)
