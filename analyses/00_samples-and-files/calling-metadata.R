# List of samples
samples <- c("APA10-N", "APA10-T", "APA12-N", "APA12-T", "APA13-N", "APA13-T", "APA15-N", "APA15-T", "APA17-N", "APA17-T", "APA18-N", "APA18-T", "APA19-N", "APA19-T", "APA1-N", "APA1-T", "APA23-N", "APA23-TS", "APA25-N", "APA25-T", "APA29-N", "APA29-T", "APA2-T", "APA31-N", "APA31-T", "APA32-N", "APA32-TS", "APA34-N", "APA34-T", "APA36-N", "APA36-T", "APA38-N", "APA38-T", "APA3-N", "APA3-T", "APA41-N", "APA41-T", "APA43-N", "APA43-T", "APA44-N", "APA44-T", "APA45-N", "APA45-T", "APA47-N", "APA47-T", "APA48-N", "APA48-T", "APA49-N", "APA49-T", "APA51-T", "APA52-N", "APA52-T", "APA53-N", "APA53-T", "APA55-T", "APA56-N", "APA56-T", "APA58-N", "APA58-T", "APA59-N", "APA59-T", "APA61-N", "APA61-T", "APA61-TS", "APA69-N", "APA69-T", "APA72-N", "APA72-T", "APA74-N", "APA74-T", "APA75-N", "APA75-T", "APA78-N", "APA78-T", "APA79-N", "APA79-T", "APA81-N", "APA81-T", "APA86-N", "APA86-T", "APA87-N", "APA87-T", "iAPA1-F", "iAPA1-NF", "iAPA1-N", "iAPA2-F", "iAPA2-NF", "iAPA2-N", "iAPA3-F", "iAPA3-NF", "iAPA3-N", "iAPA5-F", "iAPA5-NF", "iAPA5-N", "iAPA6-F", "iAPA6-NF", "iAPA6-N")

# Function to extract the individual and type
extract_info <- function(sample) {
  individual <- sub("-.*", "", sample)
  type <- sub(".*-", "", sample)
  return(c(individual, type))
}

# Create combinations
combinations <- expand.grid(sample1 = samples, sample2 = samples)
combinations <- combinations[combinations$sample1 != combinations$sample2, ]

# Extract individual and type for each sample
combinations$individual1 <- sapply(combinations$sample1, function(x) extract_info(x)[1])
combinations$individual2 <- sapply(combinations$sample2, function(x) extract_info(x)[1])
combinations$type1 <- sapply(combinations$sample1, function(x) extract_info(x)[2])
combinations$type2 <- sapply(combinations$sample2, function(x) extract_info(x)[2])

# Filter combinations where individuals match and types do not match
combinations <- combinations[combinations$individual1 == combinations$individual2 & combinations$type1 != combinations$type2, ]

# Create bam file names
combinations$bam1_file_basename <- paste0(combinations$sample1)
combinations$bam2_file_basename <- paste0(combinations$sample2)

# Create analysis column
combinations$analysis <- paste0(combinations$type1, "vs", combinations$type2)

# Select and order columns
result <- combinations[, c("sample1", "sample2", "bam1_file_basename", "bam2_file_basename", "individual1", "individual2", "analysis")]

# Identify individuals with only one sample
individuals <- sapply(samples, function(x) extract_info(x)[1])
single_sample_individuals <- names(which(table(individuals) == 1))

# Create "Tumor only" analysis for individuals with only one sample
to_analysis <- data.frame(
  sample1 = samples[sapply(samples, function(x) extract_info(x)[1]) %in% single_sample_individuals],
  sample2 = NA,
  bam1_file_basename = samples[sapply(samples, function(x) extract_info(x)[1]) %in% single_sample_individuals],
  bam2_file_basename = NA,
  individual1 = single_sample_individuals,
  individual2 = NA,
  analysis = "To"
)

# Combine with the result
final_result <- rbind(result, to_analysis)

# Print final result
print(final_result)

# Optionally, write the final result to a CSV file
write.csv(final_result, "sample_combinations_with_to.csv", row.names = FALSE)
