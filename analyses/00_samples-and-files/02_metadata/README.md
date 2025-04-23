# Instructions how to generate metadata for the samples in a project

# TODO tickets
- [ ] add instructions for genomes project part
- [ ] add instructions for exomes project part
- [ ] make the commands into a script that works with different file naming conventions

## 1) Alignment metadata for the genomes project part

# TODO: add instructions for genomes project part


## 2) Alignment metadata for the exomes project part
```
# command to list all R1 FASTQ files in a folder
# remove _R1_001.fastq.gz from file names using sed
# remove the download/ folder name using sed
find download/exomes -type f -follow -print|xargs ls | grep "_R1_001.fastq.gz$" | sed 's/_R1_001.fastq.gz//g' | sed 's/download\/exomes\///g' > fastq_files.txt

# cat the file and pipe to:
# split pipe input by "/" or "_" using awk and print all fields
cat fastq_files.txt | awk -F"[/_]" -v OFS="\t" '{for(i=1;i<=NF;i++) printf $i"\t"; print ""}' > fastq_parts.txt

# add the columns from fastq_parts.txt to the fastq_files.txt
# using paste
# and remove the subfolder column using sed
paste fastq_files.txt fastq_parts.txt | sed 's/^.*\///g' > fastq_files_with_parts.txt

# change potential splits that resulted in different number of columns using sed
# then merge the columns to generate the columns project_sample, sample_sheet_number and lane using awk
# ! Note: this is specific for each project
cat fastq_files_with_parts.txt | sed 's/T\tS\t/TS\t/g' | awk -F"\t" -v OFS="\t" '{print $1, $2, $3, $4"-"$5, $6, $7}' > fastq_files_with_parts_merged.txt

# add header using echo
echo -e "fastq_files_basename\tsubfolder\tmdc_project\tproject_sample\tsample_sheet_number\tlane" | cat - fastq_files_with_parts_merged.txt > metadata_exomes.tsv

# remove intermediate files
rm fastq_files.txt fastq_parts.txt fastq_files_with_parts.txt fastq_files_with_parts_merged.txt
```


## 3) Calling metadata for the genomes project part
```
# find all files in a directory and its subdirectories
# remove the subfolder path from the file names
# remove the file extension preserving all file names
find results/bqsr/ -type f -name "*.merged.dedup.bqsr.bam" | sed 's/results\/bqsr\///g' | sed 's/.merged.dedup.bqsr.bam//g' > final_bams.txt

# generate a list of all possible combinations of lines in a file
awk -v OFS="\t" 'NR==FNR { a[$0]; next } { for (i in a) print i, $0 }' final_bams.txt final_bams.txt > final_bams_combinations.txt

# use awk to generate the sample names for the tumor and normal samples
# by removing everything up to the last underscore ("_") from the first and second columns
cat final_bams_combinations.txt | awk -F"[\t_]" -v OFS="\t" '{print $4, $8}' > sample_names_combinations.txt

# use awk to generate the analysis from sample type
# by removing everything up to the minus ("-") from the first and second columns of sample_names_combinations.txt
# keep both sample names (columns 1 and 3)
# for columns 2 and 4 add "vs" in between
cat sample_names_combinations.txt | awk -F"[\t-]" -v OFS="\t" '{print $1, $3, $2"vs"$4}' > analysis_combinations.txt

# combine the columns from sample_names_combinations.txt, analysis_combinations.txt and final_bams_combinations.txt
# using paste
# filter to remove lines containing "NvsN", "NvsF", "NvsFN" "FvsF", "NFvsNF"
# using grep -v
# arrange the columns in the order sample_name, analysis, final_bam using sort
# filter to remove lines where the sample name is the same as the normal sample name
# using awk
paste sample_names_combinations.txt final_bams_combinations.txt analysis_combinations.txt | grep -vP "\tNvsN" | grep -vP "\tNvsF" | grep -vP "\tNvsFN" | grep -vP "\tFvsF" | grep -vP "\tNFvsNF" | sort | awk -F"\t" -v OFS="\t" '{if ($5 == $6) print $0}' > final_bams_combinations_merged.txt

# add header using echo
echo -e "sample1\tsample2\tbam1_file_basename\tbam2_file_basename\tindividual1\tindividual2\tanalysis" | cat - final_bams_combinations_merged.txt > calling_metadata.tsv

# remove intermediate files
rm final_bams.txt final_bams_combinations.txt sample_names_combinations.txt analysis_combinations.txt final_bams_combinations_merged.txt
```


## 4) Calling metadata for the exomes project part

The exomes project calling metadata is generated using an R script (`calling-metadata.R`) that creates all valid sample combinations for variant calling while handling several known sample issues.

### Setup and Usage

```bash
# Navigate to the appropriate directory
cd analyses/00_samples-and-files/

# Run the R script to generate the calling_metadata.tsv file
Rscript calling-metadata.R
```

### Script Details

The script follows a structured approach to generate accurate calling metadata:

1. **Initial Setup & Sample Status Tables**: 
   - Defines all samples and creates sample inclusion/exclusion tables
   - Documents reasons for excluding specific samples
   - Generates both a detailed sample status file (`sample_inclusion_status.tsv`) and the final calling metadata

2. **Generate Valid Sample Combinations**:
   - Creates pairwise combinations for each individual (same individual, different sample types)
   - Generates tumor-only ("To") analyses for all non-normal samples
   - Excludes unwanted analysis combinations (e.g., Normal-vs-Tumor, Normal-vs-Normal)

3. **Apply Sample-Specific Corrections**:
   - Corrects known sample identity issues
   - Maps sample names to the correct BAM file basenames

4. **Filter Out Problematic Samples**:
   - Excludes samples that failed quality control
   - Handles duplicate samples

5. **Output Generation**:
   - Writes detailed sample status information to `sample_inclusion_status.tsv`
   - Writes the final calling metadata to `calling_metadata.tsv`
   - Provides a summary of included/excluded samples and analyses

### Handling Known Sample Issues

The script addresses several known sample identity issues:

#### Sample Swaps

1. **APA12/APA13 Tumor Swap**
   - **Issue**: The tumor samples for patients APA12 and APA13 were swapped during processing.
   - **Solution**: The script updates sample IDs and BAM file basenames while maintaining the original individual identifiers to ensure correct variant calling and output naming.

2. **APA58 Tumor/Normal Swap**
   - **Issue**: The tumor and normal samples for patient APA58 were swapped.
   - **Solution**: For any analysis involving APA58 samples, the script swaps the sample IDs and BAM file basenames, and ensures all analyses are forced to "TvsN" format.

3. **iAPA5/iAPA6 Complex Sample Swap**
   - **Issue**: There's a complex reciprocal swap between iAPA5 and iAPA6 samples:
     - Sample group "iAPA5" should contain: iAPA6-F, iAPA6-N, iAPA5-NF
     - Sample group "iAPA6" should contain: iAPA5-F, iAPA5-N, iAPA6-NF
   - **Solution**: The script handles these swaps by:
     - For pairwise analyses: Reassigning individual identifiers based on which sample group they truly belong to
     - For tumor-only analyses: Maintaining the tissue type in individual identifiers for unique output names

#### Sample Exclusions

The following samples are excluded from analyses:

1. **APA2** 
   - **Reason**: MAPM (Micronodular adrenocortical hyperplasia with primary pigmentation) with unselective dissection
   - **Solution**: All samples from this individual are excluded

2. **APA3** 
   - **Reason**: Other sample contamination in tumor
   - **Solution**: All samples from this individual are excluded

3. **APA32** 
   - **Reason**: Hyperplasia - No clear tumor, micronoduli in normal tissue
   - **Solution**: All samples from this individual are excluded

4. **APA38** 
   - **Reason**: Sample switch with unknown sample
   - **Solution**: All samples from this individual are excluded

5. **APA45** 
   - **Reason**: Duplicate of iAPA5 (same individual)
   - **Solution**: All APA45 samples are excluded, analyses should use iAPA5 instead

6. **APA56** 
   - **Reason**: High contamination with other sample
   - **Solution**: All samples from this individual are excluded

7. **APA69** 
   - **Reason**: MAPM with unselective dissection
   - **Solution**: All samples from this individual are excluded

8. **APA61-TS**
   - **Reason**: Sample sequenced twice (APA61-T is from normal FFPE slides and APA61-TS is from scooped FFPE block)
   - **Solution**: Only APA61-TS is excluded; APA61-T is kept for analyses

#### Special Handling for Functional and Non-Functional Samples

For samples with Functional (F) and Non-Functional (NF) tissue types:

- **Issue**: In tumor-only ("To") analyses, multiple samples from the same individual (e.g., iAPA1-F and iAPA1-NF) would generate identical output basenames in the Snakemake pipeline (e.g., "iAPA1_To").
  
- **Solution**: For F and NF samples, the individual identifier includes the tissue type in tumor-only analyses:
  - Example: For sample "iAPA1-F", the individual identifier becomes "iAPA1-F" instead of just "iAPA1"
  - Example: For sample "iAPA1-NF", the individual identifier becomes "iAPA1-NF" instead of just "iAPA1"

- **Benefit**: This ensures unique output filenames in the Snakemake pipeline and prevents files from being overwritten.

#### Excluded Analysis Types

The following analysis types are excluded from the metadata:

1. **Normal-vs-Normal** ("NvsN"): Comparing identical tissue types from the same individual would yield no meaningful somatic variants.
2. **Normal-vs-Functional** ("NvsF"): Not relevant for the current analysis objectives.
3. **Normal-vs-Non-Functional** ("NvsNF"): Not relevant for the current analysis objectives.
4. **Functional-vs-Functional** ("FvsF"): Comparing identical tissue types would yield no meaningful variants.
5. **Non-Functional-vs-Non-Functional** ("NFvsNF"): Comparing identical tissue types would yield no meaningful variants.
6. **Normal-vs-Tumor** ("NvsT"): While standard tumor/normal comparisons are valid, for this specific project we focus on tumor-only analyses after extensive quality control.
7. **Normal-vs-Tumor-Secondary** ("NvsTS"): Similar to Normal-vs-Tumor, excluded to maintain focus on tumor-only analyses for this project.

### Sample Type Legend

The samples follow a naming convention of `{individual}-{type}` where type can be:
- **N**: Normal tissue
- **T**: Primary tumor
- **TS**: Secondary tumor
- **F**: Functional tissue
- **NF**: Non-functional tissue

### Output Files

The script generates two TSV files:

1. **sample_inclusion_status.tsv**: A comprehensive table containing:
   - Sample ID
   - Individual ID
   - Sample type
   - Inclusion/exclusion status
   - Detailed reason for exclusion or issue details

2. **calling_metadata.tsv**: The final metadata file for variant calling containing:
   - sample1, sample2: The sample identifiers
   - bam1_file_basename, bam2_file_basename: The actual BAM filenames to use (may differ from sample names due to corrections)
   - individual1, individual2: The individual identifiers
   - analysis: The type of analysis (e.g., "TvsF", "To")
