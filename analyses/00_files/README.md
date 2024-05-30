# Instructions how to generate metadata for the samples in a project

## 1) Alignment metadata for the genomes project part

# TODO: add instructions for genomes project part


## 2) Alignment metadata for the exomes project part
```
# command to list all R1 FASTQ files in a folder
# remove _R1_001.fastq.gz from file names using sed
# remove the download/ folder name using sed
find download/ -type f -follow -print|xargs ls | grep "_R1_001.fastq.gz$" | sed 's/_R1_001.fastq.gz//g' | sed 's/download\///g' > fastq_files.txt

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
echo -e "fastq_files_basename\tsubfolder\tmdc_project\tproject_sample\tsample_sheet_number\tlane" | cat - fastq_files_with_parts_merged.txt > metadata.tsv

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

# TODO: add instructions for exomes project part
