# Code for alignment and BAM reprocessing

## TODO
- [ ] update the run_alignment_v2.sh script to take a config file through the command line
- [ ] make separate config files for each step of the alignment process
- [ ] unify the 4 scripts needed for alignment into one script



## calculate all md5 checksums in input folder and subfolders
```bash
    snakemake -s md5sum_files.smk --profile=cubi-dev -j1
```

## run alignment for FASTQ files in input folder and subfolders
```bash
    sbatch run_alignment_v2.sh
```

## merge all the lane bam files into one sample bam file
```bash
    sbatch run_merge_bams_v2.sh
```

## deduplicate the merged bam files
```bash
    sbatch run_dedup_bams.sh
```

## run bqsr for the deduplicated bam files
```bash
    sbatch run_bqsr_bams.sh
```