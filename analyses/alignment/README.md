# Code for alignment and bam reprocessing

## calculate all md5 checksums in input folder and subfolders
```bash
    snakemake -s md5sum_files.smk --profile=cubi-dev -j1
```

## run alignment for FASTQ files in input folder and subfolders
```bash
    sbatch run_alignment.sh
```

## merge all the lane bam files into one sample bam file
```bash
    sbatch merge_bams.sh
```

## deduplicate the merged bam files
```bash
    sbatch run_dedup_bams.sh
```