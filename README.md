# apa-sequencing
Code for the aldosterone-producing adenoma sequencing project

## calculate all md5 checksums in input folder and subfolders
```bash
    snakemake -s md5sum_files.smk --profile=cubi-dev -j1
```

## run alignment for FASTQ files in input folder and subfolders
```bash
    sbatch run_alignment.sh
```