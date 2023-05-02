# apa-sequencing
Code for the aldosterone-producing adenoma sequencing project

## calculate all md5 checksums in input folder and subfolders
<code>snakemake -s md5sum_files.smk --profile=cubi-dev -j1<code>