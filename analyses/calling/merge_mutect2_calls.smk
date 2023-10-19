import os
import functools
import csv

# ----------------------------------------------------------------------------------- #
# Configuration and Environment Setup

# Load configuration file containing user-defined settings
configfile: "config.yaml"

# Define temporary directory using an environment variable (usually set by the cluster scheduler)
SCRATCH_DIR = os.environ.get('TMPDIR')

# Helper function to return memory based on the number of threads
def get_mem_from_threads(wildcards, threads):
    return threads * 4400
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Extract Input, Output Directories and Metadata

# Extract user-defined output directory from the configuration file
OUTPUT_DIR = config["output_folder"]

# Load metadata from the provided TSV file and store in a dictionary
metadata_dict = {}
with open('calling_metadata.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        analysis_key = f"{row['individual1']}_{row['analysis']}"
        metadata_dict[analysis_key] = row

# Define output directories for various results using functools
prefix_results = functools.partial(os.path.join, OUTPUT_DIR)
VARIANT_DIR = prefix_results('variant_calls')
MERGED_DIR = prefix_results('variant_merge')
LOG_DIR = prefix_results('logs')

# List of chromosomes for processing
chromosomes = [f"chr{i}" for i in range(1, 22)] + ["chrX", "chrY"]
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Snakemake Rules

# Main rule to define the expected final outputs
rule all:
    input:
        expand(
            [
                f"{MERGED_DIR}/{{individual1}}_{{analysis}}.vcf.gz",
                f"{MERGED_DIR}/{{individual1}}_{{analysis}}_merged.stats",
                f"{MERGED_DIR}/{{individual1}}_{{analysis}}_read-orientation-model.tar.gz"
            ],
            individual1=[row['individual1'] for row in metadata_dict.values()], 
            analysis=[row['analysis'] for row in metadata_dict.values()]
        )

# Rule to merge scattered VCF files across chromosomes
rule merge_vcfs:
    input:
        lambda wildcards: [f"{VARIANT_DIR}/{wildcards.individual1}_{wildcards.analysis}_{chrom}.vcf.gz" for chrom in chromosomes]
    output:
        f"{MERGED_DIR}/{{individual1}}_{{analysis}}.vcf.gz"
    params:
        input_files = lambda wildcards: ' '.join(['-I ' + item for item in [f"{VARIANT_DIR}/{wildcards.individual1}_{wildcards.analysis}_{chrom}.vcf.gz" for chrom in chromosomes]])
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        f"{LOG_DIR}/{{individual1}}_{{analysis}}.merge_vcfs.log"
    shell:
        """
        gatk --java-options '-Xms4000m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={resources.tmpdir}' GatherVcfs \
        {params.input_files} \
        -O {output} 2> {log}
        """

# Rule to merge stats files corresponding to VCFs
rule merge_stats:
    input:
        lambda wildcards: [f"{VARIANT_DIR}/{wildcards.individual1}_{wildcards.analysis}_{chrom}.vcf.gz.stats" for chrom in chromosomes]
    output:
        f"{MERGED_DIR}/{{individual1}}_{{analysis}}_merged.stats"
    params:
        input_files = lambda wildcards: ' '.join(['-stats ' + item for item in [f"{VARIANT_DIR}/{wildcards.individual1}_{wildcards.analysis}_{chrom}.vcf.gz.stats" for chrom in chromosomes]])
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        f"{LOG_DIR}/{{individual1}}_{{analysis}}.merge_stats.log"
    shell:
        """
        gatk --java-options '-Xms4000m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={resources.tmpdir}' MergeMutectStats \
        {params.input_files} \
        -O {output} 2> {log}
        """

# Rule to merge f1r2 orientation files across chromosomes
rule merge_f1r2:
    input:
        lambda wildcards: [f"{VARIANT_DIR}/{wildcards.individual1}_{wildcards.analysis}_{chrom}.f1r2.tar.gz" for chrom in chromosomes]
    output:
        f"{MERGED_DIR}/{{individual1}}_{{analysis}}_read-orientation-model.tar.gz"
    params:
        input_files = lambda wildcards: ' '.join(['-I ' + item for item in [f"{VARIANT_DIR}/{wildcards.individual1}_{wildcards.analysis}_{chrom}.f1r2.tar.gz" for chrom in chromosomes]])
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        f"{LOG_DIR}/{{individual1}}_{{analysis}}.merge_f1r2.log"
    shell:
        """
        gatk --java-options '-Xms4000m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={resources.tmpdir}' LearnReadOrientationModel \
        {params.input_files} \
        -O {output} 2> {log}
        """

# ----------------------------------------------------------------------------------- #
