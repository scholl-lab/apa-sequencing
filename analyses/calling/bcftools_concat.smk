import os
import functools
import csv

# ----------------------------------------------------------------------------------- #
# Load configuration file containing user-defined settings
configfile: "config.yaml"

# Define temporary directory using an environment variable
SCRATCH_DIR = os.environ.get('TMPDIR')

# Helper function to return memory based on the number of threads
def get_mem_from_threads(wildcards, threads):
    return threads * 4400  # Adjust as necessary
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Extract user-defined input and output directories from the configuration file
INPUT_DIR = config["scattered_vcf_folder"]
OUTPUT_DIR = config["output_folder"]

# Read metadata table into a dictionary from the TSV file
metadata_dict = {}
with open('calling_metadata.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        analysis_key = f"{row['individual1']}_{row['analysis']}"
        metadata_dict[analysis_key] = row
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Define result directories using functools.partial to join paths with the output folder
prefix_results = functools.partial(os.path.join, config['output_folder'])
CONCATENATED_VCF_DIR = prefix_results('concatenated_vcfs')
LOG_DIR = prefix_results('logs')

# Create output directories if they do not exist
os.makedirs(CONCATENATED_VCF_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Define the rules
rule all:
    input:
        expand(
            f"{CONCATENATED_VCF_DIR}/{{individual1}}_{{analysis}}.vcf.gz", 
            individual1=[row['individual1'] for row in metadata_dict.values()], 
            analysis=[row['analysis'] for row in metadata_dict.values()]
        )

rule concatenate_vcfs:
    input:
        lambda wildcards: expand(f"{INPUT_DIR}/{wildcards.individual1}_{wildcards.analysis}_chr{{chr}}.vcf.gz", chr=[str(i) for i in range(1, 22)] + ['X', 'Y'])
    output:
        concatenated_vcf = f"{CONCATENATED_VCF_DIR}/{{individual1}}_{{analysis}}.vcf.gz"
    threads: 2  # Adjust as necessary
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        "bcftools"
    log:
        f"{LOG_DIR}/{{individual1}}_{{analysis}}_concat.log"
    shell:
        """
        bcftools concat {input} -o {output.concatenated_vcf} -Oz 2> {log}
        """
# ----------------------------------------------------------------------------------- #