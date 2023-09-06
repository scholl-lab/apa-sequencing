import os
import functools
import csv

# ----------------------------------------------------------------------------------- #
# Load configuration file containing user-defined settings
configfile: "config.yaml"

# Define temporary directory using an environment variable (usually set by the cluster scheduler)
SCRATCH_DIR = os.environ.get('TMPDIR')

# Helper function to return memory based on the number of threads
def get_mem_from_threads(wildcards, threads):
    return threads * 4400
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Extract user-defined input and output directories and reference file from the configuration file
INPUT_DIR = config["final_bam_folder"]
OUTPUT_DIR = config["output_folder"]
REFERENCE_FILE = config["reference"]
PANEL_OF_NORMALS = config["panel_of_normals"]
AF_ONLY_GNOMAD = config["af_only_gnomad"]

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
VARIANT_DIR = prefix_results('variant_calls')
LOG_DIR = prefix_results('logs')
# ----------------------------------------------------------------------------------- #

# List of chromosomes to loop through
chromosomes = [f"chr{i}" for i in range(1, 22)] + ["chrX", "chrY"]

# ----------------------------------------------------------------------------------- #
# Rule "all": Defines the final output of the pipeline
rule all:
    input:
        expand(
            f"{VARIANT_DIR}/{{individual1}}_{{analysis}}_{{chromosome}}.vcf.gz", 
            individual1={row['individual1'] for row in metadata_dict.values()}, 
            analysis={row['analysis'] for row in metadata_dict.values()}, 
            chromosome=chromosomes
        )

# ----------------------------------------------------------------------------------- #
# Rule "mutect2_call": Performs variant calling using GATK's Mutect2
rule mutect2_call:
    input:
        bam1 = lambda wildcards: "{}/{}.merged.dedup.bqsr.bam".format(INPUT_DIR, metadata_dict["{}_{}".format(wildcards.individual1, wildcards.analysis)]['bam1_file_basename']),
        bam2 = lambda wildcards: "{}/{}.merged.dedup.bqsr.bam".format(INPUT_DIR, metadata_dict["{}_{}".format(wildcards.individual1, wildcards.analysis)]['bam2_file_basename'])
    output:
        variant_file = f"{VARIANT_DIR}/{{individual1}}_{{analysis}}_{{chromosome}}.vcf.gz"
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,      # Memory in MB based on the number of threads
        time = '72:00:00',                  # Time limit for the job
        tmpdir = SCRATCH_DIR                # Temporary directory
    conda:
        "gatk"  # This sets the Conda environment for this rule
    params:
        reference = REFERENCE_FILE,
        normal_sample = lambda wildcards: metadata_dict["{}_{}".format(wildcards.individual1, wildcards.analysis)]['sample2'],
        individual = lambda wildcards: metadata_dict["{}_{}".format(wildcards.individual1, wildcards.analysis)]['individual1'],
        analysis = lambda wildcards: metadata_dict["{}_{}".format(wildcards.individual1, wildcards.analysis)]['analysis'],
        panel_of_normals = PANEL_OF_NORMALS,
        af_only_gnomad = AF_ONLY_GNOMAD
    log:
        mutect2 = f"{LOG_DIR}/{{individual1}}_{{analysis}}_{{chromosome}}.mutect2.log"
    shell:
        """
        gatk --java-options '-Xms4000m -Xmx10g -Djava.io.tmpdir={resources.tmpdir}' Mutect2 \
            -R {params.reference} \
            -I {input.bam1} \
            -I {input.bam2} \
            -normal {params.normal_sample} \
            --germline-resource {params.af_only_gnomad} \
            --panel-of-normals {params.panel_of_normals} \
            --f1r2-tar-gz {params.individual}_{params.analysis}_{wildcards.chromosome}.f1r2.tar.gz \
            -L {wildcards.chromosome} \
            -O {output.variant_file} > {log.mutect2}
        """
# ----------------------------------------------------------------------------------- #
