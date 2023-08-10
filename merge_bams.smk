import glob
import re

# Read the configuration file
configfile: "config.yaml"

# Set a default input folder if not specified in the config file
aligned_folder = config.get("aligned_folder", "results/aligned/")

# Function to find all unique sample prefixes
def get_samples():
    samples = set()
    for file in glob.glob(aligned_folder + "*.bam"):
        # Extract just the filename without the directory
        filename = file.replace(aligned_folder, "")
        
        m = re.match(r"(.+)_DNA_(\d+)_(.+)_S\d+_L\d+_lane\d+\.bam", filename)
        if m:
            samples.add(f"{m.group(1)}_DNA_{m.group(2)}_{m.group(3)}")
    return list(samples)


# List of unique samples
samples = get_samples()

print(samples)

# Define the rule to collect all the targets
rule all:
    input:
        expand("results/merged/{sample}.merged.bam", sample=samples)

# Define the rule to merge BAM files
rule merge_bam_files:
    input:
        bam_files = lambda wildcards: glob.glob(aligned_folder + wildcards.sample + "_*.bam")
    output:
        bam = "results/merged/{sample}.merged.bam"
    params:
        list_file = "results/merged/{sample}.bamlist"
    shell:
        """
        # Write the BAM file names to merge into a list file
        echo "{input.bam_files}" | tr " " "\\n" > "{params.list_file}"
        
        # Create the merged directory if it does not exist
        mkdir -p results/merged
        
        # Use samtools to merge the BAM files
        samtools merge -@ 2 -O BAM -b "{params.list_file}" "{output.bam}"
        
        # Remove the list file
        rm "{params.list_file}"
        """