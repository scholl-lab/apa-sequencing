import os

configfile: "config.yaml"

# Define a function to obtain a list of all .fastq.gz files in a folder and its subfolders
def get_input_files(folder):
    input_files = []
    for root, _, files in os.walk(folder):
        for file in files:
            if file.endswith('.fastq.gz'):
                input_files.append(os.path.relpath(os.path.join(root, file), folder))
    return input_files

# Obtain input and output folder paths from the config file
input_folder = config["input_folder"]
output_folder = config["output_folder"]
md5sum_folder = os.path.join(output_folder, "md5sum")
input_files = get_input_files(input_folder)

# Define the Snakemake rules
rule all:
    input:
        os.path.join(output_folder, "all_md5sums.txt")

rule calculate_md5sum:
    input:
        fastq=os.path.join(input_folder, "{file}")
    output:
        md5sum=os.path.join(md5sum_folder, "{file}.md5sum")
    shell:
        "mkdir -p $(dirname {output.md5sum}); "
        "md5sum {input.fastq} > {output.md5sum}"

rule join_md5sums:
    input:
        expand(os.path.join(md5sum_folder, "{file}.md5sum"), file=input_files)
    output:
        os.path.join(output_folder, "all_md5sums.txt")
    shell:
        "cat {input} > {output}"