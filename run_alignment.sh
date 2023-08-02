#!/bin/bash
#
#SBATCH --job-name=sm_alignment_main_job
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=1200M
#SBATCH --output=slurm_logs/%x-%j.log

# First, point TMPDIR to the scratch in your home as mktemp will use thi
export TMPDIR=$HOME/scratch/tmp
# Second, create another unique temporary directory within this directory
export TMPDIR=$(mktemp -d)
# Finally, setup the cleanup trap
trap "rm -rf $TMPDIR" EXIT

mkdir -p slurm_logs
export SBATCH_DEFAULTS=" --output=slurm_logs/%x-%j.log"

date
srun snakemake -s alignment.smk --use-conda --profile=cubi-v1 -j10
date