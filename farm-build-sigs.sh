#!/bin/bash -login
#SBATCH -p bmm
#SBATCH -J sigs-gtdb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ntpierce@ucdavis.edu
#SBATCH --chdir /group/ctbrowngrp/sourmash_databases
#SBATCH -t 3-0:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=10gb

# generate sigs

# activate conda in general
. "/home/ntpierce/miniconda3/etc/profile.d/conda.sh"

conda activate sourmash_databases

set -o nounset
set -o errexit
set -x

snakemake -s index-from-csv.snakefile -j 32 --until signames_to_file --profile farm_utils/default_profile

#echo ${SLURM_JOB_NODELIST}       # Output Contents of the SLURM NODELIST

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

#sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
