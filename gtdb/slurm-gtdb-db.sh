#!/bin/bash -login
#SBATCH -p bmm
#SBATCH -J sourmash-build-gtdb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=titus@idyll.org
#SBATCH -t 3-0:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=50gb

. "/home/ctbrown/miniconda3/etc/profile.d/conda.sh"

cd /home/ctbrown/sourmash_databases/gtdb

conda activate sgc

snakemake gtdb-release89-k{31,51}.sbt.json --configfile config-release89.yml

set -o nounset
set -o errexit
set -x

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

