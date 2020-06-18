#!/bin/bash -login
#SBATCH -p bmm
#SBATCH -J sourmash-build-gtdb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=titus@idyll.org
#SBATCH -t 1-0:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=40gb

. "/home/ctbrown/miniconda3/etc/profile.d/conda.sh"

cd /home/ctbrown/sourmash_databases/gtdb

conda activate sgc

#snakemake --configfile config-release89.yml -p build/gtdb-release89-k31-2000.lca.json.gz
snakemake dayhoff-test-gtdb-release89-k57-200.lca.json.gz -p  -j 32

set -o nounset
set -o errexit
set -x

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

