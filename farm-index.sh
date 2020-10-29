#!/bin/bash -login
#SBATCH -p bmh
#SBATCH -J index-gtdb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ntpierce@ucdavis.edu
#SBATCH -t 3-0:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --chdir /home/ntpierce/sourmash_databases
#SBATCH --mem=10gb

# build sbts

# activate conda in general
. "/home/ntpierce/miniconda3/etc/profile.d/conda.sh"

conda activate sourmash_databases

set -o nounset
set -o errexit
set -x

snakemake -s index-from-csv.snakefile -j 10 --resources mem_mb=300000 --profile farm_utils/slurm_profile --cluster-config farm_utils/cluster_config.yml

#echo ${SLURM_JOB_NODELIST}       # Output Contents of the SLURM NODELIST

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

#sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
