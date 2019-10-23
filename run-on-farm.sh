#!/bin/bash -login
#SBATCH -p bmh
#SBATCH -J sgc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=titus@idyll.org
#SBATCH -t 1-0:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=130gb

. "/home/ctbrown/miniconda3/etc/profile.d/conda.sh"

cd /home/ctbrown/sourmash_databases

conda activate sgc

set -o nounset
set -o errexit
set -x

snakemake -p -j 32

#echo ${SLURM_JOB_NODELIST}       # Output Contents of the SLURM NODELIST

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

#sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
