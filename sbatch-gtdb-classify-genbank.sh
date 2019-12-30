#!/bin/bash -login
#SBATCH -p bmh
#SBATCH -J sourmash-gtdb-classify-genbank
#SBATCH --mail-type=ALL
#SBATCH --mail-user=titus@idyll.org
#SBATCH -t 1-0:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=100gb

. "/home/ctbrown/miniconda3/etc/profile.d/conda.sh"

cd /home/ctbrown/sourmash_databases

conda activate sgc

set -o nounset
set -o errexit
set -x

python ~/2019-sourmash-gtdb/bulk-classify-sbt-with-lca.py gtdb-scrub-k21-archaea gtdb/scrub-gtdb-release89-k21.lca.json.gz outputs/trees/scaled/genbank-archaea-d2-x1e5-k21.sbt.json
python ~/2019-sourmash-gtdb/bulk-classify-sbt-with-lca.py gtdb-scrub-k21-bacteria gtdb/scrub-gtdb-release89-k21.lca.json.gz outputs/trees/scaled/genbank-bacteria-d2-x1e5-k21.sbt.json

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

