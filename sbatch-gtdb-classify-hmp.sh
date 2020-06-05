#!/bin/bash -login
#SBATCH -p bmm
#SBATCH -J sourmash-gtdb-classify-hmp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=titus@idyll.org
#SBATCH -t 3-0:00:00
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

python ~/2019-sourmash-gtdb/bulk-classify-sbt-with-lca.py gtdb-scrub-k31-almeida gtdb/build/gtdb-release89-k31-lowrank.lca.json.gz HMP_mags/almeida-mags-k31.sbt.json
python ~/2019-sourmash-gtdb/bulk-classify-sbt-with-lca.py gtdb-scrub-k31-nayfach gtdb/build/gtdb-release89-k31-lowrank.lca.json.gz HMP_mags/nayfach-k31.sbt.json
python ~/2019-sourmash-gtdb/bulk-classify-sbt-with-lca.py gtdb-scrub-k31-pasolli gtdb/build/gtdb-release89-k31-lowrank.lca.json.gz HMP_mags/pasolli-mags-k31.sbt.json

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits
