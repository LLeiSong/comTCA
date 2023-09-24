#!/bin/bash
#SBATCH -c 1
#SBATCH --mem 120G
#SBATCH -p htc
#SBATCH -t 04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lsong36@asu.edu

FNAME=$1

module load mamba/latest
source activate spatial
cd ~/comTCA

python scripts/calc_tz_range_size_one.py --fname "$FNAME"