#!/bin/bash
#SBATCH -c 1
#SBATCH --mem 20G
#SBATCH -p htc
#SBATCH -t 04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lsong36@asu.edu

INDEX=$1
SPECIES=$2

module load r-4.2.1-gcc-11.2.0
module load mamba/latest
source activate spatial
cd ~/comTCA

srun Rscript scripts/rasterize_expert_ranges.R -i $INDEX -s $SPECIES
