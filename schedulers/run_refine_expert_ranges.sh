#!/bin/bash
#SBATCH -c 1
#SBATCH --mem 120G
#SBATCH -p htc
#SBATCH -t 02:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lsong36@asu.edu

SPECIES=$1
TAXON=$2
RANGE=$3

module load r-4.2.1-gcc-11.2.0
module load mamba/latest
source activate spatial
cd ~/comTCA

srun Rscript scripts/refine_expert_ranges.R -s "$SPECIES" -t $TAXON -r $RANGE
