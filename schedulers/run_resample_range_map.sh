#!/bin/bash
#SBATCH -c 1
#SBATCH --mem 50G
#SBATCH -p htc
#SBATCH -t 04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lsong36@asu.edu

FNAME=$1
FNAMETO=$2

module load r-4.2.1-gcc-11.2.0
module load mamba/latest
source activate spatial
cd ~/comTCA

gdalwarp -te 29.1796874960035 -11.8673718726578 40.429702498976 -0.878871782654164 -tr 0.000898348239477162 0.000898340425932279 $FNAME $FNAMETO