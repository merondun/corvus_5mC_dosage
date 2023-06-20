#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

RUN=$1

mkdir Bootstrapped
mkdir Bootstrapped/${RUN}
mkdir bootstrapped_results

egrep -i "Subanalysis|${RUN}" mCG-10x_LongForm_ATAC-21Dec21.txt > Bootstrapped/${RUN}/Input.txt
cp Bootstrap_FM.R Bootstrapped/${RUN}/Bootstrap_FM.R
cd Bootstrapped/${RUN}
Rscript Bootstrap_FM.R
cp Output.txt ../../bootstrapped_results/${RUN}.txt

