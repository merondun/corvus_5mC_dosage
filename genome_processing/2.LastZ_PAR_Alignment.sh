#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=12:00:00

PAR=PAR.fa
PARG=PAR-GENES.fa
Wsyn=Wscaffs_syntenyOnly.fa
Flycatch=PAR_genes_.flycatcher.fa
NCC=NCC-ZW.fa

#A#
lastz $PAR $Wsyn \
--notransition --step=20 --nogapped --format=general > lastz/PAR-Wsyn.out

#B#
lastz $PAR $NCC \
--notransition --step=20 --nogapped --format=general > lastz/PAR-NCC.out

#C#
lastz $PAR $Flycatch \
--notransition --step=20 --nogapped --format=general > lastz/PAR-FLY.out

#D#
lastz $Flycatch[multiple] $NCC \
--notransition --step=20 --nogapped --format=general > lastz/FLY-NCC.out

#E#
lastz $PARG[multiple] $Wsyn \
--notransition --step=20 --nogapped --format=general > lastz/PARG-Wsyn.out

#F#
lastz $PARG[multiple] $Flycatch \
--notransition --step=20 --nogapped --format=general > lastz/PARG-FLY.out

#G#
lastz $PARG[multiple] $NCC \
--notransition --step=20 --nogapped --format=general > lastz/PARG-NCC.out
