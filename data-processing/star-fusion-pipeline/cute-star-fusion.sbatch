#!/usr/bin/env bash

#SBATCH --job-name=star-fusion

#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=04:00:00

#SBATCH --partition=amilan
#SBATCH --qos=normal

#SBATCH --mail-type=END
#SBATCH --mail-user=dyammons@colostate.edu
#SBATCH --output=star-%j.log


##### Call bash script #####


singularity exec -e -H /scratch/alpine/dyammons@colostate.edu/star-fusion/ /projects/$USER/software/star-fusion/star-fusion.v1.12.0.simg /projects/$USER/software/star-fusion/ctat-genome-lib-builder/prep_genome_lib.pl \
--CPU $SLURM_NTASKS \
--genome_fa /projects/dyammons@colostate.edu/references/canine/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa \
--gtf /projects/dyammons@colostate.edu/references/canine/Canis_lupus_familiaris.CanFam3.1.104.gtf \
--pfam_db current \
--dfam_db human
