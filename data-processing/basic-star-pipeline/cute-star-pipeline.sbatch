#!/usr/bin/env bash

#SBATCH --job-name=star_blood

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=03:00:00

#SBATCH --partition=amilan
#SBATCH --qos=normal

#SBATCH --mail-type=END
#SBATCH --mail-user=dyammons@colostate.edu
#SBATCH --output=starTest-%j.log


##### Call bash script #####

#load in the conda env
module purge
source /curc/Asw/anaconda3/latest
conda activate dsci512

bash bulk_analysis.sh /scratch/alpine/dyammons@colostate.edu/proj02_k9_pbmc/01_input/BulkRNA/metadata.csv $SLURM_NTASKS

