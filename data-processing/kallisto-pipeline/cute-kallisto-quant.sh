#!/usr/bin/env bash

#SBATCH --job-name=kallisto_quant

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=02:00:00

#SBATCH --partition=amilan
#SBATCH --qos=normal

#SBATCH --mail-type=END
#SBATCH --mail-user=dyammons@colostate.edu
#SBATCH --output=kallisto_quant-%j.log


##### Call bash script #####

#ensure the node is clear of pre-loaded software
module purge

#excute the kallisto_quant_bulk.sh script
  #argument 1 == path to metadata in a .csv format
  #argument 2 == ntasks; will automatically populate from the above SLURM directives
bash kallisto_quant_bulk.sh metadata.csv $SLURM_NTASKS
