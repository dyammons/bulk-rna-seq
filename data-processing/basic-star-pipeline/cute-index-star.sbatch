#!/usr/bin/env bash

#SBATCH --job-name=star_mkRef

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=2:00:00

#SBATCH --partition=amilan
#SBATCH --qos=normal

#SBATCH --mail-type=END
#SBATCH --mail-user=dyammons@colostate.edu
#SBATCH --output=star-%j.log


##### Call bash script #####

bash star_ref.sh
