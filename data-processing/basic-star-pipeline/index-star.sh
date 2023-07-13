#!/usr/bin/env bash

##### Call STAR generate Genome #####
STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir /scratch/alpine/dyammons@colostate.edu/proj05_k9_tumor/02_scripts/k9_ref_genome_STAR/ \
--genomeFastaFiles /projects/dyammons@colostate.edu/references/canine/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa \
--sjdbGTFfile /projects/dyammons@colostate.edu/references/canine/Canis_lupus_familiaris.CanFam3.1.104.gtf \
--outFileNamePrefix /scratch/alpine/dyammons@colostate.edu/proj05_k9_tumor/02_scripts/k9_ref_genome_STAR/ \
--sjdbOverhang 100 \
--genomeSAindexNbases 12

# for CanFam3.1 build
# Job ID: 2318419
# Cluster: alpine
# User/Group: dyammons@colostate.edu/dyammonspgrp@colostate.edu
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 16
# CPU Utilized: 01:21:00
# CPU Efficiency: 36.16% of 03:44:00 core-walltime Job Wall-clock time: 00:14:00 Memory Utilized: 26.59 GB Memory Efficiency: 44.31% of 60.00 GB
