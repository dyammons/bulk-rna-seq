#!/usr/bin/env bash

################################################
# PROGRAM:
# kallisto_quant.sh
#
# DESCRIPTION: Completes pseudo-alignement of fastp processed data
#              The pipeline use kallisto cDNA annoations from a sepcies
#              of interest to generate TPM counts.
#
# AUTHOR:
# Dylan Ammons
#
# CREATION DATE:
# September 22, 2023
#
# LAST UPDATED DATE:
# NA
#
# DEPENDENCIES:
# 	Requires the installation of the follwing software: 
#		kallisto
#		See README.md on github for installation instructions.
#
#
# REQUIRES:
#    INPUT: .fastq files: Currently, this pseudo-alignemnt protocol is a supplemental script
#						  that is to used following a standard STAR run including fastp or 
#						  other QC/trimming program. Thus the input should be the QC filtered
#						  fastq files. For a seamless run you will need to ensure each sample
#						  is in its own directory with the fasta files having a samplename_trim_2.fastq
#						  namming structure. Using the star_align_bulk.sh script will format as required
#						  for this script.
#
#
#    INPUT: metadata file: A comma seperated (.csv) metadata file with three columns. 
#						   The first two columns are fastq file names. 
#						   The third column is a "nickname" of each sample.
#						   The metadata file can be housed anywhere, but
#						   you will have to point to it when you execute this script.
#
#
#    kallisto INDEXES: index files for the transcriptome. DO NOT generate the indexes using 
#					   the whole genome fasta file!!!!
#					   
#					   An index can be generated using kallisto index. 
#					   See README.md on github for futher instructions.
#
#
# USAGE (in a sbatch script):
# $ bash kallisto_quant.sh metadata.csv $SLURM_NTASKS
#
# USAGE (in a 1 core cmd line):
# $ bash kallisto_quant.sh metadata.csv 1
#
# OUTPUT:
#    All output will be located in the specified $outputdir path set in the MODIFY section
#
#    $outputdir/025_kallisto/: kallisto folder for each sample pseudo aligned
#    $outputdir/025_kallisto/$samplename: two files will be in these directories. They contain
#                                         the output from kallisto: abundance.tsv and run_info.json
#                                         abundance.tsv contains the counts that will be needed for
#                                         downstream analysis, while run_info.json provides details 
#                                         of the run. See kallisto manual for more information.
#
# KNOWN BUGS:
#    >>> none???
#    
# THINGS TO IMPROVE:
#    >>> add in options for alterative input and full analysis pipeline
################################################


########## MODIFY ###############

#Path to QC filtered/trimmed samples (can be relative or absolute):
inputdir="../03_output/2023-09-19_output/01_fastp/"

#this variable will populate using the 1st argument in the command used to exceute this script
#if desired (this is not reccomnded) you can place path to metadata here (again, not reccomnded, you will have to change code outside of the modify section if you choose to do this)
metadata=$1

#path to kallisto index:
kallIndex="/projects/dyammons@colostate.edu/references/canine/canFam3.1.idx"


#path to the main output_directory, subdirectories will be generated during a run:
# DATE=`date +%Y-%m-%d`
#OR
DATE='2023-09-19'
outputdir="../03_output/"$DATE"_output/"

########## END MODIFY ###############


########## BEGIN CODE ###############

echo -e ">>> INITIATING analyzer with command:\n\t$0 $@"

#Number of threads to use:
pthread=$2

# Make output directories
echo -e ">>> MAKING output directory"
echo -e "\tmkdir $outputdir"
mkdir -p $outputdir

####### META DATA #############

#These are the sample names, R1:
samples1=( $(cut -f1 -d',' --output-delimiter=' ' $metadata) )

#These are the sample names, R2:
samples2=( $(cut -f2 -d',' --output-delimiter=' ' $metadata) )

#These are the nicknames I want to give the files:
names=( $(cut -f3 -d',' --output-delimiter=' ' $metadata) )


####### PIPELINE ##############

# Report back to the user which files will be processed and which names they'll be given:
echo -e ">>> INPUT: This script will process files from the metafile:\n\t$metadata"
echo -e ">>> PLAN: This script will process the sample files into the following names: "
echo -e "\tSAMPLE1\tSAMPLE2\tNAMES"

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
    echo -e "\t${samples1[$counter]}\t${samples2[$counter]}\t${names[$counter]}"
done


echo -e "\n>>> kallisto: pseudo aligning each sample to the transcriptome"
outKall=$outputdir"025_kallisto/"
mkdir -p $outKall

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do

    samplename=${names[$counter]}

    mkdir -p $outKall/$samplename

    cmd25="/projects/$USER/software/kallisto/build/src/kallisto quant -i $kallIndex \
    -o $outKall/$samplename \
    $inputdir$samplename/$samplename"_trim_1.fastq" \
    $inputdir$samplename/$samplename"_trim_2.fastq" \
    -t $pthread"

    echo $cmd25
    echo -e "\t$ ${cmd25}"
    time eval $cmd25

done


######## VERSIONS #############

echo -e "\n>>> VERSIONS:"
echo -e "\n>>> kallisto VERSION:"
/projects/$USER/software/kallisto/build/src/kallisto version
