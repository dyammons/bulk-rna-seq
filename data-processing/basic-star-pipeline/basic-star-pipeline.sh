#!/usr/bin/env bash

################################################
# PROGRAM:
# bulk_analysis.sh
#
# DESCRIPTION:
#
# AUTHOR:
# Dylan Ammons
#
# CREATION DATE:
# July, 7, 2023
#
# LAST UPDATED DATE:
# NA
#
# DEPENDENCIES:
# 	Requires the installation of the follwing software: 
#		fastp
#		STAR
#		featureCounts
#		samtools
#		deep-tools
#
#
# REQUIRES:
#    INPUT: .fastq files.    For each sample, unzipped paired forward and reverse sequencing files
#								are required. These should be placed in an input
#								directory.
#
#    INPUT: _metadata.txt file: A metadata file with two columns. The first two columns
#								are fastq file names. The third column is a "nickname"
#								of each sample. Later columns can be included with other
#								metadata information. Metadata file should be placed
#								within the inputdir directory.
#
#
#    STAR INDEXES: index files for the genome. These are produced using STAR genomeGenerate.
#
#    GENOME SEQUENCE: .fa  or .tar.gz file for the genome. This is the sequence of the 
#                                genome.
#
#    GENOME ANNOTATION: .gtf file for the genome. This is a genome annotation file of gene
#								features. Version and coordinates must match the genome
#								sequence (.fa above).
#
# USAGE:
# $ bash bulk_analysis.sh meta.txt 16
# use awk 'NR==1 || $7>0' counts.txt | wc -l to check number of features
#
# OUTPUT:
#
# KNOWN BUGS:
#
# THINGS TO IMPROVE:
#
################################################


####### MODIFY THIS SECTION #############

#The input samples live in directory:
inputdir="../01_input/BulkRNA/"

#Metadata file. This pulls the metadata path and file from the command line
#create metadata file 

ls $inputdir | grep ".fq" | cut -d"_" -f1 | sort -u > samples.tmp
while read line
do
     echo $line"_1.fq,"$line"_2.fq,"$line >> metadata.csv
done < samples.tmp

metadata=$1

#This is where the star index files live:
starPath="/scratch/alpine/dyammons@colostate.edu/proj05_k9_tumor/02_scripts/k9_ref_genome_STAR"

#This is where the genome sequence lives:
genomefa="/projects/dyammons@colostate.edu/references/canine/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa"

#This is where the gtf file lives:
gtffile="/projects/dyammons@colostate.edu/references/canine/Canis_lupus_familiaris.CanFam3.1.104.gtf"

#This is the output_directory:
DATE=`date +%Y-%m-%d`
#OR
# DATE='2023-07-12'
outputdir="../03_output/"$DATE"_output/"


########## DONE MODIFYING ###############



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


# FASTP to remove unwanted sequences
# FASTP to determine quality
echo -e "\n>>> FASTP: Trimming excess and low-quality sequences from .fastq file; generating quality report"
mkdir -p $outputdir"01_fastp"

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
    samplename=${names[$counter]}
    sample1=${samples1[$counter]}
    sample2=${samples2[$counter]}

    ## Echo statements
    
            ##### ENTER ECHO STATEMENTS HERE #####
    
    ## Make output directories
    mkdir -p $outputdir"01_fastp/"$samplename
    
    ## execute fastp
    cmd1="fastp -i $inputdir/$sample1 \
-I $inputdir/$sample2 \
-o ${outputdir}01_fastp/${samplename}/${samplename}_trim_1.fastq \
-O ${outputdir}01_fastp/${samplename}/${samplename}_trim_2.fastq \
-h ${outputdir}01_fastp/${samplename}/${samplename}_report.html \
-j ${outputdir}01_fastp/${samplename}/${samplename}_report.json \
--detect_adapter_for_pe \
--thread $pthread \
-x -g "

     echo -e "\t$ ${cmd1}"
     time eval $cmd1

done

# STAR to align to the genome
echo -e "\n>>> STAR: aligning each sample to the genome"
outSTAR=$outputdir"02_star/"
mkdir -p $outSTAR

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
    samplename=${names[$counter]}
    sample1=${samples1[$counter]}
    sample2=${samples2[$counter]}


    ## execute star
    cmd2="STAR --genomeDir $starPath --readFilesIn $outputdir"01_fastp/"$samplename/$samplename"_trim_1.fastq" $outputdir"01_fastp/"$samplename/$samplename"_trim_2.fastq" --runThreadN $pthread --outFileNamePrefix $outSTAR/$samplename"

     echo -e "\t$ ${cmd2}"
     time eval $cmd2

done


# FEATURECOUNTS to tabulate reads per gene:
echo -e "\n>>> FEATURECOUNTS: Run featureCounts on all files to tabulate read counts per gene"
outfeature=$outputdir"03_feature/"
mkdir -p $outfeature

# Acquire a list of .sam names
samfilePath=()
for (( counter=0; counter < ${#names[@]}; counter++ ))
do
    samfile=${names[$counter]}Aligned.out.sam
    samfilePath+=(${outSTAR}${samfile})

done


# Execute featureCounts
cmd3="featureCounts -p -T $pthread -a $gtffile -o ${outfeature}counts.txt ${samfilePath[*]}"
echo -e "\t$ $cmd3"
time eval $cmd3



# SAMTOOLS and BAMCOVERAGE: to convert .sam output to uploadable .bam and .wg files
echo -e "\n>>> SAMTOOLS/BAMCOVERAGE: to convert files to uploadable _sort.bam and _sort.bam.bai files:"
samout=$outputdir"04_samtools/"
mkdir -p $samout

for seqname in ${names[@]}
do
    # echo
    echo -e "\tSamtools and BamCoverage convert: ${seqname}"
    
    # Samtools: compress .sam -> .bam
    cmd4="samtools view --threads $pthread -bS ${outSTAR}${seqname}Aligned.out.sam > ${samout}${seqname}.bam"
    echo -e "\t$ ${cmd4}"
    time eval $cmd4

    
    # Samtools: sort .bam -> _sort.bam
    cmd5="samtools sort --threads $pthread -o ${samout}${seqname}_sort.bam --reference $genomefa ${samout}${seqname}.bam"
    echo -e "\t$ ${cmd5}"
    time eval $cmd5
    
    
    # Samtools: index _sort.bam -> _sort.bam.bai
    cmd6="samtools index ${samout}${seqname}_sort.bam"
    echo -e "\t$ ${cmd6}"
    time eval $cmd6
    
    
    # bamCoverage: Create a .bw file that is normalized. This can be uploaded to IGV or UCSC
    cmd7="bamCoverage -b ${samout}${seqname}_sort.bam -o ${samout}${seqname}_sort.bw --outFileFormat bigwig -p $pthread --normalizeUsing CPM --binSize 1"
    echo -e "\t$ ${cmd7}"
    time eval $cmd7
    
done

# Clean up
rm *.tmp

######## VERSIONS #############
echo -e "\n>>> VERSIONS:"
echo -e "\n>>> FASTP VERSION:"
$fastp --version
echo -e "\n>>> STAR VERSION:"
$STAR --version
echo -e "\n>>> SAMTOOLS VERSION:"
$samtools --version
echo -e "\n>>> FEATURECOUNTS VERSION:"
$featureCounts -v
echo -e "\n>>> BAMCOVERAGE VERSION:"
$bamCoverage --version
echo -e ">>> END: Analayzer complete."
