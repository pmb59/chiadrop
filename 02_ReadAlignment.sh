#!/bin/sh

# Dependencies:
# samtools, bowtie2, deeptools, BEDtools

############################################################
# Bowtie2 alignment PE
# Trim GEMcodes from R1 with Trimmomatic (R2 unmodified)
############################################################
# Code below for mouse data (GRCm38)

#threads
T=4
data='/fastq'
cd $data
m=$(ls *_R1.fastq.gz | sed -e 's/\_R1.fastq.gz$//' )

for J in $(echo $m) ; do   

  #Trimmomatics with opcion HEADCROP=16 to TRIM the GEM barcode in R1 FASTQ
  java -jar trimmomatic-0.39.jar SE -threads $T  \
    ${J}_R1.fastq.gz headcrop_${J}_R1.fastq.gz \
    HEADCROP:16

  # bowtie2 align
  bowtie2 -q -p $T --very-sensitive-local -N 1 -x /bowtie2_indices/GRCm38 -1 headcrop_${J}_R1.fastq.gz -2 ${J}_R2.fastq.gz -S ${J}.sam
  #-I 0 -X 1000 --no-mixed --no-discordant 

  # SAM to BAM to bigWig
  samtools view -@ $T -b -o ${J}.bam ${J}.sam
  samtools sort -@ $T -o ${J}_sorted.bam  ${J}.bam
  samtools index -@ $T ${J}_sorted.bam

  #properly-paired (FLAG == 0x2) reads to BED format and MAQC>= 40
  samtools view -@ $T -f 0x2 -bh -q 40 -o q40_${J}_sorted.bam ${J}_sorted.bam
  samtools index -@ $T q40_${J}_sorted.bam

  #deeptools bam to bigWig
  bamCoverage --bam q40_${J}_sorted.bam -o q40_${J}_sorted.bw  --numberOfProcessors $T -of bigwig

  # BamToBed
  bedtools bamtobed -bed12 -i q40_${J}_sorted.bam >  q40_${J}_sorted.bed

  #https://sites.google.com/site/anshulkundaje/projects/blacklists
  #cut -c4- mm10.blacklist.bed  > mm10.blacklist_nochr.bed
  bedtools intersect -v -wa -a q40_${J}_sorted.bed -b mm10.blacklist_nochr.bed  >  q40_${J}_sorted_filtered.bed

done

# For next steps in the analysis, these files are necessary:
# q40_${J}_sorted.bam
# q40_${J}_sorted.bai
# q40_${J}_sorted_filtered.bed
