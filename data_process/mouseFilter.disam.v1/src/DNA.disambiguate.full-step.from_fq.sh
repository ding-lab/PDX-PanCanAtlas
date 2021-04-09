#!/bin/bash

# Hua Sun

# v1
# 5/6/2019, 3/23/2019 - added the remove command for running many samples in the same time

# sh disambiguate.full-step.from_fq.sh -C contig.ini -N sampleName -1 name.fq1.gz -2 name.fq2.gz -O /path/outDir
# output - the result output to outDir/sampleName folder
# memory 18 Gb # some times need to high memory

# getOptions
while getopts "C:N:1:2:O:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    N)
      NAME=$OPTARG
      ;;
    1)
      FQ1=$OPTARG
      ;;
    2)
      FQ2=$OPTARG
      ;;    
    O)
      OUTDIR=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done


source $CONFIG

if [ ! -d $OUTDIR ]; then
  echo "[ERROR] The $OUTDIR not exists!" >&2
  exit 1
fi

if [ -z "$NAME" ]; then
  echo "[ERROR] The Name is empty!" >&2
  exit 1
fi


OUT=$OUTDIR/$NAME
mkdir -p $OUT

#######################################
##           mapping to human
#######################################
# bwa 
echo "[INFO] 1: human - bwa" >&2
$BWA mem -t 8 -M -R "@RG\tID:$NAME\tPL:illumina\tLB:$NAME\tPU:$NAME\tSM:$NAME" $REF_HUMAN $FQ1 $FQ2 | $SAMTOOLS view -Shb -o $OUT/$NAME.human.bam -


# sort bam by natural name
echo "[INFO] 2: human - samtools sort" >&2
$SAMTOOLS sort -m 2G -@ 6 -o $OUT/$NAME.human.sort.bam -n $OUT/$NAME.human.bam


##-------------- remove process file for save space
rm -f $OUT/$NAME.human.bam



#######################################
##           mapping to mouse
#######################################

echo "[INFO] 3: mouse - bwa" >&2
$BWA mem -t 8 -M -R "@RG\tID:$NAME\tPL:illumina\tLB:$NAME\tPU:$NAME\tSM:$NAME" $REF_MOUSE $FQ1 $FQ2 | $SAMTOOLS view -Shb -o $OUT/$NAME.mouse.bam - 


# sort bam by natural name
echo "[INFO] 4: mouse - samtools sort" >&2
$SAMTOOLS sort -m 2G -@ 6 -o $OUT/$NAME.mouse.sort.bam -n $OUT/$NAME.mouse.bam


##-------------- remove process file for save space
rm -f $OUT/$NAME.mouse.bam


#######################################
##           Do Disambiguate 
#######################################

# Disambiguate (mouse-filter)
echo "[INFO] 5: disambiguate" >&2
$DISAMBIGUATE -s $NAME -o $OUT -a bwa $OUT/$NAME.human.sort.bam $OUT/$NAME.mouse.sort.bam


##-------------- remove all disambiguate processing files
rm -f $OUT/$NAME.human.sort.bam $OUT/$NAME.mouse.sort.bam $OUT/$NAME.disambiguatedSpeciesB.bam


#######################################
##  make new fq from filtered pdx bam 
#######################################

echo "[INFO] 6: samtools sort" >&2
$SAMTOOLS sort -m 2G -@ 6 -o $OUT/$NAME.disam.sortbyname.bam -n $OUT/$NAME.disambiguatedSpeciesA.bam

echo "[INFO] 7: samtools fq" >&2
$SAMTOOLS fastq $OUT/$NAME.disam.sortbyname.bam -1 $OUT/$NAME.disam_1.fastq.gz -2 $OUT/$NAME.disam_2.fastq.gz -0 /dev/null -s /dev/null -n -F 0x900


##-------------- remove disambiguate spaciesA and disam.sortbyname bams
rm -f $OUT/$NAME.disambiguatedSpeciesA.bam 
rm -f $OUT/$NAME.disam.sortbyname.bam


#######################################
##  mapping to human reference and create to new bam 
#######################################

echo "[INFO] 8: bwa" >&2
$BWA mem -t 8 -M -R "@RG\tID:$NAME\tPL:illumina\tLB:$NAME\tPU:$NAME\tSM:$NAME" $REF_HUMAN $OUT/$NAME.disam_1.fastq.gz $OUT/$NAME.disam_2.fastq.gz | $SAMTOOLS view -Shb -o $OUT/$NAME.disam.reAlign.pre.bam -


##-------------- remove fq files
rm -f $OUT/$NAME.disam_1.fastq.gz $OUT/$NAME.disam_2.fastq.gz


# sort
echo "[INFO] 9: picard - sortsam" >&2
$JAVA -Xmx16G -jar $PICARD SortSam \
   CREATE_INDEX=true \
   I=$OUT/$NAME.disam.reAlign.pre.bam \
   O=$OUT/$NAME.disam.reAlign.bam \
   SORT_ORDER=coordinate \
   VALIDATION_STRINGENCY=STRICT


##-------------- remove process file for save space
rm -f $OUT/$NAME.disam.reAlign.pre.bam


# remove-duplication
echo "[INFO] 10: picard - markduplicates" >&2
$JAVA -Xmx16G -jar $PICARD MarkDuplicates \
   I=$OUT/$NAME.disam.reAlign.bam \
   O=$OUT/$NAME.disam.reAlign.remDup.bam \
   REMOVE_DUPLICATES=true \
   M=$OUT/$NAME.disam.reAlign.remDup.metrics.txt

# index bam
$SAMTOOLS index $OUT/$NAME.disam.reAlign.remDup.bam


##-------------- remove process file for save space
rm -f $OUT/$NAME.disam.reAlign.bam
rm -f $OUT/$NAME.ambiguous*
rm -f $OUT/$NAME.sorted.bai
rm -f $OUT/$NAME.disam.reAlign.bai
