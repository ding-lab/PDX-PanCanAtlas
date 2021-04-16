#!/bin/bash

# Hua Sun
# 5/5/2019

# sh fq2bam_human.sh -C contig.ini -S sampleName -1 SAMPLE.fq1.gz -2 SAMPLE.fq2.gz -O /path/outDir
# -T human/mouse
# memory set 18 Gb

# getOptions
while getopts "C:S:1:2:O:" opt; do
  case $opt in
    C)
      CONFIG=$OPTARG
      ;;
    S)
      SAMPLE=$OPTARG
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


source ${CONFIG}


if [ ! -d $OUTDIR ]; then
	echo "[ERROR] The $OUTDIR not exists!" >&2
	exit 1
fi

if [ -z "$SAMPLE" ]; then
	echo "[ERROR] The SAMPLE is empty!" >&2
	exit 1
fi

if [ ! -f $GENOME ]; then
	echo "[ERROR] The $GENOME not exists!" >&2
	exit 1
fi


OUT="$OUTDIR/$SAMPLE"
mkdir -p $OUT

# mapping to reference and create to bam
$BWA mem -t 8 -M -R "@RG\tID:$SAMPLE\tPL:illumina\tLB:$SAMPLE\tPU:$SAMPLE\tSM:$SAMPLE" $GENOME $FQ1 $FQ2 | samtools view -Shb -o $OUT/$SAMPLE.bam -

# sort
$JAVA -Xmx16G -jar $PICARD SortSam \
   CREATE_INDEX=true \
   I=$OUT/$SAMPLE.bam \
   O=$OUT/$SAMPLE.sorted.bam \
   SORT_ORDER=coordinate \
   VALIDATION_STRINGENCY=STRICT

# remove sam for save space
rm -f $OUT/$SAMPLE.bam

# remove-duplication
$JAVA -Xmx16G -jar $PICARD MarkDuplicates \
   I=$OUT/$SAMPLE.sorted.bam \
   O=$OUT/$SAMPLE.remDup.bam \
   REMOVE_DUPLICATES=true \
   M=$OUT/$SAMPLE.remdup.metrics.txt

# index bam
$SAMTOOLS index $OUT/$SAMPLE.remDup.bam

# remove bam for save space
rm -f $OUT/$SAMPLE.sorted.bam $OUT/$SAMPLE.sorted.bam.bai


