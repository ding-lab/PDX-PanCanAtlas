#!/bin/bash

# Hua Sun

# v1
# 4/18/2019

# create folder within out_dir and the folder name set using -N

# sh disambiguate.full-step.from_fq.sh -C contig.ini -N test -1 name.fq1.gz -2 name.fq2.gz -O /path/outDir
# memory 32 Gb in MGI-server

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

# align sequencing reads to the genome
$STAR --runThreadN 8 --genomeDir $STAR_HUMAN_GENOME_DIR --sjdbGTFfile $GTF_HUMAN --sjdbOverhang 100 --readFilesIn $FQ1 $FQ2 --outFileNamePrefix $OUT/human. --outSAMtype BAM Unsorted --twopassMode Basic --outSAMattributes All --genomeLoad NoSharedMemory --readFilesCommand zcat
$STAR --runThreadN 8 --genomeDir $STAR_MOUSE_GENOME_DIR --sjdbGTFfile $GTF_MOUSE --sjdbOverhang 100 --readFilesIn $FQ1 $FQ2 --outFileNamePrefix $OUT/mouse. --outSAMtype BAM Unsorted --twopassMode Basic --outSAMattributes All --genomeLoad NoSharedMemory --readFilesCommand zcat

# sort bam by natural name
$SAMTOOLS sort -m 3G -@ 8 -o $OUT/human.sort.bam -n $OUT/human.Aligned.out.bam
$SAMTOOLS sort -m 3G -@ 8 -o $OUT/mouse.sort.bam -n $OUT/mouse.Aligned.out.bam

# remove process file for save space
rm -f $OUT/human.Aligned.out.bam $OUT/mouse.Aligned.out.bam

# Disambiguate (mouse-filter)
$DISAMBIGUATE -s $NAME -o $OUT -a star $OUT/human.sort.bam $OUT/mouse.sort.bam

# remove process file for save space
rm -f $OUT/human.sort.bam $OUT/mouse.sort.bam

# re-create fq
$SAMTOOLS sort -m 3G -@ 8 -o $OUT/$NAME.disam.sortbyname.bam -n $OUT/$NAME.disambiguatedSpeciesA.bam
$SAMTOOLS fastq $OUT/$NAME.disam.sortbyname.bam -1 $OUT/$NAME.disam.rna-seq.r1.fastq.gz -2 $OUT/$NAME.disam.rna-seq.r2.fastq.gz -0 /dev/null -s /dev/null -n -F 0x900

# remove process file for save space
rm -f $OUT/$NAME.disam.sortbyname.bam $OUT/$NAME.disambiguatedSpeciesA.bam $OUT/$NAME.disambiguatedSpeciesB.bam $OUT/$NAME.ambiguous*

rm -f $OUT/*.out $OUT/*.tab
rm -rf $OUT/*._STARgenome $OUT/*._STARpass1
