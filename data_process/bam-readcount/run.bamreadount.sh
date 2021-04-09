#!/bin/bash
# Hua Sun
# sh bamReadcount.sh -S sampleName -L sample_loci.txt -B sample.bam -G genome.fa -O outdir


BAMREADCOUNT=/gscmnt/gc2525/dinglab/rmashl/Software/src/bam-readcount-0.7.4/mybuild/bin/bam-readcount
BRC2VAF=~/scripts/bamreadcount/script/bamReadcount2vaf.pl

GENOME=/gscmnt/gc2737/ding/hsun/data/human_genome/gencode_GRCh38_v29/genome/GRCh38.primary_assembly.genome.fa


while getopts "S:L:B:G:O:" opt; do
  case $opt in
    S)
      SAMPLE=$OPTARG
      ;;
    L)
      LOCI=$OPTARG
      ;;
    B)
      BAM=$OPTARG
      ;;
    G)
      GENOME=$OPTARG
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

mkdir -p $OUTDIR

awk '{print $1"\t"$2"\t"$2}' $LOCI > $OUTDIR/$SAMPLE.rc.loci

$BAMREADCOUNT -q 10 -b 20 -f $GENOME -l $OUTDIR/$SAMPLE.rc.loci $BAM > $OUTDIR/$SAMPLE.rc

perl $BRC2VAF -s $SAMPLE -l $LOCI $OUTDIR/$SAMPLE.rc > $OUTDIR/$SAMPLE.rc.vaf

