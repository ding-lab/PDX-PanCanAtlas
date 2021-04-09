#!/bin/bash

# Hua Sun
# 2/27/2019

# sh trimGalor.sh -C config.ini -S test -1 name.fq1.gz -2 name.fq2.gz -L 50 -O /path/outDir
# memory set 8 Gb
# min length=50  length related accurate
# output *_val_1.fq.gz & *_val_2.fq.gz
# automatically do fastqc

# getOptions
MINLEN=50

while getopts "C:S:1:2:L:O:" opt; do
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
    L)
      MINLEN=$OPTARG
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


OUT=$OUTDIR/$SAMPLE
mkdir -p $OUT

# default
# This step removes entire read pairs if at least one of the two sequences became shorter than a certain threshold.
$TRIMGALORE --phred33 --fastqc --length $MINLEN -q 20 -o $OUT --paired $FQ1 $FQ2


# if use Bowtie(1) to mapping need to use -t
# $TRIMGALORE -t --phred33 --fastqc --length $MINLEN -q 20 -o $OUT --paired $FQ1 $FQ2

#7z x fastqc.zip

