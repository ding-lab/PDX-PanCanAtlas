#!/bin/bash

# Hua Sun
# 2/6/2020 v2.1
# star-fusion_v1.6


# USAGE
# sh starFusion.sh -S sampleName -1 /path/*.r1.fq.gz -2 /path/*.r2.fq.gz -O `pwd`/fusion_results
# result output to outDir/sampleName/starFusion 

# set memory 42GB (default 40GB)

# input *.r1.fq.gz  *.r2.fq.gz

export PATH=/gscmnt/gc3021/dinglab/hsun/software/miniconda2/envs/starFusion_v1_6/bin:$PATH

# version v1.6.0
STAR_FUSION=/gscmnt/gc3021/dinglab/hsun/software/miniconda2/envs/starFusion_v1_6/bin/STAR-Fusion
# GRCh38
CTAT_LIB_DIR=/gscmnt/gc2737/ding/hsun/data/fusionDatabase/GRCh38/GRCh38_gencode_v29_CTAT_lib_Mar272019/ctat_genome_lib_build_dir


while getopts "S:1:2:O:" opt; do
  case $opt in
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



mkdir -p $OUTDIR

OUT=$OUTDIR/$SAMPLE

mkdir -p $OUT



# auto create "StarFusionOut" folder in $OUT
starFusionOutDir=$OUT/StarFusionOut

$STAR_FUSION --genome_lib_dir $CTAT_LIB_DIR \
	--left_fq $FQ1 \
	--right_fq $FQ2 \
	--FusionInspector validate \
	--examine_coding_effect \
	--output_dir $starFusionOutDir
#	--denovo_reconstruct

# backup FusionInspector-validate for remove processing data
mkdir -p $starFusionOutDir/FusionInspector-validate2
cp $starFusionOutDir/FusionInspector-validate/finspector.FusionInspector.fusions.abridged.tsv* $starFusionOutDir/FusionInspector-validate2/
# for visualization
cp $starFusionOutDir/FusionInspector-validate/*fusion_inspector_web.* $starFusionOutDir/FusionInspector-validate2/
cp -r $starFusionOutDir/FusionInspector-validate/fi_workdir $starFusionOutDir/FusionInspector-validate2/


# remove peocessing data
rm -f $starFusionOutDir/Aligned.out.bam $starFusionOutDir/Log.* $starFusionOutDir/*.log
rm -f $starFusionOutDir/SJ.out.tab $starFusionOutDir/pipeliner.9.cmds
rm -rf $starFusionOutDir/star-fusion.preliminary $starFusionOutDir/_STAR*
rm -rf $starFusionOutDir/FusionInspector-validate


# rename to original folder name
mv $starFusionOutDir/FusionInspector-validate2 $starFusionOutDir/FusionInspector-validate


