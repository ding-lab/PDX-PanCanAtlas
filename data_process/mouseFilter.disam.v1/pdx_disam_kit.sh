#!/bin/bash

# Hua Sun
# v1 5/14/2019

# usage: 

# dna full step
# sh pdx_toolkit.sh -p dnaFull -n sample -1 fq1.gz -2 fq2.gz -o outdir

# rna full step
# sh pdx_toolkit.sh -p rnaFull -n sample -1 fq1.gz -2 fq2.gz -o outdir


# separate steps
# step-1 make human to bam
# sh pdx_toolkit.sh -p humanBam -n sample -1 fq1.gz -2 fq2.gz -o outdir

# step-1 make mouse to bam
# sh pdx_toolkit.sh -p mouseBam -n sample -1 fq1.gz -2 fq2.gz -o outdir

# step-2 disambiguate
# sh pdx_toolkit.sh -p disambiguate -n sample -d pdx_wxs_dir


# set script dir
scriptDir=/gscmnt/gc3021/dinglab/hsun/PDX-workflow/PDXToolkit/mouseFilter.disam.v1

# set config.ini
config=${scriptDir}/config/config.pdx.mgi.v2.ini


while getopts "p:n:b:1:2:d:o:" opt; do
  case $opt in
    p)
      pipeline=$OPTARG
      ;;
    n)
      name=$OPTARG
      ;;
    b)
      bam=$OPTARG
      ;;
    1)
      fq1=$OPTARG
      ;;  
    2)
      fq2=$OPTARG
      ;;
    d)
      dir=$OPTARG
      ;;
    o)
      outdir=$OPTARG
      ;;
    \?)
      echo "script usage: $(basename $0) [-t] [-n] " >&2
      exit 1
      ;;
  esac
done
#shift "$(($OPTIND -1))"


###############################
##        DNA full step
###############################

if [[ $pipeline == "dnaFull" ]]; then
  sh $scriptDir/src/DNA.disambiguate.full-step.from_fq.sh -C ${config} -N ${name} -1 ${fq1} -2 ${fq2} -O ${outdir}
fi


###############################
##        RNA full step
###############################

if [[ $pipeline == "rnaFull" ]]; then
  sh $scriptDir/src/RNA.disambiguate.full-step.from_fq.sh -C ${config} -N ${name} -1 ${fq1} -2 ${fq2} -O ${outdir}
fi


###############################
##        DNA mapping
###############################

# mapping to human
if [[ $pipeline == "humanBam" ]]; then
  sh $scriptDir/src/DNA.makeHumanBam.sh -C ${config} -N ${name} -1 ${fq1} -2 ${fq2} -O ${outdir}
fi

# mapping to mouse
if [[ $pipeline == "mouseBam" ]]; then
  sh $scriptDir/src/DNA.makeMouseBam.sh -C ${config} -N ${name} -1 ${fq1} -2 ${fq2} -O ${outdir}
fi


###############################
##       DNA disambiguate
###############################

# do disambiguate and create new bam
if [[ $pipeline == "disambiguate" ]]; then
  sh $scriptDir/src/DNA.disambiguate.sh -C ${config} -N ${name} -D ${dir}
fi


