#!/bin/bash

# Hua Sun
# hua.sun@wustl.edu

# v1.0.5  11/20/2019  added -r 
# v1.0.4  11/15/2019  replace input style
# v1.0.3  11/3/2019   run full-steps by one step
# v1.0.2  10/2/2019   change the parameter
# v1.0.1  9/30/2019   added merge function
# v1.0    9/16/2019



## USAGE:

## Note: Please set config.ini before running pipeline

## sh somaticMut.tumor-only.mutect2.sh -c <config.ini> -p <programname> -n <name> -b <bam> -o <outdir>

# -p s0      full steps as one
# -p s1      call draft mutations
# -p s2      filter germline mutations
# -p s3      annotation by VEP
# -p s4      vcf2maf and filter
# -p s5      merge all of sample mafs



## run-full steps by one (-p s0)
# sh somaticMut.tumor-only.mutect2.sh -p full -n sample -b bam -o outdir

# re-run
# sh somaticMut.tumor-only.mutect2.sh -p full -r yes -n sample -b bam -o outdir



## step 1 to step 5

# step-1 call mutations (-p s1)
# -p call (or -p s3)
# sh somaticMut.tumor-only.mutect2.sh -p call -n sample -b bam -o outdir

# step-2 filter-germline (-p s2)
# -p filter-germline (or -p s3)  it reads outdir/sample/filtered.vcf
# sh somaticMut.tumor-only.mutect2.sh -p filter-germline -n sample -o outdir

# step-3 annotation (-p s3)
# -p anno (or -p s3)    it reads outdir/sample/filtered.rem_dbSNP_noCOSMIC.vcf
# sh somaticMut.tumor-only.mutect2.sh -p anno -n sample -o outdir

# step-4 vcf2maf and filter 
# -p vcf2maf (or -p s4)
# sh somaticMut.tumor-only.mutect2.sh -p vcf2maf -n sample -o outdir

# step-5 merge maf sample files, which are in the same folder 
# -p merge (or -p s5)
# sh somaticMut.tumor-only.mutect2.sh -p merge -o outdir 



##============================================================================##

# set script dir
scriptDir=/gscuser/hua.sun/scripts/tumor-only-pipeline.mutect2

# set config.ini
config=${scriptDir}/config.mgi.ini

reRun='no'

while getopts "c:p:r:n:b:o:" opt; do
  case $opt in
    c)
      config=$OPTARG
      ;;
    p)
      program=$OPTARG
      ;;
    r)
      reRun=$OPTARG
      ;;
    n)
      name=$OPTARG
      ;;  
    b)
      bam=$OPTARG
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


source $config


###############################
##  STEP-1    Call mutations
###############################

# input name.fastq.gz
if [[ $program == "call" ]] || [[ $program == "s1" ]]; then

    if [ ! -e $bam ]; then
        echo "[ERROR] The $bam is not exists !" 1>&2
        exit
    fi
    
    if [ ! -d $outdir ]; then
        echo "[ERROR] The $outdir is not exists !" 1>&2
        exit
    fi
    
    
    if [ -d $outdir/$name ] && [ $reRun == "no" ]; then
        echo "[WARNING] The $outdir/$name is exists !" 1>&2
        exit
    fi
    
    
    # output - filtered.vcf
    sh $scriptDir/scripts/1.run.gatk4.mutect2.Tumor-OnlyVariantCallingPipeline.v1.sh -c ${config} -n ${name} -b ${bam} -o ${outdir}

fi



###############################
##  STEP-2    Filter-germline
###############################

# input name.fastq.gz
if [[ $program == "filter-germline" ]] || [[ $program == "s2" ]]; then
    
    vcf=$outdir/$name/filtered.vcf
    
    if [ ! -e $vcf ]; then
        echo "[ERROR] The $vcf is not exists !" 1>&2
        exit
    fi
    
    
    # input - filtered.vcf
    # output - filtered.rem_dbSNP_noCOSMIC.vcf
    sh $scriptDir/scripts/2.run.vcf_filter_dbSNP_noCOSMIC.sh -c ${config} -i ${vcf}

fi



###############################
##  STEP-3    Annotation
###############################

# remove duplication for bam
if [[ $program == "anno" ]] || [[ $program == "s3" ]]; then
    
    vcf=$outdir/$name/filtered.rem_dbSNP_noCOSMIC.vcf
    
    if [ ! -e $vcf ]; then
        echo "[ERROR] The $vcf is not exists !" 1>&2
        exit
    fi
    
    # input - filtered.rem_dbSNP_noCOSMIC.vcf
    # output - filtered.rem_dbSNP_noCOSMIC.vep.vcf
    sh $scriptDir/scripts/3.run.vcf_annotation4pass.vep.sh -c ${config} -i ${vcf}

fi




###############################
##  STEP-4   vcf2maf and filter
###############################

# vcf2maf
if [[ $program == "vcf2maf" ]] || [[ $program == "s4" ]]; then

    vcf=$outdir/$name/filtered.rem_dbSNP_noCOSMIC.vcf

    if [ ! -e $vcf ]; then
        echo "[ERROR] The $vcf is not exists !" 1>&2
        exit
    fi
    
    
    fileName=$outdir/$name/filtered.rem_dbSNP_noCOSMIC

    # format vcf form as NORMAL TUMOR
    echo "[INFO] Make vcf ... !" 1>&2
    # col-8 set=mutect2
    perl -wlni.bak -e 'if(/^#/){if(/^#CHROM/){s/FORMAT\t.*$/FORMAT\tNORMAL\tTUMOR/;print}else{print}}else{@F=split("\t"); if(scalar @F==10){$F[7]=~s/$/\;set=mutect2/; $F[8]=~s/$/\t\.\/\./}; print join("\t",@F)}' ${vcf}
    
    
    ## vcf2maf
    echo "[INFO] The vcf2maf ... !" 1>&2
    # input - filtered.rem_dbSNP_noCOSMIC.vcf
    # output - filtered.rem_dbSNP_noCOSMIC.maf    
    fileName=${vcf%%.vcf}
    perl $scriptDir/scripts/4.vcf2maf.vep.pl --input-vcf ${vcf} --output-maf ${fileName}.maf --tumor-id ${name} --normal-id 'NORMAL' --ref-fasta ${VEP_GENOME} --file-tsl ${FILE_TSL}
    
    
    ## filter maf and output coding_region variants
    echo "[INFO] Making filtered.rem_dbSNP_noCOSMIC.coding_filtered.maf ... !" 1>&2
    # input - filtered.rem_dbSNP_noCOSMIC.vcf
    # output - filtered.rem_dbSNP_noCOSMIC.coding_filtered.maf
    perl $scriptDir/scripts/5.maf_filter.vep.tumor-only.pl --vcf ${vcf} --maf ${fileName}.maf --min-vaf ${MIN_VAF} --min-dep ${MIN_DEP} --min-mut ${MIN_MUT} -o ${fileName}.coding_filtered.maf
    
    
    ## remove snv that nearby indel
    echo "[INFO] Making filtered.rem_dbSNP_noCOSMIC.coding_filtered.rem_nearbyIndel_snv.maf ... !" 1>&2
    # input - filtered.rem_dbSNP_noCOSMIC.coding_filtered.maf
    # output-1 - filtered.rem_dbSNP_noCOSMIC.coding_filtered.rem_nearbyIndel_snv.maf
    # output-2 - filtered.rem_dbSNP_noCOSMIC.coding_filtered.rem_nearbyIndel_snv.maf.removed
    perl $scriptDir/scripts/6.remove_nearby_snv.pl ${fileName}.coding_filtered.maf
    
    # final output is '*.rem_nearbyIndel_snv.maf'
    
fi



###############################
##   Merge mafs for all of samples
###############################

# merge all of mafs
if [[ $program == "merge" ]] || [[ $program == "s5" ]]; then
    
    outName=`echo $outdir | sed 's/\/$//' | sed 's/.*\///'`
    dir=`echo $outdir | sed 's/\/$//'`
    
    # output - *.merged.maf
    cat $dir/*/filtered.rem_dbSNP_noCOSMIC.coding_filtered.rem_nearbyIndel_snv.maf | perl -ne 'print if (!/^Hugo_Symbol/ && $.>1 || /^Hugo_Symbol/ && $.==1)' > $dir/$outName.merged.maf    
    
    # output status
    cut -f 16 $dir/$outName.merged.maf | sed '1d' | sort | uniq -c | perl -pe 's/ +/\t/g; s/^\t//' | sort -k1,1n > $dir/$outName.merged.maf.status

fi



###############################
##   Run - Full pipelines by one step
###############################

if [[ $program == "full" ]] || [[ $program == "s0" ]]; then

    if [ ! -e $bam ]; then
        echo "[ERROR] The $bam is not exists !" 1>&2
        exit
    fi
    
    if [ ! -d $outdir ]; then
        echo "[ERROR] The $outdir is not exists !" 1>&2
        exit
    fi
    
    if [ -d $outdir/$name ] && [ $reRun == "no" ]; then
        echo "[WARNING] The $outdir/$name is exists !" 1>&2
        exit
    fi
    
    # output - filtered.vcf
    sh $scriptDir/scripts/run.full-pipelines.perSample.sh -c ${config} -d ${scriptDir} -n ${name} -b ${bam} -o ${outdir}

fi



