
# Hua Sun
# Version 1.3
# 11/20/2019


# Mutect 2 - pipeline
# https://software.broadinstitute.org/gatk/documentation/article?id=24057


# USAGE
# sh run.sh -c config.ini -n <sampleName> -b <bam> -o <outdir>

# set 8Gb memory for MGI-server

# Set in MGI-Server
export PERL_PATH=/gsc
export PATH=$PERL_PATH/bin:$PATH
export PERL_BIN=$PERL_PATH/bin/perl
export PERL5LIB=$PERL_PATH/lib/perl5/5.8.7/:$PERL5LIB


# default
name=''
outdir=`pwd`/mutect2_result

while getopts "c:d:n:b:r:o:" opt; do
  case $opt in
    c)
      CONFIG=$OPTARG
      ;;
    d)
      scriptDir=$OPTARG
      ;;
    n)
      name=$OPTARG
      ;;  
    b)
      tumor_bam=$OPTARG
      ;;
    r)
      ref_fa=$OPTARG
      ;;
    o)
      outdir=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done


source $CONFIG


# check input
if [[ $name == '' ]];then
  echo "[ERROR] Please input -n sampleName ..." >&2
  exit
fi

if [ ! -f "$tumor_bam" ];then
  echo "[ERROR] The $tumor_bam is not existing ..." >&2
  exit
fi



# make outdir folder
mkdir -p $outdir

OUT=$outdir/$name
mkdir -p $OUT


##========================================##
##               STEP - 1
##========================================##


##==================== Make softlink for bam ====================##
# Create softlink bam
if [ ! -f "$OUT/$name.bam" ];then
  ln -s $tumor_bam $OUT/$name.bam
fi

if [ -f "$tumor_bam.bai" ];then
  ln -s $tumor_bam.bai $OUT/$name.bam.bai
else
  ${SAMTOOLS} index $OUT/$name.bam
fi




##==================== 1. Make a unfilter VCF ====================##
## Make unfilter file
# *f1r2.tar.gz  *unfiltered.vcf/idx/status
for chromosome in {'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'}; do
  
  ${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
        -Dsamjdk.use_async_io_write_samtools=true \
        -Dsamjdk.use_async_io_write_tribble=false \
        -Dsamjdk.compression_level=2 \
        -Xmx4g -jar ${GATK} Mutect2 \
        -I $OUT/$name.bam \
        -R ${REF_FA} \
        -L ${chromosome} \
        --germline-resource ${GNOMAD_VCF} \
        -pon ${PANEL_OF_NORMALS_VCF} \
        --f1r2-tar-gz ${OUT}/${chromosome}-f1r2.tar.gz \
        -O ${OUT}/${chromosome}-unfiltered.vcf
 
done



## Merge unfiltered-vcf
all_unfiltered_input=`for chrom in {1..22}; do printf -- "I=${OUT}/chr${chrom}-unfiltered.vcf "; done; for chrom in {'X','Y'}; do printf -- "I=${OUT}/chr${chrom}-unfiltered.vcf "; done`

${JAVA} -Xmx4G -jar ${PICARD} GatherVcfs \
    $all_unfiltered_input \
    O=${OUT}/merged-unfiltered.vcf


# Merged vcf.stats
all_unfiltered_stats_input=`for chrom in {1..22}; do printf -- "-stats ${OUT}/chr${chrom}-unfiltered.vcf.stats "; done; for chrom in {'X','Y'}; do printf -- "-stats ${OUT}/chr${chrom}-unfiltered.vcf.stats "; done`

${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx4g -jar ${GATK} MergeMutectStats \
      $all_unfiltered_stats_input \
      -O ${OUT}/merged-unfiltered.vcf.stats





##==================== 2. Make contamination table ====================##
## Make read-orientation-model
# make 'read-orientation-model.tar.gz'
all_f1r2_input=`for chrom in {1..22}; do printf -- "-I ${OUT}/chr${chrom}-f1r2.tar.gz "; done; for chrom in {'X','Y'}; do printf -- "-I ${OUT}/chr${chrom}-f1r2.tar.gz "; done`
  
# it must write to like this
${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
  -Dsamjdk.use_async_io_write_samtools=true \
  -Dsamjdk.use_async_io_write_tribble=false \
  -Dsamjdk.compression_level=2 \
  -Xmx4g -jar ${GATK} LearnReadOrientationModel \
  ${all_f1r2_input} \
  -O ${OUT}/read-orientation-model.tar.gz


## Make 'getpileupsummaries.table'
# set 8Gb for solve OutOfMemoryError 
${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx8g -jar ${GATK} GetPileupSummaries \
      -I $OUT/$name.bam \
      -V ${COMMON_BIALLELIC} \
      -L ${COMMON_BIALLELIC} \
      -O ${OUT}/getpileupsummaries.table


## Make 'contamination.table' & 'segments.table'
${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx4g -jar ${GATK} CalculateContamination \
      -I ${OUT}/getpileupsummaries.table \
      --tumor-segmentation ${OUT}/segments.table \
      -O ${OUT}/contamination.table



# 8/25/2019
# NOTE: the FilterMutectCalls will not run when the contamination.table vale is 'NaN'
# Solution is change the value to 0 0
perl -i -pe 's/\tNaN\t1\.0/\t0\t0/ if /\tNaN\t1\.0/' ${OUT}/contamination.table


##==================== 3. Make a final filter VCF ====================##
# Make a finial filtered vcf
${JAVA} -Dsamjdk.use_async_io_read_samtools=false \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=false \
      -Dsamjdk.compression_level=2 \
      -Xmx8g -jar $GATK FilterMutectCalls \
      -V ${OUT}/merged-unfiltered.vcf \
      -R ${REF_FA} \
      --tumor-segmentation ${OUT}/segments.table \
      --contamination-table ${OUT}/contamination.table \
      --ob-priors ${OUT}/read-orientation-model.tar.gz \
      -O ${OUT}/filtered.vcf




##==================== remove temporary file ====================##
rm -f ${OUT}/chr*


# 9/26
# check result file 'fail or pass'
mut=`grep -v ^# ${OUT}/filtered.vcf | wc -l`; if [ $mut == 0 ]; then touch $OUT/s1.failed; else touch $OUT/s1.passed; fi





##========================================##
##               STEP - 2
##========================================##

VCF=${OUT}/filtered.vcf

# pre-filter for saving running time of SnpSift
grep -v '^#' $VCF | awk -F['\t'] '$7=="PASS"' > $VCF.pass.tmp
grep '^#' $VCF > $VCF.header.tmp
cat $VCF.header.tmp $VCF.pass.tmp > $VCF.pass.h.tmp


vcf_input=$VCF.pass.h.tmp


# filter dbSNP
$JAVA -Xmx8G -jar ${SNPSIFT} annotate -id ${dbSNP_noCOSMIC} $vcf_input | $JAVA -Xmx8G -jar ${SNPSIFT} filter -n "(exists ID) & (ID =~ 'rs' )" > $OUT/filtered.rem_dbSNP_noCOSMIC.vcf

# output: filtered.rem_dbSNP_noCOSMIC.vcf

# remove tempfile
rm -f $VCF.pass.tmp $VCF.header.tmp $VCF.pass.h.tmp



##========================================##
##               STEP - 3
##========================================##

VCF=${OUT}/filtered.rem_dbSNP_noCOSMIC.vcf
fileName="${OUT}/filtered.rem_dbSNP_noCOSMIC"

# format VCF for VEP
grep ^# ${VCF} > ${VCF}.header

# extract 'PASS' variant sites
# For Mutect2 DP >= 5
grep -v ^# ${VCF} | awk -F['\t'] '$7=="PASS"' | perl -ne '$DP=$1 if /DP=(\d+)/; print if $DP>=5' > ${VCF}.fmt4vep
# with 'chr' or without 'chr' both okay


# if use conda version, no need to perl 
perl ${VEP} --buffer_size 300 --offline --cache \
    --dir ${VEP_CACHE_DIR} \
    --assembly ${ASSEMBLY} \
    --fasta ${VEP_GENOME} \
    --fork 8 --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol \
    --numbers --domains --canonical --protein --biotype --uniprot --tsl --pubmed \
    --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number \
    --no_escape --xref_refseq --failed 1 --minimal --flag_pick_allele --force_overwrite \
    --pick_order canonical,tsl,biotype,rank,ccds,length \
    --format vcf --vcf -i ${VCF}.fmt4vep -o ${VCF}.fmt4vep.vep_anno


mv ${VCF}.fmt4vep.vep_anno ${fileName}.vep.vcf

# remove
rm -f ${VCF}.header ${VCF}.fmt4vep



##========================================##
##               STEP - 4
##========================================##

VCF=${OUT}/filtered.rem_dbSNP_noCOSMIC.vcf
fileName="${OUT}/filtered.rem_dbSNP_noCOSMIC"


# format vcf form as NORMAL TUMOR
echo "[INFO] Make vcf ... !" 1>&2
perl -wlni.bak -e 'if(/^#/){if(/^#CHROM/){s/FORMAT\t.*$/FORMAT\tNORMAL\tTUMOR/;print}else{print}}else{@F=split("\t"); if(scalar @F==10){$F[7]=~s/$/\;set=mutect2/; $F[8]=~s/$/\t\.\/\./}; print join("\t",@F)}' ${vcf}

## vcf2maf
echo "[INFO] The vcf2maf ... !" 1>&2
# input - filtered.rem_dbSNP_noCOSMIC.vcf
# output - filtered.rem_dbSNP_noCOSMIC.maf
perl $scriptDir/scripts/4.vcf2maf.vep.pl --input-vcf ${VCF} --output-maf ${fileName}.maf --tumor-id ${name} --normal-id 'NORMAL' --ref-fasta ${VEP_GENOME} --file-tsl ${FILE_TSL}


## filter maf and output coding_region variants
echo "[INFO] Making filtered.rem_dbSNP_noCOSMIC.coding_filtered.maf ... !" 1>&2
perl $scriptDir/scripts/5.maf_filter.vep.tumor-only.pl --vcf ${VCF} --maf ${fileName}.maf --min-vaf ${MIN_VAF} --min-dep ${MIN_DEP} --min-mut ${MIN_MUT} -o ${fileName}.coding_filtered.maf
  
  
## remove snv that nearby indel
echo "[INFO] Making filtered.rem_dbSNP_noCOSMIC.coding_filtered.rem_nearbyIndel_snv.maf ... !" 1>&2
perl $scriptDir/scripts/6.remove_nearby_snv.pl ${fileName}.coding_filtered.maf
  
# final output is '*.rem_nearbyIndel_snv.maf'




