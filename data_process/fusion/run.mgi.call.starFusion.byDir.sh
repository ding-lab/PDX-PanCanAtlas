
# 2/5/2020;2/11/2020

fq_dir=$1


lsf_submit=~/pdxNet/hsun/Toolkit/LSF_mgi/lsf_submit.sh

bashScript=./call.starFusion.full.v2.sh

OUTDIR=`pwd`/starFusion
mkdir -p $OUTDIR



ls ${fq_dir} | grep -v '^\.' | sed 's/\///' | while read sampleID
do

  FQ1=${fq_dir}/${sampleID}/*.r1.fastq.gz
  FQ2=${fq_dir}/${sampleID}/*.r2.fastq.gz
  
  sh ${lsf_submit} 48 1 ${sampleID} bash ${bashScript} -S ${sampleID} -1 ${FQ1} -2 ${FQ2} -O ${OUTDIR}
  
done
