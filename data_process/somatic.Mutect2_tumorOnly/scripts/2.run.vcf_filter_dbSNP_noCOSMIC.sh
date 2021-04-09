
# Hua Sun
# 8/13/2019

# USAGE
# sh run.sh -i <vcf> -d <dbSNP_noCOSMIC>
# *.vcf (not *.vcf.gz)

# http://snpeff.sourceforge.net/SnpSift.html#filter


# set memory 8GB

while getopts "c:i:" opt; do
  case $opt in
  	c)
  		CONFIG=$OPTARG
      ;;
    i)
    	vcf=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done


source $CONFIG


if [ ! -f ${dbSNP_noCOSMIC} ]; then
	echo "[ERROR] No dbSNP file or direction ..." >&2
	exit 1
fi

echo "[INFO] Filter dbSNP location from vcf ..."  >&2

filename=${vcf%.vcf}
OUT=$filename.rem_dbSNP_noCOSMIC.vcf


# pre-filter for saving running time of SnpSift
grep -v '^#' $vcf | awk -F['\t'] '$7=="PASS"' > $vcf.pass.tmp
grep '^#' $vcf > $vcf.header.tmp
cat $vcf.header.tmp $vcf.pass.tmp > $vcf.pass.h.tmp

vcf_input=$vcf.pass.h.tmp


# filter dbSNP
$JAVA -Xmx8G -jar ${SNPSIFT} annotate -id ${dbSNP_noCOSMIC} $vcf_input | $JAVA -Xmx8G -jar ${SNPSIFT} filter -n "(exists ID) & (ID =~ 'rs' )" > $OUT


# remove tempfile
rm -f $vcf.pass.tmp $vcf.header.tmp $vcf.pass.h.tmp

