
# 8/16/2020

##--------- Somatic Mut
somaticMut=../2.filtered_ncimatch_exclusions.v2/somaticMut/1.somaticMut.with_arm.filtered.out

outdir=./1.ncimaitch_alt
mkdir -p $outdir

sed '1d' $somaticMut | cut -f 2-3 | sort -u > $outdir/ncimatch.somaticMut.alt.out


##--------- CNV
cnv_amp=../2.filtered_ncimatch_exclusions.v2/cnv/1.cn_amp.filtered.out

sed '1d' $cnv_amp | cut -f 2 | sort -u | perl -pe 's/$/\tAmp/' > $outdir/ncimatch.cn_amp.out


cnv_del=../2.filtered_ncimatch_exclusions.v2/cnv/1.cn_del.filtered.out

sed '1d' $cnv_del | cut -f 2 | sort -u | perl -pe 's/$/\tDel/' > $outdir/ncimatch.cn_del.out


##--------- Fusion
fusion=../2.filtered_ncimatch_exclusions.v2/fusion/1.fusion.with_arm.filtered.out

sed '1d' $fusion | cut -f 3 | sort -u > $outdir/ncimatch.fusion.out



wc -l $outdir/*


