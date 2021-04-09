
# 8/15/2020

##============================= Extract from db
depo_civic=/Users/huasun/GoogleDrive/Database/DrugDB/depo_civic/depo_civic.merged.edited.out

# manual check using /Users/huasun/GoogleDrive/Database/DrugDB/depo_civic/depo_civic.merged.out


outdir=./2.drugdb
mkdir -p $outdir


#--- only for checking purpose
sed '1d' /Users/huasun/GoogleDrive/Database/DrugDB/depo_civic/depo_civic.merged.out | cut -f 1 | sort -u > $outdir/depo_civic.gene
#--- \


# extract only alt from drugDB
sed '1d' $depo_civic | cut -f 1-2 | sort -u > $outdir/depo_civic.alt.out

# any
awk -F['\t'] '$2~"^any" || $2~"^Any"' $outdir/depo_civic.alt.out | sort -u > $outdir/alt_any.gene

# cnv-amp
grep AMPLIFICATION $outdir/depo_civic.alt.out | sort -u > $outdir/cnv_amp.out

# cnv-del
grep DELETION $outdir/depo_civic.alt.out | sort -u > $outdir/cnv_del.out




##============================= Compare alt
outdir=./3.compared
mkdir -p $outdir

#------ mutation
cat 1.ncimaitch_alt/ncimatch.somaticMut.alt.out 2.drugdb/depo_civic.alt.out | sort | uniq -d > $outdir/ncimatch.known.mut.out

cat 1.ncimaitch_alt/ncimatch.somaticMut.alt.out $outdir/ncimatch.known.mut.out | sort | uniq -u > $outdir/ncimatch.unknown.mut.out



#------ cnv-amp
cat 1.ncimaitch_alt/ncimatch.cn_amp.out 2.drugdb/cnv_amp.out | cut -f 1 | sort | uniq -d > $outdir/ncimatch.known.amp.out

cat 1.ncimaitch_alt/ncimatch.cn_amp.out $outdir/ncimatch.known.amp.out | cut -f 1 | sort | uniq -u > $outdir/ncimatch.unknown.amp.out



#------ cnv-del
cat 1.ncimaitch_alt/ncimatch.cn_del.out 2.drugdb/cnv_del.out | cut -f 1 | sort | uniq -d > $outdir/ncimatch.known.del.out

cat 1.ncimaitch_alt/ncimatch.cn_del.out $outdir/ncimatch.known.del.out | cut -f 1 | sort | uniq -u > $outdir/ncimatch.unknown.del.out



wc -l $outdir/*
