# 8/15/2020

# Gene	Variant	Disease	Drug	EvidenceLevel	Effect	PubMedID



##============== format CiViC
input=../CiViC-20200730/01-Jul-2020-ClinicalEvidenceSummaries.tsv


# 1) evidence_status = accepted
# 2) evidence_direction = Supports
# 3) clinical_significance = Positive
cut -d$'\t' -f 1,3-4,7-12,14,20 $input > civic.tmp
awk -F['\t'] '$4!="" && $7=="Supports"' civic.tmp | cut -f 1-4,8-10 | perl -pe 's/$/\tCiViC/' > civic.formal.form.out



##============== format DEPO
input=../depo/DEPO_final_20180328.rjm.tsv

perl -pe 's/ +/\,/g' $input | awk -F['\t'] '{print $1,$3,$2,$6,$8,$4,$9}' | perl -pe 's/ /\t/g; s/\,/ /g' | perl -pe 's/$/\tDEPO/' | grep -v 'PubMed' > depo.formal_form.out



##============== Merge

cat depo.formal_form.out civic.formal.form.out | perl -pe 's/^/Gene\tVariant\tDisease\tDrug\tEvidenceLevel\tEffect\tPubMedID\n/ if $.==1' > depo_civic.merged.out



rm *.tmp civic.formal.form.out depo.formal_form.out



# manually further edited 'depo_civic.merged.out' and made 'depo_civic.merged.edited.out'


