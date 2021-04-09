
input=./01-Jul-2020-ClinicalEvidenceSummaries.tsv
output=./civic20200730.ClinicalEvidenceSummaries.Mut.tsv


# 2) evidence_status = accepted
# 3) evidence_direction = Supports
# 3) clinical_significance = Positive
#perl -pe 's/ +/_/g' $input | awk -F['\t'] '$7!=""' | awk -F['\t'] '$10=="Supports" && $20=="accepted"' | awk -F['\t'] '{print $1,$3,$4,$7,$8,$9,$10,$11,$12}' | perl -pe 's/ +/\t/g' > $output.tmp
cut -d$'\t' -f 1,3-4,7-12,14,20 $input > $output.tmp
head -n 1 $output.tmp > head.tmp
awk -F['\t'] '$4!="" && $7=="Supports"' $output.tmp | cat head.tmp - > $output


# for merging with depo
cut -f 1-4,8-10 $output | perl -pe 's/$/\tCiViC/' > $output.formal_form.out

rm *.tmp

