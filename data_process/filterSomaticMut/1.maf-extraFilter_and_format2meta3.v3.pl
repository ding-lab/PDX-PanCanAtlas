=head1
    
    Hua Sun
    1/16/2020 

    Filter maf mutations and format to meta3 form

        Filter-tumorNormal
            General trustable mutation : t-total-depth>=14, n-total-depth>=8, t-vaf>=0.05, n-vaf<=0.01, min-mut-reads>=4, t-mut-reads/n-mut-reads >=10
            Extra_re-call_mutation: (re-called mutation sites (driver genes))
                t-total-depth>=8, n-total-depth>=8, t-vaf>=0.05, n-mut < 2, min-mut-reads 4, remove the tumor-alt/normal-alt < 20 when normal-alt >0
        
        Filter-tumorOnly
            General trustable mutation : t-total-depth>=20, t-vaf>=0.1, min-mut-reads>=4
            Extra_re-call_mutation: (re-called mutation sites (driver genes))
                t-total-depth>=8, t-vaf>=0.1, min-mut-reads>=4


    perl run.maf-extraFilter_and_format2meta.pl --maf maf --tag tumorNormal --delNameTag yes --outdir outdir

    -maf         input.maf
    -outdir      ./outdir
    -mutType     non-silent(default) / all
    -tag         tumorNormal/tumorOnly
    -delNameTag  'yes'     # '_T'  remove end of '_T' from sample name
    -refGenome   'hg38'

    outdir/tumorNormal.xxx

=cut


use strict;
use Getopt::Long;

my $maf='';
my $outdir;
my $mutType='non-silent';
my $tag='tumorNormal';
my $delNameTagFromMaf='no';
my $refGenome='GRCh38';
my $lenIndel=50;
my $help;
GetOptions(
    "maf:s" => \$maf,
    "outdir:s" => \$outdir,
    "tag:s" => \$tag,
    "delNameTag:s" => \$delNameTagFromMaf,
    "mutType:s" => \$mutType,
    "refGenome:s" => \$refGenome,
    "lenIndel:i" => \$lenIndel,
    "h|help" => \$help
);

die `pod2text $0` if ($help || $maf eq '');
die `echo [Error] Please set --mutType non-silent or all ...` if ($mutType ne 'non-silent' && $mutType ne 'all');



`mkdir $outdir` if ( ! -d $outdir );

`head -n 1 $maf > $outdir/head_maf`;


## 1. separate general maf and extra-driverGene maf

`echo "[INFO] 1. filter the standared calling mut and re-call mut ..." >&2`;
&extra_filter_maf_tumorNormal($maf, $tag, $lenIndel, $outdir) if ($tag eq 'tumorNormal');

&extra_filter_maf_tumorOnly($maf, $tag, $lenIndel, $outdir) if ($tag eq 'tumorOnly');



## 2. merge general & extra-called maf & make non-silent mutation maf
`echo "[INFO] 2. merged.filtered.$mutType.maf ..." >&2`;
if ($mutType eq 'non-silent') {

    `cat $outdir/$tag.general.filter.maf $outdir/$tag.extra.filter.maf | awk -F['\t'] '\$9!="Silent"' | grep -v Tumor_Sample_Barcode | cat $outdir/head_maf - > $outdir/$tag.filtered.$mutType.maf`;
}

if ($mutType eq 'all') {

    `cat $outdir/$tag.general.filter.maf $outdir/$tag.extra.filter.maf | grep -v Tumor_Sample_Barcode | cat $outdir/head_maf - > $outdir/$tag.filtered.$mutType.maf`;
}



## 3. format maf to meta3tsv

`echo "[INFO] 3. merged.filtered.$mutType.meta3.tsv ..." >&2`;
`perl -e 'print "Tumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tChromosome\tStart_position\tEnd_position\tReference_Allele\tTumor_Seq_Allele2\tTranscript_ID\tHugo_Symbol\tHGVSp_Short\tVariant_Classification\tVariant_Type\tt_depth\tt_ref_count\tt_alt_count\tt_vaf\tn_depth\tn_ref_count\tn_alt_count\tn_vaf\tNCBI_Build\tMutation_Group\n"' > $outdir/$tag.header.m3.tmp`;

if ($delNameTagFromMaf eq 'yes'){
    `sed '1d' $outdir/$tag.filtered.$mutType.maf | perl -pe 's/ +//g' | awk -F['\t'] '{tumor_vaf="."; normal_vaf="."; mutTag="somaticMut"; if(\$3=="Extra"){mutTag="somaticMut:Extra"}; if(\$40!="" && \$40>0){tumor_vaf=\$42/\$40}; if(\$43!="" && \$43>0){normal_vaf=\$45/\$43}; {tumorID=substr(\$16, 1,length(\$16)-2); normalID=substr(\$17, 1,length(\$17)-2); print tumorID,normalID,\$5,\$6,\$7,\$11,\$13,\$38,\$1,\$37,\$9,\$10,\$40,\$41,\$42,tumor_vaf,\$43,\$44,\$45,normal_vaf,"$refGenome",mutTag}}' | perl -pe 's/ /\t/g' | cat $outdir/$tag.header.m3.tmp - > $outdir/$tag.filtered.$mutType.meta3.tsv`;
} else {
    `sed '1d' $outdir/$tag.filtered.$mutType.maf | perl -pe 's/ +//g' | awk -F['\t'] '{tumor_vaf="."; normal_vaf="."; mutTag="somaticMut"; if(\$3=="Extra"){mutTag="somaticMut:Extra"}; if(\$40!="" && \$40>0){tumor_vaf=\$42/\$40}; if(\$43!="" && \$43>0){normal_vaf=\$45/\$43}; {print \$16,\$17,\$5,\$6,\$7,\$11,\$13,\$38,\$1,\$37,\$9,\$10,\$40,\$41,\$42,tumor_vaf,\$43,\$44,\$45,normal_vaf,"$refGenome",mutTag}}' | perl -pe 's/ /\t/g' | cat $outdir/$tag.header.m3.tmp - > $outdir/$tag.filtered.$mutType.meta3.tsv`;
}



## 4: only chr1-22,X,Y
`echo "[INFO] 4. extract chr1-22,X,Y ..." >&2`;
`perl -ne 'print if (/\tchr[XY1234567890]+\t/)' $outdir/$tag.filtered.$mutType.meta3.tsv | cat $outdir/$tag.header.m3.tmp - > $outdir/$tag.filtered.$mutType.meta3.chr1-22XY.tsv`;

if ($tag eq 'tumorOnly'){
    `perl -i -ne '\@F=split("\t"); if(\$F[1] eq "NORMAL"){\$F[1]="PoolNormal"}; print join("\t",\@F)' $outdir/$tag.filtered.$mutType.meta3.chr1-22XY.tsv`;
}



## 5: mut status
`echo "[INFO] 5. mutation status ..." >&2`;
`sed '1d' $outdir/$tag.filtered.$mutType.meta3.chr1-22XY.tsv | cut -f 1 | sort | uniq -c | perl -pe 's/ +/\t/g; s/^\t//' | sort -k1,1n > $outdir/$tag.filtered.$mutType.meta3.chr1-22XY.tsv.status`;


# remove before filtered file
`rm -f $outdir/$tag.header.m3.tmp $outdir/$tag.extra.filter.maf $outdir/$tag.general.filter.maf $outdir/$tag.filtered.$mutType.meta3.tsv`;





## extra_filter_maf
#------------------ tumorNormal
sub extra_filter_maf_tumorNormal
{
    my ($maf, $tag, $lenIndel, $outdir) = @_;
        
    my @data = `grep -v Tumor_Sample_Barcode $maf`;
        
    open my $outExtra, '>', "$outdir/$tag.extra.filter.maf";
    open my $outGeneral, '>', "$outdir/$tag.general.filter.maf";

        
    my $start;
    my $end;
    my $len;
    my $t_vaf;
    my $n_vaf;
    my $t_depth;
    my $t_alt_count;
    my $n_depth;
    my $n_alt_count;

    foreach (@data){
        my @arr = split/\t/;

        $start = $arr[5];
        $end = $arr[6];
        $len = $end - $start + 1;
        next if ($len > $lenIndel);  # filter large indel

        $t_depth = $arr[39];
        $t_alt_count = $arr[41];

        $n_depth = $arr[42];
        $n_alt_count = $arr[44];

        $t_vaf = $t_alt_count/$t_depth;
        $n_vaf = $n_alt_count/$n_depth;

                
        # t_depth
        next if ($t_depth == 0 || $t_depth eq '' || $t_depth eq '.');
                
                
        # Center
        if ($arr[2] eq 'Extra'){

            if ($t_vaf >= 0.05 && $t_depth >=8 && $t_alt_count >3){                                
                    # n_depth
                if ($n_depth>=8 && $n_vaf <= 0.01){
                        # tumorMut reads >= normalMut reads *20
                    print $outExtra $_ if ($n_alt_count < 2 && $n_alt_count>0 && $t_alt_count/$n_alt_count>=20);
                    print $outExtra $_ if ($n_alt_count==0);
                }
            }
        } else {
            # tumor-normal -  t-total-depth>=14, n-total-depth>=8, t-vaf>=0.05, n-vaf<=0.01, min-mut-reads>=4, t-mut-reads/n-mut-reads >=10
            if ($t_vaf >= 0.05 && $t_depth >=14 && $t_alt_count >3){
                if ($n_depth>=8 && $n_vaf <= 0.01){
                    print $outGeneral $_ if ($n_vaf<=0.01 && $n_alt_count>0 && $t_alt_count/$n_alt_count>=10);
                    print $outGeneral $_ if ($n_alt_count==0);
                }
            }
        }
    }
}




#------------------ tumorOnly
sub extra_filter_maf_tumorOnly
{
    my ($maf, $tag, $lenIndel, $outdir) = @_;
        
    my @data = `grep -v Tumor_Sample_Barcode $maf`;
        
    open my $outExtra, '>', "$outdir/$tag.extra.filter.maf";
    open my $outGeneral, '>', "$outdir/$tag.general.filter.maf";

    my $start;
    my $end;
    my $len;        
    my $t_vaf;
    my $t_depth;
    my $t_alt_count;

    foreach (@data){
        my @arr = split/\t/;

        $start = $arr[5];
        $end = $arr[6];
        $len = $end - $start + 1;
        next if ($len > $lenIndel);  # filter large indel

        $t_depth = $arr[39];
        $t_alt_count = $arr[41];

        $t_vaf = $t_alt_count/$t_depth;
                
        next if ($t_depth == 0 || $t_depth eq '' || $t_depth eq '.');
                
                
        if ($arr[2] eq 'Extra'){

            print $outExtra $_ if ($t_vaf >= 0.1 && $t_depth >=8 && $t_alt_count >3);

        } else {
            # tumor-only
            # t-total-depth>=20, t-vaf>=0.1, min-mut-reads>=4
            print $outGeneral $_ if ($t_vaf >= 0.1 && $t_depth >=20 && $t_alt_count >3);
        }
    }
}


