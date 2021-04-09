=head1
    
    Hua Sun
    3/19/2020 make meta3
    
    Filtered make non-silent meta3-maf form

    perl maf2meta3.v2.pl --maf maf --prefix tumorNormal --outdir outdir

    -maf          input.maf
    -outdir       ./outdir
    -tag          tumorNormal/tumorOnly
    -refGenome    'hg38'

    outdir/*.maf2meta3Form.tsv

=cut

use strict;
use Getopt::Long;

my $maf='';
my $outdir;
my $prefix='tumorOnly';
my $refGenome='GRCh38';

my $min_depth = 20;
my $min_mut = 4;
my $min_vaf = 0.1;
my $max_indel = 50;

my $non_silent;

my $help;
GetOptions(
    "maf:s" => \$maf,
    "outdir:s" => \$outdir,
    "prefix:s" => \$prefix,
    "refGenome:s" => \$refGenome,
    "min_vaf:f" => \$min_vaf,
    "min_depth:i" => \$min_depth,
    "min_mut:i" => \$min_mut,
    "max_indel:i" => \$max_indel,
    "non_silent" => \$non_silent,
    "h|help" => \$help
);

die `pod2text $0` if ($help || $maf eq '');

`mkdir $outdir` if ( ! -d $outdir );

&Maf2meta3Form($prefix, $maf, $refGenome, $outdir, $min_vaf, $min_depth, $min_mut, $max_indel, $non_silent);



sub Maf2meta3Form
{
    my ($prefix, $maf, $refGenome, $outdir, $min_vaf, $min_depth, $min_mut, $max_indel, $non_silent) = @_;

    ## 3. format maf to meta3tsv
    # output 'somaticMut.filtered.non-silent.meta3.tsv'
    `echo "[INFO] make meta3 form from maf file ..." >&2`;
    # header
    open my $OUT, '>', "$outdir/$prefix.maf2meta3Form.tsv";
    print $OUT "Tumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2\tTranscript_ID\tHugo_Symbol\tHGVSp_Short\tVariant_Classification\tVariant_Type\tt_depth\tt_ref_count\tt_alt_count\tt_vaf\tn_depth\tn_ref_count\tn_alt_count\tn_vaf\tNCBI_Build\tMutation_Group\n";
    # make meta3tsv
    # add mutation-tag as somaitcMut/somaticMut:Extra

    my @data = `sed '1d' $maf`;

    my $opt;
    foreach (@data){
        chomp;
        my @arr = split("\t");

        $opt->{'t_vaf'} = "";
        $opt->{'n_vaf'} = "";
        $opt->{'mutationGroup'} = "somaticMut";

        $opt->{'tumor_sample'} = $arr[15];
        $opt->{'normal_sample'} = "PoolNormal";
        $opt->{'chr'} = $arr[4];
        $opt->{'start'} = $arr[5];
        $opt->{'end'} = $arr[6];
        $opt->{'ref'} = $arr[10];
        $opt->{'alt'} = $arr[12];
        $opt->{'transcriptID'} = $arr[37];
        $opt->{'gene'} = $arr[0];
        $opt->{'HGVSp'} = $arr[36];
        $opt->{'variantClass'} = $arr[8];
        $opt->{'variantType'} = $arr[9];
        $opt->{'t_depth'} = $arr[39];
        $opt->{'t_ref_count'} = $arr[40];
        $opt->{'t_alt_count'} = $arr[41];
        $opt->{'t_vaf'} = sprintf("%.4f", $arr[41]/$arr[39]) if ($arr[39]!="" && $arr[39]>0);
        $opt->{'n_depth'} = $arr[42];
        $opt->{'n_ref_count'} = $arr[43];
        $opt->{'n_alt_count'} = $arr[44];
        $opt->{'n_vaf'} = sprintf("%.4f", $arr[44]/$arr[42]) if ($arr[42]!="" && $arr[42]>0);
        $opt->{'refGenome'} = $refGenome;
        $opt->{'mutationGroup'} = "somaticMut:Extra" if ($arr[2] eq 'Extra');

        # indel length
        my $ref_len = length($opt->{'ref'});
        my $alt_len = length($opt->{'alt'});
        my $indel_len = $ref_len;
        $indel_len = $alt_len if ($alt_len > $ref_len);

        # output
        if ($non_silent){
            next if ($opt->{'variantClass'} eq 'Silent');
        }

        if ($indel_len <= $max_indel &&  $opt->{'t_vaf'} >= $min_vaf && $opt->{'t_depth'} >= $min_depth &&  $opt->{'t_alt_count'} >= $min_mut){
            print $OUT "$opt->{'tumor_sample'}\t$opt->{'normal_sample'}\t$opt->{'chr'}\t$opt->{'start'}\t$opt->{'end'}\t",
                    "$opt->{'ref'}\t$opt->{'alt'}\t$opt->{'transcriptID'}\t$opt->{'gene'}\t$opt->{'HGVSp'}\t$opt->{'variantClass'}\t",
                    "$opt->{'variantType'}\t$opt->{'t_depth'}\t$opt->{'t_ref_count'}\t$opt->{'t_alt_count'}\t$opt->{'t_vaf'}\t",
                    "$opt->{'n_depth'}\t$opt->{'n_ref_count'}\t$opt->{'n_alt_count'}\t$opt->{'n_vaf'}\t",
                    "$opt->{'refGenome'}\t$opt->{'mutationGroup'}\n";
        }
    }

}




