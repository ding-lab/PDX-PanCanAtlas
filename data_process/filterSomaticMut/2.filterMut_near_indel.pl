
=head1

    Hua Sun
    1/14/2020

    filter near INDEL mutations

    perl filterMut_near_indel.pl --form meta3 input.meta3.tsv -o outdir
    
    --form meta3 (default)
    -w 20 (default)
    -o results (default)
    -outName outname (default input file name add *.remMutNearIndel.tsv)
    

=cut


use strict;
use warnings;
use Getopt::Long;

my $form = 'meta3';
my $outdir = "results";
my $window=20;
GetOptions(
    "form:s" => \$form,
    "w|window:s" => \$window,
    "o|outdir:s" => \$outdir
);


my $file = shift;
my $fileName = $file;
$fileName =~ s/.*\///;

my $outfile = "$fileName.remMutNearIndel.tsv";



# main
if ($form eq 'meta3'){

    # separate sample
    `echo "[INFO] Separate file ..." >&2`;
    &Separate_samples($file, 1, $outdir);
    
    # filter mut per file
    `echo "[INFO] Filter mutations per file ..." >&2`;
    my @sampleFile = readpipe("ls $outdir/temp_dir/*.org.tsv");
    foreach (@sampleFile){
        chomp;
        &Meta3_FilterMutations_per_sample($_, $window);
    }
    
    # merge filtered files
    `echo "[INFO] Merge files ..." >&2`;
    system("head -n 1 $file > $outdir/header");
    system("cat $outdir/header $outdir/temp_dir/*.filtered > $outdir/$outfile");
    system("sed '1d' $outdir/$outfile | cut -f 1 | sort | uniq -c | perl -pe 's/ +/\\t/g; s/^\\t//' | sort -k1,1n > $outdir/$outfile.status");
    system("rm -f $outdir/header");
    system("rm -rf $outdir/temp_dir");

}



## separate_samples
sub Separate_samples
{
    my ($file, $sampleCol, $outdir) = @_;

    `mkdir -p $outdir`;
    `mkdir -p $outdir/temp_dir`;

    system("sed '1d' $file | cut -f $sampleCol | sort -u | while read sample; do awk -F['\\t'] '\$$sampleCol==\"'\$sample'\"' $file > $outdir/temp_dir/\$sample.org.tsv; done");
}



sub Meta3_FilterMutations_per_sample
{
    my ($file, $window) = @_;

    my @data = readpipe("cat $file");

    open my $OUT_removed, '>', "$file.removed";
    open my $OUT_passed, '>', "$file.filtered";

    # add indel mutation sites
    my %hashInDel;
    my $Variant_Type;
    my ($gene, $startLoci, $endLoci);
    foreach (@data){
        
        my @arr = split("\t");
        
        $Variant_Type = $arr[11];
        $gene = $arr[8];
        $startLoci = $arr[3];
        $endLoci = $arr[4];

        # solve strange issue
        if ($startLoci=~/\D/ || $endLoci=~/\D/){
            `echo [WARNING] the loci is not integer $startLoci $endLoci ... >&2`;
            next;
        }

        if ($Variant_Type eq 'INS' || $Variant_Type eq 'DEL'){
            if (exists $hashInDel{$gene}){
                my $newLoci = "$startLoci $endLoci";
                my $orgLoci = $hashInDel{$gene};
                if (&Check_addIndel_or_not($orgLoci, $newLoci, $window) == 0){
                    $hashInDel{$gene} .= ",$startLoci $endLoci";
                    print $OUT_passed $_;
                } else {
                    print $OUT_removed $_;
                }
            } else {
                $hashInDel{$gene} = "$startLoci $endLoci";
                print $OUT_passed $_;
            }
        }
    }
    

    # remove snp based on indel regeion
    my $flag_remove;
    foreach (@data){
        my @arr = split("\t");

        $Variant_Type = $arr[11];
        $gene = $arr[8];
        $startLoci = $arr[3];
        $endLoci = $arr[4];

        $flag_remove = 0;

        if ($Variant_Type ne 'INS' && $Variant_Type ne 'DEL'){
            if (exists $hashInDel{$gene}){
                my $str = $hashInDel{$gene};

                # compare with multiple indel per gene
                if ($str =~ /\,/){
                    my @subArr = split(',', $str);
                    foreach my $s (@subArr){
                        my ($start, $end) = split(' ', $s);
                        if ($startLoci>=($start-$window) && $endLoci<=($end+$window)){
                            print $OUT_removed $_;
                            $flag_remove = 1;
                            last;
                        }
                    }
                    print $OUT_passed $_ if ($flag_remove == 0);
                } else {
                    # compare with single indel per gene
                    my ($start, $end) = split(' ', $str);
                    if ($startLoci>=($start-$window) && $endLoci<=($end+$window)){
                        print $OUT_removed $_;
                    } else {
                        print $OUT_passed $_;
                    }
                }

            } else {
                print $OUT_passed $_;
            }
        }
    }
}



sub Check_addIndel_or_not
{
    my ($orgLoci, $newLoci, $window) = @_;

    my ($new_start, $new_end) = split(' ', $newLoci);

    if ($orgLoci =~ /\,/){
        my @arr = split(',');
        my ($start, $end);
        foreach (@arr){
            ($start, $end) = split(' ');

            # solve strange issue
            if ($start=~/\D/ || $end=~/\D/){
                `echo [WARNING] the loci is not integer $start $end ... >&2`;
                return 0;
            }

            my $upstream = $start - $window;
            my $downstream = $end + $window;
            return 1 if ($new_start >= $upstream && $new_start <= $downstream);
            return 1 if ($new_end >= $upstream && $new_end <= $downstream);
        }
    } else {
        my ($start, $end) = split(' ', $orgLoci);
        my $upstream = $start - $window;
        my $downstream = $end + $window;

        return 1 if ($new_start >= $upstream && $new_start <= $downstream);
        return 1 if ($new_end >= $upstream && $new_end <= $downstream);
    }

    return 0;
}



