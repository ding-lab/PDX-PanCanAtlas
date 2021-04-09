=head1

    Hua Sun
    9/15/2019
    v0.1
    

    ToDo:
        Filter MAF
        
    Input:
        VCF & MAF, which get from Mutect2 tumor-only pipeline
        # .vcf
        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR
    
    Output:
      Filtered coding-region variants of maf, which add reads depth info from vcf
    
    Usage:
        perl maf_filter.vep.tumor-only.pl --vcf mutect2.vcf --maf mutation.maf --min-vaf 0.05 -o output
        
        # The best setting in PDX samples
        perl maf_filter.vep.tumor-only.pl --vcf mutect2.vcf --maf mutation.maf --min-vaf 0.1 --min-dep 20 --min-mut 5 -o output
        
    
    Options:
        --vcf         <file.vcf>    # add vcf for call
        --maf         <file.maf>    # original maf no count info
        --min-vaf       [0.001]      # minimum tumor vaf
        --min-dep       [8]        # minimum total depth
        --min-mut       [3]         # minimum alternate reads
        --max-indel-len [100]        # minimum indel length
        -o           <output.maf>  # output name

=cut


use strict;
use Getopt::Long;

my $vcf = '';
my $maf = '';

my $min_vaf = 0.04;
my $min_depth = 8;
my $min_mut = 4;
my $max_indel_len = 50;  
# PMID: 20858594 (2-16bp indel 96%)

my $output = './output.filtered.maf';
my $help;

GetOptions(
        "vcf:s" => \$vcf,
        "maf:s" => \$maf,
        "min-vaf:f" => \$min_vaf,
        "min-dep:i" => \$min_depth,
        "min-mut:i" => \$min_mut,
        "max-indel-len:i" => \$max_indel_len,
        "o|output:s" => \$output,
        "help" => \$help
);


die `pod2text $0 !` if ( $help || $vcf eq '' || $maf eq '' );


# maf file
my @vcf_data = `grep -v ^# $vcf`;
my @maf_data = `grep -v ^# $maf | grep -v Hugo_Symbol`;
my $head = `grep Hugo_Symbol $maf`;


# make index for loci 
my %vcf_hash;
my $key;

# FORMAT                      TUMOR
# GT:AD:AF:DP:F1R2:F2R1:SB    0/1:24,3:0.135:27:7,0:17,3:20,4,2,1
my ($chr, $pos, $ref, $alt, $info, $key);
my ($GT, $AD, $AF, $DP, $F1R2, $F2R1, $SB);
my ($ref_count, $alt_count, $vaf);
my ($len_ref_seq, $len_alt_seq);

foreach (@vcf_data){
        chomp;
        my @arr = split(/\t/);
    
        $chr = $arr[0];
        $pos = $arr[1];
        $ref = $arr[3];
        $alt = $arr[4];
        $info = $arr[9];
        ($GT, $AD, $AF, $DP, $F1R2, $F2R1, $SB) = split(/:/, $info);
        ($ref_count, $alt_count) = split(",", $AD);
        next if ($DP==0);
        $vaf = $alt_count/$DP;
        $key = "$chr:$pos:$ref:$alt";
        $vcf_hash{$key} = "$GT:$DP:$ref_count:$alt_count:$vaf";
}


# extract filtered maf
# output header
open my ($OUT), '>', "$output";
print $OUT $head;

my $maf;
my $anno;
my ($ref_seq, $alt_seq);
foreach (@maf_data){
        my @arr = split(/\t/);

        $chr = $arr[4];
        $pos = $arr[5];
        $ref = $arr[10];
        $alt = $arr[12];
    
        $key = "$chr:$pos:$ref:$alt";
    
        $anno = $arr[8];
        
        # filter variants
        if ( exists $vcf_hash{$key} ){
                ($GT, $DP, $ref_count, $alt_count, $vaf) = split(':', $vcf_hash{$key});
                        
                # filter variations by read counts
                if ( $DP >= $min_depth && $alt_count >= $min_mut && $vaf >= $min_vaf ){
                        $arr[39] = $DP;
                        $arr[40] = $ref_count;
                        $arr[41] = $alt_count;
                        
                        # filer indel > 10 bp
                        $len_ref_seq = length($arr[10]);
                        $len_alt_seq = length($arr[12]);
                        
                        next if ( $len_ref_seq > $max_indel_len || $len_alt_seq > $max_indel_len);
                        
                        # extract coding region
                        # ref: somaticwrapper
                        if ($anno=~/Frame_/ || $anno=~/Missense_/ || $anno=~/Nonsense_/ ||  $anno=~/Nonstop_/ || $anno=~/Silent/ || $anno=~/Splice_Site/) {
                        
                                print $OUT join("\t", @arr);
                        }
                }
        }
        
    
}



