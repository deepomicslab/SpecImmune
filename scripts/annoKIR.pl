#!/usr/bin/perl -w
use FindBin qw($Bin);
use Getopt::Long;

my ($sample, $dir, $pop, $wxs, $g_nom, $db, $rna, $help);
$g_nom = 0;
$rna='0';
GetOptions(
           "s=s"     =>      \$sample,
           "i=s"     =>      \$dir,
           "p=s"     =>      \$pop,
           "r=s"     =>      \$wxs,
           "g=s"     =>      \$g_nom,
           "d=s"     =>      \$db,
           "t=s"     =>      \$rna,
           "h"       =>      \$help
);

my $usage = <<USE;
Usage:
description: HLAtyping annotation of SpecHLA phased sequence
usage: perl $0 [options]
        Common options:
        -s       <tr>    sample name
        -i       <tr>    the directory of phased sequence
        -p       <tr>    population information "Asian|Black|Caucasian|Unknown|nonuse"
        -r       <tr>    focus region "exon|whole|tgs" ("exon" is suitable for WES or RNAseq; "whole" is suitable for WGS; "tgs" is suitable for TGS )
        -g       <tr>     G-translate 1|0"
        -d       <tr>    the directory of database
        -t       <tr>    Is RNAseq "0|1"
        -help|?           print help information
e.g.:
        perl $0 -s samplename -i indir -p Unknown -r exon -g 1
USE
die $usage unless ($sample && $dir && $pop && $wxs && $db);

print "parameter:\tsample:$sample\tdir:$dir\tpop:$pop\twxs:$wxs\tG_nom:$g_nom\n";

my $cutoff = 100;
my $bin = "$Bin/../bin";

my @genes = (
    "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DP1", "KIR2DS1", 
    "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", 
    "KIR3DL3", "KIR3DP1", "KIR3DS1", "KIR2DL5A", "KIR2DL5B"
);

my (%hash, %hashs, $ref);
my $fadir = $dir;
my $workdir = "$dir/tmp";
open BOUT, ">$dir/KIR.blast.summary.txt";
`mkdir -p $workdir`;

foreach my $class (@genes) {
    my $ref = "$db/KIR/$class/$class";
    $cutoff = 800 if ($class =~ /KIR3DL3/);
    
    for (my $i = 1; $i <= 2; $i++) {
        my $fa = "$fadir/KIR.allele.$i.$class.fasta";
        # print file path
              print "$fa\n";
        # if fa file not exist or empty, skip
        next if (!-e $fa || -z $fa);

        my $tag = "$class\_$i";
        if ($class eq "KIR2DL5") {
            `blastn -query $fa -out $workdir/$tag.blast.out1 -db $ref -outfmt 7 -num_threads 4 -max_target_seqs 10000000`;
            # print command
                print "blastn -query $fa -out $workdir/$tag.blast.out1 -db $ref -outfmt 7 -num_threads 4 -max_target_seqs 10000000\n";
       #      `blastn -query $ref -out $fadir/$tag.blast.out2 -db $fa -outfmt 7 -num_threads 4 -max_target_seqs 10000000`;
        } 
        elsif ($class eq "KIR3DL2") {
            `blastn -query $fa -out $workdir/$tag.blast.out1 -db $ref -outfmt 7 -gapopen 1 -gapextend 1 -num_threads 4 -max_hsps 1 -max_target_seqs 10000000`;
        } 
        
        else {
            `blastn -query $fa -out $workdir/$tag.blast.out1 -db $ref -outfmt 7 -num_threads 4 -max_target_seqs 10000000 -perc_identity 95 -qcov_hsp_perc 10`;
       #      `blastn -query $ref -out $fadir/$tag.blast.out2 -db $fa -outfmt 7 -num_threads 4 -max_target_seqs 10000000 -perc_identity 95 -qcov_hsp_perc 10`;
        }

        my (%hash11, %hash12, %hasha, $gene, $score);
        open IN1, "$workdir/$tag.blast.out1" or die "$!\n";
        while (<IN1>) {
            chomp;
            next if (/^#/);
            my ($hla, $iden, $t, $m, $d, $a) = (split)[1, 2, 3, 4, 5, 11];
            next if ($t < $cutoff);
            $hash11{$hla} += $t;
            
            $hasha{$hla} += $a;
            if ($class eq "KIR3DL2") {
                $hash12{$hla} += (100-$iden) * $t / 100;

            }
            else{
                $hash12{$hla} += $m + $d;
            }
        }
        close IN1;

        $score = 85;
        my $size = $cutoff;
        foreach my $key (sort keys %hash11) {
       #      next if (!exists $hash21{$key});
            my @tt = split /:/, $key;
       #      my $kid = "$tt[0]".":"."$tt[1]";
       #      $kid = "$tt[0]:$tt[1]" if ($#tt > 0);
            my $s1 = 100 * (1 - $hash12{$key} / $hash11{$key});
            print BOUT "$tag\t$key\t$hash12{$key}\t$hash11{$key}\t$s1\n";
            next if ($hash11{$key} < $cutoff);
       #      next if ($hash21{$key} < $cutoff);
            my $scorel = $s1;
            if ($scorel > $score) {
                $score = $scorel;
                $gene = $key;
                $size = $hash11{$key};
            } elsif ($scorel == $score) {
                if ($hash11{$key} > $size) {
                    $gene = $key; 
                    $size = $hash11{$key};
                }
            }
        }
        $hash{$tag} = $gene;
        $hashs{$tag} = $score;
        `rm -rf $fadir/$class.temp*`;
    }
}

# open OUT, ">$fadir/allele.phase.result.exon.txt";
# print OUT "Sample\tType\tAllele\tScore\n";
# foreach my $kk (sort keys %hash) {
#     my $score = $hashs{$kk};
#     my $hla = $hash{$kk};
#     print OUT "$sample\t$kk\t$hla\t$score\n";
# }
# close OUT;