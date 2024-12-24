#!/usr/bin/perl -w
use FindBin qw($Bin);
use Getopt::Long;

my ($sample, $dir, $gene_class, $pop, $wxs, $g_nom, $help);
$g_nom=0;
GetOptions(
           "s=s"     =>      \$sample,
           "i=s"     =>      \$dir,
	   "c=s"     =>      \$gene_class,
           "p=s"     =>      \$pop,
           "r=s"     =>      \$wxs,
           "g=s"     =>      \$g_nom,
           "h"       =>      \$help
);
my $usage = <<USE;
Usage:
description:Typing annotation of SpecComplex phased sequence
usage: perl $0 [options]
        Common options:
        -s       <tr>    sample name
        -i       <tr>    the directory of phased sequence
	-c       <tr>    gene class "HLA|KIR|CYP"
        -p       <tr>    population information "Asian|Black|Caucasian|Unknown|nonuse"
        -r       <tr>    focus region "exon|whole" ("exon" is suitable for WES or RNAseq; "whole" is suitable for WGS )
        -g       <tr>     G-translate 1|0"
        -help|?           print help information
e.g.:
        perl $0 -s samplename -i indir -p Unknown -r exon -g 1 -c HLA
USE
die $usage unless ($sample && $dir && $pop && $wxs && $gene_class );
print "parameter:\tsample:$sample\tdir:$dir\tpop:$pop\twxs:$wxs\tGeneclass:$gene_class\tG_nom:$g_nom\n";

my $cutoff = 100;
my $db="$Bin/../db";
my $bin="$Bin/../bin";

my @genes;
if($gene_class eq "HLA"){
         @genes = ( 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DPB2', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-E', 'HLA-F', 'HLA-G', 'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L', 'HLA-P', 'HLA-V', 'HLA-DQA2', 'HLA-DPA2', 'HLA-N', 'HLA-S', 'HLA-T', 'HLA-U', 'HLA-W', 'MICA', 'MICB', 'TAP1', 'TAP2', 'HFE' );
}elsif($gene_class eq "KIR"){
	@genes = ("KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DP1", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DP1", "KIR3DS1", "KIR2DL5");
}elsif($gene_class eq "CYP"){
        @genes=("CYP19A1", "CYP1A1", "CYP1B1", "CYP26A1", "CYP2A13", "CYP2A6", "CYP2B6", "CYP2C19", "CYP2C8", "CYP2C9", "CYP2D6", "CYP2F1", "CYP2J2", "CYP2R1", "CYP2S1", "CYP2W1", "CYP3A4", "CYP3A43", "CYP4A22", "CYP4B1", "CYP4F2", "CYP3A5", "CYP3A7", "CYP8A1");
}else{print "$gene_class\terror\n"}

my (%hash,%hashs,$ref);
my $fadir = $dir;
foreach my $class(@genes){
	if($wxs eq "exon"){$ref="$db/$gene_class/ref/split/$class.exon.fasta";}
	else{$ref="$db/whole/$class.fasta"}
	
        if($class =~ /KIR3DL3/){$cutoff = 800}	
	#my $tfa = "$fadir/result.$class.fasta";
	#`less $tfa |head -2 >$fadir/allele.1.$class.fasta`;
	#`less $tfa |tail -2 >$fadir/allele.2.$class.fasta`;
	for(my $i=1;$i<=2;$i++){
		my $fa = "$fadir/allele.$i.$class.fasta";
		my $tag = "$class"."_"."$i";
		#my $ref="$db/split/$class.fasta";
               `$bin/makeblastdb -in $fa -dbtype nucl -parse_seqids -out $fa`;
                if($class eq "KIR2DL5"){
			`$bin/blastn -query $fa -out $fadir/$tag.blast.out1 -db $ref -outfmt 7 -num_threads 4 -max_target_seqs 1000 `;
                        `$bin/blastn -query $ref -out $fadir/$tag.blast.out2 -db $fa -outfmt 7 -num_threads 4 -max_target_seqs 1000 `;}
               else{
	       `$bin/blastn -query $fa -out $fadir/$tag.blast.out1 -db $ref -outfmt 7 -num_threads 4 -max_target_seqs 200 -perc_identity 95 -qcov_hsp_perc 10`;
               `$bin/blastn -query $ref -out $fadir/$tag.blast.out2 -db $fa -outfmt 7 -num_threads 4 -max_target_seqs 200 -perc_identity 95 -qcov_hsp_perc 10`;}

	       my (%hash11,%hash12, %hash21, %hash22,$gene,$score,%hasha);
	       open IN1, "$fadir/$tag.blast.out1" or die "$!\n";
               while(<IN1>){
                       chomp;
                       next if(/^#/);
                       my ($hla, $t, $m,$d, $a) = (split)[1,3,4,5,11];
                       #next if($hla =~ /[N|Q]$/);
                       next if($t <$cutoff);
                       $hash11{$hla} += $t;
                       $hash12{$hla} += $m + $d;
		       $hasha{$hla} += $a;
              }
              close IN1;

              open IN2, "$fadir/$tag.blast.out2" or die "$!\n";
              while(<IN2>){
                      chomp;
                      next if(/^#/);
                      my ($hla, $t, $m,$d, $a) = (split)[0,3,4,5,11];
                      #next if($hla =~ /[N|Q]$/);
                      next if($t <$cutoff);
                      $hash21{$hla} += $t;
                      $hash22{$hla} += $m + $d;
              }
              close IN2;
              $score=85;
	      my $size = $cutoff;
              foreach my $key(sort keys %hash11){
                      next if(!exists $hash21{$key});
                      my @tt = (split /:/, $key);
                      my $kid = $tt[0];
		      if($#tt>0){$kid = "$tt[0]".":"."$tt[1]";}
		      #my $s1 = $hasha{$key};
		      my $s1 = 100 * (1 - $hash12{$key}/$hash11{$key});
                      my $s2 = 100 * (1 - $hash22{$key}/$hash21{$key});
                      next if($hash11{$key} < $cutoff);
                      next if($hash21{$key} < $cutoff);
                      my $scorel=$s1;
		      if($scorel > $score){
                             $score = $scorel;
                             $gene = $key;
			     $size = $hash11{$key}
                      }elsif($scorel == $score){
			      #if($score == 100){print "$key\t$size\t$hash11{$key}\n"}
			      if($hash11{$key} > $size){$gene = $key; $size=$hash11{$key}}
		      }
             }
	     $hash{$tag} = $gene;
	     $hashs{$tag} = $score;
	     `rm -rf $fa.*`;
     }
}

open OUT, ">$fadir/allele.phase.result.exon.txt";
#my $line1 = "Sample";
#my $line2 = $sample;
#my $line3 = "Score";
#foreach my $kk(sort keys %hash){
#        $line1 .= "\t$kk";
#        my $hla = $hash{$kk};
#	my $score = $hashs{$kk};
#        $line2 .= "\t$hla";
#	$line3 .= "\t$score";
#}
#print OUT "$line1\n$line2\n$line3\n";
print OUT "Sample\tType\tAllele\tScore\n";
foreach my $kk(sort keys %hash){
	my $score = $hashs{$kk};
	my $hla = $hash{$kk};
	print OUT "$sample\t$kk\t$hla\t$score\n";
}
close OUT;
