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
        -r       <tr>    focus region "exon|whole|tgs" ("exon" is suitable for WES or RNAseq; "whole" is suitable for WGS, "tgs" is suitable for TGS, pacbio )
        -g       <tr>    G-translate 1|0"
        -help|?          print help information
e.g.:
        perl $0 -s samplename -i indir -p Unknown -r exon -g 1 -c HLA 
USE
die $usage unless ($sample && $dir && $pop && $wxs && $gene_class);
print "parameter:\tsample:$sample\tdir:$dir\tpop:$pop\twxs:$wxs\tGeneclass:$gene_class\tG_nom:$g_nom\n";

my $cutoff = 100;
my $db="$Bin/../db";
my $bin="$Bin/../bin";
my $bias = 1;
if($wxs eq "tgs"){$bias = 0.5} # gap score in TGS
my @genes;
if($gene_class eq "HLA"){
         @genes = ( 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DPB2', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-E', 'HLA-F', 'HLA-G', 'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L', 'HLA-P', 'HLA-V', 'HLA-DQA2', 'HLA-DPA2', 'HLA-N', 'HLA-S', 'HLA-T', 'HLA-U', 'HLA-W', 'MICA', 'MICB', 'TAP1', 'TAP2', 'HFE' );
}elsif($gene_class eq "KIR"){
	@genes = ("KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DP1", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DP1", "KIR3DS1", "KIR2DL5");
}elsif($gene_class eq "CYP"){
        @genes=("CYP19A1", "CYP1A1", "CYP1B1", "CYP26A1", "CYP2A13", "CYP2A6", "CYP2B6", "CYP2C19", "CYP2C8", "CYP2C9", "CYP2D6", "CYP2F1", "CYP2J2", "CYP2R1", "CYP2S1", "CYP2W1", "CYP3A4", "CYP3A43", "CYP4A22", "CYP4B1", "CYP4F2", "CYP3A5", "CYP3A7", "CYP8A1");
}else{print "$gene_class\terror\n"}

my (%hash,%hashs,$ref,%hashp,%hashpp,%hashp2);
my %hashg = ("HLA-A"=>1, "HLA-B"=>1, "HLA-C" => 1, "HLA-DPA1"=>1, "HLA-DPB1"=>1, "HLA-DQA1"=>1, "HLA-DQB1"=>1, "HLA-DRB1"=>1);
my $fadir = $dir;
if( ! -d "$fadir/tmp" ){`mkdir $fadir/tmp`};
#population frequency
open FIN, "$db/HLA/HLA_FREQ_HLA_I_II.txt" or die "$!\n";
<FIN>;
while(<FIN>){
    chomp;
    my ($gene,$c,$b,$a) = (split);
    $gene = "HLA-"."$gene";
    $a = sprintf "%.5f",$a;
    $b = sprintf "%.5f",$b;
    $c = sprintf "%.5f",$c;
    $hashpp{$gene} = "$c;$b;$a";
    $hashp2{$gene} = ($a+$b+$c)/3;
    if($pop eq "Unknown"){$hashp{$gene} = ($a+$b+$c)/3}
    if($pop eq "Asian"){$hashp{$gene} = $a}
    if($pop eq "Black"){$hashp{$gene} = $b}
    if($pop eq "Caucasian"){$hashp{$gene} = $c}
    if($pop eq "nonuse"){$hashp{$gene} = 0}
}
close FIN;

foreach my $class(@genes){
	if($wxs eq "exon"){$ref="$db/$gene_class/ref/split/$class.exon.fasta";}
	elsif($class eq "HLA-DRB1"){$ref="$db/$gene_class/ref/split/$class.exon.fasta";}
	elsif($wxs ne "tgs" && ($class eq "CYP3A5" || $class eq "CYP2D6" || $class eq "CYP2A13"|| $class eq "CYP2C8" || $class eq "CYP2C9")) {$ref="$db/$gene_class/ref/split/$class.exon.fasta";}
	else{$ref="$db/whole/$class.fasta"}
	$cutoff = 100;
	if($class =~ /KIR/){$cutoff = 600}
        if($class =~ /KIR3DL3/){$cutoff = 800}	
	if($class eq "HLA-DOB"){$cutoff = 500}
	if($wxs eq "exon"){
		my $fa = "$fadir/allele.1.$class.fasta";
                `echo ">$class.1" >$fadir/tmp/allele.1.$class.fasta`;
                `less $fa |grep -v ">" >>$fadir/tmp/allele.1.$class.fasta`;
                $fa = "$fadir/allele.2.$class.fasta";
                `echo ">$class.2" >$fadir/tmp/allele.2.$class.fasta`;
                `less $fa |grep -v ">" >>$fadir/tmp/allele.2.$class.fasta`;
        }else{
		if($class eq "HLA-B"){
			`$bin/samtools faidx $fadir/allele.1.$class.fasta HLA-B_0:200-3000  > $fadir/tmp/allele.1.$class.fasta`;
			`$bin/samtools faidx $fadir/allele.2.$class.fasta HLA-B_1:200-3000  > $fadir/tmp/allele.2.$class.fasta`;}
                elsif($class eq "HLA-C"){
                        `$bin/samtools faidx $fadir/allele.1.$class.fasta HLA-C_0:200-3500  > $fadir/tmp/allele.1.$class.fasta`;
                        `$bin/samtools faidx $fadir/allele.2.$class.fasta HLA-C_1:200-3500  > $fadir/tmp/allele.2.$class.fasta`;
		}
		 elsif($class eq "MICB"){
                        `$bin/samtools faidx $fadir/allele.1.$class.fasta MICB_0:7000-10000  > $fadir/tmp/allele.1.$class.fasta`;
                        `$bin/samtools faidx $fadir/allele.2.$class.fasta MICB_1:7000-10000  > $fadir/tmp/allele.2.$class.fasta`; 
		}
		elsif($class eq "CYP4B1"){
		        `$bin/samtools faidx $fadir/allele.1.$class.fasta CYP4B1_0:1-18500  > $fadir/tmp/allele.1.$class.fasta`;
		        `$bin/samtools faidx $fadir/allele.2.$class.fasta CYP4B1_1:1-18500  > $fadir/tmp/allele.2.$class.fasta`; 
		}
		elsif($wxs ne "tgs" && ($class eq "CYP3A5"|| $class eq "CYP2D6" || $class eq "CYP2A13"|| $class eq "CYP2C8" || $class eq "CYP2C9")){
			my $exon_region = `less $db/CYP/CYP.ref.exon.bed |grep $class |xargs`;
			chomp $exon_region;
			my $vcf = "$fadir/$class.rephase.vcf.gz";
			if($wxs eq "tgs"){$vcf = "$fadir/$sample.$class.phased.vcf.gz"}
			`echo ">"$class"_0" >$fadir/tmp/allele.1.$class.fasta`;
		        `$bin/samtools faidx $db/CYP/ref/CYP.ref.fasta $exon_region | $bin/bcftools consensus -H 1 $fadir/$class.rephase.vcf.gz |grep -v ">CYP" >> $fadir/tmp/allele.1.$class.fasta`;
			`echo ">"$class"_1" >$fadir/tmp/allele.2.$class.fasta`;
                        `$bin/samtools faidx $db/CYP/ref/CYP.ref.fasta $exon_region | $bin/bcftools consensus -H 2 $fadir/$class.rephase.vcf.gz |grep -v ">CYP" >> $fadir/tmp/allele.2.$class.fasta`;
		       
		}
		else{
		    `cp $fadir/allele.1.$class.fasta $fadir/tmp`; 
		    `cp $fadir/allele.2.$class.fasta $fadir/tmp`;
		    }
        }
	for(my $i=1;$i<=2;$i++){
		my $fa = "$fadir/tmp/allele.$i.$class.fasta";
		my $tag = "$class"."_"."$i";
		#my $ref="$db/split/$class.fasta";
               `$bin/makeblastdb -in $fa -dbtype nucl -parse_seqids -out $fa`;
	       if($class eq "KIR2DL5" || $class eq "HLA-DRB1"){
                        `$bin/blastn -query $fa -out $fadir/tmp/$tag.blast.out1 -db $ref -outfmt 7 -num_threads 4 -max_target_seqs 1000 -gapopen 2 -gapextend 1 `;
			# `$bin/blastn -query $ref -out $fadir/tmp/$tag.blast.out2 -db $fa -outfmt 7 -num_threads 4 -max_target_seqs 1000 `;
			}
               elsif($wxs eq "tgs"){
                        `$bin/blastn -query $fa -out $fadir/tmp/$tag.blast.out1 -db $ref -outfmt 7 -num_threads 4 -max_target_seqs 2000 -perc_identity 90 -qcov_hsp_perc 10 -gapopen 2 -gapextend 1`;
			#`$bin/blastn -query $ref -out $fadir/tmp/$tag.blast.out2 -db $fa -outfmt 7 -num_threads 4 -max_target_seqs 2000 -perc_identity 90 -qcov_hsp_perc 10 -gapopen 2 -gapextend 1`;
			} 
               else{
	       `$bin/blastn -query $fa -out $fadir/tmp/$tag.blast.out1 -db $ref -outfmt 7 -num_threads 4 -max_target_seqs 2500 -perc_identity 95 -qcov_hsp_perc 10`;
	       #`$bin/blastn -query $ref -out $fadir/tmp/$tag.blast.out2 -db $fa -outfmt 7 -num_threads 4 -max_target_seqs 20 -perc_identity 95 -qcov_hsp_perc 10`;
	       }

	       my (%hash11,%hash12, %hash21, %hash22,$gene,$score,%hasha);
	       open IN1, "$fadir/tmp/$tag.blast.out1" or die "$!\n";
               while(<IN1>){
                       chomp;
                       next if(/^#/);
                       my ($hla, $t, $m,$d, $a) = (split)[1,3,4,5,11];
                       #next if($hla =~ /[N|Q]$/);
                       next if($t <$cutoff);
		       next if($hla eq "CYP1B1*REFGRCH38P7");
		       #next if($hla eq "HLA-A*24:23");
		       #next if($hla eq "HLA-A*26:26");
                       $hash11{$hla} += $t;
                       $hash12{$hla} += $m + $d * $bias;
		       $hasha{$hla} += $a;
              }
              close IN1;

	      # open IN2, "$fadir/tmp/$tag.blast.out2" or die "$!\n";
	      # while(<IN2>){
	      #        chomp;
	      #        next if(/^#/);
	      #        my ($hla, $t, $m,$d, $a) = (split)[0,3,4,5,11];
	      #        #next if($hla =~ /[N|Q]$/);
	      #        next if($t <$cutoff);
	      #        $hash21{$hla} += $t;
	      #        $hash22{$hla} += $m + $d;
	      #}
	      #close IN2;
              $score=85;
	      my $size = $cutoff;
              foreach my $key(sort keys %hash11){
		      #       next if(!exists $hash21{$key});
                      my @tt = (split /:/, $key);
                      my $kid = $tt[0];
		      if($#tt>0){$kid = "$tt[0]".":"."$tt[1]";}
		      #my $s1 = $hasha{$key};
		      my $s1 = 100 * (1 - $hash12{$key}/$hash11{$key});
		      #my $s2 = 100 * (1 - $hash22{$key}/$hash21{$key});
                      next if($hash11{$key} < $cutoff);
		      #next if($hash21{$key} < $cutoff);
		      my @ds = split /:/,$key;
		      my $d4;
		      if(exists $hashg{$class}){$d4 = "$ds[0]".":"."$ds[1]";}
		      next if(exists $hashg{$class} && $pop ne "nonuse" && !exists $hashp{$d4});
		      next if(exists $hashg{$class} && $pop ne "nonuse" && $hashp{$d4} ==0);
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

open OUT, ">$fadir/allele.phase.result.txt";
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
