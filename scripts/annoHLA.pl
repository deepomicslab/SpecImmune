#!/usr/bin/perl 
use FindBin qw($Bin);
use Getopt::Long;

my ($sample, $dir, $pop, $wxs, $g_nom, $help);
$g_nom=0;
$rna=1;
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
die $usage unless ($sample && $dir && $pop && $wxs );
print "parameter:\tsample:$sample\tdir:$dir\tpop:$pop\twxs:$wxs\tG_nom:$g_nom\tdb:$db\n";

my $version="";
my $k = 2;
my $bias = 1;
if($wxs eq "tgs"){$bias = 1} # gap score in TGS
my (%hashp,%hashp2, %hashpp, %hashg, %hashc, %hash,%hashdd);
# my $db="$Bin/../db/HLA";
my $bin="$Bin/../bin";
my @hlas = (
    'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DPB2', 
    'HLA-DQA1', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-E', 'HLA-F', 'HLA-G', 
    'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L', 'HLA-P', 'HLA-V', 'HLA-DQA2', 'HLA-DPA2', 'HLA-N', 'HLA-S', 'HLA-T', 'HLA-U', 
    'HLA-W', 'MICA', 'MICB', 'TAP1', 'TAP2', 'HFE'
);
my $fadir=$dir;
my $workdir = "$dir/tmp";
`mkdir  -p $workdir`;
# if is rna , subdir=HLA_CDS, else subdir=HLA
my $subdir = "HLA";
if($rna == 1){$subdir = "HLA_CDS"}

sub tgs_blast{
        open BOUT, ">$dir/HLA.blast.summary.txt";  # new 
	foreach my $class(@hlas){
		my $ref = "$db/$subdir/$class/$class";
              print "blast $ref\n";
		for(my $i=1;$i<=$k;$i++){
			my $tag = "$class"."_"."$i";
			my $fa = "$fadir/HLA.allele.$i.$class.fasta";
                     # if fa file not exist or empty, skip
                     if(!-e $fa || -z $fa){next}
			system("blastn -query $fa -out $fadir/tmp/$tag.blast.out1 -db $ref -outfmt 7 -num_threads 4 -max_target_seqs 10000000 -perc_identity 95 -qcov_hsp_perc 10");
			my (%hash_max,%hash11,%hash12, %hash21, %hash22,$gene,$score);
			open IN1, "$workdir/$tag.blast.out1" or die "$!\n";
			while(<IN1>){
				chomp;
				next if(/^#/);
				my ($hla, $t, $m,$d,$s) = (split)[1,3,4,5,8];
				$hash11{$hla} += $t;
				$hash12{$hla} += $m + $d * $bias;
			}
			close IN1;
			$score=85;
                     my $ff=0;
                     foreach my $key(sort keys %hash11){
                            my @tt = (split /:/, $key);
                            my $kid = "$tt[0]".":"."$tt[1]";
                     my ($fre,$fre2)=(0,0);
                     if(exists $hashp{$kid}){$fre=$hashp{$kid};$fre2=$hashp2{$kid}} ## population frequency of 4 digit hla allele
                     my $s = 100 * (1 - $hash12{$key}/$hash11{$key}); #blast score
                     print BOUT "$tag\t$key\t$hash12{$key}\t$hash11{$key}\t$s\n";  # new 
                     next if($pop ne "nonuse" && $fre2 == 0);
                     my $scorel=$s;

                     if($scorel >= $score){
                            $score = $scorel;
                            $gene = $key;
                            $ff=$fre;
                            $hash_max{$scorel} .= "$gene;$s\t";
                     }
             		}
             		$hash{$tag} = $hash_max{$score};
             		`rm -rf $fadir/$class.temp*`;
        }
    }
}

open COUT, ">$dir/hla.result.txt";
print COUT "$version\n";
print COUT "Sample\tHLA_A_1\tHLA_A_2\tHLA_B_1\tHLA_B_2\tHLA_C_1\tHLA_C_2\tHLA_DPA1_1\tHLA_DPA1_2\tHLA_DPB1_1\tHLA_DPB1_2\tHLA_DQA1_1\tHLA_DQA1_2\tHLA_DQB1_1\tHLA_DQB1_2\tHLA_DRB1_1\tHLA_DRB1_2\n";

open OUT, ">$dir/hla.result.details.txt";
print OUT "$version\n";
print OUT "Gene\tG_best\tallele\tdetails:allele;Score;Caucasian;Black;Asian\n";
if($wxs eq "exon"){
       &exon_blast; 
}
if($wxs eq "whole"){
       &whole_blast;
}
if($wxs eq "tgs"){
       &tgs_blast;
}
sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}

#print the alleles annotation of best score for each HLA gene
my $hout = $sample;
foreach my $hla(@hlas){
       for(my $i=1;$i<=$k;$i++){
             my $id = "$hla"."_"."$i";
             my ($line1,$line2,$line3,$pout,$out) = ("","","","","-");
             my (%ggs, %hashmm);
             if(!exists $hash{$id}){ $hout .= "\t-"; print OUT "$id\t-\t-\t-\n"}
             next if(!exists $hash{$id});
             my @arrs = (split /\t/,$hash{$id});
             my ($max,$pfre) = (0,0);
             foreach my $oo(@arrs){
                 my ($allele,$score) = (split /;/,$oo)[0,1];
                 $score = sprintf "%.3f", $score;
                 #DRB1*14:01 and DRB1*14:54 differ in HLA_DRB1:9519
                 if($allele =~ /DRB1\*14:01/){
                          system("samtools  mpileup -r HLA-DRB1:9519-9519 -t DP -t SP -uvf $db/hla.ref.extend.fa $dir/$sample.realign.sort.bam --output $workdir/snp.vcf");
                          open TE, "$workdir/snp.vcf" or die "$!\n";
                          while(<TE>){
                                 chomp;
                                 next if(/^#/);
                                 my $alt = (split)[4];
                                 if($alt =~ /T/){print "$allele\n"} else{$allele = "DRB1*14:54";}
                          }
                          close TE;
                 }
                 #C*07:01 and C*07:18 differ in HLA_C:4061
                 if($allele =~ /C\*07:01/ && $wxs eq "exon"){
                          system("samtools  mpileup -r HLA-C:4061-4061 -t DP -t SP -uvf $db/hla.ref.extend.fa $dir/$sample.realign.sort.bam --output $workdir/snp.vcf");
                          open TE, "$workdir/snp.vcf" or die "$!\n";
                          while(<TE>){
                                 chomp;
                                 next if(/^#/);
                                 my $alt = (split)[4];
                                 if($alt =~ /T/){$allele = "C*07:18:01:01";}
                          }
                          close TE;
                 }

                 my @tt = (split /:/, $allele);
                 my $kid = "$tt[0]".":"."$tt[1]";
                 #$hashmm{$kid} .= "$allele\t";
                 $line2 .= "$allele;";
                 $oo = "$kid".";"."$score";
                 if(exists $hashg{$allele}){
                      my @ttt = (split /:/,$hashg{$allele});
                      my $tid = "$ttt[0]".":"."$ttt[1]";
                      $ggs{$tid} += 1;
                      $hashmm{$tid} .= "$allele\t";
                 }
                 else{$ggs{$kid} += 1; $hashmm{$kid} .= "$allele\t"}
                 if(exists $hashpp{$kid}){$line3 .= "$oo;$hashpp{$kid}\t"; 
                     if($pfre < $hashp{$kid}){$pfre = $hashp{$kid};$pout = $allele;
                     }
                 }
                 else{$line3 .= "$oo;-;-;-\t";} 
             }
             my @lines3 = (split /\t/,$line3); @lines3 = uniq(@lines3);
             my @lines2 = (split /;/,$line2); $line3 = join("\t",@lines3);
             if($pop eq "nonuse" || !$pout){
               foreach my $gg(sort {$ggs{$b} <=> $ggs{$a}} keys %ggs ){
                   my $agg = (split /\t/, $hashmm{$gg})[0];
                   if($ggs{$gg} >= $max){
                       $max = $ggs{$gg};$line1 .= "$agg;";
                        if($out eq "-"){$out = $agg;}
                   }
                   if($g_nom == 1 && exists $hashg{$out}){$out = $hashg{$out}}
                }
             }else{ 
                   if(exists $hashg{$pout}){$pout = $hashg{$pout};}
                   my $temp1 = "$pout".":01";
		   my $temp2 = "$pout".":01:01";
		   if(exists $hashg{$temp1}){$pout = $hashg{$temp1};}
	           if(exists $hashg{$temp2}){$pout = $hashg{$temp2};}
             #      else{$pout = $lines2[0]}
                   $out = $pout;   $line1 = $pout; 
             }
             if(!$out){$out="-"}
             $hout .= "\t$out";
             print OUT "$id\t$line1\t$line2\t$line3\n";
       }    
}
close OUT;
print COUT "$hout\n";
close COUT;
