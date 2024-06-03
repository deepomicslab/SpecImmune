#!/usr/bin/perl -w
my ($sample,$outdir) = @ARGV;
my (%hash, %hashc, %hashb,%hashr, %hash1,%hash2);
my $db_fasta = "../db/whole/HLA.full.fasta";
my $hap_txt = "./HLA.ref.hap.txt";
my $hap_fa ="/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/down/HLA.ref.hap.fa";
open IN, $db_fasta or die "$!\n";
# while(<IN>){
# 	chomp;
# 	s/^>//;
# 	my $id = $_;
# 	my $fa = <IN>;
# 	chomp $fa;
# 	### used to extract the middle region
# 	# my $ll = length($fa);
# 	# my $e = $ll - 600;
# 	# my $ffa = substr($fa,299,$e);
# 	my $ffa = $fa;  ## no need to extract the middle region
# 	$hash{$id} = $ffa;
# 	my $gene = (split /\*/,$id)[0];
# 	if($gene =~ /KIR2DL5/){$gene = "KIR2DL5"}
# 	$hashc{$gene} .= "$id\t";
# }
my $id = '';
while(<IN>){
    chomp;
    if (/^>/) {
        $id = substr($_, 1);  # remove '>'
        $hash{$id} = '';  # initialize sequence
        my $gene = (split /\*/, $id)[0];
        if($gene =~ /KIR2DL5/){$gene = "KIR2DL5"}
        $hashc{$gene} .= "$id\t";
    } else {
        $hash{$id} .= $_;  # append to sequence
    }
}

close IN;
open COUT, ">$outdir/$sample.HLA.hap.alleles.txt";
print COUT "Gene\tHap1\tHap2\n";
foreach my $gene(sort keys %hashc){
	my @alls = (split /\t/,$hashc{$gene});
	my $s = $#alls;
        my $i = int(rand($s));
        my $j = int(rand($s));
	my $a1 = $alls[$i];
	my $a2 = $alls[$j];
	$hash1{$gene} = $a1;
	$hash2{$gene} = $a2;
	print COUT "$gene\t$a1\t$a2\n";
}
close COUT;

`less $hap_txt |sort -k 3 -n >HLA.ref.hap.sort.txt`;
my @regions;
open CI, "HLA.ref.hap.sort.txt" or die "$!\n";
while(<CI>){
	chomp;
	my ($gene,$chr,$s,$e,$flag) = (split)[0,1,2,3,4];
	$s += 1;
	$e += 1;
	my $region = "$chr".":"."$s"."-"."$e";
	$hashb{$region} = "$gene\t$flag";
	push @regions, $region;
}
close CI;
open FA, $hap_fa or die "$!\n";
open SEP, ">$outdir/$sample.HLA.sep.fa";
<FA>;
my $ref=<FA>;
close FA;
chomp $ref;
my ($nref1,$nref2);
my ($ss,$end) = (1,1);
foreach my $re(@regions){
	my ($gene,$flag) = (split /\t/,$hashb{$re})[0,1];
	my ($chr,$se) = (split /:/,$re)[0,1];
	my ($s,$e) = (split /-/,$se)[0,1];
	my $ee = $s - 1;
	my $pre = "$chr".":"."$ss"."-"."$ee";
	my $pfa = `samtools faidx $hap_fa $pre |grep -v ">" `;
	chomp $pfa;
	$pfa =~ s/\s//g;
	my $allele1 = $hash1{$gene};
	my $allele2 = $hash2{$gene};
	my $fa1 = $hash{$allele1};
	my $fa2 = $hash{$allele2};

	print SEP ">$allele1\n";
	print SEP "$fa1\n";
	print SEP ">$allele2\n";
	print SEP "$fa2\n";

	if($flag == 16){
		$fa1 = uc $fa1;
                $fa1 = reverse $fa1;
                $fa1 =~ tr/ATCG/TAGC/;
		$fa2 = uc $fa2;
                $fa2 = reverse $fa2;
                $fa2 =~ tr/ATCG/TAGC/;
	}
	
	$nref1 .= $pfa; $nref1 .= $fa1;
	$nref2 .= $pfa; $nref2 .= $fa2;
	$end = $e;
	$ss = $e + 1;
	#my $l1 = length($pfa); my $l2 = length($fa1); my $l3 = length($fa2);
	#print "$pre\t$l1\t$re\t$l2\t$l3\n";
}
close SEP;


my $len = length($ref);
$end += 1;
my $region0 = "HLA.hap:"."$end"."-"."$len";
my $pe = `samtools faidx $hap_fa $region0 | grep -v ">"`;
$pe =~ s/\s//g;
$nref1 .= $pe;
$nref2 .= $pe;
#my $ll = length($pe);
#my $la1 = length($nref1);
#my $la2 = length($nref2);
#print "$region0\t$ll\nhap1\t$la1\thap2\t$la2\n";
open O1, ">$outdir/$sample.HLA.hap1.fa";
open O2, ">$outdir/$sample.HLA.hap2.fa";
print O1 ">HLA.hap1\n$nref1\n";
print O2 ">HLA.hap2\n$nref2\n";
close O1;
close O2;


