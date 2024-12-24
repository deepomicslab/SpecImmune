#!/usr/bin/perl -w
use FindBin qw($Bin);
my ($ifile,$ofile) = @ARGV;
my $db="$Bin/../db";
my $snpfile = "$db/snp/merge.allele.snp.vcf";
my %hash;
open IN, "$snpfile" or die "$!\n";
while(<IN>){
	chomp;
        my ($allele,$pos,$ref,$alt) = (split /\t/,$_)[0,1,3,4];
	next if($allele eq "CYP2D6*82");
	my $gene = (split /\*/,$allele)[0];
	my $key = "$gene\t$pos";
        if(!$ref){$ref = "-"}
	if(!$alt){$alt = "-"}
	$hash{$key} = "$ref\t$alt";
}
close IN;

open OUT, ">$ofile";
open CIN, "$ifile" or die "$!\n";
while(<CIN>){
	chomp;
	if(/^#/){print OUT "$_\n"}
	else{   
		my $tag = "N";
		my ($gene,$pos,$ref,$alt) = (split)[0,1,3,4];
		my $len = length($ref);
		for(my $i=0;$i<$len;$i++){
			my $r = substr($ref, $i, 1);
			#my $a = substr($alt, $i, 1);
			my $p = $pos + $i;
		        my $k = "$gene\t$p";
		        if(exists $hash{$k}){$tag = "T";}
		}
		if($tag eq "T"){print OUT "$_\n"}
	}
}
close OUT;
close CIN;

