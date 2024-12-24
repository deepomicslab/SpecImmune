#!/usr/bin/perl -w
use FindBin qw($Bin);
my ($ifile,$ofile) = @ARGV;
my $db="$Bin/../db";
my $snpfile = "$db/KIR/merge.KIR.snp.txt";
my %hash;
open IN, "$snpfile" or die "$!\n";
while(<IN>){
	chomp;
        my ($allele,$gene,$pos,$ref,$alt) = (split);
	my $key = "$gene\t$pos";
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

