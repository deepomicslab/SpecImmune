#!/usr/bin/perl -w
use FindBin qw($Bin);
my $db="$Bin/../db/IG_TR";
my $bin="$Bin/../bin";

my ($vcf, $fa, $outdir) = @ARGV;
my (%hash, %hashf1, %hashf2);
`$bin/bcftools consensus -H 2 -f $fa $vcf > $outdir/phased2.fa`;
`$bin/bcftools consensus -H 1 -f $fa $vcf > $outdir/phased1.fa` ;

open FA, "$outdir/phased1.fa" or die "$!\n";
my $key;
while(<FA>){
	chomp;
	if(/^>/){s/^>//;$key=$_}
	else{$hashf1{$key} .= $_;}
}
close FA;
open FA, "$outdir/phased2.fa" or die "$!\n";
while(<FA>){
	chomp;
	if(/^>/){s/^>//;$key=$_}
	else{$hashf2{$key} .= $_;}
}
close FA;

open (IN, "gunzip -c  $vcf| ") or die "$!\n";
while(<IN>){
	chomp;
	next if(/^#/);
	my ($gene,$pos,$ref,$alt,$line) = (split /\t/, $_)[0,1,3,4,9];
	my $genotype = (split /:/, $line)[0];
	my $base = "$ref".","."$alt";
        my $value = "$pos\t$base\t$genotype\n";
	$genotype =~ m/(\d)[\||\/](\d)/;
	my ($g1,$g2) = ($1,$2);
	next if($g1 == $g2);
	$hash{$gene} += 1;
}
close IN;

open OUT, ">$outdir/VDJ.phased.fasta";
foreach my $gene(sort keys %hashf1){
	if(exists $hash{$gene}){
		my $fa1 = $hashf1{$gene};
		my $fa2 = $hashf2{$gene};
		my $len = length($fa1) - 200;
		#$fa1 = substr($fa1,100,$len);
		#$fa2 = substr($fa2,100,$len);
		my $id1 = "$gene".".1";
		my $id2 = "$gene".".2";
		print OUT ">$id1\n$fa1\n";
		print OUT ">$id2\n$fa2\n";
	}else{
		my $fa = $hashf1{$gene};
		my $len = length($fa) - 200;
		#$fa = substr($fa,100,$len);
		print OUT ">$gene\n$fa\n";
	}
}
close OUT;


#foreach my $key(keys %hash){
#	my @lines = (split /\n/, $hash{$key});
#	my $fa = $hashf{$key};
#	my @fas = split //, $fa;
#	my ($fa1, $fa2) = ($fa, $fa);
#	my (@fas1, @fas2) = (@fas, @fas);
#	foreach my $line(@lines){
#		my ($pos,$base,$genotype) = (split /\t/, $line)[0,1,2];
#		$genotype =~ m/(\d)[\||\/](\d)/;
#		my ($g1,$g2) = ($1,$2);
#		my @bases = (split /,/,$alt);
#		my $len = length($bases[0]);
#		my $prefix = substr($fa1, 0, $pos-2);
#		my $ref = substr($fa1, $pos-1, $len);
#		my $suffix = substr($fa1

		


