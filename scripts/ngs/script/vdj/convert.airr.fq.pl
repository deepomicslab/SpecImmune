#!/usr/bin/perl -w

my ($infile,$outfile) = @ARGV;
open OUT, ">$outfile";
open IN, "$infile" or die "$!\n";
<IN>;
while(<IN>){
	chomp;
	my ($id,$seq) = (split)[2,3];
	my $qual = $seq;
	$qual =~ s/[A|T|C|G|N]/F/g;
	print OUT "\@$id\n$seq\n+\n$qual\n";
}
close IN;
close OUT;

