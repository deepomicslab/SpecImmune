#!/usr/bin/perl -w
my ($file,$fa) = @ARGV;
open OUT, ">$fa";
open IN, "$file" or die "$!\n";
<IN>;
my $i = 1;
#my (%hash,%hashi);
while(<IN>){
	chomp;
	my ($re,$frame,$id,$seq) = (split /\t/, $_)[52,54,94,99];
	next if($frame ne "In");
	my $line = "Contig_"."$i";
	$i += 1;
	if($seq){print OUT ">$line\t$id\n$seq\n";}
	else{print OUT ">$line\t$id\n$re\n"}
	#$hash{$line} = $seq;
	#$hashi{$line} = $id;
}
close IN;
close OUT;


