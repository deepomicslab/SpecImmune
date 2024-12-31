#!/usr/bin/perl -w
my ($ifile,$ofile) = @ARGV;
my (%hashs,%hash);
open IN, "$ifile" or die "$!\n";
while(<IN>){
	chomp;
	next if(/^#/);
	if(/^[V|D|J]/){
		my ($type,$id,$allele,$score) = (split)[0,1,2,-1];
		#my $gene = (split /\*/, $allele)[0];
		my $k = "$type\t$id";
		if(!exists $hash{$k}){$hash{$k} = $allele; $hashs{$k} = $score}
		elsif($score > $hashs{$k}){$hash{$k} = $allele; $hashs{$k} = $score}
		elsif($score == $hashs{$k}){$hash{$k} .= ",$allele";}
	}
}
close IN;

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
my (%hashv,%hashd,%hashj);
open OUT, ">$ofile";
foreach my $key(sort keys %hash){
	#print OUT "$key\t$hash{$key}\n";
	my ($type, $id) = (split /\t/,$key)[0,1];
	my @arr = (split /,/, $hash{$key});
	my @all = uniq(@arr);
	my $genes = join(";",@all);
        if($type eq "V"){$hashv{$id} = $genes;}
	if($type eq "D"){$hashd{$id} = $genes;}
	if($type eq "J"){$hashj{$id} = $genes;}
}

foreach my $id(sort keys %hashv){
	print OUT "$id\t$hashv{$id}\t$hashj{$id}";
	if(exists $hashd{$id}){
		print OUT "\t$hashd{$id}\n";}
	else{print OUT "\n";}
}
close OUT;
