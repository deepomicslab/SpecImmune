#!/usr/bin/perl -w
use FindBin qw($Bin);
my $db="$Bin/../db/IG_TR";
my $bin="$Bin/../bin";
my ($infile, $outpre, $tag) = @ARGV;
my $ref="$db/imgt_human_tr_ig_allele.fasta";
my (%hash,%hashh,%hashc,%hashn);
open (IN, "gunzip -c $infile |") or die "$!\n";
while(<IN>){
	chomp;
	if(/^#/){
		if(/^##contig/){
			s/##contig=<ID=(\S+),length=\d+>//;
			$gene = $1;
                        if($gene =~ /$tag/){
			        $hashc{$gene} = 1;}
		}
	}
	else{
		my ($gene, $qual, $info) = (split)[0,5,9];
	        next if($qual < 1);
		my $class = (split /\*/,$gene)[0];
		$hashn{$class}=1;
	        my $gt = (split /:/, $info)[0];
	        if($gt ne "1/1"){$hashh{$gene} += 1}
	        $hash{$gene} += 1;
        }
}
close IN;

open OUT, ">$outpre.fa";
foreach my $gene(sort keys %hash){
	my $class = substr($gene,3,1);
	my $fa = `$bin/samtools faidx $ref $gene | $bin/bcftools consensus -H 1 $infile |grep -v ">" `;
	$fa =~ s/\s//g;
        if($class eq "V"){$fa = substr($fa,80)}
	my $id = (split /\*/, $gene)[0];
	print OUT ">$id.1\n$fa\n";
	if(exists $hashh{$gene}){
		my $fa2 = `$bin/samtools faidx $ref $gene | $bin/bcftools consensus -H 2 $infile |grep -v ">" `;
		$fa2 =~ s/\s//g;
		if($class eq "V"){$fa2 = substr($fa2,80)}
		print OUT ">$id.2\n$fa2\n";
	}
}
close OUT;

`$bin/blastn -db db/imgt_human_tr_ig_uniq -query $outpre.fa -out $outpre.blast.out -outfmt 7`;

my (%hashs,%hashi);
open CI, "$outpre.blast.out" or die "$!\n";
while(<CI>){
      chomp;
      next if(/^#/);
      my ($id1,$id2,$score,$mis) = (split)[0,1,2,4];
      if(!exists $hashs{$id1}){$hashs{$id1} = $score; $hashi{$id1} = "$id2\t$score\t$mis";}
      elsif($score >= $hashs{$id1}){$hashs{$id1} = $score; $hashi{$id1} = "$id2\t$score\t$mis";}
}
close CI;
open OO, ">$outpre.blast.txt";
print OO "ID\tGene\tScore\tmismatches\n";
foreach my $ii(sort keys %hashi){
	print OO "$ii\t$hashi{$ii}\n";
}
foreach my $k(sort keys %hashc){
	my $kk = (split /\*/,$k)[0];
	next if(exists $hashn{$kk});
	print OO "$kk\t$k\t100\t0\n";
}
close OO;
