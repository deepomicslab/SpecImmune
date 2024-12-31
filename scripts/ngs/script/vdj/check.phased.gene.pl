#!/usr/bin/perl -w
#perl $0 ESCC001N/read.vdj.file.txt ESCC001N/VDJ.phased.fasta ESCC001N/check.phase.txt
use FindBin qw($Bin);
my $db="$Bin/../db/IG_TR";
my $bin="$Bin/../bin";

my ($vfile, $sfile, $out) = @ARGV;
`rm -rf $out`;
my (%hash, %hashs, %hashv);
my $key;
open FA, "$db/IMGT.VDJ.fasta" or die "$!\n";
while(<FA>){
	chomp;
	if(/^>/){
	      s/^>//;
	      $key = $_;
	}else{
	      $hash{$key} .= $_
	}
}
close FA;

open IN, "$vfile" or die "$!\n";
while(<IN>){
	chomp;
	my ($class, $genes) = split;
	#if($genes =~ m/,/){
              $hashv{$class} = "$genes";
	      #}
}
close IN;

open SI, "$sfile" or die "$!\n";
while(<SI>){
	chomp;
	if(/^>/){
		s/^>//;
		$key = $_;
	}else{
		$hashs{$key} .= $_;
	}
}
close SI;

foreach my $class(sort keys %hashv){
	next if($class =~ /\//);
	my ($gene1, $gene2) = ($hashv{$class}, $hashv{$class});
	if($hashv{$class} =~ /,/){($gene1,$gene2) = (split /,/, $hashv{$class})[0,1];}
	my $phase1 = "$gene1".".1";
	my $phase2 = "$gene1".".2";
	
	my $phase_out = "tmp/$class.phased.fa";
	my $select_out = "tmp/$class.select.fa";
	open O1, ">$phase_out";
	open O2, ">$select_out";
	if(exists $hashs{$phase1}){
                print O1 ">$phase1\n$hashs{$phase1}\n";
		print O1 ">$phase2\n$hashs{$phase2}\n";
	}else{
		print O1 ">$phase1\n$hashs{$gene1}\n";
		print O1 ">$phase2\n$hashs{$gene1}\n";
	}
	if($gene1 eq $gene2){
		print O2 ">$gene1\n$hash{$gene1}\n";
	}else{

		print O2 ">$gene1\n$hash{$gene1}\n";
		print O2 ">$gene2\n$hash{$gene2}\n";
	}
	        close O1;
		close O2;
		`$bin/makeblastdb -in $select_out -dbtype nucl -parse_seqids -out $select_out`;
                `$bin/blastn -query $phase_out -out $phase_out.out -db $select_out -outfmt 7`;
		`less $phase_out.out|grep -v "#"|sort -k 3 -n -r  >> $out`;
		#open BL, "$phase_out.out" or die "$!\n";
		#while(<BL>){
		#	chomp;
		#	next if(/^#/);
		#	my ($a1,$a2,$score,$num) = (split)[0,1,2,3];
		#	next if($a1 eq $a2);
			
		#}
		#close BL;

		#}else{
	
		#	print "$class\t$hashv{$class}\n";
		#}

}

open CI, "$out" or die "$!\n";
my (%hash1,%hash2,%hashs1,%hashs2);
open OUT, ">$out.report";
while(<CI>){
	chomp;
	my ($phase,$allele,$score) = (split)[0,1,2];
	my ($gene,$type) = (split /\./,$phase)[0,1];
	$gene = (split /\*/,$gene)[0];
	if($type == 1){
		if(!exists $hash1{$gene}){$hash1{$gene} = $allele; $hashs1{$gene} = $score}
		elsif($hashs1{$gene} < $score){$hash1{$gene} = $allele; $hashs1{$gene} = $score}
	}
	if($type == 2){
		if(!exists $hash2{$gene}){$hash2{$gene} = $allele; $hashs2{$gene} = $score}
		elsif($hashs2{$gene} < $score){$hash2{$gene} = $allele; $hashs2{$gene} = $score}
	}
}
close CI;
print OUT "Gene\tAllele1\tAllele2\tScore1\tScore2\n";
foreach my $gene(sort keys %hash1){
	print OUT "$gene\t$hash1{$gene}\t$hash2{$gene}\t$hashs1{$gene}\t$hashs2{$gene}\n";
}
close OUT;
