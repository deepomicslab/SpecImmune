#!/usr/bin/perl -w
my ($sample,$outdir) = @ARGV;
my (%hash, %hashc, %hashb,%hashr, %hash1,%hash2, %hashf, %hasho);
open IN, "/home/wangmengyao/SpecComplex/db/IG_TR/imgtrefseq.human.fasta" or die "$!\n";
while(<IN>){
	chomp;
	s/^>//;
	my $id = $_;
	my $fa = <IN>;
	chomp $fa;
	$hash{$id} = $fa;
	my $gene = (split /\*/,$id)[0];
	$hashc{$gene} .= "$id\t";
}
close IN;
open COUT, ">$outdir/$sample.IG.TR.hap.alleles.txt";
print COUT "Gene\tChr\tHap1\tHap2\n";
foreach my $gene(sort keys %hashc){
	my @alls = (split /\t/,$hashc{$gene});
	my $s = $#alls;
        my $i = int(rand($s));
        my $j = int(rand($s));
	my $a1 = $alls[$i];
	my $a2 = $alls[$j];
	$hash1{$gene} = $a1;
	$hash2{$gene} = $a2;
	$hasho{$gene} = "$a1\t$a2";
}


`less IG.TR.ref.hap.txt |sort -k 3 -n >IG.TR.ref.hap.sort.txt`;
my @regions;
open CI, "IG.TR.ref.hap.sort.txt" or die "$!\n";
while(<CI>){
	chomp;
	my ($gene,$chr,$s,$e,$flag) = (split)[0,1,2,3,4];
	print COUT "$gene\t$chr\t$hasho{$gene}\n";
	$s += 1;
	$e += 1;
	my $region = "$chr".":"."$s"."-"."$e";
	$hashb{$region} = "$gene\t$flag";
	push @regions, $region;
}
close CI;
close COUT;
my $id;
open FA, "/mnt/delta_WS_1/wangmengyao/Complex/ig.tr.merge.hg38.fa" or die "$!\n";
while(<FA>){
	chomp;
	if(/^>/){
	    s/^>//;
	    $id = $_;
        }else{$hashf{$id} .= $_}
}
close FA;

my (%nref1,%nref2,%hashe);
my ($ss,$end) = (1,1);
foreach my $re(@regions){
	my ($gene,$flag) = (split /\t/,$hashb{$re})[0,1];
	my ($chr,$se) = (split /:/,$re)[0,1];
	my ($s,$e) = (split /-/,$se)[0,1];
	my $ee = $s - 1;
	my $pre = "$chr".":"."$ss"."-"."$ee";
	my $pfa = `samtools faidx ig.tr.merge.hg38.fa $pre |grep -v ">" `;
	chomp $pfa;
	$pfa =~ s/\s//g;
	my $allele1 = $hash1{$gene};
	my $allele2 = $hash2{$gene};
	my $fa1 = $hash{$allele1};
	my $fa2 = $hash{$allele2};
	if($flag == 16){
		$fa1 = uc $fa1;
                $fa1 = reverse $fa1;
                $fa1 =~ tr/ATCG/TAGC/;
		$fa2 = uc $fa2;
                $fa2 = reverse $fa2;
                $fa2 =~ tr/ATCG/TAGC/;
	}
	$nref1{$chr} .= $pfa; $nref1{$chr} .= $fa1;
	$nref2{$chr} .= $pfa; $nref2{$chr} .= $fa2;
	$end = $e;
	$ss = $e + 1;
	$hashe{$chr} = $end;
	#my $l1 = length($pfa); my $l2 = length($fa1); my $l3 = length($fa2);
	#print "$pre\t$l1\t$re\t$l2\t$l3\n";
}
open O1, ">$outdir/$sample.IG.TR.hap1.fa";
open O2, ">$outdir/$sample.IG.TR.hap2.fa";
foreach my $cc(sort keys %hashe){
	my $end0 = $hashe{$cc};
	my $ref = $hashf{$cc};
        my $len = length($ref);
        $end0 += 1;
        my $region0 = "Gene.hap:"."$end"."-"."$len";
        my $pe = `samtools faidx ig.tr.merge.hg38.fa $region0 | grep -v ">"`;
        $pe =~ s/\s//g;
        $nref1{$cc} .= $pe;
        $nref2{$cc} .= $pe;

        print O1 ">$cc.hap1\n$nref1{$cc}\n";
        print O2 ">$cc.hap2\n$nref2{$cc}\n";
}
close O1;
close O2;


