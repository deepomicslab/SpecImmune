
use FindBin qw($Bin);
my ($snpfile,$svfile,$outfile) = @ARGV;
# my @chrs = ("chr7_trb","chr9_trb","chr14_tra","chr7_trg","chr7_alt_trb","chr14_igh","chr15_igh","chr16_igh","chr2_igk","chr22_igl");
my @chrs = ("chr1","chr2","chr7","chr14","chr22");
my (%hash,%hash2);
open SV, "$svfile" or die "$!\n";
while(<SV>){
	chomp;
	next if(/^#/);
	my ($chr,$start,$id, $ref,$alt,$filter,$info,$format,$aa) = (split /\t/,$_)[0,1,2,3,4,6,7,8,9];
        my $tag = 1; 
        next if($filter ne "PASS");
	$info =~ s/IMPRECISE;//;
	next if($chr eq "chr2_igk" && $start >2042500 && $start < 2042700);
	my @arrs = split /;/,$info;
        my $type = (split /=/,$arrs[0])[1];
        my $end = (split /=/,$arrs[1])[1];
        my $len = (split /=/,$arrs[2])[1];
        next if($type ne "DEL");
        $len =~ s/-//;
        next if($len < 500);
        print "$chr\t$start\t$end\t$len\n";
	#my $key = "$chr\t$start";
        my $format1 = "GT:DP:AD:RO:QR:AO:QA";
        my @tes = (split /:/, $aa);
	my $gt = $tes[0];
	my $ad = $tes[1];
	my $dp = $tes[2];
	my $ro = (split /,/,$ad)[0];
	my $ao = (split /,/,$ad)[1];
	my $qr = $ro * 5;
	my $qa = $ao * 5;
	my $outf = "$gt".":"."$dp".":"."$ad".":"."$ro".":"."$qr".":"."$ao".":"."$qa";
        my $info1 = $info;
	my $line = "$chr\t$start\t.\t$ref\t$alt\t100\tPASS\t$info1\t$format1\t$outf";
	$hash{$chr}{$start} = $line;
	if($gt =~ /0/){$tag = 0}
	for(my $i=$start;$i<=$end;$i++){
		my $key = "$chr\t$i";
		$hash2{$key} = $tag;
	}
        
}
close SV;

open IN, "$snpfile" or die "$!\n";
`cat $Bin/head.vcf >$outfile`;
open OUT, ">>$outfile";
while(<IN>){
	chomp;
	next if(/^##/);
	if(/^#/){
		print OUT "$_\n"}
	else{
		my @arrs = split;
		my $outf = $arrs[-1];
		my $format = $arrs[-2];
		my $info = $arrs[-3];
		my @infos = split /;/,$info;
		my $dp = (split /=/,$infos[0])[1];
		my $ad = (split /=/,$infos[1])[1];
		my $format1 = "GT:DP:AD:RO:QR:AO:QA";
        	my $gt = (split /:/, $outf)[0];
                my $ro = (split /,/,$ad)[0];
       	        my $ao = (split /,/,$ad)[1];
        
		my $qr = $ro * 5;
                my $qa = $ao * 5;
		#my $outf1 = "$gt".":"."$dp".":"."$ad".":"."$ro".":"."$qr".":"."$ao".":"."$qa";
                my $info1 = $info;
		#$arrs[-3]=$info1; $arrs[-2] = $format1; $arrs[-1] = $outf1;
		#my $line = join("\t",@arrs);
		my $key = "$arrs[0]\t$arrs[1]";
		next if(exists $hash2{$key} && $hash2{$key} == 1);
		if(exists $hash2{$key} && $hash2{$key} == 0){
			if($ao>$ro){$gt = "1/1"}else{$gt="0/0"}
		}
		my $outf1 = "$gt".":"."$dp".":"."$ad".":"."$ro".":"."$qr".":"."$ao".":"."$qa";
		$arrs[-3]=$info1; $arrs[-2] = $format1; $arrs[-1] = $outf1;
		my $line = join("\t",@arrs);
		$hash{$arrs[0]}{$arrs[1]} = $line;
                
	}
}
close IN;

foreach my $chr(@chrs){
	next if(!exists $hash{$chr});
	my $hashs = $hash{$chr};
	foreach my $pos(sort {$a <=> $b} keys %$hashs){
		print OUT "$hash{$chr}{$pos}\n";
	}
}
close OUT;
