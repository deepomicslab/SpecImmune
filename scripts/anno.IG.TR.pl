#!/usr/bin/perl -w
my ($sample,$hap1,$hap2,$outdir,$db) = @ARGV;
# my $db ="$db_dir/IG_TR/IG.TR.merge.allele.fasta";
print $db;
# my @chrs = ("chr14_igh", "chr15_igh", "chr16_igh", "chr2_igk", "chr22_igl",  "chr7_alt_trb", "chr7_trb", "chr9_trb", "chr14_tra", "chr7_trg");
# my @chr7_alt_trbs = ("TRBV5-3", "TRBV5-7", "TRBV5-8", "TRBV6-4", "TRBV6-8", "TRBV6-9", "TRBV7-3", "TRBV7-7", "TRBV7-8", "TRBV7-9", "TRBV9");
# my @chrs = ("chr14");
my @chrs = ("chr1","chr2","chr7","chr14","chr22");
my @chr7_alt_trbs = ();
my %hasht;
foreach my $gg(@chr7_alt_trbs){
     $hasht{$gg} = $gg;
}
open OUT, ">$outdir/$sample.IG.TR.allele.txt";
print OUT "sample\tgene\tallele\tscore\tlen\tstart\tend\tchr\thap\n";
foreach my $chr(@chrs){
	`samtools faidx $hap1`;
	`samtools faidx $hap2`;
	`samtools faidx $hap1 $chr > $outdir/tmp/$sample.$chr.hap1.fasta`;
	`samtools faidx $hap2 $chr > $outdir/tmp/$sample.$chr.hap2.fasta`;
        for(my $i=1;$i<=2;$i++){
	    `blastn -query $outdir/tmp/$sample.$chr.hap$i.fasta -out $outdir/tmp/$sample.$chr.hap$i.blast.txt -db $db -outfmt 7 -max_target_seqs 3000 -num_threads 4 `;
            my (%hashg, %hashp, %hashs, %hashl);
	    open IN, "$outdir/tmp/$sample.$chr.hap$i.blast.txt" or die "$!\n";
	    while(<IN>){
                    chomp;
                    next if(/^#/);
                    my ($allele,$identity,$len,$mis,$gap,$start,$end,$score) = (split)[1,2,3,4,5,6,7,11];
		    #next if($mis >= 5);
		    my $gene = (split /\*/, $allele)[0];
		    next if($identity < 97 && $gene =~ /TRB/);
                    next if($identity < 98 && $gene =~ /TRA/);
		    next if($identity < 97 && $gene =~ /IGK/);
		    next if($identity < 97 && $gene =~ /IGL/);
		    next if($identity < 98 && $gene =~ /IGH/ && $gene ne "IGHV1-2");
		    next if($allele =~ /IGKV1\/ORY-1/);
		    next if($allele =~ /IGKV1-NL1/);
		    next if($chr eq "chr22_igk" && $start <1000);

		    next if($chr eq "chr7_alt_trb" && !exists $hasht{$gene});
		    next if($chr eq "chr7_trb" && exists $hasht{$gene});
		    next if($allele =~ /IGHV\S+D/);
		    next if($allele =~ /IGHV4-30-4/);
		    next if($allele =~ /IGKV1\/OR15-118/);
		    next if($allele =~ /IGKV1\/OR1-1/);
		    next if($allele =~ /IGKV1\/OR2-118/);
		    next if($allele =~ /IGKV1\/OR10-1/);
		    next if($allele =~ /IGKV1\/OR-2/);
		    next if($allele =~ /IGKV1\/OR9-2/);
		    next if($allele =~ /IGKV1\/OR-4/);
		    next if($allele =~ /IGKV1\/OR22-5/);
		    next if($allele =~ /IGKV1\/OR9-1/);
		    # next if($allele =~ /TRD/ && $chr eq "chr14_tra");
		    #next if($allele =~ /TRA/ && $chr eq "chr14_trd");
		    next if($allele =~ /IGHV/ && $len < 150);
		    next if($allele =~ /TR[A|B|D|G]V/ && $len < 100);
		    next if($allele =~ /IGHV\d+\/OR/ && $chr eq "chr14_igh");
                    next if($chr eq "chr15_igh" && ! ($allele =~ /IGHV\d+\/OR15/));
		    next if($chr eq "chr16_igh" && ! ($allele =~ /IGHV\d+\/OR16/));
		    
                    if(!exists $hashg{$gene}){
                            $hashg{$gene} = $allele;
                            $hashs{$gene} = $identity;
                            $hashl{$gene} = $len;
                            $hashp{$gene} = "$start\t$end";
                    }else{
                            if($len >= $hashl{$gene} - 5 && $identity > $hashs{$gene}){
                                    $hashg{$gene} = $allele;
                                    $hashs{$gene} = $identity;
                                    $hashl{$gene} = $len;
                                    $hashp{$gene} = "$start\t$end";
                            }
                    }
            }
            close IN;
            foreach my $g(sort keys %hashg){
                    print OUT "$sample\t$g\t$hashg{$g}\t$hashs{$g}\t$hashl{$g}\t$hashp{$g}\t$chr\thap$i\n"
            }

    }

}
close OUT;