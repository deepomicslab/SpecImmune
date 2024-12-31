#!/usr/bin/perl -w
use FindBin qw($Bin);
my $db="$Bin/../db/IG_TR";
my $bin="$Bin/../bin";
my ($vcf,$sample) = @ARGV;
my $db ="$db/imgtrefseq.human";
my $ref="$db/ig.tr.merge.hg38.fa";
`$bin/bgzip -f $vcf`;
`$bin/tabix -f $vcf.gz`;
`$bin/samtools faidx $ref`;
my %hash;
open LI, "$ref.fai" or die "$!\n";
open OUT, ">$sample.allele.report.out.txt";
print OUT "gene\tallele\tscore\tlen\tstart\tend\tchr\n";
while(<LI>){
	chomp;
	my ($chr,$len) = (split)[0,1];
	my $region = "$chr".":1-"."$len";
        for(my $i=1;$i<=2;$i++){
	    `$bin/samtools faidx $ref $region | $bin/bcftools consensus -H $i $vcf.gz > $chr.$i.fa`;
	    `$bin/blastn -query $chr.$i.fa -out $chr.$i.blast.txt -db $db -outfmt 7 -max_target_seqs 3000 -num_threads 4 `;
	    my (%hashg, %hashp, %hashs, %hashl);
	    open IN, "$chr.$i.blast.txt" or die "$!\n";
	    while(<IN>){
		    chomp;
		    next if(/^#/);
		    my ($allele,$identity,$len,$mis,$gap,$start,$end,$score) = (split)[1,2,3,4,5,6,7,11];
		    next if($mis > 5);
		    next if($allele =~ /IGHV/ && $len < 150);
		    my $gene = (split /\*/, $allele)[0];
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
		    print OUT "$g\t$hashg{$g}\t$hashs{$g}\t$hashl{$g}\t$hashp{$g}\t$chr\tP$i\n"
	    }

    }          
	
}
close LI;

