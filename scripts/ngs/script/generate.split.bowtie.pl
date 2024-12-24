#!/usr/bin/perl
use FindBin qw($Bin);
#perl $0 KIR fq1.list fq2.list output
my ($class,$fq1,$fq2,$outdir) = @ARGV;
my (%hashr, %hash1, %hash2);

my $db="$Bin/../db";
my $bin="$Bin/../bin";
#my $bowtie2 = "$bin/bowtie2";
my $samtools = "$bin/samtools";
my $freebayes = "$bin/freebayes";
my $group = "\@RG\tID:$class\tSM:$class";
`mkdir $outdir/tmp`;
open RE, "$db/KIR/ref/KIR.ref.extend.fasta" or die "$!\n";
while(<RE>){
	chomp;
	s/^>//;
        my $id = $_;
        my $fa = <RE>;
        chomp $fa;	
	next if($id eq "KIR2DL5B");
        if($id =~ /KIR2DL5/){$id = "KIR2DL5"}
	$hashr{$id} = $fa;	
}
close RE;

open IN1, "$fq1" or die "$!\n";
while(<IN1>){
	chomp;
	s/\*$//;
	my $gene = (split /\./, (split /\//, $_)[-1])[0];
	$hash1{$gene} = $_;
}
close IN1;

open IN2, "$fq2" or die "$!\n";
while(<IN2>){
	chomp;
	s/\*$//;
	my $gene = (split /\./, (split /\//, $_)[-1])[0];
	$hash2{$gene} = $_;
}
close IN2;

open OUT, ">$outdir/$class.ref.merge.fa";
foreach my $gene(sort keys %hash1){
	my $fq1 = $hash1{$gene};
	my $fq2 = $hash2{$gene};
	my $ref = $hashr{$gene};
	$id = $gene;
	open TE, ">$outdir/tmp/$gene.ref.fa";
        print TE ">$gene\n$ref\n";
	print OUT ">$gene\n$ref\n";
	close TE;
	my $option1 = "-N 1 --local -X 1000  ";
        my $option2 = "--score-min L,1,0.6 -N 1 --local -X 1000 ";
	`bowtie2-build $outdir/tmp/$gene.ref.fa $outdir/tmp/$gene.ref.bt2`;
	`bowtie2 $option1 -p 8 -x $outdir/tmp/$gene.ref.bt2 -1 $fq1 -2 $fq2 -S $outdir/tmp/$gene.split.sam`;
       	`$samtools view -bS $outdir/tmp/$gene.split.sam |$samtools sort - >$outdir/tmp/$gene.split.bam`;
	`$samtools index $outdir/tmp/$gene.split.bam`;
	`rm -rf $outdir/tmp/$gene.ref.fa* $outdir/tmp/$gene.split.sam`;
}
close OUT;
`bowtie2-build $outdir/$class.ref.merge.fa $outdir/$class.ref.merge.bt2`;
`$samtools faidx $outdir/$class.ref.merge.fa`;
my ($fq1,$fq2) = ($hash1{$id}, $hash2{$id});
`bowtie2 -p 8 -x $outdir/$class.ref.merge.bt2 -1 $fq1 -2 $fq2 -S $outdir/$class.tmp.sam`;
`$samtools view -H $outdir/$class.tmp.sam  >$outdir/header.sam`;
`$samtools merge -f -h $outdir/header.sam $outdir/$class.tmp.bam $outdir/tmp/*.split.bam `;
`$samtools addreplacerg -o $outdir/merge.$class.bam -O BAM -r "\@RG\tID:$class\tLB:$class\tPL:ILLUMINA\tSM:$class" $outdir/$class.tmp.bam `;
`$samtools index $outdir/merge.$class.bam`;
`rm -rf $outdir/$class.tmp.bam $outdir/$class.tmp.sam`;
`rm -rf $outdir/tmp/*.split.bam $outdir/tmp/*.split.bam.bai`;
`$freebayes -f $outdir/$class.ref.merge.fa -p 3 -m 0 -V $outdir/merge.$class.bam > $outdir/$class.freebayes.snp.vcf`;

