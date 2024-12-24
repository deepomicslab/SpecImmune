#!/usr/bin/perl
use FindBin qw($Bin);
#perl $0 KIR fq1.list fq2.list output
my ($class,$fq1,$fq2,$outdir) = @ARGV;
my (%hashr, %hash1, %hash2);

my $db="$Bin/../db";
my $bin="$Bin/../bin";

my $bwa = "$bin/bwa";
my $samtools = "$bin/samtools";
my $freebayes = "$bin/freebayes";
my $group = "\@RG\tID:$class\tSM:$class";
`mkdir $outdir/tmp`;

my $option = "-Y -T 20 -k 10 ";
if($class eq "HLA"){$option = " -U 10000 -L 10000,10000 ";}
open RE, "$db/immune.complex.gene.ref.extend.fasta" or die "$!\n";
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
	`$bwa index $outdir/tmp/$gene.ref.fa`;
	`$bwa mem -t 10 $option $outdir/tmp/$gene.ref.fa $fq1 $fq2 |$samtools view -bS -F 0x800 -|$samtools sort - >$outdir/tmp/$gene.split.bam`;
	`$samtools index $outdir/tmp/$gene.split.bam`;
	`rm -rf $outdir/tmp/$gene.ref.fa*`;
}
close OUT;
`$bwa index $outdir/$class.ref.merge.fa`;
`$samtools faidx $outdir/$class.ref.merge.fa`;
my ($fq1,$fq2) = ($hash1{$id}, $hash2{$id});
`$bwa mem $outdir/$class.ref.merge.fa $fq1 $fq2 | $samtools view -H >$outdir/header.sam`;
`$samtools merge -f -h $outdir/header.sam $outdir/$class.tmp.bam $outdir/tmp/*.split.bam `;
`$samtools addreplacerg -o $outdir/merge.$class.bam -O BAM -r "\@RG\tID:$class\tLB:$class\tPL:ILLUMINA\tSM:$class" $outdir/$class.tmp.bam `;
`$samtools index $outdir/merge.$class.bam`;
`rm -rf $outdir/$class.tmp.bam`;
`rm -rf $outdir/tmp/*.split.bam $outdir/tmp/*.split.bam.bai`;
`$freebayes -f $outdir/$class.ref.merge.fa -p 3 $outdir/merge.$class.bam > $outdir/$class.freebayes.snp.vcf`;

