#!/usr/bin/perl -w
use FindBin qw($Bin);
my $db="$Bin/../db/IG_TR";
my $bin="$Bin/../bin";
#perl $0 read.vdj.file.txt read1.assign_file.txt read2.assign_file.txt ./test escc001n ESCC001N.sort.bam
#my $nn = "N"x20;
my $nn = "";

my ($infile, $afile1,$afile2, $outdir, $sample, $bam) = @ARGV;
my $bwa = "$bin/bwa";
my $samtools = "$bin/samtools";
my $freebayes = "$bin/freebayes";
my $group = "\@RG\tID:$sample\tSM:$sample";
`mkdir -p $outdir`;
`$samtools fastq -1 $outdir/$sample.r1.fq -2 $outdir/$sample.r2.fq $bam`;
my %hash;
my $key;
open FA, "$db/IMGT.VDJ.fasta" or die "$!\n";
while(<FA>){
    chomp;
    if(/^>/){
         s/^>//;
         $key = $_;
    }
    else{
         $hash{$key} .= $_;
    }
}
close FA;

my (%hashr,%hasha,%hashs);
open AI, "$afile1" or die "$!\n";
while(<AI>){
    chomp;
     my ($id,$gene) = (split)[0,1];
     $id .= "_1";
     $hashr{$gene} .= "$id\t";
     $hasha{$gene} += 1;
}
close AI;
open AI, "$afile2" or die "$!\n";
while(<AI>){
    chomp;
    my ($id,$gene) = (split)[0,1];
    $id .= "_2";
    $hashr{$gene} .= "$id\t";
    $hasha{$gene} += 1;
}
close AI;

open OUT, ">$outdir/$sample.VDJ.ref.fasta";
open IN, "$infile" or die "$!\n";
while(<IN>){
    chomp;
    my ($class, $genes) = split;
    my @gens = split /,/, $genes;
    my $id = $gens[0];
    my $seq = $hash{$id};
    #my $sseq = "$nn"."$seq"."$nn";
    #my $class = (split /\*/, $id)[0];
    $hashs{$class} = "$id";
    print OUT ">$id\n$seq\n";
    
}
close IN;
close OUT;

my (%hash1,%hash2);
open F1, "$outdir/$sample.r1.fq" or die "$!\tfq1\n";
while(<F1>){
     chomp;
     my $id = $_;
     $id =~ s/^@//;
     my $line1 = $_;
     my $line2 = <F1>;
     my $line3 = <F1>;
     my $line4 = <F1>;
     $id .= "_1";
     $line1 .= "_1";
     chomp $line2;chomp $line3;chomp $line4;
     my $out = "$line1\n$line2\n$line3\n$line4";
     $hash1{$id} = $out;
} 
close F1;

open F2, "$outdir/$sample.r2.fq" or die "$!\tfq2\n";
while(<F2>){
    chomp;
    my $id = $_;
    $id =~ s/^@//;
    my $line1 = $_;
    my $line2 = <F2>;
    my $line3 = <F2>;
    my $line4 = <F2>;
    $id .= "_2";
    $line1 .= "_2";
    chomp $line2; chomp $line3; chomp $line4;
    my $out = "$line1\n$line2\n$line3\n$line4";
    $hash1{$id} = $out;
}
close F2;

foreach my $class(sort keys %hasha){
     $class =~ s/\//_/;
     next if(!exists $hashs{$class});
     next if(!exists $hashr{$class});
     my $gene =  $hashs{$class};
     my @reads = (split /\t/, $hashr{$class});
     open O1, ">$outdir/$sample.$class.read.fq";
     #open O2, ">$outdir/$sample.$class.R2.fq";
     open OO, ">$outdir/$sample.$class.split.fa";
     #print "@genes\n";
     
     next if(!exists $hash{$gene});
     my $fa = "$hash{$gene}";
     print OO ">$gene\n$fa\n";
     
     my $sum = 0;
     foreach my $read(@reads){
           next if(!exists $hash1{$read});
	   # next if(!exists $hash2{$read});
           my $f1 = $hash1{$read};
	   #my $f2 = $hash2{$read};
           $sum += 1;
           print O1 "$f1\n";
	   #print O2 "$f2\n";
     }
     close OO;
     close O1;
     #close O2;
     next if($sum<10);

    `$bwa index $outdir/$sample.$class.split.fa`;
    `$bwa mem $outdir/$sample.$class.split.fa $outdir/$sample.$class.read.fq |$samtools view -bS -F 0x800 - |$samtools sort - >$outdir/$sample.$class.split.bam`;
    `$samtools index $outdir/$sample.$class.split.bam`;
}
#print "$group\n";
`$bwa index $outdir/$sample.VDJ.ref.fasta`;
`$samtools faidx $outdir/$sample.VDJ.ref.fasta`;
`$bwa mem $outdir/$sample.VDJ.ref.fasta $outdir/$sample.r1.fq $outdir/$sample.r1.fq |$samtools view -H >$outdir/header.sam`;
`$samtools merge -f -h $outdir/header.sam $outdir/merge.$sample.bam $outdir/$sample.*.split.bam `;
`$samtools index $outdir/merge.$sample.bam`;
`rm -rf $outdir/*.read.fq $outdir/header.sam`;
`rm -rf $outdir/*.split.bam* $outdir/*.split.fa*`;
`$freebayes -f $outdir/$sample.VDJ.ref.fasta -m 30 -q 10 -R 0 -S 0 -p 2 -C 2 $outdir/merge.$sample.bam >$outdir/$sample.vcf`;

