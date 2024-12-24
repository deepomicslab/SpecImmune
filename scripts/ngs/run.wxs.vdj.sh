sample=$1
fq1=$2
fq2=$3

dir=/home/wangmengyao/SpecComplex/
db=$dir/db
bin=$dir/bin
sdir=$dir/script/vdj
outdir=$(pwd)/$sample
mkdir -p $outdir
mkdir $outdir
ref=$db/imgt_human_tr_ig_uniq.fasta
$bin/bwa mem -t 8 -a -R "@RG\\tID:${sample}\\tPL:ILLUMINA\\tSM:${sample}" $ref $fq1 $fq2 > $outdir/$sample.sam
$bin/samtools view --threads 8 -Sbh -F 4 $outdir/$sample.sam > $outdir/$sample.bam
rm -rf $outdir/$sample.sam
$bin/samtools sort --threads 8 $outdir/$sample.bam > $outdir/$sample.sort.bam
rm -rf $outdir/$sample.bam
$bin/samtools index $outdir/$sample.sort.bam

$bin/python3  $sdir/VDJ_reads_assign.wgs.nm.v2.py -b $outdir/$sample.sort.bam -o $outdir -nm 2 -p 3 -r 0.05  >  $outdir/log2.txt

perl $sdir/selected.vdj.db.PE.pl $outdir/read.vdj.file.txt $outdir/read1.assign_file.txt $outdir/read2.assign_file.txt $outdir $sample $outdir/$sample.sort.bam

vcf=$outdir/$sample.vcf
less $vcf | grep "#" > $outdir/$sample.filter.vcf
less $vcf |grep -v "#"|awk '$6>1' >> $outdir/$sample.filter.vcf
ref=$outdir/$sample.VDJ.ref.fasta
bam=$outdir/merge.$sample.bam
$bin/extractHAIRS_spec --indels 1 --ref $ref --bam $bam --VCF $outdir/$sample.filter.vcf --out $outdir/$sample.frament.file
$bin/bgzip -f $outdir/$sample.filter.vcf
$bin/tabix $outdir/$sample.filter.vcf.gz

$bin/SpecHap  --vcf $outdir/$sample.filter.vcf.gz --frag $outdir/$sample.frament.file --out $outdir/$sample.specHap.phased.vcf
$bin/bgzip -f $outdir/$sample.specHap.phased.vcf
$bin/tabix $outdir/$sample.specHap.phased.vcf.gz

perl $sdir/vdj.gene.add.snp2.pl $outdir/$sample.specHap.phased.vcf.gz $ref $outdir

perl $sdir/check.phased.gene.pl $outdir/read.vdj.file.txt $outdir/VDJ.phased.fasta $outdir/check.phase.txt

