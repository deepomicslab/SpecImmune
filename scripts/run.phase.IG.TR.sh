sample=$1
reads=$2
outdir=$3
db=$4
# ref=$db_dir/IG_TR/merge.IG.TR.ref.fasta
# ref=$db_dir/IG_TR/ig.tr.merge.hg38.fa
# ref=/mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/chr14.fa
threads=$5
low_depth=$6
ref=$7



dir=$(cd `dirname $0`; pwd)

###### minimap2 -t $threads -R "@RG\tID:$sample\tSM:$sample" -k 19 -w 19 -g 10k -A 1 -B 4 -O 6,26 -E 2,1 -s 200 -a $ref $reads  | samtools view -bS -F 0x800 -| samtools sort - >$outdir/$sample.bam
bwa mem -t $threads -R "@RG\tID:$sample\tSM:$sample" $ref $reads | samtools view -bS -F 0x800 -| samtools sort - >$outdir/$sample.bam

samtools index $outdir/$sample.bam
longshot -F -S --sample_id $sample  --bam $outdir/$sample.bam --ref $ref --out $outdir/$sample.longshot.vcf

pbsv discover $outdir/$sample.bam $outdir/$sample.svsig.gz
pbsv call $ref $outdir/$sample.svsig.gz $outdir/$sample.sv.vcf

perl $dir/add.deletion.pl $outdir/$sample.longshot.vcf $outdir/$sample.sv.vcf $outdir/$sample.merge.vcf
bgzip -f $outdir/$sample.merge.vcf
tabix -f $outdir/$sample.merge.vcf.gz

whatshap phase -o $outdir/$sample.phase.vcf.gz -r $ref --indels $outdir/$sample.merge.vcf.gz $outdir/$sample.bam
tabix -f $outdir/$sample.phase.vcf.gz


samtools depth -a $outdir/$sample.bam > $outdir/$sample.depth.txt

python $dir/mask_low_depth_region.py -c $outdir/$sample.depth.txt -o $outdir -w 100 -d $low_depth 

bcftools norm -f $ref -O z -o $outdir/$sample.phase.norm.vcf.gz $outdir/$sample.phase.vcf.gz
tabix -f $outdir/$sample.phase.norm.vcf.gz

bcftools consensus  -f $ref -H 1 $outdir/$sample.phase.norm.vcf.gz >$outdir/$sample.hap1.raw.fasta
bcftools consensus  -f $ref -H 2 $outdir/$sample.phase.norm.vcf.gz >$outdir/$sample.hap2.raw.fasta

perl $dir/anno.IG.TR.py $sample $outdir/$sample.hap1.raw.fasta $outdir/$sample.hap2.raw.fasta $outdir $db $threads