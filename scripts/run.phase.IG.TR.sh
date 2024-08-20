sample=$1
reads=$2
outdir=$3
db=$4
threads=$5
low_depth=$6
ref=$7



dir=$(cd `dirname $0`; pwd)

###### minimap2 -t $threads -R "@RG\tID:$sample\tSM:$sample" -k 19 -w 19 -g 10k -A 1 -B 4 -O 6,26 -E 2,1 -s 200 -a $ref $reads  | samtools view -bS -F 0x800 -| samtools sort - >$outdir/$sample.bam
bwa mem -t $threads -R "@RG\tID:$sample\tSM:$sample" $ref $reads | samtools view -bS -F 0x800 -| samtools sort - >$outdir/$sample.bam
##### minimap2 -t $threads -R "@RG\tID:$sample\tSM:$sample" -a $ref $reads  | samtools view -bS -F 0x800 -| samtools sort - >$outdir/$sample.bam

samtools index $outdir/$sample.bam
longshot -F -S --sample_id $sample  --bam $outdir/$sample.bam --ref $ref --out $outdir/$sample.longshot.vcf # --min_alt_frac 0.3

pbsv discover $outdir/$sample.bam $outdir/$sample.svsig.gz
pbsv call $ref $outdir/$sample.svsig.gz $outdir/$sample.sv.vcf


bgzip -f $outdir/$sample.sv.vcf
bgzip -f $outdir/$sample.longshot.vcf
tabix -f $outdir/$sample.sv.vcf.gz
tabix -f $outdir/$sample.longshot.vcf.gz

python3 $dir/filter_vcf.py $outdir/$sample.longshot.vcf.gz $outdir/$sample.longshot.2.vcf.gz 0.3
tabix -f $outdir/$sample.longshot.2.vcf.gz
bcftools concat -a $outdir/$sample.sv.vcf.gz $outdir/$sample.longshot.2.vcf.gz -Oz -o $outdir/$sample.merge.unsort.vcf.gz

# bcftools concat -a $outdir/$sample.sv.vcf.gz $outdir/$sample.longshot.vcf.gz -Oz -o $outdir/$sample.merge.unsort.vcf.gz
bcftools sort $outdir/$sample.merge.unsort.vcf.gz -Oz -o $outdir/$sample.merge.vcf.gz
tabix -f $outdir/$sample.merge.vcf.gz

whatshap phase -o $outdir/$sample.phase.vcf.gz -r $ref --indels $outdir/$sample.merge.vcf.gz $outdir/$sample.bam
tabix -f $outdir/$sample.phase.vcf.gz


samtools depth -b $dir/../gene_dist/IG_TR_chrs.bed -a $outdir/$sample.bam > $outdir/$sample.depth.txt

# python $dir/mask_low_depth_region.py -c $outdir/$sample.depth.txt -o $outdir  -d $low_depth -w 3 #-r ../gene_dist/IG_TR.gene.bed -i IG_TR  #-w 100
touch $outdir/low_depth.bed

bcftools norm -f $ref -O z -o $outdir/$sample.phase.norm.vcf.gz $outdir/$sample.phase.vcf.gz
tabix -f $outdir/$sample.phase.norm.vcf.gz