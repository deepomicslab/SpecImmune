sample=$1
reads=$2
outdir=$3
db=$4
threads=$5
low_depth=$6
ref=$7
gene_work_dir=$8


dir=$(cd `dirname $0`; pwd)

###### minimap2 -t $threads -R "@RG\tID:$sample\tSM:$sample" -k 19 -w 19 -g 10k -A 1 -B 4 -O 6,26 -E 2,1 -s 200 -a $ref $reads  | samtools view -bS -F 0x800 -| samtools sort - >$outdir/$sample.bam
# bwa mem -t $threads -R "@RG\tID:$sample\tSM:$sample" $ref $reads | samtools view -bS -F 0x800 -| samtools sort - >$outdir/$sample.bam
minimap2 -t $threads -a $ref $reads  | samtools view -bS -F 0x800 -| samtools sort - >$outdir/$sample.bam

# call variant

bash $dir/run_dv.sh $ref $outdir/$sample.bam $outdir/$sample.dv.vcf $outdir/$sample.dv.gvcf $threads chr22

small_vars_vcf=$outdir/$sample.dv.vcf
refined_sv=$gene_work_dir/$sample.snisv.vcf
filtered_sv=$gene_work_dir/$sample.snisv.filtered.vcf

run_sniffles() {
    local genotype_error=$1
    local minsupport=$2
    local mapq=$3
    local cluster_binsize=$4
    local cluster_r=$5
    local cluster_merge_pos=$6
    sniffles -i $bam \
        -v $refined_sv \
        --phase \
        --minsvlen 50 \
        --minsupport $minsupport \
        --allow-overwrite \
        --reference $ref \
        --mapq $mapq \
        --no-qc \
        --genotype-error $genotype_error \
        --cluster-binsize $cluster_binsize \
        --cluster-r $cluster_r \
        --cluster-merge-pos $cluster_merge_pos \
        --sample-id 'SAMPLE' \
        --output-rnames
    echo "sv calling done!"
}

run_sniffles $genotype_error $minsupport $mapq $cluster_binsize $cluster_r $cluster_merge_pos

#whatshap phase
phased_snv_vcf=$outdir/$sample.dv.small.phased.vcf.gz
whatshap phase -o $phased_snv_vcf --reference=$ref $small_vars_vcf $bam --ignore-read-group
whatshap haplotag \
    --ignore-read-groups  \
    -o $gene_work_dir/haplotagged.bam \
    --reference $ref \
    $phased_snv_vcf \
    $bam \
    --output-haplotag-list $gene_work_dir/hap.tsv

whatshap split --output-h1 $outdir/h0.bam --output-h2 $gene_work_dir/h1.bam $gene_work_dir/haplotagged.bam $gene_work_dir/hap.tsv --output-untagged $gene_work_dir/untag.bam


$longphase phase -s $phased_snv_vcf \
    -b $gene_work_dir/haplotagged.bam \
    -r $ref \
    --sv-file $refined_sv \
    -o $gene_work_dir/longphase \
    --pb

bcftools view -i 'GT!="0/0" && GT!="." && GT!="./." && INFO/PRECISE=1' $gene_work_dir/longphase_SV.vcf -o $filtered_sv


snv_sv_merged=$gene_work_dir/$sample.snv_sv.merged.vcf.gz
bcftools concat $filtered_sv $gene_work_dir/longphase.vcf -Oz -o $snv_sv_merged
sorted_snv_sv_merged=$gene_work_dir/$sample.snv_sv.merged.sorted.vcf.gz
bcftools sort $snv_sv_merged -Oz -o $sorted_snv_sv_merged
# whatshap split --output-h1 $outdir/h0.bam --output-h2 $gene_work_dir/h1.bam $gene_work_dir/haplotagged.bam $gene_work_dir/hap.tsv --output-untagged $gene_work_dir/untag.bam
