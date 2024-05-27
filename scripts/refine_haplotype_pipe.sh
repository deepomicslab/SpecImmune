#!/bin/bash

# 参数列表
bam=$1
hla_ref=$2
hla=$3
interval=$4
mask_bed=$5
gene_work_dir=$6
threads=$7
i=$8
sample=$9
scripts_dir=$(dirname $0)

# Sniffles 
refined_sv=$gene_work_dir/HLA_$hla.snisv.vcf
filtered_sv=$gene_work_dir/HLA_$hla.snisv.filtered.vcf
snv_vcf=$gene_work_dir/../$sample.$hla.phased.vcf.gz

# def
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
        --reference $hla_ref \
        --mapq $mapq \
        --no-qc \
        --genotype-error $genotype_error \
        --cluster-binsize $cluster_binsize \
        --cluster-r $cluster_r \
        --cluster-merge-pos $cluster_merge_pos \
        --sample-id 'SAMPLE' \
        --output-rnames
}

# parameters
genotype_error=0.05
minsupport="auto"
mapq=20
cluster-binsize=100
cluster-r=2.5
cluster-merge-pos=150


case "$hla" in
    "DQB1")
        genotype_error=0.00
        minsupport=1
        mapq=0
        cluster-binsize=10
        cluster-r=0.1
        cluster-merge-pos=1000
        ;;
    "DRB1")
        genotype_error=0.01
        minsupport=1
        mapq=0
        cluster-binsize=10
        cluster-r=0.1
        cluster-merge-pos=1000
        ;;
esac


run_sniffles $genotype_error $minsupport $mapq $cluster_binsize $cluster_r $cluster_merge_pos

# Whatshap haplotag
whatshap haplotag \
    --ignore-read-groups \
    -o $gene_work_dir/haplotagged.$hla.bam \
    --reference $hla_ref \
    $snv_vcf \
    $bam


whatshap split --output-h1 $outdir/h0.bam --output-h2 $gene_work_dir/h1.bam $gene_work_dir/haplotagged.bam $gene_work_dir/hap.tsv --output-untagged $gene_work_dir/untag.bam
whatshap split --output-h1 $outdir/h0_untag.bam --output-h2 $gene_work_dir/h1_untag.bam $gene_work_dir/haplotagged.bam $gene_work_dir/hap.untag.tsv --add-untagged

samtools index $gene_work_dir/h0.bam
samtools index $gene_work_dir/h1.bam
samtools index $gene_work_dir/h0_untag.bam
samtools index $gene_work_dir/h1_untag.bam


samtools index $gene_work_dir/haplotagged.$hla.bam

# Longphase
# todo:: change to spechap here 
longphase phase -s $snv_vcf \
    -b $gene_work_dir/haplotagged.$hla.bam \
    -r $hla_ref \
    --sv-file $refined_sv \
    -o $gene_work_dir/longphase \
    --ont


bcftools view -i 'GT!="0/0" && GT!="." && GT!="./." && INFO/PRECISE=1' $gene_work_dir/longphase_SV.vcf -o $filtered_sv

fmt_sv=$gene_work_dir/HLA_$hla.snisv.filtered.fmt.vcf
python $scripts_dir/vcf2seq.py $filtered_sv $hla_ref $fmt_sv
snv_sv_merged=$gene_work_dir/$sample.$hla.snv_sv.merged.vcf.gz
bcftools concat $fmt_sv $gene_work_dir/longphase.vcf -Oz -o $snv_sv_merged
sorted_snv_sv_merged=$gene_work_dir/$sample.$hla.snv_sv.merged.sorted.vcf.gz
bcftools sort $snv_sv_merged -Oz -o $sorted_snv_sv_merged
tabix -f $sorted_snv_sv_merged

if [[ "$hla" =~ ^(A|B|C)$ ]]; then
    sorted_snv_sv_merged=$gene_work_dir/../$sample.$hla.phased.vcf.gz
else
    sorted_snv_sv_merged=$gene_work_dir/$sample.$hla.snv_sv.merged.sorted.vcf.gz
fi

# matching step

# python $current_dir/parsesv2seg_new.py \
#     -sv $outdir/all_phased_sv.s.f.merge.vcf \
#     -p $hap1_dir/hap1 \
#     -hid 1 \
#     -ob $raw_bam \
#     -h0b $h0_bam \
#     -h1b $h1_bam \
#     -r $ref \
#     -snp $outdir/all_snp.s.vcf.gz

# ############################################################################ matching
# echo "matching for $sample !"
# h0_g=$hap0_dir/hap0_graph.txt
# h0_path=$hap0_dir/hap0_path.txt
# h0_c_path=$hap0_dir/hap0_c_path.txt
# h0_new_g=$hap0_dir/hap0_new_graph.txt

# # check
# sv_flag_file=$hap0_dir/hap0_nosv.flag
# if [ -s "$sv_flag_file" ]; then
#     echo "$sv_flag_file exists and is not empty."
#     cp $hap0_dir/hap0_seg.fa $hap0_dir/hap0_asm.fa
#     sed -i '1c\>mhc_hap_0' $hap0_dir/hap0_asm.fa
# else
#     echo "$sv_flag_file does not exist or is empty."
#     /home/gzpan2/scratch/virus/wxd/app/seqGraph/build2/matching \
#     -b \
#     --model 1 \
#     -v 1 \
#     -g $h0_g \
#     -r $h0_path \
#     -c $h0_c_path \
#     -m $h0_new_g \
#     --break_c
#     # make fasta for hap0
#     samtools faidx $hap0_dir/hap0_seg.fa
#     python $current_dir/make_fa_from_path.py \
#         $h0_path \
#         $hap0_dir/hap0_seg.txt \
#         $hap0_dir/hap0_seg.fa \
#         $hap0_dir/hap0_asm.fa \
#         0 
# fi



# h1_g=$hap1_dir/hap1_graph.txt
# h1_path=$hap1_dir/hap1_path.txt
# h1_c_path=$hap1_dir/hap1_c_path.txt
# h1_new_g=$hap1_dir/hap1_new_graph.txt

# # check
# sv_flag_file=$hap1_dir/hap1_nosv.flag
# if [ -s "$sv_flag_file" ]; then
#     echo "$sv_flag_file exists and is not empty."
#     cp $hap1_dir/hap1_seg.fa $hap1_dir/hap1_asm.fa
#     sed -i '1c\>mhc_hap_1' $hap1_dir/hap1_asm.fa
# else
#     echo "$sv_flag_file does not exist or is empty."
#     /home/gzpan2/scratch/virus/wxd/app/seqGraph/build2/matching \
#     -b \
#     --model 1 \
#     -v 1 \
#     -g $h1_g \
#     -r $h1_path \
#     -c $h1_c_path \
#     -m $h1_new_g \
#     --break_c
#     # make fasta for hap1
#     samtools faidx $hap1_dir/hap1_seg.fa
#     python $current_dir/make_fa_from_path.py \
#         $h1_path \
#         $hap1_dir/hap1_seg.txt \
#         $hap1_dir/hap1_seg.fa \
#         $hap1_dir/hap1_asm.fa \
#         1 
# fi

echo "gene : $hla"
echo "file : $sorted_snv_sv_merged"
samtools faidx $hla_ref $interval | bcftools consensus -H 1 --mask $mask_bed $sorted_snv_sv_merged > $gene_work_dir/$hla.1.raw.fa
samtools faidx $hla_ref $interval | bcftools consensus -H 2 --mask $mask_bed $sorted_snv_sv_merged > $gene_work_dir/$hla.2.raw.fa