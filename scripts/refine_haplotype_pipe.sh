#!/bin/bash

# 参数列表
bam=$1
hla_ref=$2
hla=$3
interval=$4
mask_bed=$5
gene_work_dir=$6
threads=$7
sample=$8
scripts_dir=$(dirname $0)
longphase=$scripts_dir/../bin/longphase
data_type=$9

# mask_low_script=$scripts_dir/mask_low_depth_region.py
get_intron_script=$scripts_dir/get_intron.py


# seq_type dict:
# pacbio :pacbio
# nanopore: ont
# traditional: pacbio
# 2D: ont
# Direct: ont
# SIRV: pacbio
seq_type="pb"
if [[ "$data_type" =~ ^(nanopore|2D|Direct)$ ]]; then
    seq_type="ont"
fi

# Sniffles 
refined_sv=$gene_work_dir/HLA_$hla.snisv.vcf
filtered_sv=$gene_work_dir/HLA_$hla.snisv.filtered.vcf
snv_vcf=$gene_work_dir/../$sample.$hla.phased.vcf.gz
samtools index $bam

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
    echo "sv calling done!"
}

# parameters
genotype_error=0.05
minsupport="auto"
mapq=20
cluster_binsize=100
cluster_r=2.5
cluster_merge_pos=150


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

echo "sniffles for $sample !"
run_sniffles $genotype_error $minsupport $mapq $cluster_binsize $cluster_r $cluster_merge_pos

# Whatshap haplotag
echo "haplotag for $sample !"
# fix dup HD in header
if [ $(samtools view -H $bam | grep -c '^@HD') -gt 1 ]; then
    fixed_bam=$gene_work_dir/fixed.bam
    samtools view -H $bam | awk 'BEGIN {last_hd_line=-1} /^@HD/ {last_hd_line=NR} {lines[NR]=$0} END {for (i=1; i<=NR; i++) if (i != last_hd_line) print lines[i]}' | samtools reheader - $bam > $fixed_bam
else
  echo "The BAM file does not contain multiple @HD lines."
  fixed_bam=$bam
fi

samtools index $fixed_bam



echo """whatshap haplotag --ignore-read-groups -o $gene_work_dir/haplotagged.$hla.bam --reference $hla_ref $snv_vcf $fixed_bam --output-haplotag-list $gene_work_dir/hap.tsv"""
whatshap haplotag \
    --ignore-read-groups  \
    -o $gene_work_dir/haplotagged.$hla.bam \
    --reference $hla_ref \
    $snv_vcf \
    $fixed_bam \
    --output-haplotag-list $gene_work_dir/hap.tsv

samtools index $gene_work_dir/haplotagged.$hla.bam
whatshap split --output-h1 $gene_work_dir/h0.bam --output-h2 $gene_work_dir/h1.bam $gene_work_dir/haplotagged.$hla.bam $gene_work_dir/hap.tsv
samtools index $gene_work_dir/h0.bam
samtools index $gene_work_dir/h1.bam
samtools depth -aa $gene_work_dir/h0.bam > $gene_work_dir/h0.depth
samtools depth -aa $gene_work_dir/h1.bam > $gene_work_dir/h1.depth
# whatshap split --output-h1 $outdir/h0_untag.bam --output-h2 $gene_work_dir/h1_untag.bam $gene_work_dir/haplotagged.bam $gene_work_dir/hap.untag.tsv --add-untagged

# samtools index $gene_work_dir/h0.bam
# samtools index $gene_work_dir/h1.bam
# samtools index $gene_work_dir/h0_untag.bam
# samtools index $gene_work_dir/h1_untag.bam


samtools index $gene_work_dir/haplotagged.$hla.bam

# Longphase
# todo:: change to spechap here 
echo "longphase for $sample !"
$longphase phase -s $snv_vcf \
    -b $gene_work_dir/haplotagged.$hla.bam \
    -r $hla_ref \
    --sv-file $refined_sv \
    -o $gene_work_dir/longphase \
    --${seq_type}

echo "filtering for $sample !"
bcftools view -i 'GT!="0/0" && GT!="." && GT!="./." && INFO/PRECISE=1' $gene_work_dir/longphase_SV.vcf -o $filtered_sv

fmt_sv=$gene_work_dir/HLA_$hla.snisv.filtered.fmt.vcf
python $scripts_dir/vcf2seq.py $filtered_sv $hla_ref $fmt_sv
snv_sv_merged=$gene_work_dir/$sample.$hla.snv_sv.merged.vcf.gz
bcftools concat $fmt_sv $gene_work_dir/longphase.vcf -Oz -o $snv_sv_merged
sorted_snv_sv_merged=$gene_work_dir/$sample.$hla.snv_sv.merged.sorted.vcf.gz
bcftools sort $snv_sv_merged -Oz -o $sorted_snv_sv_merged
tabix -f $sorted_snv_sv_merged

if [[ "$hla" =~ ^(HLA-A|HLA-B|HLA-C)$ ]]; then
    sorted_snv_sv_merged=$snv_vcf
else
    sorted_snv_sv_merged=$gene_work_dir/$sample.$hla.snv_sv.merged.sorted.vcf.gz
    # sorted_snv_sv_merged=$snv_vcf
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

# if seq_type == traditional or 2D or Direct or SIRV, generate low depth file for hap0 and hap1
if [[ "$data_type" =~ ^(traditional|2D|Direct|SIRV)$ ]]; then

    # assemble hap0 and hap1
    stringtie $gene_work_dir/h0.bam -o $gene_work_dir/h0.gtf -c 10
    stringtie $gene_work_dir/h1.bam -o $gene_work_dir/h1.gtf -c 10

    python $get_intron_script $gene_work_dir/h0.gtf $hla_ref $gene_work_dir/h0.intron.bed
    python $get_intron_script $gene_work_dir/h1.gtf $hla_ref $gene_work_dir/h1.intron.bed

    # echo "masking low depth region for $sample !"
    # avg_depth=$(samtools depth  $gene_work_dir/h0.bam | awk '{sum+=$3} END {print sum/NR}')
    # echo "average depth: $avg_depth for $gene_work_dir/h0.bam"
    # # assign int value of 0.5*avg_depth to set_dp
    # set_dp=$(awk -v avg=$avg_depth 'BEGIN {print int(0.8*avg)}')
    # python $mask_low_script -f False -c $gene_work_dir/h0.depth -o $gene_work_dir -w 1 -d $set_dp
    # mv $gene_work_dir/low_depth.bed $gene_work_dir/h0.low_depth.bed

    # avg_depth=$(samtools depth  $gene_work_dir/h1.bam | awk '{sum+=$3} END {print sum/NR}')
    # echo "average depth: $avg_depth for $gene_work_dir/h1.bam"
    # set_dp=$(awk -v avg=$avg_depth 'BEGIN {print int(0.8*avg)}')
    # python $mask_low_script -f False -c $gene_work_dir/h1.depth -o $gene_work_dir -w 1 -d $set_dp
    # mv $gene_work_dir/low_depth.bed $gene_work_dir/h1.low_depth.bed
    

    # samtools faidx $hla_ref $interval | bcftools consensus -H 1 --mask $gene_work_dir/h0.low_depth.bed $sorted_snv_sv_merged > $gene_work_dir/$hla.1.raw.fa
    # samtools faidx $hla_ref $interval | bcftools consensus -H 2 --mask $gene_work_dir/h1.low_depth.bed $sorted_snv_sv_merged > $gene_work_dir/$hla.2.raw.fa

    samtools faidx $hla_ref $interval | bcftools consensus -H 1 --mask $gene_work_dir/h0.intron.bed $sorted_snv_sv_merged > $gene_work_dir/$hla.1.raw.fa
    samtools faidx $hla_ref $interval | bcftools consensus -H 2 --mask $gene_work_dir/h1.intron.bed $sorted_snv_sv_merged > $gene_work_dir/$hla.2.raw.fa
else
    echo "gene : $hla"
    echo "file : $sorted_snv_sv_merged"
    samtools faidx $hla_ref $interval | bcftools consensus -H 1 --mask $mask_bed $sorted_snv_sv_merged > $gene_work_dir/$hla.1.raw.fa
    samtools faidx $hla_ref $interval | bcftools consensus -H 2 --mask $mask_bed $sorted_snv_sv_merged > $gene_work_dir/$hla.2.raw.fa
fi



