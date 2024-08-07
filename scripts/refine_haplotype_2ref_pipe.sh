#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 9 ]; then
    echo "Usage: $0 <bam> <hla_ref> <hla> <interval> <mask_bed> <gene_work_dir> <threads> <sample> <allele_idx>"
    exit 1
fi

# Assign command-line arguments to variables
bam=$1
hla_ref=$2
hla=$3
interval=$4
mask_bed=$5
gene_work_dir=$6
threads=$7
sample=$8
allele_idx=$9

# Define paths
scripts_dir=$(dirname "$0")
longphase=$scripts_dir/../bin/longphase

# Ensure the Python script exists
python_script="$scripts_dir/check_ad.py"
if [ ! -f "$python_script" ]; then
    echo "Python script $python_script not found."
    exit 1
fi

# Define output files
refined_sv=$gene_work_dir/HLA_$hla.snisv.vcf
filtered_sv=$gene_work_dir/HLA_$hla.snisv.filtered.vcf
snv_vcf=$gene_work_dir/../$sample.$hla.$allele_idx.phased.vcf.gz

# Index the BAM file
samtools index "$bam"

# Run the Python script and capture the most supported allele index
most_supported_allele=$(python "$python_script" "$snv_vcf")

# Check if the Python script ran successfully
if [ $? -ne 0 ]; then
    echo "Error running Python script."
    exit 1
fi

# Check if the output indicates no heterozygous sites were found
if [ "$most_supported_allele" == "No heterozygous sites found." ]; then
    echo "No heterozygous sites found."
    most_supported_allele=1
else
    # Print the captured allele index
    echo "Most supported allele index: $most_supported_allele"
fi

# Define Sniffles run function
run_sniffles() {
    local genotype_error=$1
    local minsupport=$2
    local mapq=$3
    local cluster_binsize=$4
    local cluster_r=$5
    local cluster_merge_pos=$6
    sniffles -i "$bam" \
        -v "$refined_sv" \
        --phase \
        --minsvlen 50 \
        --minsupport "$minsupport" \
        --allow-overwrite \
        --reference "$hla_ref" \
        --mapq "$mapq" \
        --no-qc \
        --genotype-error "$genotype_error" \
        --cluster-binsize "$cluster_binsize" \
        --cluster-r "$cluster_r" \
        --cluster-merge-pos "$cluster_merge_pos" \
        --sample-id 'SAMPLE' \
        --output-rnames
    echo "SV calling done!"
}

# Parameters
genotype_error=0.05
minsupport="auto"
mapq=20
cluster_binsize=100
cluster_r=2.5
cluster_merge_pos=150

# Set parameters based on HLA type
case "$hla" in
    "DQB1")
        genotype_error=0.00
        minsupport=1
        mapq=0
        cluster_binsize=10
        cluster_r=0.1
        cluster_merge_pos=1000
        ;;
    "DRB1")
        genotype_error=0.01
        minsupport=1
        mapq=0
        cluster_binsize=10
        cluster_r=0.1
        cluster_merge_pos=1000
        ;;
esac

echo "Running Sniffles for $sample!"
run_sniffles "$genotype_error" "$minsupport" "$mapq" "$cluster_binsize" "$cluster_r" "$cluster_merge_pos"

echo "Filtering for $sample!"
bcftools view -i 'INFO/PRECISE=1' -g hom "$refined_sv" -o "$filtered_sv"

fmt_sv=$gene_work_dir/HLA_$hla.$allele_idx.snisv.filtered.fmt.vcf
python "$scripts_dir/vcf2seq.py" "$filtered_sv" "$hla_ref" "$fmt_sv"
snv_sv_merged=$gene_work_dir/$sample.$hla.$allele_idx.snv_sv.merged.vcf.gz
bcftools concat "$fmt_sv" "$snv_vcf" -Oz -o "$snv_sv_merged"
sorted_snv_sv_merged=$gene_work_dir/$sample.$hla.$allele_idx.snv_sv.merged.sorted.vcf.gz
bcftools sort "$snv_sv_merged" -Oz -o "$sorted_snv_sv_merged"
tabix -f "$sorted_snv_sv_merged"

if [[ "$hla" =~ ^(HLA-A|HLA-B|HLA-C)$ ]]; then
    sorted_snv_sv_merged=$snv_vcf
else
    sorted_snv_sv_merged=$gene_work_dir/$sample.$hla.$allele_idx.snv_sv.merged.sorted.vcf.gz
fi

echo "Gene: $hla"
echo "File: $sorted_snv_sv_merged"
echo "Mask bed:"
cat "$mask_bed"
hap_idx=$(expr $allele_idx + 1)
# Generate the consensus sequence
echo """
samtools faidx "$hla_ref" "$interval" | bcftools consensus -H "$most_supported_allele" --mask "$mask_bed" "$sorted_snv_sv_merged" > "$gene_work_dir/$hla.$hap_idx.raw.fa"
"""

samtools faidx "$hla_ref" "$interval" | bcftools consensus -H "$most_supported_allele" --mask "$mask_bed" "$sorted_snv_sv_merged" > "$gene_work_dir/$hla.$hap_idx.raw.fa"

# Additional command example (ensure the path is correct for your environment)
# bash /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/SpecLong/scripts/refine_haplotype_2ref_pipe.sh test_fh14_2/FH14//HLA-A.0.bam test_fh14_2/FH14/individual_ref/HLA-A/HLA-A.1.fasta HLA-A HLA-A_ref1 test_fh14_2/FH14//low_depth.bed test_fh14_2/FH14//HLA-A_0_work 20