#!/bin/bash

# ==============================================================================
# Setup and Parameters
# ==============================================================================

# Ensure pipe failures are caught, but do NOT use 'set -e' to allow caller failures
set -o pipefail

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
matching=$scripts_dir/../bin/matching
data_type=$9

# Helper scripts
get_intron_script=$scripts_dir/get_intron.py

# Determine sequencing type
seq_type="pb"
if [[ "$data_type" =~ ^(nanopore|2D|Direct)$ ]]; then
    seq_type="ont"
fi

# Define Output Files
refined_sv=$gene_work_dir/HLA_$hla.snisv.vcf
refined_sv2=$gene_work_dir/HLA_$hla.cutesv.vcf
refined_sv3=$gene_work_dir/HLA_$hla.dysgu.vcf
merged_sv=$gene_work_dir/HLA_$hla.snisv.cutesv.dysgu.vcf
snv_vcf=$gene_work_dir/../$sample.$hla.phased.vcf.gz

# Index input BAM if needed
if [ ! -f "${bam}.bai" ]; then
    samtools index $bam
fi

# ==============================================================================
# SV Caller 1: Sniffles
# ==============================================================================

# Default parameters
genotype_error=0.05
minsupport="auto"
mapq=20
cluster_binsize=100
cluster_r=2.5
cluster_merge_pos=150

# Gene-specific parameter adjustments
# FIX: Changed hyphens to underscores for variable assignment
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

run_sniffles() {
    echo "Running Sniffles for $sample..."
    # Added '|| true' to prevent script exit on failure
    sniffles -i $bam \
        -v $refined_sv \
        --phase \
        --minsvlen 50 \
        --minsupport $2 \
        --allow-overwrite \
        --reference $hla_ref \
        --mapq $3 \
        --no-qc \
        --genotype-error $1 \
        --cluster-binsize $4 \
        --cluster-r $5 \
        --cluster-merge-pos $6 \
        --sample-id 'SAMPLE' \
        --output-rnames || echo "WARNING: Sniffles encountered an error."
}

run_sniffles $genotype_error $minsupport $mapq $cluster_binsize $cluster_r $cluster_merge_pos

# ==============================================================================
# SV Caller 2: CuteSV
# ==============================================================================

cute_work_dir=$gene_work_dir/cutesv
mkdir -p $cute_work_dir

# Set parameters based on seq_type
if [[ "$seq_type" =~ ^(pb|pacbio)$ ]]; then
    max_cluster_bias_INS=1000
    diff_ratio_merging_INS=0.9
    max_cluster_bias_DEL=1000
    diff_ratio_merging_DEL=0.5
elif [[ "$seq_type" =~ ^(ont|nanopore)$ ]]; then
    max_cluster_bias_INS=100
    diff_ratio_merging_INS=0.3
    max_cluster_bias_DEL=100
    diff_ratio_merging_DEL=0.3
else
    echo "Unknown seq_type: $seq_type, using defaults."
    max_cluster_bias_INS=100
    diff_ratio_merging_INS=0.3
    max_cluster_bias_DEL=100
    diff_ratio_merging_DEL=0.3
fi

echo "Running CuteSV for $sample..."
cuteSV $bam $hla_ref $refined_sv2 $cute_work_dir \
    --max_cluster_bias_INS $max_cluster_bias_INS \
    --diff_ratio_merging_INS $diff_ratio_merging_INS \
    --max_cluster_bias_DEL $max_cluster_bias_DEL \
    --diff_ratio_merging_DEL $diff_ratio_merging_DEL \
    --genotype -S "cutesv" || echo "WARNING: CuteSV encountered an error."

# ==============================================================================
# SV Caller 3: Dysgu
# ==============================================================================

dysgu_work_dir=$gene_work_dir/dysgu
# FIX: Corrected variable name from cute_work_dir to dysgu_work_dir
mkdir -p $dysgu_work_dir 

echo "Running Dysgu for $sample..."
if [[ "$seq_type" =~ ^(pb|pacbio)$ ]]; then
    mode="pacbio"
elif [[ "$seq_type" =~ ^(ont|nanopore)$ ]]; then
    mode="nanopore"
else
    mode="nanopore" # Default fallback
fi

dysgu call --divergence auto --mode $mode -x $hla_ref $dysgu_work_dir $bam > $refined_sv3 || echo "WARNING: Dysgu encountered an error."

# ==============================================================================
# Process and Merge SVs (Fault Tolerant)
# ==============================================================================

VCF_LIST=""

# Check Sniffles
if [[ -s $refined_sv ]]; then
    bgzip -f $refined_sv
    tabix -f $refined_sv.gz
    VCF_LIST="$VCF_LIST $refined_sv.gz"
else
    echo "Sniffles output missing or empty."
fi

# Check CuteSV
if [[ -s $refined_sv2 ]]; then
    bgzip -f $refined_sv2
    tabix -f $refined_sv2.gz
    VCF_LIST="$VCF_LIST $refined_sv2.gz"
else
    echo "CuteSV output missing or empty."
fi

# Check Dysgu
if [[ -s $refined_sv3 ]]; then
    bgzip -f $refined_sv3
    tabix -f $refined_sv3.gz
    VCF_LIST="$VCF_LIST $refined_sv3.gz"
else
    echo "Dysgu output missing or empty."
fi

echo "Merging SVs from: $VCF_LIST"

if [[ -n "$VCF_LIST" ]]; then
    # Merge available VCFs
    bcftools merge -m none $VCF_LIST --force-samples -Ov -o $merged_sv
else
    echo "WARNING: All SV callers failed or produced empty output."
    # Create a dummy header-only VCF to allow downstream tools to proceed without crashing
    echo "##fileformat=VCFv4.2" > $merged_sv
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> $merged_sv
fi

bgzip -f $merged_sv
tabix -f $merged_sv.gz

# ==============================================================================
# Truvari Collapse
# ==============================================================================

truvari_merge=$gene_work_dir/HLA_$hla.snisv.cutesv.dysgu.truvari.vcf
truvari_collapsed_ref=$gene_work_dir/HLA_$hla.snisv.cutesv.dysgu.truvari.collapsed.vcf

# Only run Truvari if the merged VCF has content (more than just header)
if [[ $(zgrep -v "^#" $merged_sv.gz | wc -l) -gt 0 ]]; then
    echo "Running Truvari collapse..."
    truvari collapse -i $merged_sv.gz -o $truvari_merge -c $truvari_collapsed_ref -k first --gt het --intra || echo "Truvari failed."
else
    echo "No SVs to collapse. Creating empty Truvari output."
    cp $merged_sv.gz ${truvari_merge}.gz # Copy empty/dummy VCF
    gunzip -f ${truvari_merge}.gz
fi

# ==============================================================================
# Whatshap Haplotagging
# ==============================================================================

echo "Haplotagging for $sample !"

# fix dup HD in header
if [ $(samtools view -H $bam | grep -c '^@HD') -gt 1 ]; then
    fixed_bam=$gene_work_dir/fixed.bam
    samtools view -H $bam | awk 'BEGIN {last_hd_line=-1} /^@HD/ {last_hd_line=NR} {lines[NR]=$0} END {for (i=1; i<=NR; i++) if (i != last_hd_line) print lines[i]}' | samtools reheader - $bam > $fixed_bam
else
    echo "The BAM file does not contain multiple @HD lines."
    fixed_bam=$bam
fi

samtools index $fixed_bam

echo "Running Whatshap haplotag..."
whatshap haplotag \
    --ignore-read-groups  \
    -o $gene_work_dir/haplotagged.$hla.bam \
    --reference $hla_ref \
    $snv_vcf \
    $fixed_bam \
    --output-haplotag-list $gene_work_dir/hap.tsv

samtools index $gene_work_dir/haplotagged.$hla.bam

# Split BAMs
whatshap split --output-h1 $gene_work_dir/h0.bam --output-h2 $gene_work_dir/h1.bam $gene_work_dir/haplotagged.$hla.bam $gene_work_dir/hap.tsv
samtools index $gene_work_dir/h0.bam
samtools index $gene_work_dir/h1.bam
samtools depth -aa $gene_work_dir/h0.bam > $gene_work_dir/h0.depth
samtools depth -aa $gene_work_dir/h1.bam > $gene_work_dir/h1.depth

h0_bam=$gene_work_dir/h0.bam
h1_bam=$gene_work_dir/h1.bam

# ==============================================================================
# Longphase
# ==============================================================================

echo "Running Longphase for $sample !"

# Check if Truvari output exists and is valid for Longphase
SV_ARG=""
if [[ -s $truvari_merge ]] && [[ $(grep -v "^#" $truvari_merge | wc -l) -gt 0 ]]; then
    SV_ARG="--sv-file $truvari_merge"
else
    echo "No valid SV file for Longphase, running without SVs."
fi

$longphase phase -s $snv_vcf \
    -b $gene_work_dir/haplotagged.$hla.bam \
    -r $hla_ref \
    $SV_ARG \
    -o $gene_work_dir/longphase \
    --${seq_type}

# ==============================================================================
# Filtering and Merging SNV/SV
# ==============================================================================

echo "Filtering SVs..."
filtered_sv=$gene_work_dir/HLA_$hla.snisv.filtered.vcf
longphase_sv_vcf=$gene_work_dir/longphase_SV.vcf

# Check if Longphase produced SV output
if [[ -s $longphase_sv_vcf ]]; then
    bcftools view -i 'GT!="0/0" && GT!="." && GT!="./." && INFO/PRECISE=1' $longphase_sv_vcf -o $filtered_sv
else
    echo "Longphase SV output missing. Creating empty filtered SV file."
    echo "##fileformat=VCFv4.2" > $filtered_sv
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> $filtered_sv
fi

# Phased SV formatting
fmt_sv=$gene_work_dir/HLA_$hla.snisv.filtered.fmt.vcf
python $scripts_dir/vcf2seq.py $filtered_sv $hla_ref $fmt_sv

snv_sv_merged=$gene_work_dir/$sample.$hla.snv_sv.merged.vcf.gz

# Concat SNVs and SVs (handle missing SVs gracefully)
if [[ -s $fmt_sv ]] && [[ $(grep -v "^#" $fmt_sv | wc -l) -gt 0 ]]; then
    bcftools concat $fmt_sv $gene_work_dir/longphase.vcf -Oz -o $snv_sv_merged
else
    echo "No SVs to merge, using only Longphase SNVs."
    bcftools view $gene_work_dir/longphase.vcf -Oz -o $snv_sv_merged
fi

sorted_snv_sv_merged=$gene_work_dir/$sample.$hla.snv_sv.merged.sorted.vcf.gz
bcftools sort $snv_sv_merged -Oz -o $sorted_snv_sv_merged
tabix -f $sorted_snv_sv_merged

# ==============================================================================
# Consensus Sequence Generation
# ==============================================================================

if [[ "$data_type" =~ ^(traditional|2D|Direct|SIRV)$ ]]; then
    echo "Processing Traditional/2D/Direct/SIRV data..."
    
    # Assemble hap0 and hap1
    if [ $(bcftools view -H $snv_vcf -g het | wc -l) -eq 0 ]; then
         echo "No het variants in $snv_vcf"
         samtools index $fixed_bam
         stringtie $fixed_bam -o $gene_work_dir/h0_g1_mixed.gtf -c 10
         python $get_intron_script $gene_work_dir/h0_g1_mixed.gtf $hla_ref $gene_work_dir/h0_g1_mixed.intron.bed
         samtools faidx $hla_ref $interval | bcftools consensus -H 1 --mask $gene_work_dir/h0_g1_mixed.intron.bed $snv_vcf > $gene_work_dir/$hla.1.raw.fa
         samtools faidx $hla_ref $interval | bcftools consensus -H 2 --mask $gene_work_dir/h0_g1_mixed.intron.bed $snv_vcf > $gene_work_dir/$hla.2.raw.fa
    else
        stringtie $gene_work_dir/h0.bam -o $gene_work_dir/h0.gtf -c 10
        stringtie $gene_work_dir/h1.bam -o $gene_work_dir/h1.gtf -c 10

        python $get_intron_script $gene_work_dir/h0.gtf $hla_ref $gene_work_dir/h0.intron.bed
        python $get_intron_script $gene_work_dir/h1.gtf $hla_ref $gene_work_dir/h1.intron.bed
        samtools faidx $hla_ref $interval | bcftools consensus -H 1 --mask $gene_work_dir/h0.intron.bed $snv_vcf > $gene_work_dir/$hla.1.raw.fa
        samtools faidx $hla_ref $interval | bcftools consensus -H 2 --mask $gene_work_dir/h1.intron.bed $snv_vcf > $gene_work_dir/$hla.2.raw.fa
    fi

else
    # Standard Consensus Generation
    echo "Generating consensus for gene: $hla"
    echo "Using file: $sorted_snv_sv_merged"
    
    samtools faidx $hla_ref $interval | bcftools consensus -H 1 --mask $mask_bed $sorted_snv_sv_merged > $gene_work_dir/$hla.1.raw.fa
    samtools faidx $hla_ref $interval | bcftools consensus -H 2 --mask $mask_bed $sorted_snv_sv_merged > $gene_work_dir/$hla.2.raw.fa

    if [[ "$hla" =~ ^(HLA-A|HLA-B|HLA-C)$ ]]; then
        sorted_snv_sv_merged=$snv_vcf
        echo "HLA Class I detected, exiting early as per logic."
        exit 0
    else
        sorted_snv_sv_merged=$gene_work_dir/$sample.$hla.snv_sv.merged.sorted.vcf.gz
    fi

    # ==============================================================================
    # Matching Step (Graph Assembly)
    # ==============================================================================
    
    # --- Hap0 ---
    hap0_dir=$gene_work_dir/hap0
    mkdir -p $hap0_dir

    python $scripts_dir/parsesv2seg_speclong.py \
        -sv $fmt_sv \
        -p $hap0_dir/hap0 \
        -hid 0 \
        -ob $fixed_bam \
        -h0b $h0_bam \
        -h1b $h1_bam \
        -r $gene_work_dir/$hla.1.raw.fa || echo "Warning: parsesv2seg failed for hap0"

    # --- Hap1 ---
    hap1_dir=$gene_work_dir/hap1
    mkdir -p $hap1_dir

    python $scripts_dir/parsesv2seg_speclong.py \
        -sv $fmt_sv \
        -p $hap1_dir/hap1 \
        -hid 1 \
        -ob $fixed_bam \
        -h0b $h0_bam \
        -h1b $h1_bam \
        -r $gene_work_dir/$hla.2.raw.fa || echo "Warning: parsesv2seg failed for hap1"

    echo "Matching for $sample !"
    
    # Process Hap0 Graph
    h0_g=$hap0_dir/hap0_graph.txt
    h0_path=$hap0_dir/hap0_path.txt
    h0_c_path=$hap0_dir/hap0_c_path.txt
    h0_new_g=$hap0_dir/hap0_new_graph.txt
    sv_flag_file=$hap0_dir/hap0_nosv.flag

    if [ -s "$sv_flag_file" ]; then
        echo "$sv_flag_file exists and is not empty (No SVs)."
    else
        echo "$sv_flag_file does not exist or is empty. Running matching..."
        # Only run matching if graph file exists
        if [ -f "$h0_g" ]; then
            $matching \
            -b \
            --model 1 \
            -v 1 \
            -g $h0_g \
            -r $h0_path \
            -c $h0_c_path \
            -m $h0_new_g \
            --break_c || echo "Matching binary failed for hap0"

            # make fasta for hap0
            if [ -f "$h0_path" ]; then
                python $scripts_dir/make_fa_from_path.py \
                    $h0_path \
                    $hap0_dir/hap0_seg.txt \
                    $hap0_dir/hap0_seg.fa \
                    $hap0_dir/hap0_asm.fa \
                    0 
                
                if [ -f "$hap0_dir/hap0_asm.fa" ]; then
                    mv -f $hap0_dir/hap0_asm.fa $gene_work_dir/$hla.1.raw.fa
                fi
            fi
        fi
    fi

    # Process Hap1 Graph
    h1_g=$hap1_dir/hap1_graph.txt
    h1_path=$hap1_dir/hap1_path.txt
    h1_c_path=$hap1_dir/hap1_c_path.txt
    h1_new_g=$hap1_dir/hap1_new_graph.txt
    sv_flag_file=$hap1_dir/hap1_nosv.flag

    if [ -s "$sv_flag_file" ]; then
        echo "$sv_flag_file exists and is not empty (No SVs)."
    else
        echo "$sv_flag_file does not exist or is empty. Running matching..."
        if [ -f "$h1_g" ]; then
            $matching \
            -b \
            --model 1 \
            -v 1 \
            -g $h1_g \
            -r $h1_path \
            -c $h1_c_path \
            -m $h1_new_g \
            --break_c || echo "Matching binary failed for hap1"

            # make fasta for hap1
            if [ -f "$h1_path" ]; then
                python $scripts_dir/make_fa_from_path.py \
                    $h1_path \
                    $hap1_dir/hap1_seg.txt \
                    $hap1_dir/hap1_seg.fa \
                    $hap1_dir/hap1_asm.fa \
                    1 
                
                if [ -f "$hap1_dir/hap1_asm.fa" ]; then
                    mv -f $hap1_dir/hap1_asm.fa $gene_work_dir/$hla.2.raw.fa
                fi
            fi
        fi
    fi
fi

echo "Pipeline finished for $sample"