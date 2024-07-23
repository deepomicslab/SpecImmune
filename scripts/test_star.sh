#!/bin/bash
#SBATCH --partition=tiny
#SBATCH --job-name="star"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=log.speclong.mg.log
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00
SAMPLE_LIST=$1
eval "$(conda shell.bash hook)"
# conda activate /scratch/project/cs_shuaicli/wxd/app/spec_env
threads=1
src=..//
reads_dir=//mnt/d/HLAPro_backup/Nanopore_optimize/data/reads_cyp_hpc/
outdir=//mnt/d/HLAPro_backup/Nanopore_optimize/data/reads_cyp_hpc/
REFERENCE_GENOME="/mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa"
db=..//db/
if [[ ! -d "${outdir}" ]]; then
    mkdir -p "${outdir}"
fi




# Function to process a single sample for all specified genes
process_sample() {
    local sample_id=$1
    echo "Processing sample ${sample_id} ..."
    sample_outdir="${outdir}/${sample_id}"
    if [[ ! -d "${sample_outdir}" ]]; then
        mkdir -p "${sample_outdir}"
    fi

    fq=$reads_dir/$sample_id.CYP.fastq.gz

    minimap2 -t $threads -ax map-hifi $REFERENCE_GENOME $fq | samtools view -bS - | samtools sort -@ $threads -o $sample_outdir/$sample_id.sorted.bam
    samtools index $sample_outdir/$sample_id.sorted.bam
    bam=$sample_outdir/$sample_id.sorted.bam
    extract_reads=$sample_outdir/$sample_id.extracted_clip.fa
    python extract_clip_seq.py $bam $extract_reads 200 false

    clipbam=$sample_outdir/$sample_id.sorted.clip.bam
    minimap2 -t $threads -ax map-hifi $REFERENCE_GENOME $extract_reads | samtools view -bS - | samtools sort -@ $threads -o $clipbam
    samtools index $clipbam
    
    mergedbam=$sample_outdir/$sample_id.merge_clip.bam
    samtools merge $mergedbam $bam $clipbam -f
    samtools index $mergedbam

    bash run_dv.sh $REFERENCE_GENOME $mergedbam $sample_outdir/$sample_id.merge.dv.vcf $sample_outdir/$sample_id.merge.dv.gvcf $threads chr22
    vcf=$sample_outdir/$sample_id.merge.dv.vcf
    bam=$mergedbam



    # bash run_dv.sh $REFERENCE_GENOME $bam $sample_outdir/$sample_id.dv.vcf $sample_outdir/$sample_id.dv.gvcf $threads chr22
    # vcf=$sample_outdir/$sample_id.dv.vcf
    # call vcf 

    echo """
    python $src/packages/PStrain/scripts/single_species_long.py -b $bam \
    -v $vcf \
    -o $sample_outdir \
    -n $sample_id\_merge
    
    """

    python $src/packages/PStrain/scripts/single_species_long.py -b $bam \
    -v $vcf \
    -o $sample_outdir \
    -n $sample_id\_merge

}

# Read the sample list file and process each sample
while IFS=' ' read -r sample_id; do
    echo "Processing sample ${sample_id} ..."
    # Ensure both sample_id and sample_cram are not empty
    fq=$reads_dir/$sample_id.CYP.fastq.gz
    if [[ -z "${sample_id}" || ! -f "${fq}" ]]; then
        echo $fq
        echo "ERROR: Invalid sample_id or sample_cram: ${sample_id}"
        continue
    fi
    process_sample "${sample_id}"
    
done < "${SAMPLE_LIST}"

