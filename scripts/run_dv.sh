ref=$1
bam=$2
out_vcf=$3
out_gvcf=$4
threads=$5
region=$6
sif=$7
data_type='ONT_R104'
BIN_VERSION="1.5.0"
# data_type='PACBIO'


bam_dir="$(dirname "$bam")"
out_dir="$(dirname "$out_vcf")"
ref_dir="$(dirname "$ref")"
bam_file="$(basename "$bam")"
out_vcf_file="$(basename "$out_vcf")"
out_gvcf_file="$(basename "$out_gvcf")"
ref_file="$(basename "$ref")"

## if $sif is not given then
if [ -z "$sif" ]; then
    SIF_FILE="/home/xuedowang2/scratch/test_mhc_indiv_genome/deepvariant-1.5.0.sif"
else
    SIF_FILE=$sif
fi

# activate conda environment (optional)

echo ref_dir $ref_dir
# Run DeepVariant.

echo """
bam: $bam
out_vcf: $out_vcf
out_gvcf: $out_gvcf
threads: $threads
region: $region
data_type: $data_type
ref: $ref
sif: $sif

"""

singularity run -B /usr/lib/locale/:/usr/lib/locale/ \
  -B "$bam_dir":"/bam_dir" \
  -B "$out_dir":"/out_dir" \
  -B "$ref_dir":"/ref_dir" \
  "${SIF_FILE}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=$data_type \
  --ref=/ref_dir/$ref_file \
  --reads=/bam_dir/$bam_file \
  --regions=$region \
  --output_vcf=/out_dir/$out_vcf_file \
  --output_gvcf=/out_dir/$out_gvcf_file \
  --intermediate_results_dir /out_dir/intermediate_results_dir \
  --num_shards=$threads \
  --sample_name='SAMPLE' \
  --make_examples_extra_args="--allow_empty_examples=true"

