#!/usr/bin/bash

# Automatically set SKIRT_WD to the parent directory of this script,
# unless it is already set by the caller.
if [ -z "$SKIRT_WD" ]; then
  SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
  SKIRT_WD="$(realpath "${SCRIPT_DIR}/../")"
fi

# Check for required arguments
if [ $# -ne 2 ]; then
  echo "Usage: $0 <input_fasta> <output_directory>"
  exit 1
fi

# Input FASTA and output path
fasta_path="$1"
output_path="$2"
python_bin="$3"

# Reference files
nuc_query_path="${SKIRT_WD}/IPDKIR/kir_nuc.fasta"
gen_query_path="${SKIRT_WD}/IPDKIR/kir_gen.fasta"
Eds1_query_path="${SKIRT_WD}/IPDKIR/fasta/KIR3DS1_nuc.fasta"
zds4_fusion_path="${SKIRT_WD}/IPDKIR/fasta/KIR2DS4-00101e124567_3DL1-03501e89_nuc.fasta"
zdl3_fusion_path="${SKIRT_WD}/IPDKIR/fasta/KIR2DL3-00101e1245_2DP1-00201e6789_nuc.fasta"
novel_allele_path="${SKIRT_WD}/novel.fa"

# Output filename prefix
fastafile="$(basename -- "$fasta_path")"
outputname="${fastafile%.fa*}"
output_prefix="${output_path}/${outputname}"
paf_path="${output_prefix}.paf"

# Create output directory
mkdir -p "${output_path}"

# Merge known novel allele sequences into a single novel.fa file
cat "$zds4_fusion_path" > "$novel_allele_path"
cat "$zdl3_fusion_path" >> "$novel_allele_path"
cat "$nuc_query_path" >> "$novel_allele_path"

# Run minimap2 alignments
minimap2 -cx splice:hq -G16k -y --cs -t32 -2 "$fasta_path" "$novel_allele_path" > "$paf_path" 2>"${paf_path}.err"
minimap2 -cx splice:hq -G16k -k8 --end-seed-pen 5 -y --cs -t32 -2 "$fasta_path" "$Eds1_query_path" >> "$paf_path" 2>>"${paf_path}.err"
minimap2 -cx splice:hq -G16k -y --cs -t32 -2 "$fasta_path" "$gen_query_path" > "${output_prefix}.gen.paf" 2>"${output_prefix}.gen.paf.err"

# Run the main Python script to process the alignment results
echo "Using Python: $python_bin"
$python_bin "${SKIRT_WD}/scripts/skirt.py" -asm "$fasta_path" -g "${output_prefix}.gen.paf" "$paf_path" "$output_prefix" > "${output_prefix}.skirt.err"
# python3 "${SKIRT_WD}/scripts/skirt.py" -asm "$fasta_path" -g "${output_prefix}.gen.paf" "$paf_path" "$output_prefix" > "${output_prefix}.skirt.err"