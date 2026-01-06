#!/bin/bash
# Batch run SpecImmune for simulated reads

# Configuration
SPECIMMUNE="/data4/wangxuedong/test_specimmune/SpecImmune/scripts/main.py"
DB="/data4/wangxuedong/test_specimmune/SpecImmune/db"
THREADS=15

# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir> [platforms] [depths]"
    echo ""
    echo "Examples:"
    echo "  $0 simulated_reads results"
    echo "  $0 simulated_reads results 'pacbio_hifi ont_r10' '10 20 30'"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
PLATFORMS=${3:-"pacbio_clr pacbio_hifi ont"}
DEPTHS=${4:-"2 4 6 8 10 20 30 40"}

# Platform to seqtype mapping
declare -A SEQTYPE=(
    ["pacbio_clr"]="pacbio"
    ["pacbio_hifi"]="pacbio-hifi"
    ["ont"]="nanopore"
    ["ont_r10"]="nanopore"
)

echo "=========================================="
echo "Batch SpecImmune HLA Typing"
echo "=========================================="
echo "Input:  $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "Platforms: $PLATFORMS"
echo "Depths: $DEPTHS"
echo "=========================================="
echo ""

total=0
success=0
failed=0

# Process each platform and depth
for platform in $PLATFORMS; do
    seqtype=${SEQTYPE[$platform]}
    
    for depth in $DEPTHS; do
        depth_dir="$INPUT_DIR/$platform/depth_${depth}x"
        
        if [ ! -d "$depth_dir" ]; then
            continue
        fi
        
        echo ""
        echo ">>> Processing: $platform / depth_${depth}x (type: $seqtype)"
        
        # Process each FASTQ file
        for fastq in "$depth_dir"/*.{fastq,fq,fastq.gz,fq.gz}; do
            if [ ! -f "$fastq" ]; then
                continue
            fi
            
            ((total++))
            
            filename=$(basename "$fastq")
            sample="${filename%.*}"
            sample="${sample%.*}"
            
            output="$OUTPUT_DIR/$platform/depth_${depth}x/$sample"
            
            echo "  [$total] $filename"
            
            # Run SpecImmune
            python $SPECIMMUNE \
                -r "$fastq" \
                -j $THREADS \
                -i HLA \
                -n "${platform}_${depth}_${sample}" \
                -o "$output" \
                --align_method_1 minimap2 \
                -y "$seqtype" \
                --db "$DB" 
            
            if [ $? -eq 0 ]; then
                echo "      ✓ Success"
                ((success++))
            else
                echo "      ✗ Failed"
                ((failed++))
            fi
        done
    done
done

echo ""
echo "=========================================="
echo "Summary"
echo "=========================================="
echo "Total:   $total"
echo "Success: $success"
echo "Failed:  $failed"
echo "Output:  $OUTPUT_DIR"
echo "=========================================="