#!/bin/bash

# Specify the file list file
file_list_file="/mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/pacbio_hifi_hla.list"
outdir="/mnt/d/HLAPro_backup/Nanopore_optimize/out_pac_kir"



tool="/home/wangshuai/softwares/SpecLong/scripts/main.py"

# Read the file list file line by line
while IFS= read -r file_path; do
    # Extract the base name of the file
    file_basename=$(basename "$file_path" ".hla.extract.fastq.gz")



    echo $file_basename
    if [ "1" == "1" ]; then
    # if [ "$file_basename" == "fredhutch-hla-GO85" ]; then

        python $tool -n $file_basename -o $outdir -j 15 -y pacbio -i HLA -r $file_path --db ../db/ --hete_p 0.3 --max_read_num 500 --candidate_allele_num 200
    fi

done < "$file_list_file"

### A -b 0.0009 100%   -b 0.0009 82%
