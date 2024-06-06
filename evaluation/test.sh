#!/bin/bash

# Specify the file list file
file_list_file="/mnt/d/HLAPro_backup/Nanopore_optimize/data/fastq.list"
outdir="/mnt/d/HLAPro_backup/Nanopore_optimize/output5"

# file_list_file="/mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/pacbio_hifi_hla.list"
# outdir="/mnt/d/HLAPro_backup/Nanopore_optimize/output_pac"



tool="/home/wangshuai/softwares/SpecLong/scripts/main_test.py"

# Read the file list file line by line
while IFS= read -r file_path; do
    # Extract the base name of the file
    file_basename=$(basename "$file_path" ".fastq.gz")
    # file_basename=$(basename "$file_path" ".hla.extract.fastq.gz")


    echo $file_basename
    if [ "1" == "1" ]; then
    # if [ "$file_basename" == "fredhutch-hla-RSH" ]; then

        # python $tool -n $file_basename -o $outdir -j 15 -y pacbio -i HLA -r $file_path --db /mnt/d/HLAPro_backup/Nanopore_optimize/SpecComplex/db/ 

        # python $tool --mode 3 -n $file_basename -o $outdir -j 15 -y nanopore -i HLA -r $file_path --db /mnt/d/HLAPro_backup/Nanopore_optimize/SpecComplex/db/ 

        python $tool --mode 3 -n $file_basename -o $outdir -j 15 -y nanopore -i HLA -r $file_path --db ../db/ 
        # python ~/softwares/SpecLong/scripts/select_best_reference_alleleV2.py -b 0.0007 -n $file_basename  -o output0 -j 15 -y nanopore
        #python /mnt/d/HLAPro_backup/Nanopore_optimize/SpecHLA/script//refine_typing.py -n $file_basename  -o output/$file_basename/ 
        # python ~/softwares/SpecHLA/script/long_read_typing.py -n $file_basename -r $file_path  -j 15 -o output0 #-m 10086  #--strand_bias_pvalue_cutoff 0 #   --max_depth 500
        #python SpecHLA/script/long_read_typing.py -k 5 -g 1 -n $file_basename -r $file_path  -j 15
        #perl /mnt/d/HLAPro_backup/Nanopore_optimize/SpecHLA/script/whole/annoHLA.pl -p nonuse   -s $file_basename -i output/$file_basename -r tgs 
        #python ~/softwares//SpecHLA/script//refine_typing.py -n $file_basename  -o output/$file_basename/
    fi
    #cat output/$file_basename/hla.result.txt
    #break;
        # Add your code here to perform actions on the file

done < "$file_list_file"

### A -b 0.0009 100%   -b 0.0009 82%
