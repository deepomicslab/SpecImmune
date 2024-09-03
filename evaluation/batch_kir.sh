threads=15
samples=(HG00512 HG00513 HG00514 HG00731 HG00732 HG00733 NA19238 NA19239 NA19240)
# samples=(NA19239)
src=/home/wangshuai/softwares//SpecLong/scripts/


# reads_dir=/mnt/d/HLAPro_backup/Nanopore_optimize/data/reads_kir_hpc/
# outdir=/mnt/d/HLAPro_backup/Nanopore_optimize/KIR_HGSCV2_hifi2
# for sample in ${samples[@]}; do
#     echo "Processing $sample"
#     fq=$reads_dir/$sample.KIR.fastq.gz
#     python3 $src/main.py -r $fq -n $sample -i KIR -j $threads -o $outdir --hete_p 0.2 #--align_method_1 minimap2 #--align_method_2 bwa #--align_method_1 minimap2
#     # break
# done



reads_dir=/mnt/d/HLAPro_backup/Nanopore_optimize/data/hprc_hifi_kir/
outdir=/mnt/d/HLAPro_backup/Nanopore_optimize/KIR_HPRC_hifi
file_list_file=/mnt/d/HLAPro_backup/Nanopore_optimize/data/hprc_hifi_kir/hprc_hifi_kir.list
while IFS= read -r file_path; do
    # Extract the base name of the file
    sample=$(basename "$file_path" ".KIR.fastq.gz")
    # file_basename=$(basename "$file_path" ".hla.extract.fastq.gz")


    echo $sample
    # if [ "1" == "1" ]; then
    if [ "$sample" == "HG00621" ]; then
        echo "Processing $sample"
        
        
        fq=$reads_dir/$sample.KIR.fastq.gz
        # sample=HG02559_2
        python3 $src/main.py -r $fq -n $sample -i KIR -j $threads -o $outdir --hete_p 0.2 #--align_method_1 minimap2 #--align_method_2 bwa #--align_method_1 minimap2
        # break
        
    fi

done < "$file_list_file"