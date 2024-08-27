threads=15
# samples=(HG00512 HG00513 HG00514 HG00731 HG00732 HG00733 NA19238 NA19239 NA19240)
samples=(HG00733)
src=/home/wangshuai/softwares//SpecLong/scripts/
reads_dir=/mnt/d/HLAPro_backup/Nanopore_optimize/data/reads_kir_hpc/
outdir=/mnt/d/HLAPro_backup/Nanopore_optimize/KIR_HGSCV2_hifi2
for sample in ${samples[@]}; do

    echo "Processing $sample"
    
    fq=$reads_dir/$sample.KIR.fastq.gz
    python3 $src/main.py -r $fq -n $sample -i KIR -j $threads -o $outdir #--align_method_1 minimap2 #--align_method_2 bwa #--align_method_1 minimap2
    # break


done


# python /home/wangshuai/softwares/SpecLong/scripts/main_test.py  -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/reads_hla_hpc/NA19239.HLA.fastq.gz -n NA19239 -i HLA -j 10 -o pacbio_hla
#python /home/wangshuai/softwares/SpecLong/scripts/main.py  -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/reads_kir_hpc/HG00512.KIR.fastq.gz -n HG00512 -i KIR -j 10
#python /home/wangshuai/softwares/SpecLong/scripts/main.py  -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/reads_kir_hpc/NA19239.KIR.fastq.gz -n NA19239 -i KIR -j 10
#python /home/wangshuai/softwares/SpecLong/scripts/main.py  -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/NA19238.kir.extract.fastq.gz -n NA19238 -i KIR -j 10
#python /home/wangshuai/softwares/SpecLong/scripts/main.py  -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/HG00514.kir.extract.fastq.gz -n HG00514 -i KIR
