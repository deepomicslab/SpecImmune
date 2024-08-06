# python3 ../scripts/main.py -n HG00514 -o /mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results/ -j 15 -y pacbio -i IG_TR --hg38 /mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/chr14.fa -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/HG00514.IG.TR.extract.fastq.gz
# python3 ../scripts/main.py -n HG00514_all -o /mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results/ -j 15 -y pacbio -i IG_TR -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/HG00514.IG.TR.extract.fastq.gz
# python3 ../scripts/main.py -n test -o /mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results/ -k 1 -j 15 -y pacbio -i IG_TR -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/VDJ/one.fq
# python3 ../scripts/main.py -n large -o /mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results/ -k 10 -j 15 -y nanopore -i IG_TR -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/VDJ/nanoporelymphotrack.1000.fastq
# python3 ../scripts/main.py -n large2 -o /mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results/ -k 10 -j 15 -y nanopore -i IG_TR -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/VDJ/nanoporelymphotrack.filtered.fastq
# python3 ../scripts/main.py -n SRR24993720 -o /mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results/\
#  -j 15 -y nanopore -i IG_TR \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/VDJ/SRR24993720.fastq.gz

# python3 ../scripts/main.py -n NA18506_new -o /mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results/\
#  -j 15 -y nanopore -i IG_TR \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/VDJ/SRR19372746.fastq.gz

# python3 ../scripts/main.py -n NA19240 -o /mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results/\
#  -j 15 -y nanopore -i IG_TR \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/NA19240.IG.TR.extract.fastq.gz

# python3 ../scripts/main.py --hg38 ../CYP_ref/CYP.segment.fa -n NA19240_4 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/reads_cyp_hpc/NA19240.CYP.fastq.gz 

# python3 ../scripts/main.py --hg38 ../CYP_ref/CYP.segment.fa -n NA19239_4 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/reads_cyp_hpc/NA19239.CYP.fastq.gz 

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n NA19239_5 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/reads_cyp_hpc/NA19239.CYP.fastq.gz 

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n HG00733_5 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/reads_cyp_hpc/HG00733.CYP.fastq.gz 

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n NA19819 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/1000G_ont/downloads/NA19819/NA19819.CYP.fastq.gz

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n NA18565 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/CYP2D6/shuai/NA18565.CYP.fastq.gz

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n NA18632   -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y nanopore -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/CYP2D6/shuai/NA18632.CYP.fastq.gz --align_method_1 minimap2

# while IFS= read -r line; do

# ## sample name is the first column
# sample=$(echo $line | cut -d ' ' -f 1)
# echo $sample

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n $sample\
#    -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results2/\
#  -j 15 -y nanopore -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/CYP2D6/shuai/$sample.CYP.fastq.gz --align_method_1 minimap2

# done < "cyp/ont_truth.csv"

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n NA18642_3   -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y nanopore -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/NA18642_4/NA18642_4.noSV_1.fastq.gz --align_method_1 minimap2


# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n NA18980_2 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y nanopore -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/NA18980/NA18980.haplotype_2_39.fastq.gz --align_method_1 minimap2

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n NA07439 \
#  -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/amplicon/\
#  -j 15 -y nanopore -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/CYP2D6/SRR15476234.fastq.gz  --align_method_1 minimap2

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n test -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y nanopore -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/CYP2D6/shuai/HG00436.CYP.fastq.gz --mode 1 --align_method_1 minimap2

# file_list_file="/mnt/d/HLAPro_backup/Nanopore_optimize/data/reads_cyp_hpc/hgscv_hifi_cyp.list"
# # Read the file list file line by line
# while IFS= read -r file_path; do
#     # Extract the base name of the file
#     file_basename=$(basename "$file_path" ".fastq.gz")
#     echo $file_basename
#     if [ "1" == "1" ]; then

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n $file_basename\
#  -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/hgscv_hifi\
#  -j 15 -y pacbio -i CYP \
# -r $file_path --align_method_1 minimap2
#         # python $tool --mode 3 -n $file_basename -o $outdir -j 15 -y nanopore -i HLA -r $file_path --db ../db/ 
#     break
#     fi

# done < "$file_list_file"


# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n HG02572\
#  -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/hprc_hifi/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/CYP2D6/CYP/HG02572/HG02572.CYP.fastq.gz\
#  --mode 1 --align_method_1 minimap2



# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n NA07439 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y nanopore -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/amplicon/NA17300/NA07439.haplotype_3_13.fastq.gz.gz

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n HG00732_5 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/HG00732_5/HG00732_5.noSV.fastq.gz.gz

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n HG03781 \
#  -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/1000G_ont/downloads/HG03781/HG03781.CYP.fastq.gz

python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n HG01108 \
 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/1KGP\
 -j 15 -y nanopore -i CYP \
-r /mnt/d/HLAPro_backup/Nanopore_optimize/data/CYP2D6/shuai/HG01108.CYP.fastq.gz --align_method_1 minimap2


# python3 ../scripts/main.py --hg38 ../CYP_ref/CYP.segment.fa -n NA19239 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/NA19239.cyp.extract.fastq.gz 

# python3 ../scripts/main.py --hg38 ../CYP_ref/CYP.segment.fa -n SRR15476229 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/CYP2D6/SRR15476229.fastq.gz