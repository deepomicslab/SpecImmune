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

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n NA19207 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/CYP2D6/shuai/NA19207.CYP.fastq.gz

python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n NA12336 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
 -j 15 -y pacbio -i CYP \
-r /mnt/d/HLAPro_backup/Nanopore_optimize/data/CYP2D6/shuai/NA12336.CYP.fastq.gz

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n NA17300_test -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y nanopore -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/amplicon/NA17300/NA17300.haplotype_3_13.fastq.gz.gz

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n HG00732_5 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/HG00732_5/HG00732_5.noSV.fastq.gz.gz

# python3 ../scripts/main.py --hg38 //mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n HG03781 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/1000G_ont/downloads/HG03781/HG03781.CYP.fastq.gz


# python3 ../scripts/main.py --hg38 ../CYP_ref/CYP.segment.fa -n NA19239 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/NA19239.cyp.extract.fastq.gz 

# python3 ../scripts/main.py --hg38 ../CYP_ref/CYP.segment.fa -n SRR15476229 -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
#  -j 15 -y pacbio -i CYP \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/CYP2D6/SRR15476229.fastq.gz