# python3 ../scripts/main.py -n HG00514 -o /mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results/ -j 15 -y pacbio -i IG_TR --hg38 /mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/chr14.fa -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/HG00514.IG.TR.extract.fastq.gz
# python3 ../scripts/main.py -n HG00514_all -o /mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results/ -j 15 -y pacbio -i IG_TR -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/HG00514.IG.TR.extract.fastq.gz
# python3 ../scripts/main.py -n test -o /mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results/ -k 1 -j 15 -y pacbio -i IG_TR -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/VDJ/one.fq
# python3 ../scripts/main.py -n large -o /mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results/ -k 10 -j 15 -y nanopore -i IG_TR -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/VDJ/nanoporelymphotrack.1000.fastq
# python3 ../scripts/main.py -n large2 -o /mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results/ -k 10 -j 15 -y nanopore -i IG_TR -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/VDJ/nanoporelymphotrack.filtered.fastq
# python3 ../scripts/main.py -n SRR24993720 -o /mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results/\
#  -j 15 -y nanopore -i IG_TR \
# -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/VDJ/SRR24993720.fastq.gz

python3 ../scripts/main.py -n NA18506 -o /mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results/\
 -j 15 -y nanopore -i IG_TR \
-r /mnt/d/HLAPro_backup/Nanopore_optimize/data/VDJ/SRR19372746.fastq.gz