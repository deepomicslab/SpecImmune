db=/mnt/d/HLAPro_backup/Nanopore_optimize/SpecComplex/db/

# python ../scripts/long_read_typing.py -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/hifi_fq/HG00514_pacbio.fastq.gz -j 15 -i HLA -m 0 -n HG00514 -o test --db $db

# python ../scripts/select_best_reference_alleleV2.py -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/HG00514.1.hla.extract.fastq.gz -j 15 -i HLA -n HG00514_HLA -o test --db $db
# python ../scripts/main.py -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/HG00514.1.hla.extract.fastq.gz -j 15 -i HLA -n HG00514_HLA -o test --db $db

python ../scripts/main.py -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/HG00514.1.kir.extract.fastq.gz -j 15 -i KIR -n HG00514_KIR -o test --db $db

# python ../scripts/main.py -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/HG00514.1.cyp.extract.fastq.gz -j 15 -i CYP -n HG00514_CYP -o test --db $db

# python ../scripts/main.py -r /mnt/d/HLAPro_backup/Nanopore_optimize/data/complex_reads/HG00514.1.IG.TR.extract.fastq.gz -i IG_TR -n HG00514_IG -o test -j 15 --db $db