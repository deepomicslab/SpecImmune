db=../db/
#### please change the db path if your db is not in the default path

python ../scripts/main.py -r HLA/test_HLA_lite.fastq.gz -j 15 -i HLA -n test_HLA -o test --align_method_1 minimap2 -y pacbio --db $db
python ../scripts/main.py -r KIR/KIR_dp50_acc98_1.fastq.gz -j 15 -i KIR -n test_KIR -o test --hete_p 0.2 --align_method_1 minimap2 -y pacbio --db $db
python ../scripts/main.py -r CYP/HG03579.CYP.fastq.0.1.fq.gz -j 15 -i CYP -n test_CYP -o test -y pacbio --db $db
