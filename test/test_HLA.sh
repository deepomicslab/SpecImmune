db=../db/
#### please change the db path if your db is not in the default path

python ../scripts/main.py -r HLA/test_HLA_lite.fastq.gz -j 15 -i HLA -n test_HLA -o test --align_method_1 minimap2 -y pacbio --db $db
