db=../db/
#### please change the db path if your db is not in the default path

python ../scripts/main.py -r CYP/HG03579.CYP.fastq.0.1.fq.gz -j 15 -i CYP -n test_CYP -o test -y pacbio --db $db --hg38 $hg38
