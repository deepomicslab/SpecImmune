hg38=/path-to-hg38/hg38.fa
db=../db/
#### please change the db path if your db is not in the default path
python3 ../scripts/main.py -n test_IG_TR -o test -j 10 -y nanopore -i IG_TR -r IG_TR/vdj.fq.gz --db $db --hg38 $hg38
