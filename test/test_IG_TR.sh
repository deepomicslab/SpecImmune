hg38=/path-to-hg38/hg38.fa
db=../db/
#### please change the db path if your db is not in the default path
echo "only ten gene loci are included in the test data"
python3 ../scripts/main.py -n test_IG_TR -o test -j 10 -y nanopore -i IG_TR -r IG_TR/NA20300.first_10_igtr.fq.gz --db $db --hg38 $hg38
### python3 ../scripts/main.py --hg38 /home/shuaiw/borg/hg38/hg38_no_alt.fa  --db /home/shuaiw/methylation/data/spec/db -n test_IG_TR -o test -j 10 -y nanopore -i IG_TR -r IG_TR/NA20300.first_10_igtr.fq.gz
