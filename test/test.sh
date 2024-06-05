db=../db/

python ../scripts/main.py -r test.fastq.gz -j 15 -i HLA -n test_HLA -o test --db $db


python3 ../evaluation/assess_typing.py -i HLA --true test.HLA.hap.alleles.txt --infer test/test_HLA/test_HLA.HLA.type.result.txt 
