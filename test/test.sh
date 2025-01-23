db=../db/
#### please change the db path if your db is not in the default path
ref=hg38_no_alt.fa



python ../scripts/main.py -r HLA/test_HLA_lite.fastq.gz -j 15 -i HLA -n test_HLA -o test --align_method_1 minimap2 -y pacbio --db $db
python ../scripts/main.py -r KIR/HG00131.KIR.fastq.gz -j 15 -i KIR -n test_KIR -o test --hete_p 0.2 --align_method_1 minimap2 -y pacbio --db $db
python3 main.py --hg38 $ref -n test_CYP -o $outdir -j 10 -y pacbio -i CYP -r CYP/HG03579.CYP.fastq.0.1.fq.gz --align_method_1 minimap2
python3 main.py --hg38 $ref -n test_IG_TR -o $outdir -j 10 -y pacbio -i IG_TR -r IG_TR/vdj.fq.gz
