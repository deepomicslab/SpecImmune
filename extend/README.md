# Showcase how to extend SpecImmune to other genes

For easy extend, SpecImmune supports a built in gene family named extend, so just prepare the database, and update the gene list.
Let's take NHKIR genes for example.

### step 1: prepare the database
Download the NHKIR database
```
wget https://raw.githubusercontent.com/ANHIG/IPDNHKIR/refs/heads/Latest/NHKIR_gen.fasta
```

Make the db, all the extended genes should have a family name `extend`, so set `-i extend`
```
python make_db.py -i extend --extend_fa /home/shuaiw/methylation/data/hla/NHKIR_gen.fasta  -o /home/shuaiw/methylation/data/hla/db
```

### Test
Prepare the test data and run like
```
 python main.py -i extend --db /home/shuaiw/methylation/data/hla/db -r /home/shuaiw/methylation/data/hla/Mamu-KIR1D_test.fq.gz -n test -o /home/shuaiw/methylation/data/hla/extend_test

```



