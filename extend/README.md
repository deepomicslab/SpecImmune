# Showcase how to extend SpecImmune to other genes

Let's take NHKIR genes for example

### step 1: prepare the database
download the NHKIR database
```
wget https://raw.githubusercontent.com/ANHIG/IPDNHKIR/refs/heads/Latest/NHKIR_gen.fasta
```
The database looks like
```
>IPD-NHKIR:NHP00048 Mamu-KIR1D*002 14567 bp
ATGTCGCTCATGGTCGTTAGCGTGGCGTGTGTTGGTGAGTACTGGAAGGGAATCGACGGA
GGAAGCACTGGGTGGAGATCTGGGCCCGGAGGTGGAGATATGGGCCTAGAGGTGGAGATA
CGGGTCTAGAGGTGGCGTTACGGGCCTGGAGTGGAGATCTGGGCCTAGAGGTGGAGATAT
GGGTCTAGAGATGGAGTTATGGGCCTGAAGTGGAGATGTGGGCCTGGAGAGGAGATCTCG
```

Next, change the contig name to standard allele name in the fasta file. Use the given script like:
```
python change_db_name.py NHKIR_gen.fasta 
```  
The renamed file is like:
```
>Mamu-KIR1D*002 IPD-NHKIR:NHP00048 Mamu-KIR1D*002 14567 bp
ATGTCGCTCATGGTCGTTAGCGTGGCGTGTGTTGGTGAGTACTGGAAGGGAATCGACGGA
GGAAGCACTGGGTGGAGATCTGGGCCCGGAGGTGGAGATATGGGCCTAGAGGTGGAGATA
CGGGTCTAGAGGTGGCGTTACGGGCCTGGAGTGGAGATCTGGGCCTAGAGGTGGAGATAT
GGGTCTAGAGATGGAGTTATGGGCCTGAAGTGGAGATGTGGGCCTGGAGAGGAGATCTCG
```

All the extended genes should have a family name `extend`
```
python make_db.py -i extend --extend_fa /home/shuaiw/methylation/data/hla/NHKIR_gen.fasta  -o /home/shuaiw/methylation/data/hla/db
```

### step 2: get the focused gene list
Take `['Mamu-KIR1D']` for example, let's call the extended gene with a name of `extend` :
Add the genes to `get_focus_gene` function of `determine_gene.py` 
```
    elif gene_class == "extend":
        gene_list = ['Mamu-KIR1D']
```



