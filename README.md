# SpecImmune: A Tool for Accurate Typing of Diverse Immune-Related Gene Families from Long-Read Sequencing Data

**SpecImmune** is a bioinformatics software tool designed to accurately type five key immune-related gene families—**HLA, KIR, IG, TCR, and CYP**—from long-read sequencing data. These genes are critical for human immune functions and drug metabolism, but their genetic complexity makes them difficult to decode using traditional short-read sequencing methods. SpecImmune leverages the advantages of long-read sequencing technologies, such as **Nanopore** and **PacBio**, to provide highly accurate genotyping of these gene families.

### Key features of **SpecImmune** include:

1. **Accurate Typing of Immune-Related Genes**  
   SpecImmune can type **HLA, KIR, IG, TCR, and CYP** genes with high accuracy by categorizing long reads to specific loci and selecting the best-matching alleles from a reference database.

2. **Broad Compatibility**  
   It supports whole-genome sequencing (WGS) and targeted amplicon sequencing data from various long-read sequencing platforms like ONT and PacBio.

3. **Superior Performance**  
   SpecImmune outperforms existing tools such as **SpecImmune**, **HLA*LA**, and **Pangu** in typing accuracy, particularly for **HLA** and **CYP2D6** genes. It is also the only tool capable of typing **KIR** and **IG/TCR** from long-read data.

4. **Consensus Sequence Reconstruction**  
   It bins reads to alleles and reconstructs consensus sequences, ensuring high-quality haplotype sequences for each typed gene.

5. **Visualization of Results**  
   SpecImmune provides visual reports in an **IGV-like report**, allowing users to observe novel variants and the confidence of typing results, making it easier to interpret and validate findings.

6. **Efficient and User-Friendly**  
   SpecImmune is computationally efficient, making it suitable for use on personal computers, enabling convenient use in clinical settings.

## Quick start
### Install  
First, create the env with conda or mamba, and activate the env. 

**Use conda**
```
git clone git@github.com:deepomicslab/SpecImmune.git
cd SpecImmune/
conda env create -n SpecImmune -f environment.yml
conda activate SpecImmune
```

**For faster installation and environment management, we recommend using mamba.**

```
mamba env create -n SpecImmune -f environment.yml
mamba activate SpecImmune
```


Second, make the softwares in bin/ executable.
```
chmod +x -R bin/*
```

### Database construction
Third, build the allele database. You can build a database for all gene families, or just the ones you need.
For HLA, KIR, and IG, TCR:
```
python scripts/make_db.py -o SpecImmune/db  -i HLA

python scripts/make_db.py -o SpecImmune/db  -i KIR

python scripts/make_db.py -o SpecImmune/db  -i IG_TR
```
For CYP, download the complete pharmvar database at [Pharmvar](https://www.pharmvar.org/download), unzip it, merge the alleles of all CYP loci into a single `fasta` file, and afford the path to the `fasta` file to SpecImmune:
```
cat pharmvar-6.1.2.1/*/*.haplotypes.fasta >CYP.all.fasta
python scripts/make_db.py -o SpecImmune/db  -i CYP --CYP_fa  CYP.all.fasta
``` 


### Run & test
Perform SpecImmune with
```
python3 script/main.py -h
```
Please go to the `test/` folder, run SpecImmune with given scripts, and check results.

Note:
- SpecImmune now supports Linux and Windows WSL systems.



## Basic Usage  

### Main functions
| Scripts | Description |
| --- | --- |
|scripts/ExtractReads.sh| Extract gene-region-related reads from enrichment-free data.|
|scripts/main.py| Typing with naopore or pacbio data  |



### Extract gene-region-related reads
First extract gene reads with enrichment-free data. Otherwise, Typing would be slow. Map reads onto the `hg38`, then use `ExtracReads.sh` to extract reads by
```
usage() {
  echo "Usage: $0 -s <sample_id> -i <input_bam_or_cram> -g <gene_class> -o <output_directory> [-r <reference>]"
  echo "  -s  Sample ID or gene ID (required)"
  echo "  -i  Input BAM or CRAM file mapped to hg38 (required)"
  echo "  -g  Gene class, one of: HLA, KIR, CYP, IG_TR (required)"
  echo "  -o  Output directory (required)"
  echo "  -r  Reference file (required if input is CRAM)"
  exit 1
}
```

## Typing 

### HLA Typing

Perform four-field HLA typing by
```
python3 SpecImmune/main.py \
        -r <fastq> \
        -j <threads> \
        -i HLA \
        -n <sample_id> \
        -o <outdir> \
        --db SpecImmune/db \
        -y nanopore 
```
Example cmd:
```
xx
```

Perform full-resolution HLA typing with long-read RNA data
```
python3 /scratch/project/cs_shuaicli/wxd/app/SpecImmune/scripts/main.py \
        -r <fastq> \
        -j <threads> \
        -i HLA \
        -n <sample_id> \
        -o <outdir> \
        --db SpecImmune/db  \
        --seq_tech rna \
        --RNA_type traditional
```
### KIR Typing
```
python3 SpecImmune/main.py \
        -r <fastq> \
        -j <threads> \
        -i KIR \
        -n <sample_id> \
        -o <outdir> \
        --db SpecImmune/db \
        -y <datatype> 

```
Example cmd:
```
python3 SpecImmune/scripts/main.py -n $sample -o $outdir -j 10 -y pacbio -i KIR -r $fq --hete_p 0.2
```

### CYP Typing
```
python3 SpecImmune/main.py \
        -r <fastq> \
        -j <threads> \
        -i CYP \
        -n <sample_id> \
        -o <outdir> \
        -hg38 <no_alt_ref> \
        --db SpecImmune/db \
        -y <datatype>
```
Example cmd:
```
python3 SpecImmune/scripts/main.py --hg38 $ref -n $sample -o $outdir -j 10 -y nanopore -i CYP -r $fq --align_method_1 minimap2
```

### IG&TCR Typing
```
python3 SpecImmune/main.py \
        -r <fastq> \
        -j <threads> \
        -i IG_TCR \
        -n <sample_id> \
        -o <outdir> \
        --db SpecImmune/db \
        -y <datatype> \
        -hg38 <no_alt_ref>
```
Example cmd:
```
python3 SpecImmune/scripts/main.py --hg38 $ref -n $sample -o $outdir -j 10 -y pacbio -i IG_TR -r $fq
```


Full arguments can be seen in
```
usage: python3 main.py -h

Typing with only long-read data.

Required arguments:
  -r                  Long-read fastq file. PacBio or Nanopore. (default: None)
  -n                  Sample ID (default: None)
  -o                  The output folder to store the typing results. (default: ./output)
  -i                  HLA,KIR,CYP,IG_TR (default: HLA)

Optional arguments:
  -j                  Number of threads. (default: 5)
  --mode              4 represents all steps, 3 skip first, 2 skip two, 3, skipt three (default: 4)
  --analyze_method    phase/assembly (default: phase)
  -k                  The mean depth in a window lower than this value will be masked by N, set 0 to avoid masking (default: 5)
  -y                  Read type, [nanopore|pacbio|pacbio-hifi]. (default: pacbio)
  --db                db dir. (default: /run/media/wangxuedong/One Touch/SpecImmune_latest/SpecImmune/scripts/../db/)
  --hg38              referece fasta file, used by IG_TR typing, generated by extract_VDJ_segments_from_hg38.py (default: /run/media/wangxuedong/One Touch/SpecImmune_latest/SpecImmune/scripts/../VDJ_ref/IG_TR.segment.fa)
  -f, --first_run   set False for rerun (default: True)
  --min_identity      Minimum identity to assign a read. (default: 0.85)
  --hete_p            Hete pvalue. (default: 0.3)
  --candidate_allele_num 
                        Maintain this number of alleles for ILP step. (default: 200)
  --min_read_num      min support read number for each locus. (default: 2)
  --max_read_num      max support read number for each locus. (default: 500)
  -rt, --RNA_type   traditional,2D,Direct,SIRV (default: traditional)
  --seq_tech          Amplicon sequencing or WGS sequencing [wgs|amplicon]. (default: wgs)
  --align_method_1    align method in read binning, bwa or minimap2 (default: bwa)
  --align_method_2    align method in typing, bwa or minimap2 (default: minimap2)
  -v, --version         Display the version number (default: False)
  -h, --help
```




## Interpret output
In the denoted outdir, the results of each sample are saved in a folder named as the sample ID.  

In the directory of one specific sample, you will find the below files:
| Output | Description |
| --- | --- |
| sample_id.GENE.final.type.result.formatted.txt | GENE-typing results for all alleles |
| sample_id.pdf | Visualization report of the sample |
| Sequences/*fasta | Reconstructed allele sequences (the low-depth region is masked by N) |
| Genes_step2/*phased.vcf.gz | Phased vcf file for each gene  |
|Reads/*| Locus-specific binned reads|
|sample_id.GENE.type.result.txt| Best-matched allele pair from the database|
| Genes/* | alignment summary of binned reads mapped to the database  |


If you performed RNA-Seq typing, GENE-typing result file is formated as follow:
| Output | Description |
| --- | --- |
| sample_id.GENE.final.rna.type.result.txt | Typing results for RNA-Seq reads |
| sample_id.GENE.final.rna.type.result.g.txt | Typing results at G group resolution for RNA-Seq reads|


1. **An example for `sample_id.GENE.final.type.result.formatted.txt` (HLA/KIR) is as below:** 
```
# version: IPD-IMGT/HLA 3.56.0
Locus	Chromosome	Genotype	Match_info	Reads_num	Step1_type	One_guess
HLA-A	1	HLA-A*02:151	HLA-A*02:151|3516|1.0	25	HLA-A*02:151	HLA-A*02:151
HLA-A	2	HLA-A*03:01:01:01 HLA-A*03:01:01:01|3516|1.0	HLA-A*03:01:01:01	HLA-A*03:01:01:01
HLA-U	1	HLA-U*01:04	HLA-U*01:04|730|1.0	27	HLA-U*01:04	HLA-U*01:04
HLA-U	2	HLA-U*01:03	HLA-U*01:03|732|1.0	27	HLA-U*01:03	HLA-U*01:03
HLA-DMA	1	HLA-DMA*01:01:01:04	HLA-DMA*01:01:01:04|5013|1.0	37	HLA-DMA*01:01:01:04	HLA-DMA*01:01:01:04
HLA-DMA	2	HLA-DMA*01:01:01:02	HLA-DMA*01:01:01:02|5013|1.0	37	HLA-DMA*01:01:01:02 HLA-DMA*01:01:01:02
HLA-J	1	HLA-J*01:01:01:05	HLA-J*01:01:01:05|3544|1.0	51	HLA-J*01:01:01:05	HLA-J*01:01:01:05
HLA-J	2	HLA-J*01:01:01:04	HLA-J*01:01:01:04|3544|1.0	51	HLA-J*01:01:01:04	HLA-J*01:01:01:04
HLA-DPA1	1	HLA-DPA1*01:03:01:02	HLA-DPA1*01:03:01:02|9775|1.0	56	HLA-DPA1*01:03:01:02;HLA-DPA1*01:03:01:30;HLA-DPA1*01:03:17;HLA-DPA1*01:03:01:03;HLA-DPA1*01:03:01:32	HLA-DPA1*01:03:01:02
HLA-DPA1	2	HLA-DPA1*01:03:01:05;HLA-DPA1*01:03:01:15;HLA-DPA1*01:03:01:74	HLA-DPA1*01:03:01:05|9757|1.0;HLA-DPA1*01:03:01:15|9756|1.0;HLA-DPA1*01:03:01:74|9718|1.0	56	HLA-DPA1*01:03:01:05;HLA-DPA1*01:03:01:15	HLA-DPA1*01:03:01:05
HLA-DQA1	1	HLA-DQA1*01:03:01:02	HLA-DQA1*01:03:01:02|6492|1.0	58	HLA-DQA1*01:03:01:02 HLA-DQA1*01:03:01:02
HLA-DQA1	2	HLA-DQA1*01:05:01:01	HLA-DQA1*01:05:01:01|6485|1.0	58	HLA-DQA1*01:05:01:01	HLA-DQA1*01:05:01:01
HLA-DPA2	1	HLA-DPA2*01:01:01:01	HLA-DPA2*01:01:01:01|6743|1.0	46	HLA-DPA2*01:01:01:01	HLA-DPA2*01:01:01:01
HLA-DPA2	2	HLA-DPA2*01:01:01:02	HLA-DPA2*01:01:01:02|6743|1.0	46	HLA-DPA2*01:01:01:02	HLA-DPA2*01:01:01:02
HLA-G	1	HLA-G*01:01:01:01	HLA-G*01:01:01:01|3138|1.0	45	HLA-G*01:01:01:01	HLA-G*01:01:01:01
HLA-G	2	HLA-G*01:01:01:05	HLA-G*01:01:01:05|3138|1.0	45	HLA-G*01:01:01:05	HLA-G*01:01:01:05
HLA-P	1	HLA-P*02:01:01:01	HLA-P*02:01:01:01|2931|1.0	39	HLA-P*02:01:01:01	HLA-P*02:01:01:01
...
```

Interpret each column in the annotation line as
|Column| Description |
| --- | --- |
|1st| gene name|
|  2nd |  1 for haplotype 1, 2 for haplotype 2 |
| 3rd  | typing result, may contain ambiguity  |
|  4th | Matched alleles with their length and align ratio, separate by `\|` |
|  5th | Support reads count on the locus  |
|  6th | Best-matched allele in the database |
|  7th | One guess allele result, to handle ambiguity  |

2. **An example for `sample_id.CYP.merge.type.result.txt` (CYP) is as below:** 
```
#       version:        N/A
#       CYP2D6 Diplotype:       *101/*4
#       CYP2D6 Phenotype:       poor_metabolizer
#       CYP2D6 Activity_score:  0.0
#       CYP2D6 Detailed_diplotype       *101.001/*4.013
Locus   Chromosome      Genotype        Match_info      Reads_num       Step1_type      One_guess
CYP2A6  1       CYP2A6*15.002   CYP2A6*15.002|13902|0.999425    38      CYP2A6*15.002   CYP2A6*15.002
CYP2A6  2       CYP2A6*54.001   CYP2A6*54.001|13482|0.974203    38      CYP2A6*54.001   CYP2A6*54.001
NAT2    1       NAT2*7.002      NAT2*7.002|16952|0.999057       25      NAT2*7.002;NAT2*4;NAT2*4.001;NAT2*14;NAT2*14.001        NAT2*7.002
NAT2    2       NAT2*4.005      NAT2*4.005|16962|0.999646       25      NAT2*4.005;NAT2*7.002;NAT2*7.003        NAT2*4.005
...
```
Top lines with `#` start are the inferred Diplotype and phenotype for **CYP2D6**. Below lines are the typing results of other CYP loci, with format same as the above file.

3. **An example for `sample_id.IG_TR_typing_result.txt` (IG,TCR) is as below:** 
```
gene    depth   phase_set       allele_1        score_1 length_1        hap_1   allele_2        score_2 length_2        hap_2   hg38_chrom      hg38_len        variant_num     hete_variant_num
TRAV1-1 35.98   20983386        TRAV1-1*01      100.0   275     hap1    TRAV1-1*02      100.0   269     hap2    chr14   729     3       2
TRAV1-2 29.74   NA      TRAV1-2*01      100.0   267     hap1    TRAV1-2*01      100.0   267     hap2    chr14   689     1       0
TRAV2   31.39   NA      TRAV2*01        100.0   263     hap1    TRAV2*01        100.0   263     hap2    chr14   522     0       0
TRAV3   22.87   NA      TRAV3*01        100.0   286     hap1    TRAV3*01        100.0   286     hap2    chr14   608     1       1
TRAV4   36.98   20983386        TRAV4*01        100.0   277     hap1    TRAV4*01        100.0   277     hap2    chr14   830     1       1
```
Interpret each column in the annotation line as
|Column| Description |
| --- | --- |
|gene| gene locus name|
|  depth | average sequencing depth at this locus |
| phase_set  | phase block ID of this locus, NA if this locus is homozygous  |
|  allele_1/2 | first/second typed allele |
|  score_1/2 | alignment identity between the typed allele and the personlized haplotype   |
|  length_1/2 | alignment length between the typed allele and the personlized haplotype  |
|  hap_1 | locate on the first/second personlized haplotype  |
|  hg38_chrom | which chromosome the locus locate on the hg38 reference |
|  hg38_len | length of this locus on the hg38 reference  |
|  variant_num | Number of variants detected on this locus |
|  hete_variant_num | Number of heterozygous variants detected on this locus  |



## Dependencies 

### Systematic requirement
SpecImmune requires `conda 4.12.0+`, `cmake 3.16.3+`, and `GCC 9.4.0+` for environment construction and software installation.

### Programming 
* Python 3.8.15 or above  

### Third party packages
SpecImmune enables automatic installation of these third party packages using `conda` or `mamba`. 


## Getting help
Should you have any queries, please feel free to contact us, we will reply as soon as possible (swang66-c@my.cityu.edu.hk).
