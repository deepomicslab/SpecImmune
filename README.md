# SpecLong

## install

```
conda env create --name speclong -f environment.yml
```

## extract reads from WGS data
run.extract.reads.sh

## run
python3 main.py

```
usage: python3 main.py -h

HLA Typing with only long-read data.

Required arguments:
  -r                 Long-read fastq file. PacBio or Nanopore. (default: None)
  -n                 Sample ID (default: None)
  -o                 The output folder to store the typing results. (default: ./output)
  -i                 HLA,KIR,CYP,IG_TR (default: HLA)

Optional arguments:
  -j                 Number of threads. (default: 5)
  -k                 The mean depth in a window lower than this value will be masked by N, set 0 to avoid
                       masking (default: 5)
  -y                 Read type, [nanopore|pacbio]. (default: pacbio)
  --db               db dir. (default: /home/wangshuai/softwares/SpecLong/scripts/../db/)
  -f, --first_run  set False for rerun (default: True)
  -h, --help
```

## output
HLA,KIR,CYP: `HG00514.HLA.type.result.txt`, `HG00514.KIR.type.result.txt`, `HG00514.CYP.type.result.txt`
```
Locus   Chromosome      Allele  Reads_num
HLA-A   1       HLA-A*02:01:01:01;HLA-A*01:01:01:01;HLA-A*01:01:01:02N;HLA-A*01:01:01:03;HLA-A*01:01:01:04;HLA-A*01:01:01:05;HLA-A*01:01:01:06;HLA-A*01:01:01:07;HLA-A*01:01:01:08;HLA-A*01:01:01:09;HLA-A*01:01:01:10;HLA-A*01:01:01:11;HLA-A*01:01:01:12;HLA-A*01:01:01:13;HLA-A*01:01:01:14;HLA-A*01:01:01:16;HLA-A*01:01:01:17;HLA-A*01:01:01:18;HLA-A*01:01:01:19;HLA-A*01:01:01:20;HLA-A*01:01:01:21;HLA-A*01:01:01:22;HLA-A*01:01:01:23;HLA-A*01:01:01:24;HLA-A*01:01:01:25;HLA-A*01:01:01:26;HLA-A*01:01:01:27;HLA-A*01:01:01:28;HLA-A*01:01:01:29;HLA-A*01:01:01:30    15
HLA-A   2       HLA-A*02:01:01:01;HLA-A*01:01:01:01;HLA-A*01:01:01:02N;HLA-A*01:01:01:03;HLA-A*01:01:01:04;HLA-A*01:01:01:05;HLA-A*01:01:01:06;HLA-A*01:01:01:07;HLA-A*01:01:01:08;HLA-A*01:01:01:09;HLA-A*01:01:01:10;HLA-A*01:01:01:11;HLA-A*01:01:01:12;HLA-A*01:01:01:13;HLA-A*01:01:01:14;HLA-A*01:01:01:16;HLA-A*01:01:01:17;HLA-A*01:01:01:18;HLA-A*01:01:01:19;HLA-A*01:01:01:20;HLA-A*01:01:01:21;HLA-A*01:01:01:22;HLA-A*01:01:01:23;HLA-A*01:01:01:24;HLA-A*01:01:01:25;HLA-A*01:01:01:26;HLA-A*01:01:01:27;HLA-A*01:01:01:28;HLA-A*01:01:01:29;HLA-A*01:01:01:30    15
HLA-B   1       HLA-B*40:01:02:04;HLA-B*40:01:02:31     27
HLA-B   2       HLA-B*46:01:01:01       27
HLA-C   1       HLA-C*01:02:01:01       30
HLA-C   2       HLA-C*03:04:01:12;HLA-C*03:04:01:02     30
HLA-DMA 1       HLA-DMA*01:01:01:02     56
```

IG_TR: `HG00514_IG.IG.TR.allele.txt`:
```
sample  gene    allele  score   len     start   end     chr     hap
HG00514_IG      IGHD2-15        IGHD2-15*01     100.000 31      39027   39057   chr14_igh       hap1
HG00514_IG      IGHD2-21        IGHD2-21*01     100.000 28      29648   29675   chr14_igh       hap1
HG00514_IG      IGHD2-8 IGHD2-8*02      100.000 31      48281   48311   chr14_igh       hap1
HG00514_IG      IGHD2/OR15-2a   IGHD2/OR15-2a*01        100.000 31      57895   57925   chr14_igh       hap1
HG00514_IG      IGHD3-22        IGHD3-22*01     100.000 31      27176   27206   chr14_igh       hap1
HG00514_IG      IGHD3-9 IGHD3-9*01      100.000 31      45751   45781   chr14_igh       hap1
HG00514_IG      IGHD3/OR15-3a   IGHD3/OR15-3a*01        100.000 31      55428   55458   chr14_igh       hap1
HG00514_IG      IGHJ1   IGHJ1*01        100.000 52      6964    7015    chr14_igh       hap1
HG00514_IG      IGHJ2   IGHJ2*01        100.000 53      6756    6808    chr14_igh       hap1
HG00514_IG      IGHJ3   IGHJ3*02        100.000 50      6144    6193    chr14_igh       hap1
HG00514_IG      IGHJ4   IGHJ4*02        100.000 48      5772    5819    chr14_igh       hap1
HG00514_IG      IGHJ5   IGHJ5*02        100.000 51      5371    5421    chr14_igh       hap1
HG00514_IG      IGHJ6   IGHJ6*02        100.000 62      4756    4817    chr14_igh       hap1
HG00514_IG      IGHV1-18        IGHV1-18*01     100.000 296     311366  311661  chr14_igh       hap1
HG00514_IG      IGHV1-2 IGHV1-2*04      100.000 296     112598  112893  chr14_igh       hap1
HG00514_IG      IGHV1-24        IGHV1-24*01     100.000 296     402977  403272  chr14_igh       hap1
HG00514_IG      IGHV1-3 IGHV1-3*01      100.000 296     131126  131421  chr14_igh       hap1
```
