# PStrain: An Iterative Microbial Strains Profiling Algorithm for Shotgun Metagenomic Sequencing Data

## Quick start

### PStrain with Metaphlan4
First, get the source code and construct the env
```
git clone https://github.com/wshuai294/PStrain.git --depth 1
cd PStrain/
conda env create --name pstrain -f pstrain_metaphlan4_env.yml
conda activate pstrain
```
Then get the metaphlan database
```
bash scripts/collect_metaphlan_datbase.sh -x mpa_vOct22_CHOCOPhlAnSGB_202403 -m 4 -d ./   #constructing the reference for Metaphlan4 
python3 scripts/PStrain.py --help
```
test PStrain by:
```
cd test/
bash test_PStrain_V40.sh  # remember edit the --bowtie2db and -x if you have downloaded the database yourself.
```
Also, you can install it using Bioconda directly by:
```
conda create --name pstrain -c bioconda -c conda-forge pstrain
conda activate pstrain
```
or install it in an existing environment:
```
conda install bioconda::pstrain
```

### PStrain with Metaphlan3
First, get the source code and construct the env
```
git clone https://github.com/wshuai294/PStrain.git --depth 1
cd PStrain/
conda env create --name pstrain3 -f pstrain_metaphlan3_env.yml
conda activate pstrain3
```
Then get the metaphlan database
```
bash scripts/collect_metaphlan_datbase.sh -x mpa_v31_CHOCOPhlAn_201901 -m 3 -d ./   #constructing the reference for Metaphlan3 
python3 scripts/PStrain.py --help
```
test PStrain by:
```
cd test/
bash test_PStrain_V30.sh  # remember edit the --bowtie2db and -x if you have downloaded the database yourself.
```


### PStrain with Metaphlan2
```
git clone https://github.com/wshuai294/PStrain.git --depth 1
cd PStrain/

conda env create --name pstrain2 -f pstrain_metaphlan2_env.yml
conda activate pstrain2
python3 scripts/PStrain_m2.py -h
```
Download `db/` and `packages/` from [Google Drive](https://drive.google.com/drive/folders/1NgBNlfduLwG61P3JMBZts6pVUJn3UR8n?usp=sharing) or [Zenodo](https://zenodo.org/doi/10.5281/zenodo.10457543), and put them in the same path with the scripts/ folder. Also, you can appoint the path of the dependencies by `--dbdir`, `--prior` and `--metaphlan2`. 

test PStrain (with Metaphlan2) by:
```
cd test/
bash test_run.sh
```

Note:
- PStrain supports both single-end and paired-end reads.

## Basic Usage  

### Main functions
| Scripts | Description |
| --- | --- |
|scripts/PStrain.py| Strain profiling based on the Metaphlan4 or Metaphlan3|
|scripts/PStrain_m2.py| Strain profiling based on the Metaphlan2|
|scripts/collect_metaphlan_datbase.sh| collect the Metaphlan4/3 database and build the bowtie2 index |
|scripts/single_species.py| phase variants for a single genome, based on the BAM and VCF file|
|scripts/merge.py| Merge the strain results for downstream analyses  |


### Run PStrain based on Metaphlan4 or Metaphlan3
```
usage: python3 PStrain.py [options] -c/--conf <config file> -o/--outdir <output directory> --bowtie2db <metaphlan database> -x <metaphlan db index>

    Example: python3 PStrain.py -c config.txt -o out --bowtie2db ../mpa_vOct22_CHOCOPhlAnSGB_202403/ -x mpa_vOct22_CHOCOPhlAnSGB_202403 # Metaphlan 4
             python3 PStrain.py --metaphlan_version 3 -c config.txt -o outdir/ --bowtie2db ../mpa_v31_CHOCOPhlAn_201901/ -x mpa_v31_CHOCOPhlAn_201901 # Metaphlan 3

    The config file should follow the format:

    //
    sample: [sample1_ID]
    fq1: [forward reads fastq]
    fq2: [reverse/mate reads fastq]
    //
    sample: [sample2_ID]
    fq1: [forward reads fastq]
    fq2: [reverse/mate reads fastq]
    ...

    Help information can be found by python3 PStrain.py -h/--help, config file format for single end reads , and additional information can be found    in README.MD or https://github.com/wshuai294/PStrain.


PStrain: profile strains in shotgun metagenomic sequencing reads.

required arguments:
  -c , --conf           The configuration file of the input samples. (default: None)
  -o , --outdir         The output directory. (default: None)
  --bowtie2db           Path to metaphlan bowtie2db. (default:
                        /home/wangshuai/softwares/PStrain/scripts/../mpa_vOct22_CHOCOPhlAnSGB_202403)
  -x , --metaphlan_index
                        metaphlan bowtie2db index. (default: mpa_vOct22_CHOCOPhlAnSGB_202403)

options:
  -h, --help            show this help message and exit
  --metaphlan_version   Metaphlan version [3 or 4]. Used to download the corresponding database. If you build the
                        metaphlan database yourself, just ignore this. (default: 4)
  --bowtie2             Path to bowtie2 binary. (default: bowtie2)
  --bowtie2-build       Path to bowtie2 binary. (default: bowtie2-build)
  --samtools            Path to samtools binary. (default: samtools)
  --metaphlan           Path to metaphlan. (default: metaphlan)
  --metaphlan_output_files
                        If you have run metaphlan already, give metaphlan result file in each line, the order should be
                        the same with config file. In particular, '--tax_lev s' should be added while running
                        metaphlan. (default: )
  -p , --proc           The number of process to use for parallelizing the whole pipeline, run a sample in each
                        process. (default: 1)
  -n , --nproc          The number of CPUs to use for parallelizing the mapping with bowtie2. (default: 1)
  -w , --weight         The weight of genotype frequencies while computing loss, then the weight of linked read type
                        frequencies is 1-w. The value is between 0~1. (default: 0.3)
  --lambda1             The weight of prior knowledge while rectifying genotype frequencies. The value is between 0~1.
                        (default: 0.0)
  --lambda2             The weight of prior estimation while rectifying second order genotype frequencies. The value is
                        between 0~1. (default: 0.0)
  --species_dp          The minimum depth of species to be considered in strain profiling step (default is 5).
                        (default: 5)
  --snp_ratio           The SNVs of which the depth are less than the specific ratio of the species's mean depth would
                        be removed. (default: 0.45)
  --qual                The minimum quality score of SNVs to be considered in strain profiling step. (default: 20)
  --similarity          The similarity cutoff of hierachical clustering in merge step. (default: 0.95)
  --elbow               The cutoff of elbow method while identifying strains number. If the loss reduction ratio is
                        less than the cutoff, then the strains number is determined. (default: 0.24)
  --gatk                Path to gatk binary. (default: /home/wangshuai/softwares/PStrain/scripts//GenomeAnalysisTK.jar)
  --picard              Path to picard binary. (default: /home/wangshuai/softwares/PStrain/scripts//picard.jar)
  --prior               The file of prior knowledge of genotype frequencies in the population. Not providing this is
                        also ok. (default: None)
```

If you have run metaphlan3 already for each file, to save runtime, you can build a file to record 
the metaphlan3 resultfile of each sample in each line (In particular,  `--tax_lev s` should be added while running metaphlan3.). 
Please afford this file to `PStrain.py` with the parameter `--metaphlan_output_files`.
The file should be like
```
sample1_metaphlan_output.txt
sample2_metaphlan_output.txt
sample3_metaphlan_output.txt
...
``` 
Then PStrain will skip running metaphlan, and use the existing metaphlan result files.
```
python3 PStrain.py --metaphlan_output_files exisiting_metaphlan_results.txt -c test.config -o test_dir 
```



For metaphlan3/4, we have not built prior genotype frequencies database which will not impact the performance very much. So just ignore --prior parameter. Also, you can build a prior genotype frequencies database by yourself. I will show how to do this in the following section.
In this way, PStrain will run Metaphlan to find the species profile in the sample and profile the strains in each species.





### Run PStrain based on Metaphlan2
Example:
```
Example: python3 PStrain_m2.py -p 10 -n 8 --similarity 0.95 -o ./outdir -c config.txt \
--dbdir /home/wangshuai/strain/00.simulation/15.needle/db \
--prior /home/wangshuai/strain/00.simulation/15.needle/db/prior_beta.pickle\
--metaphlan2  /path-to/metaphlan2
```
Command:
```
PStrain: profile strains in shotgun metagenomic sequencing reads.

required arguments:
  -c , --conf       The configuration file of the input samples. (default:
                    None)
  -o , --outdir     The output directory. (default: None)
  --dbdir           Path to marker gene database, can be downloaded from
                    https://zenodo.org/doi/10.5281/zenodo.10457543. (default:
                    /home/wangshuai/softwares/PStrain/scripts/../db/)

optional arguments:
  -h, --help        show this help message and exit
  -p , --proc       The number of process to use for parallelizing the whole
                    pipeline, run a sample in each process. (default: 1)
  -n , --nproc      The number of CPUs to use for parallelizing the mapping
                    with bowtie2. (default: 1)
  -w , --weight     The weight of genotype frequencies while computing loss,
                    then the weight of linked read type frequencies is 1-w.
                    The value is between 0~1. (default: 0.3)
  --lambda1         The weight of prior knowledge while rectifying genotype
                    frequencies. The value is between 0~1. (default: 0.0)
  --lambda2         The weight of prior estimation while rectifying second
                    order genotype frequencies. The value is between 0~1.
                    (default: 0.0)
  --species_dp      The minimum depth of species to be considered in strain
                    profiling step. (default: 5)
  --snp_ratio       The SNVs of which the depth are less than the specific
                    ratio of the species's mean depth would be removed.
                    (default: 0.45)
  --qual            The minimum quality score of SNVs to be considered in
                    strain profiling step. (default: 20)
  --similarity      The similarity cutoff of hierachical clustering in merge
                    step. (default: 0.95)
  --elbow           The cutoff of elbow method while identifying strains
                    number. If the loss reduction ratio is less than the
                    cutoff, then the strains number is determined. (default:
                    0.24)
  --bowtie2         Path to bowtie2 binary. (default: bowtie2)
  --bowtie2-build   Path to bowtie2 binary. (default: bowtie2-build)
  --samtools        Path to samtools binary. (default: samtools)
  --metaphlan       Path to metaphlan. (default: /home/wangshuai/softwares/PSt
                    rain/scripts/../packages/metaphlan2/metaphlan2.py)
  --gatk            Path to gatk binary. (default: /home/wangshuai/softwares/P
                    Strain/scripts//GenomeAnalysisTK.jar)
  --picard          Path to picard binary. (default:
                    /home/wangshuai/softwares/PStrain/scripts//picard.jar)
  --prior           The file of prior knowledge of genotype frequencies in the
                    population. (default: /home/wangshuai/softwares/PStrain/sc
                    ripts/../db/prior_beta.pickle)
  --skip_prior      Whether use prior knowledge of genotype frequencies.
                    (default: False)
```

### Run PStrain with single-end data
To use single-end data, the only difference is in the config file, just leave the `fq2:` line empty. An config example for single-end data is:
```
//
sample: [sample1_ID]
fq1: [single-end reads fastq]
fq2: 
//
sample: [sample2_ID]
fq1: [single-end reads fastq]
fq2: 
```

## Interpret output
In the output directory specified by `-o`, there are directories for each sample named by sample ID and one directory named merge.  

In the directory of one specific sample, you will find a folder called `result`, in which you can find the below files and directories:
```
strain_RA.txt  
strain_merged_RA.txt 
species1_seq.txt
species2_seq.txt
...

```
The formats for `strain_RA.txt` and `strain_merged_RA.txt` are similar. An example is as below:
```
# Species	Species_RA	Strain_ID	Cluster_ID	Strain_RA
Escherichia_coli 88.00 str-1 clu-1 0.1
Escherichia_coli 88.00 str-2 clu-2 0.2
Escherichia_coli 88.00 str-3 clu-3 0.7
Leifsonia_rubra 12.00 str-1 clu-1 0.4
Leifsonia_rubra 12.00 str-2 clu-2 0.6
```
This means in this sample, there are two species: Escherichia coli with three strains and Leifsonia rubra with two strains.   

* Species_RA column means the species relative abundance in the sample. 

* Strain_ID column means the strains ID of each species in the sample.  

* Cluster_ID column means the strains' corresponding cluster in all of the samples.

* Strain_RA column means the strains abundance of each species in the sample. 

In the merge directory, you can find below files.
```
strain_number.txt
species_1_seq.txt
species_1_clu.txt
species_2_seq.txt
species_2_clu.txt
...
```
* `strain_number.txt` contains the strain number of each species in  all of the samples.

* `species_seq.txt` contains the strain sequence of in all of the samples.

* `species_clu.txt` contains the consensus sequence of each strains cluster in all of the samples.

The file listed below is one `*_clu.txt` example, the formate of `*_clu.txt` and `*_seq.txt` are similar.
```
# Gene	Locus	Ref	Alt	clu-1	clu-2	clu-3	clu-4	clu-5	
gi|545276179|ref|NZ_KE701962.1|:c548151-546829	90	T	G	1	0	0	1	1	
gi|545276179|ref|NZ_KE701962.1|:c548151-546829	303	G	A	0	1	0	1	1	
gi|545276179|ref|NZ_KE701962.1|:c548151-546829	310	A	T	0	0	1	0	1	
gi|545276179|ref|NZ_KE701962.1|:c548151-546829	312	C	A	1	1	1	0	0
...
```
Gene, Locus, Ref, and Alt columns represent the locus, following columns mean the composition of the consensus sequence of each cluster.

## Analyzing the results of PStrain
First, we can use merge.py to cluster the strains obtained from all the 
samples with user-defined parameters.
```
usage: python3 merge.py [options] -c/--conf <config file> -o/--outdir <output directory>

Merge the results from all of the samples.

required arguments:
  -c , --conf     The configuration file of the input samples, which is same as mentioned above.
  -o , --outdir   The output directory.

optional arguments:
  -h, --help      show this help message and exit
  --similarity    The similarity cutoff of hierachical clustering in merge step (default is 0.95).
```
Furthermore, we can adopt [**PStrain-tracer**](https://github.com/deepomicslab/PStrain-tracer) for downstream analysis and visualization.

## Assign the GenBank ID for identified strains
Please refer to [**PStrain-tracer**](https://github.com/deepomicslab/PStrain-tracer).

## Single Species Mode
If there is only one species, you can map reads and call SNPs by yourself. After that, you can use single_species.py from scripts/ folder to profile strains.

If you have the prior knowledge of strains number, you can add this info to argument -k/--strainsNum. Once the -k/--strainsNum argument is set to 1, the consensus sequence with major alleles will be generated.

Detailed usage information can be found as follows:
```
Usage: python3 single_species.py [options] -b/--bam <bam file> -v/--vcf <vcf file> -o/--outdir <output directory>
Example: python3 single_species.py -b mapped.bam -v mapped.vcf.gz --prior /home/wangshuai/strain/00.simulation/15.needle/db/prior_beta.pickle -o phase_result -k 2

Before using this script, you should align the reads to your reference, and call SNPs. Then you can provide the bam and vcf file. Once you have the prior knowledge of the strains number, you can provie it the to script by -k/--strainsNum parameter. Otherwise,it will choose the strains number by elbow method. 
Help information can be found by python3 single_species.py -h/--help, additional information can be found in README.MD or https://github.com/wshuai294/PStrain.

PStrain: phase strains in samples with single species.

required arguments:
  -b , --bam          The bam file of the input samples.
  -v , --vcf          The vcf file of the input samples.
  -o , --outdir       The output directory.

optional arguments:
  -h, --help          show this help message and exit
  -k , --strainsNum   The number of strains in the sample (default is chosen
                      by elbow method).
  --snp_dp            The minimum depth of SNPs to be considered in strain
                      profiling step (default is 5).
  --prior             The file of prior knowledge of genotype frequencies in
                      the population.
  --elbow             The cutoff of elbow method while identifying strains
                      number. If the loss reduction ratio is less than the
                      cutoff, then the strains number is determined.
  --qual              The minimum quality score of SNVs to be considered in
                      strain profiling step (default is 20).
  -w , --weight       The weight of genotype frequencies while computing loss,
                      then the weight of linked read type frequencies is 1-w.
                      The value is between 0~1. (default is 0.0)
  --lambda1           The weight of prior knowledge while rectifying genotype
                      frequencies. The value is between 0~1. (default is 0.0)
  --lambda2           The weight of prior estimation while rectifying second
                      order genotype frequencies. The value is between 0~1.
                      (default is 0.0)
```

## Run PStrain on your concerned species
You can build the species-specific pipeline easily by yourself, the code single_species.py in the scripts/ folder can be very useful.

You can build the pipeline as follow:
```
Step 1. Given the raw fastq file of a sample, run Metaphlan3 to get the species 
        abundance.
Step 2. For each species in the sample, extract the species' corresponding marker 
        genes from Metaphlan3 as the reference.
Step 3. For each species, map sample's reads to the reference of the species with 
        bowtie2 and call SNV with GATK, you will get the BAM and VCF file.
Step 4. For each species, run python3 single_species.py [options] -b/--bam <bam file>
         -v/--vcf <vcf file> -o/--outdir <output directory>.
Step 5. Then you get the strain profiles of each species in each sample, you can now 
        analyze them.   
```
You can also use other marker gene database in the same way.

## Build your own prior database
Although PStrain's performance is good enough without the prior database except for the gene with too low depth. You can build your own prior database to improve the     performance. PStrain will assume the prior beta is 0.5 if not finding prior knowledge.

To generate the prior database, please follow these steps:

1. Run PStrain.py (Metaphlan2) or PStrain_V30.py (Metaphlan3) directly to get the bam file in each sample, and call SNP with GATK in gvcf formate which will indicate the depth for each base.  
2. Merge the SNP by computing mean frequency of each base in all samples. And generate a file named "prior_beta.txt" in this formate:
```
locus   A       T       C       G
gi|378760316|ref|NZ_AGDG01000038.1|:8694-9629|Bacteroides_faecis_517    1       0       0       0
gi|410936685|ref|NZ_AJRF02000012.1|:54814-55320|Enterococcus_faecium_299        0       1       0       0
gi|139439993|ref|NZ_AAVN02000016.1|:7311-8633|Collinsella_aerofaciens_783       0.0102827763496144      0.0205655526992288      0       0.969151670951157
gi|239633764|ref|NZ_GG670286.1|:165060-166319|Enterococcus_gallinarum_626       0       0       0       1
gi|345651640|ref|NZ_JH114383.1|:276898-279936|Bacteroides_vulgatus_2725 0       0       0       1
GeneID:16385027|Lactococcus_phage_jm2_95        1       0       0       0
gi|488673303|ref|NZ_KB850950.1|:c879018-878179|Clostridium_hathewayi_511        0       0       0       1
```
3. Then you can convert it to pickle formate by
```
import pickle
popu={}
f=open('prior_beta.txt','r')
for line in f:
    break
for line in f:
    line=line.strip()
    array=line.split()
    atcg=array[1:]
    popu[array[0]]=atcg  
 with open('prior_beta.pickle', 'wb') as f:
     pickle.dump(popu, f, pickle.HIGHEST_PROTOCOL) 
```
Then please afford the file to PStrain with parameter --prior.



## Dependencies 
* Python 3.7.0 or above  
* Python libraries:  
    * Numpy  
    * Scipy  
    * Pandas  
    * Pulp
    * Pysam  
    * pickle

PStrain needs two database files, one is the metaphlan2 marker gene database, and the other is the prior genotype frequencies database. We have uploaded the two database in the db/ folder to the [Google Drive](https://drive.google.com/drive/folders/1NgBNlfduLwG61P3JMBZts6pVUJn3UR8n?usp=sharing), you can directly download the db/ folder, and put it in the same path with the scripts/ folder. PStrain could use them by default. Also, you can appoint the path of these two databases by corresponding arguments --dbdir and --prior.


PStrain relies on these software packages. PStrain could use them by default. Also, you can install them by yourself, and appoint them by --samtools, --bowtie2-build, --bowtie2, --metaphlan2, --gatk, and --picard arguments.   
* Third-party pipelines:
    * SamTools-1.3.1
    * bowtie2-2.3.1
    * metaphlan2/3
    * GATK-3.5 (java-1.8.0_241)
    * picard-tools-2.1.0

Note:
- dependencies can also be obtained at [Zenodo](https://zenodo.org/doi/10.5281/zenodo.10457543).

## Citation
Shuai Wang, Yiqi Jiang, Shuaicheng Li, PStrain: An Iterative Microbial Strains Profiling Algorithm for Shotgun Metagenomic Sequencing Data, Bioinformatics, , btaa1056, https://doi.org/10.1093/bioinformatics/btaa1056

## Getting help
Should you have any queries, please feel free to contact me, I will reply as soon as possible (swang66-c@my.cityu.edu.hk).
