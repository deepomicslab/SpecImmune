#!/bin/bash

###
### Gene typing with paired-end reads. This script can use PacBio, Nanopore,
### Hi-C, and 10X sequencing data to improve the phasing performance if provided.
### 
###
### Usage:
###   sh SpecComplex.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -t <sample.pacbio.fq.gz> -p <Asian> -i <HLA>
###
### Options:
###   -n        Sample ID. <required>
###   -1        The first fastq file. <required>
###   -2        The second fastq file. <required>
###   -o        The output folder to store the typing results.
###   -u        Choose full-length or exon typing. 0 indicates full-length, 1 means exon, 
###             default is to perform full-length typing.
###   -p        The population of the sample [Asian, Black, Caucasian, Unknown, nonuse]. 
###             Default is Unknown, meaning use mean frequency. nonuse indicates only adopting 
###             mapping score and considering alleles with frequency as zero.
###   -i        gene class [HLA, KIR, CYP].
###   -t        Pacbio fastq file.
###   -e        Nanopore fastq file.
###   -c        fwd hi-c fastq file.
###   -d        rev hi-c fastq file.
###   -x        Path of folder created by 10x demultiplexing. Prefix of the filenames of FASTQs
###             should be the same as Sample ID. Please install Longranger in the system env.
###   -w        How to use linkage info from allele imbalance [0, 0.5, 1], default is 0 that means 
###             not use, 0.5 means use both reads and imbalance info, 1 means only use imbalance info.
###   -j        Number of threads [5]
###   -m        The maximum mismatch number tolerated in assigning gene-specific reads. Deault
###             is 2. It should be set larger to infer novel alleles.
###   -y        The minimum different mapping score between the best- and second-best aligned gene. 
###             Discard the read if the score is lower than this value. Deault is 0.1. 
###   -v        True or False. Consider long InDels if True, else only consider short variants. 
###             Default is False. 
###   -q        Minimum variant quality. Default is 0.01. Set it larger in high quality samples.
###   -s        Minimum variant depth. Default is 5.
###   -a        Use this long InDel file if provided.
###   -r        The minimum Minor Allele Frequency (MAF), default is 0.05 for full length and
###             0.1 for exon typing.
###   -g        Whether use G-translate in annotation [1|0], default is 0.
###   -k        The mean depth in a window lower than this value will be masked by N, default is 5.
###             Set 0 to avoid masking.
###   -z        Whether only mask exon region, True or False, default is False.
###   -f        The trio infromation; child:parent_1:parent_2 [Example: NA12878:NA12891:NA12892]. 
###             Note: this parameter should be used after performing SpecHLA once.
###   -b        Whether use database for phasing [1|0], default is 1.
###   -h        Show this message.

help() {
    sed -rn 's/^### ?//;T;p' "$0"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
    help
    exit 1
fi

while getopts ":n:1:2:p:f:m:i:v:q:t:a:e:x:c:d:r:y:o:j:w:u:s:g:k:z:y:f:b:" opt; do
  case $opt in
    n) sample="$OPTARG"
    ;;
    1) fq1="$OPTARG"
    ;;
    2) fq2="$OPTARG"
    ;;
    p) pop="$OPTARG"
    ;;
    m) nm="$OPTARG"
    ;;
    i) gene_class="$OPTARG"
    ;;
    v) long_indel="$OPTARG"
    ;;
    q) snp_quality="$OPTARG"
    ;;
    t) tgs="$OPTARG"
    ;;
    a) sv="$OPTARG"
    ;;
    e) nanopore_data="$OPTARG"
    ;;
    x) tenx_data="$OPTARG"
    ;;
    c) hic_data_fwd="$OPTARG"
    ;;
    d) hic_data_rev="$OPTARG"
    ;;
    r) maf="$OPTARG"
    ;;
    o) given_outdir="$OPTARG"
    ;;
    j) num_threads="$OPTARG"
    ;;
    w) weight_imb="$OPTARG"
    ;;
    u) focus_exon="$OPTARG"
    ;;
    s) snp_dp="$OPTARG"
    ;;
    g) trans="$OPTARG"
    ;;
    k) mask_depth="$OPTARG"
    ;;
    z) mask_exon="$OPTARG"
    ;;
    y) mini_score="$OPTARG"
    ;;
    f) trio="$OPTARG"
    ;;
    b) use_database="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done


dir=$(cd `dirname $0`; pwd)
# export LD_LIBRARY_PATH=path/to/SpecHLA/spechla_env/lib/:$LD_LIBRARY_PATH
bin=$dir/bin
db=$dir/db
#hlaref=$db/immune.complex.gene.ref.extend.fasta

if [ ${given_outdir:-NA} == NA ]
  then
    outdir=$(pwd)/output/$sample
  else
    outdir=$given_outdir/$sample   
fi
focus_exon_flag=${focus_exon:-0} # default is perform full-length typing

echo Start profiling ${gene_class:-HLA} for $sample. 
mkdir -p $outdir
# exec >$outdir/$sample.log 2>&1 #redirect log info to the outdir
group='@RG\tID:'$sample'\tSM:'$sample
echo use ${num_threads:-5} threads.

:<<!
# ################ remove the repeat read name #################
python3 $dir/script/uniq_read_name.py $fq1 $outdir/$sample.uniq.name.R1.gz
python3 $dir/script/uniq_read_name.py $fq2 $outdir/$sample.uniq.name.R2.gz
fq1=$outdir/$sample.uniq.name.R1.gz
fq2=$outdir/$sample.uniq.name.R2.gz
# ###############################################################



# ################### assign the reads to original gene######################################################
echo map the reads to database to assign reads to corresponding genes.
license=$dir/bin/novoalign.lic
if [ $gene_class == "HLA" ];then
  database_prefix=$db/HLA/ref/HLA.extend
  genes=(HLA-A HLA-B HLA-C HLA-DMA HLA-DMB HLA-DOA HLA-DOB HLA-DPA1 HLA-DPB1 HLA-DPB2 HLA-DQA1 HLA-DQB1 HLA-DRA HLA-DRB1 HLA-DRB3 HLA-DRB4 HLA-DRB5 HLA-E HLA-F HLA-G HLA-H HLA-J HLA-K HLA-L HLA-P HLA-V HLA-DQA2 HLA-DPA2 HLA-N HLA-S HLA-T HLA-U HLA-W MICA MICB TAP1 TAP2 HFE)
fi
if [ $gene_class == "KIR" ];then
  database_prefix=$db/KIR/ref/KIR.extend.select
  genes=(KIR2DL1 KIR2DL2 KIR2DL3 KIR2DL4 KIR2DL5 KIR2DP1 KIR2DS1 KIR2DS2 KIR2DS3 KIR2DS4 KIR2DS5 KIR3DL1 KIR3DL2 KIR3DL3 KIR3DP1 KIR3DS1)
fi
if [ $gene_class == "CYP" ];then
  database_prefix=$db/CYP/ref/CYP.merge
  genes=(CYP19A1 CYP1A1 CYP1B1 CYP26A1 CYP2A13 CYP2A6 CYP2B6 CYP2C19 CYP2C8 CYP2C9 CYP2D6 CYP2F1 CYP2J2 CYP2R1 CYP2S1 CYP2W1 CYP3A4 CYP3A43 CYP4A22 CYP4B1 CYP4F2 CYP8A1 CYP3A5 CYP3A7)
fi

#:<<!
if [ -f "$license" ];then
    $bin/novoalign -d $database_prefix.ndx -f $fq1 $fq2 -F STDFQ -o SAM \
    -o FullNW -r All 100000 --mCPU ${num_threads:-5} -c 10  -g 20 -x 3  | $bin/samtools view \
    -Sb - | $bin/samtools sort -  > $outdir/$sample.map_database.bam
else
    bowtie2 --very-sensitive -p ${num_threads:-5} -k 30 -x $database_prefix.fasta -1 $fq1 -2 $fq2| \
    $bin/samtools view -bS -| $bin/samtools sort - >$outdir/$sample.map_database.bam
fi
$bin/samtools index $outdir/$sample.map_database.bam
#!

if [ $gene_class == "CYP" ]; then
        python3 $dir/script/assign_reads_to_genes.py -a $gene_class -1 $fq1 -2 $fq2 -n $bin -o $outdir -d ${mini_score:-0.1} -b ${outdir}/${sample}.map_database.bam -nm ${nm:-2} -s 0
else
        python3 $dir/script/assign_reads_to_genes.py -a $gene_class -1 $fq1 -2 $fq2 -n $bin -o $outdir -d ${mini_score:-0.1} -b ${outdir}/${sample}.map_database.bam -nm ${nm:-1} -s 1
fi

###############################################################################################################


############# align the gene-specific reads to the corresponding gene reference################################
ls $outdir/*.R1.fq.gz >$outdir/fqlist1
ls $outdir/*.R2.fq.gz >$outdir/fqlist2
if [ $gene_class == "KIR" ]; then

       perl $dir/script/generate.split.bowtie.pl $gene_class $outdir/fqlist1 $outdir/fqlist2 $outdir
       python $dir/script/requal.py -i $outdir/merge.$gene_class.bam -o $outdir/merge.$gene_class.requal.bam
       $bin/samtools index $outdir/merge.$gene_class.requal.bam
       $bin/freebayes -f $outdir/$gene_class.ref.merge.fa -p 3 $outdir/merge.$gene_class.requal.bam > $outdir/$sample.$gene_class.freebayes.snp.re.vcf
       perl $dir/script/filter.kir.snp.pl $outdir/$sample.$gene_class.freebayes.snp.re.vcf $outdir/$sample.$gene_class.realign.vcf
else
       perl $dir/script/generate.split.bwa.pl $gene_class $outdir/fqlist1 $outdir/fqlist2 $outdir
fi

#################################################################################################### 
#:<<!
generef=$outdir/$gene_class.ref.merge.fa

#################################### local assembly and realignment #################################
echo start realignment...
if [ $focus_exon_flag == 1 ];then #exon
  assemble_region=$db/$gene_class/select.region.exon.txt
else # full length
  assemble_region=$db/$gene_class/select.region.txt
fi
if [ $gene_class == "KIR" ];then
	bam=$outdir/merge.$gene_class.requal.bam
	vcf=$outdir/$sample.$gene_class.realign.vcf
else
   sh $dir/script/run.assembly.realign.sh $sample $outdir/merge.$gene_class.bam $outdir 70 $gene_class $assemble_region ${num_threads:-5}
   $bin/freebayes -a -f $generef -p 3 $outdir/$sample.$gene_class.realign.sort.bam > $outdir/$sample.$gene_class.realign.vcf 
   bam=$outdir/$sample.$gene_class.realign.sort.bam
fi

#rm -rf $outdir/$sample.$gene_class.realign.vcf.gz 
bgzip -f $outdir/$sample.$gene_class.realign.vcf
tabix -f $outdir/$sample.$gene_class.realign.vcf.gz
zless $outdir/$sample.$gene_class.realign.vcf.gz |grep "#" > $outdir/$sample.$gene_class.realign.filter.vcf

echo BAM and VCF are ready.
if [ $focus_exon_flag == 1 ];then #exon
    $bin/bcftools filter -R $db/$gene_class/$gene_class.ref.exon.extend.txt $outdir/$sample.$gene_class.realign.vcf.gz |grep -v "#" |awk '$6>1'  >> $outdir/$sample.$gene_class.realign.filter.vcf  
else # full length
    $bin/bcftools filter -R $db/$gene_class/gene_extend.bed $outdir/$sample.$gene_class.realign.vcf.gz |grep -v "#"|awk '$6>1' >> $outdir/$sample.$gene_class.realign.filter.vcf  
fi
vcf=$outdir/$sample.$gene_class.realign.filter.vcf
####################################################################################################


#################### assign long reads to gene ###################
if [ ${tgs:-NA} != NA ];then
    python3 $dir/script/long_read_typing.py -r ${tgs} -n $sample -i $gene_class -m 0 -o $outdir -j ${num_threads:-5} -a pacbio
fi
if [ ${nanopore_data:-NA} != NA ];then
    python3 $dir/script/long_read_typing.py -r ${nanopore_data} -n $sample -i $gene_class -m 0 -o $outdir -j ${num_threads:-5} -a nanopore
fi

# ###################### mask low-depth region #############################################
$bin/samtools depth -a $bam >$bam.depth  

if [ $focus_exon_flag == 1 ];then my_mask_exon=True; else my_mask_exon=${mask_exon:-False}; fi
bam=$outdir/$sample.$gene_class.realign.sort.bam
python3 $dir/script/mask_low_depth_region.py -c $bam.depth -o $outdir -w 20 -d ${mask_depth:-5} -f ${mask_exon:-False} -g $gene_class


# ###################### call long indel #############################################
if [ ${long_indel:-False} == True ] && [ $focus_exon_flag != 1 ]; #don't call long indel for exon typing
    then
    port=$(date +%N|cut -c5-9)
    bfile=$outdir/$sample.long.InDel.breakpoint.txt
    #generef=$db/$gene_class/ref/$gene_class.ref.extend.fasta
    if [ ${tgs:-NA} != NA ] # detect long Indel with pacbio
        then
        $bin/pbmm2 align -j ${num_threads:-5} $generef ${tgs:-NA} $outdir/$sample.movie1.bam --sort --sample $sample --rg '@RG\tID:movie1'
        $bin/samtools view -H $outdir/$sample.movie1.bam >$outdir/header.sam
 
        for gene in ${genes[@]}; do
                gene_ref=$db/split_ref/$gene.fasta
                $bin/pbmm2 align -j ${num_threads:-5} $gene_ref $outdir/$sample/$gene.pacbio.fq.gz $outdir/$gene.splitgene.bam --sort --sample $sample --rg '@RG\tID:movie1'
                $bin/samtools index $outdir/$gene.splitgene.bam
        done
        $bin/samtools merge -f -h $outdir/header.sam $outdir/$sample.pacbio.bam $outdir/*.splitgene.bam
        $bin/samtools index $outdir/$sample.pacbio.bam
	rm -rf $outdir/*.splitgene.bam


        $bin/pbsv discover -l 100 $outdir/$sample.pacbio.bam $outdir/$sample.svsig.gz
        $bin/pbsv call -t DEL,INS -m 150 -j ${num_threads:-5} $generef $outdir/$sample.svsig.gz $outdir/$sample.var.vcf
        python3 $dir/vcf2bp.py $outdir/$sample.var.vcf $outdir/$sample.tgs.breakpoint.txt
        cat $outdir/$sample.tgs.breakpoint.txt >$bfile
    else # detect long Indel with pair end data.
        sh $dir/../ScanIndel/run_scanindel_sample.sh $sample $bam $outdir $port
        cat $outdir/Scanindel/$sample.breakpoint.txt >$bfile
    fi
else
    bfile=nothing
fi
if [ ${sv:-NA} != NA ]
    then
    bfile=$sv
fi
# #############################################################################################

bam=$outdir/$sample.$gene_class.realign.sort.bam
#bam=$outdir/merge.$gene_class.requal.bam
vcf=$outdir/$sample.$gene_class.realign.filter.vcf

# ########### phase, link blocks, calculate haplotype ratio, give typing results ##############
if [ "$maf" == "" ];then
    if [ $focus_exon_flag != 1 ]; then
        my_maf=0.05
    else
        my_maf=0.1
    fi
else
    my_maf=$maf
fi

echo Minimum Minor Allele Frequency is $my_maf.

#genes=(MICA)
for gene in ${genes[@]}; do
gene_ref=$db/split_ref/$gene.fasta
python3 $dir/script/phase_variants.py \
  -o $outdir \
  -b $bam \
  -s $bfile \
  -v $vcf \
  --fq1 $outdir/$gene.R1.fq.gz \
  --fq2 $outdir/$gene.R2.fq.gz \
  --gene $gene \
  --freq_bias $my_maf \
  --snp_qual ${snp_quality:-0.01} \
  --snp_dp ${snp_dp:-5} \
  --ref $gene_ref \
  --tgs ${tgs:-NA} \
  --nanopore ${nanopore_data:-NA} \
  --hic_fwd ${hic_data_fwd:-NA} \
  --hic_rev ${hic_data_rev:-NA} \
  --tenx ${tenx_data:-NA} \
  --sa $sample \
  --weight_imb ${weight_imb:-0} \
  --exon $focus_exon_flag \
  --thread_num ${num_threads:-5} \
  --use_database ${use_database:-1} \
  --trio ${trio:-None}
done
# ##################################################################################################
!
# ############################ annotation ####################################
echo start annotation...
# perl $dir/annoHLApop.pl $sample $outdir $outdir 2 $pop
if [ $focus_exon_flag == 1 ];then #exon
    perl $dir/script/annoHLA.pl -s $sample -i $outdir -p ${pop:-Unknown} -g ${trans:-0} -r exon -c $gene_class
else
    perl $dir/script/annoHLA.pl -s $sample -i $outdir -p ${pop:-Unknown} -g ${trans:-0} -r whole -c $gene_class 
fi
# #############################################################################


rm -rf $outdir/*_hap2.fasta $outdir/*_hap1.fasta $outdir/*_hap2.fasta.out $outdir/*_hap1.fasta.out
rm -rf $outdir/fragment* $outdir/newref_insertion*
# sh $dir/../clear_output.sh $outdir/
cat $outdir/allele.phase.result.txt
echo $sample is done.

