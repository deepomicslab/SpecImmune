import re, os


reported_alleles_dict = {
    'NA07439': '*4xN,*41',
    'NA10005': '*17,*29',
    'NA12244': '*35,*41',
    'NA17052': '*1,*1',
    'NA17058': '*10,*10',
    'NA17203': '*4,*35',
    'NA17246': '*4,*35',
    'NA17252': '*4,*5',
    'NA17280': '*2,*3',
    'NA17300': '*1,*6'
}


def read_sra(sra_file):
    with open(sra_file, 'r') as f:
        for line in f:
            if line.startswith('Run'):
                continue
            if not re.search("CYP2D6", line):
                continue
            line = line.strip().split(',')
            sra_id = line[0]
            sample_id = line[-4]
            # if sample_id != 'NA17280':
            #     continue
            # if sra_id != 'SRR15476227':
            #     # print(f"{sra_id} {sample_id}")
            #     continue
            if sample_id in reported_alleles_dict:
                print(f"{sra_id} {sample_id} {reported_alleles_dict[sample_id]}")
            # else:
            #     print(f"{sra_id} {sample_id} NA"
                cmd = f"""
                ID={sra_id}
                sample={sample_id}
                # prefetch -O {fq_dir} $ID
                # fastq-dump -O {fq_dir} --split-3 {fq_dir}/$ID/$ID.sra  --gzip
                # rm -r {fq_dir}/$ID

                python3 ../scripts/main.py --hg38 /mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n $sample \
                -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/amplicon\
                -j 15 -y pacbio -i CYP \
                -r {fq_dir}/$ID.fastq.gz --seq_tech amplicon
                """
                os.system(cmd)

                cmd = f"""
                    sample={sample_id}
                    ID={sra_id}
                    # fq=/mnt/d/HLAPro_backup/Nanopore_optimize/data/reads_cyp_hpc/$sample.CYP.fastq.gz
                    fq={fq_dir}/$ID.fastq.gz
                    outdir=/mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/star/
                    ref=/mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa
                    
                    # bwa mem -R "@RG\\tID:$sample\\tSM:$sample" -t 15 $ref $fq | samtools view -bS -F 0x800 -| samtools sort - >$outdir/$sample.bam
                    # samtools index $outdir/$sample.bam

                    # longshot -F -S --sample_id $sample  --bam $outdir/$sample.bam --ref $ref --out $outdir/$sample.longshot.vcf
                    # python ../scripts/CYP_star_caller.py cyp2d6 $outdir/$sample.longshot.vcf chr22:42126499-42130865  $outdir/$sample.star.tab


                    gatk3 -T HaplotypeCaller -R $ref -I $outdir/$sample.bam -o $outdir/$sample.vcf -L chr22:42126499-42130865
                    whatshap phase -o $outdir/$sample.phase.vcf.gz -r $ref --indels $outdir/$sample.vcf $outdir/$sample.bam

                    python ../scripts/CYP_star_caller.py cyp2d6 $outdir/$sample.phase.vcf.gz chr22:42126499-42130865 $outdir/$sample.star.tab

                """
                # print (cmd)
                # os.system(cmd)
                # break


sra_file = "./cyp/PRJNA754842.csv"
fq_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/CYP2D6"
read_sra(sra_file)

