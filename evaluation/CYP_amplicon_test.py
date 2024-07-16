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
            if sample_id != 'NA17280':
                continue
            if sample_id in reported_alleles_dict:
                print(f"{sra_id} {sample_id} {reported_alleles_dict[sample_id]}")
            # else:
            #     print(f"{sra_id} {sample_id} NA"
                cmd = f"""
                ID={sra_id}
                sample={sample_id}
                prefetch -O {fq_dir} $ID
                fastq-dump -O {fq_dir} --split-3 {fq_dir}/$ID/$ID.sra  --gzip

                python3 ../scripts/main.py --hg38 ../CYP_ref/CYP.segment.fa -n $sample -o /mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/\
                -j 15 -y pacbio -i CYP \
                -r {fq_dir}/$ID.fastq.gz
                """
                os.system(cmd)


sra_file = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/CYP2D6/PRJNA754842.csv"
fq_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/CYP2D6"
read_sra(sra_file)

