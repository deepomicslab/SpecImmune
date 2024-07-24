import sys
import pandas as pd

from read_binning import filter_fq

def get_read_meta(meta_file):
    read_type_dict = {}
    sv_type_list = set()
    
    with open(meta_file, 'r') as f:
        df = pd.read_csv(f)
        for index, row in df.iterrows():
            read_name = row["hifi_read"]
            SV = row["SV"]
            SV_label = row["SVlabel"]
            HP = row["HP"]
            if HP == '-1':
                HP='noHap'
            read_type_dict[read_name] = HP
            sv_type_list.add(HP)
            # skip if SV_label is empty
            # print (f"{read_name} {SV_label} {SV}")
            # if pd.isna(SV_label):
            #     read_type_dict[read_name] = ['noSV']
            #     sv_type_list.add('noSV')
            # else:
            #     # print (f"{read_name} {SV_label}", row["SV"])
            #     # pure_sv_lable = SV_label.split("_")[0]
            #     # print(f"{read_name} {pure_sv_lable}")
            #     read_type_dict[read_name] = [SV]
            #     sv_type_list.add(SV)
    return read_type_dict, sv_type_list

def get_read_bam(labeled_bam):
    read_type_dict = {}
    sv_type_list = set()
    ## read bam with pysam
    import pysam
    bamfile = pysam.AlignmentFile(labeled_bam, "rb")
    for read in bamfile:
        read_name = read.query_name
        ## check if the read has the HP tag
        if not read.has_tag('HP'):
            HP = 'noHap'
        HP = read.get_tag('HP')
        ## replace space and * to _
        HP = HP.replace(' ', '_')
        HP = HP.replace('*', '')
        HP = HP.replace(':', '')

        read_type_dict[read_name] = HP
        sv_type_list.add(HP)

    return read_type_dict, sv_type_list

if __name__ == "__main__":   
    prefix = sys.argv[1]
    raw_fq = sys.argv[2]

    meta_file = f"{prefix}.readMeta.csv"
    labeled_bam = f"{prefix}_labeledReads.bam"
    # raw_fq = '/mnt/d/HLAPro_backup/Nanopore_optimize/data/reads_cyp_hpc/HG00733.CYP.fastq.gz'

    # read_type_dict, sv_type_list = get_read_meta(meta_file)
    read_type_dict, sv_type_list = get_read_bam(labeled_bam)
    print (sv_type_list)

    for sv_label in sv_type_list:
        outfile = f"{prefix}.{sv_label}.fastq.gz"
        filter_fq(sv_label, read_type_dict, raw_fq, outfile)
