"""
Function 1 : assign long reads to the gene
"""

import sys
import os
import pysam
import gzip
import argparse
import pickle
from collections import defaultdict


from read_objects import My_read, My_locus, Read_bin
from determine_gene import get_focus_gene, get_folder_list
from db_objects import My_db
from alignment_modules import Read_Type, read_bin_map2db
from folder_objects import My_folder

def load_gene_distance():
    # read the dict stored by pickle in ../gene_dist/HLA.distance_matrix.pkl, check if the file exists
    matrix_file = f"{sys.path[0]}/../gene_dist/{args['i']}.distance_matrix.pkl"
    if os.path.isfile(matrix_file):
        with open(matrix_file, 'rb') as f:
            distance_matrix = pickle.load(f)
    else:
        print(matrix_file, "File does not exist.")
        distance_matrix = {}
    return distance_matrix


class Score_Obj():
    # determine which gene to assign
    def __init__(self):
        self.loci_score = {}
        self.primary_dict = defaultdict(set)
        self.read_loci = {}
    
    def add_read(self, read_obj):

        if read_obj.gap_ends_flag:  ## read with clip-match-clip, if the match is too short, this read could be noise
            if read_obj.match_num < 2000 and read_obj.loci_name == "HLA-DRB1":
                return
            if read_obj.loci_name != 'CYP2D7':
                if read_obj.match_num < my_db.gene_min_len[read_obj.loci_name] * 0.5 and read_obj.loci_name != "HLA-U":
                    # print ("too short",read_obj.read_name, read_obj.match_num, my_db.gene_min_len[read_obj.loci_name] * 0.5, read_obj.loci_name)
                    return
        if args['i'] == "CYP" and read_obj.match_num < 500:
            return
        # if read_obj.loci_name == "CYP2D6" and read_obj.match_num < 500:
        #     return

        if read_obj.identity < args["min_identity"]:
            return

        if read_obj.primary:
            self.primary_dict[read_obj.read_name].add(read_obj.allele_name)
        # else:
        #     return 0
        # if read_obj.match_num < 100:
        #     return 0

        
        if read_obj.read_name not in self.loci_score:
            self.loci_score[read_obj.read_name] = {}
        
        if read_obj.loci_name not in self.loci_score[read_obj.read_name]:
            locus_obj = My_locus()
            locus_obj.add_record(read_obj)

            self.loci_score[read_obj.read_name][read_obj.loci_name] = locus_obj

        else:
            self.loci_score[read_obj.read_name][read_obj.loci_name].add_record(read_obj)
        

    
    def assign(self, assign_file):
        f = open(assign_file, 'w')
        # print (len(self.loci_score))
        for read_name in self.loci_score: # for each read

            read_bin = Read_bin(self.loci_score[read_name])
            # identity_cutoff = 0.85, identity_diff = 0.01, dist_cutoff = 500
            assigned_locus = read_bin.assign_multiple(distance_matrix, read_name, args["min_identity"], args["max_identity_diff"], args["max_distance"]) 
            # print ("\n\n")
            # if "d51d293b-36fc-4a12-9978-46fbad8a7c18" == read_name:  ## fredhutch-hla-KOSE
            #     print (read_name, self.primary_dict[read_name], assigned_locus)            

            print (read_name, ",".join(assigned_locus), sep = "\t", file = f)
            self.read_loci[read_name] = assigned_locus
        f.close()
        return self.read_loci


class Binning():

    def __init__(self):
        self.db = my_db.full_db
        if args["seq_tech"] == "rna":
            self.cds_db = my_db.full_cds_db
        self.sam = f"""{my_folder.sample_prefix}.db.bam"""
        
        if args["m"] != 2:
            read_bin_map2db(args, my_db, args['align_method'] )

        self.bamfile = pysam.AlignmentFile(self.sam, 'rb')   
        self.assign_file = f"{my_folder.sample_prefix}.read_binning.txt"

    def read_bam(self):
        # observe each read, assign it to gene based on alignment records
        scor = Score_Obj()
        for read in self.bamfile:
            if read.is_unmapped:
                continue
            # print (read)
            # read_obj = Read_Obj(read)
            read_obj = My_read()
            read_obj.load_bam(read)
            scor.add_read(read_obj)
            # if read_obj.read_name == "m64076_200603_055852/5440181/ccs":
            #     print (read_obj.read_name, read_obj.mismatch_rate, read_obj.allele_name, read_obj.gap_ends_flag, read_obj.match_num )
        read_loci = scor.assign(self.assign_file)
        for gene in gene_list:
            outfile = my_folder.reads_dir + '/%s.%s.fq'%(gene, args["a"])
            filter_fq(gene, read_loci, args["r"], outfile)
        print ("reads-binning done.")

def filter_fq(gene, dict, raw_fq, outfile):
    # output the assigned reads to the fastq file of each gene
    i = 0
    #gene = 'A'
    
    out = open(outfile, 'w')
    flag = False
    if raw_fq.split(".")[-1] == "gz":
        f = gzip.open(raw_fq,'rt')
    else:
        f = open(raw_fq)
    for line in f:
        line = line.strip()
        if i % 4 == 0:
            read_name = line.split()[0][1:]
            if read_name in dict.keys() and gene in dict[read_name]:
                flag = True
                num = 1
                print (line, file = out)
        elif flag:
            print (line, file = out)
            num += 1
            if num == 4:
                flag = False
        i += 1
    f.close()
    out.close()
    os.system('gzip -f %s'%(outfile))


if __name__ == "__main__":   

    parser = argparse.ArgumentParser(description="Read binning.", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    required.add_argument("-r", type=str, help="Long-read fastq file. PacBio or Nanopore.", metavar="\b")
    required.add_argument("-n", type=str, help="Sample ID", metavar="\b")
    required.add_argument("-o", type=str, help="The output folder to store the typing results.", metavar="\b", default="./output")
    required.add_argument("-i", type=str, help="HLA,KIR,CYP",metavar="\b", default="HLA")
    optional.add_argument("-p", type=str, help="The population of the sample [Asian, Black, Caucasian, Unknown, nonuse] for annotation. Unknown means use mean allele frequency in all populations. nonuse indicates only adopting mapping score and considering zero-frequency alleles.", metavar="\b", default="Unknown")
    optional.add_argument("-j", type=int, help="Number of threads.", metavar="\b", default=5)
    optional.add_argument("-d", type=float, help="Minimum score difference to assign a read to a gene.", metavar="\b", default=0.001)
    optional.add_argument("--min_identity", type=float, help="Minimum identity to assign a read.", metavar="\b", default=0.85)
    optional.add_argument("--max_distance", type=int, help="max distance diff between read and hg38.", metavar="\b", default=2000)
    optional.add_argument("--max_identity_diff", type=float, help="A read assigned to two loci if the identity difference is lower than this.", metavar="\b", default=0.1)
    optional.add_argument("-g", type=int, help="Whether use G group resolution annotation [0|1].", metavar="\b", default=0)
    optional.add_argument("-m", type=int, help="1 represents typing, 0 means only read assignment", metavar="\b", default=1)
    optional.add_argument("-k", type=int, help="The mean depth in a window lower than this value will be masked by N, set 0 to avoid masking", metavar="\b", default=5)
    optional.add_argument("-a", type=str, help="Prefix of filtered fastq file.", metavar="\b", default="long_read")
    optional.add_argument("-y", type=str, help="Read type, [nanopore|pacbio|pacbio-hifi].", metavar="\b", default="pacbio")
    optional.add_argument("--minimap_index", type=int, help="Whether build Minimap2 index for the reference [0|1]. Using index can reduce memory usage.", metavar="\b", default=1)
    optional.add_argument("--db", type=str, help="db dir.", metavar="\b", default=sys.path[0] + "/../db/")
    optional.add_argument("--strand_bias_pvalue_cutoff", type=float, help="Remove a variant if the allele observations are biased toward one strand (forward or reverse). Recommand setting 0 to high-depth data.", metavar="\b", default=0.01)
    # optional.add_argument("-u", type=str, help="Choose full-length or exon typing. 0 indicates full-length, 1 means exon.", metavar="\b", default="0")
    optional.add_argument("--align_method", type=str, help="bwa or minimap2", metavar="\b", default='bwa')
    optional.add_argument("--seed", type=int, help="seed to generate random numbers", metavar="\b", default=8)
    optional.add_argument("--max_depth", type=int, help="maximum depth for each HLA locus. Downsample if exceed this value.", metavar="\b", default=2000)
    optional.add_argument("-rt", "--RNA_type", type=str, help="traditional,2D,Direct,SIRV",metavar="\b", default="traditional")
    optional.add_argument("--seq_tech", type=str, help="Amplicon sequencing or WGS sequencing [wgs|amplicon|rna].", metavar="\b", default="wgs")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    Min_score = 0  #the read is too long, so the score can be very low.
    Min_diff = args["d"]  #0.001

    print ("start read binning...")

    # gene_list, interval_dict =  get_focus_gene(args)

    my_db = My_db(args)
    my_folder = My_folder(args)
    

    # db_folder=os.path.dirname(my_db.full_cds_db) if args["seq_tech"] == "rna" else os.path.dirname(my_db.full_db)
    db_folder = os.path.dirname(my_db.full_db)
    gene_list = get_folder_list(db_folder)

    
    read_type = Read_Type(args["seq_tech"], args["y"], args["RNA_type"])
    if args["seq_tech"] == "rna":
        bwa_para = read_type.get_bwa_param()
    else:
        minimap_para = read_type.get_minimap2_param()


    distance_matrix = load_gene_distance()
    ###assign reads
    if args["m"] == 10086:
        print ("skip assignment, just for testing")
    else:
        pbin = Binning()
        pbin.read_bam()        

    print ("Finished.")




