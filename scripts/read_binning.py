"""
Function 1 : assign long reads to the gene
"""

import sys
import os
import pysam
import gzip
import argparse
from collections import defaultdict


from downsample_bam import main
from read_objects import My_read, My_locus, Read_bin
from determine_gene import get_focus_gene
from db_objects import My_db



class Score_Obj():
    # determine which gene to assign
    def __init__(self):
        self.loci_score = {}
        self.primary_dict = defaultdict(set)
        self.read_loci = {}
    
    def add_read(self, read_obj):
        # score = round(read_obj.match_rate * (1 - read_obj.mismatch_rate), 6)
        if read_obj.primary:
            self.primary_dict[read_obj.read_name].add(read_obj.allele_name)
        # else:
        #     return 0
        # if read_obj.match_num < 100:
        #     return 0
        score = read_obj.match_rate
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
            # if "508262af-89b0-42c3-8d64-07cbcf21d8b0" != read_name:
            #     continue
            # print (read_name, self.primary_dict[read_name])
            read_bin = Read_bin(self.loci_score[read_name])
            assigned_locus = read_bin.assign_multiple()
            # print ("\n\n")

            print (read_name, ",".join(assigned_locus), sep = "\t", file = f)
            self.read_loci[read_name] = assigned_locus
        f.close()
        return self.read_loci


class Score_Obj_bk():
    # determine which gene to assign
    def __init__(self):
        self.loci_score = {}
        self.loci_mismatch_score = {}
        self.read_loci = {}
    
    def add_read(self, read_obj):
        # score = round(read_obj.match_rate * (1 - read_obj.mismatch_rate), 6)
        score = read_obj.match_rate
        if read_obj.read_name not in self.loci_score:
            self.loci_score[read_obj.read_name] = {}
            self.loci_score[read_obj.read_name][read_obj.loci_name] = [score, read_obj.match_num]
        elif read_obj.loci_name not in self.loci_score[read_obj.read_name]:
            self.loci_score[read_obj.read_name][read_obj.loci_name] = [score, read_obj.match_num]
        else:
            if score > self.loci_score[read_obj.read_name][read_obj.loci_name][0]:
                self.loci_score[read_obj.read_name][read_obj.loci_name] = [score, read_obj.match_num]
    
    def assign(self, assign_file):
        f = open(assign_file, 'w')
        # print (len(self.loci_score))
        for read_name in self.loci_score: # for each read
            assigned_locus = []
            gene_score = sorted(self.loci_score[read_name].items(), key=lambda item: item[1][0], reverse = True)
            gene_match_len = sorted(self.loci_score[read_name].items(), key=lambda  x: x[1][1], reverse = True)
            # if len(gene_score) > 1 and (gene_score[0][0] == "DQB1"):
            #     print (read_name, gene_score[:2])
            if gene_score[0][1][0] <= Min_score:
                continue
            if len(gene_score) == 1: # mapped to only one gene, directly assign to that gene
                assigned_locus = [gene_score[0][0]]
            else:
                # real-data based adjustment
                

                # if gene_score[0][0] == "DRB1" or gene_score[1][0] == "DRB1":
                #     print (read_name, gene_score[0][0], gene_score[:5], gene_match_len[:5])
                if gene_score[0][0] in ["HLA-U"] and gene_score[1][0] == "HLA-A" :
                    assigned_locus = ["HLA-A"]                
                elif gene_score[0][0] == "HLA-DRB1" and gene_score[0][1][0] - gene_score[1][1][0] < 0.05:  # 0.02 0.05
                    continue
                elif gene_score[0][0] == "HLA-DQB1" and gene_score[0][1][0] < 0.9:
                    continue
                elif gene_score[0][0] == 'HLA-DPB2' and gene_score[1][0] == "HLA-DPA1":
                    assigned_locus = ["HLA-DPA1"]
                elif gene_score[0][0] in ['HLA-DPB1', "HLA-DPA1"] and gene_score[1][0] in ['HLA-DPB1', "HLA-DPA1"]:
                    assigned_locus = ['HLA-DPB1', "HLA-DPA1"]
                # map to more than one gene, check the score difference
                elif gene_score[0][1][0] - gene_score[1][1][0] >= Min_diff:
                    assigned_locus = [gene_score[0][0]]
                # score diff too small, can not determine which gene to assign
                # discard this read
                else:
                    continue
            # print ("assigned locus", assigned_locus)
            print (read_name, assigned_locus, file = f)
            self.read_loci[read_name] = assigned_locus
        f.close()
        return self.read_loci

class Pacbio_Binning():

    def __init__(self):
        self.db = my_db.lite_db

        self.sam = f"""{parameter.outdir}/{parameter.sample}.db.bam"""
        
        if args["m"] != 2:
            self.map2db()

        self.bamfile = pysam.AlignmentFile(self.sam, 'rb')   
        self.assign_file = f"{parameter.outdir}/{parameter.sample}.assign.txt"

    def index_db(self):
        ref_index = self.db[:-5] + args["y"] + ".mmi"
        # print ("search the reference index:", ref_index)
        if not os.path.isfile(ref_index):
            print ("start build Minimap2 index for the reference...")
            os.system(f"minimap2 {minimap_para} -d {ref_index} {self.db} ")
        else:
            print (f"Detect Minimap2 index for the reference: {ref_index}")
        self.db = ref_index

    def map2db(self):
        if args["minimap_index"] == 1:
            self.index_db()
        # map raw reads to database
        alignDB_order = f"""
        fq={parameter.raw_fq}
        ref={self.db}
        outdir={parameter.outdir}
        bin={sys.path[0]}/../bin
        sample={parameter.sample}
        # minimap2 -t {parameter.threads} {minimap_para} -a $ref $fq |samtools view -bS -o {self.sam}
        bwa mem -R '@RG\\tID:foo\\tSM:bar' -t {parameter.threads} {my_db.lite_db} $fq |samtools view -bS -o {self.sam}
        echo alignment done.
        """
        # print (alignDB_order)
        os.system(alignDB_order)

    def read_bam(self):
        # observe each read, assign it to gene based on alignment records
        scor = Score_Obj()
        for read in self.bamfile:
            if read.is_unmapped:
                continue
            # print (read)
            # read_obj = Read_Obj(read)
            read_obj = My_read(read)
            scor.add_read(read_obj)
            # print (read_obj.read_name, read_obj.mismatch_rate, read_obj.allele_name )
        read_loci = scor.assign(self.assign_file)
        for gene in gene_list:
            outfile = parameter.outdir + '/%s.%s.fq'%(gene, args["a"])
            filter_fq(gene, read_loci, parameter.raw_fq, outfile)
        print ("reads-binning done.\n\n\n")

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


class Parameters():

    def __init__(self):

        self.sample = args["n"]
        self.raw_fq = args["r"]
        outdir = args["o"]
        self.population = args["p"]
        self.threads = args["j"]
        self.bin = "%s/../bin/"%(sys.path[0])      
        self.outdir = "%s/%s/"%(outdir, self.sample)
        self.whole_dir = "%s/whole/"%(sys.path[0])

        if not os.path.exists(args["o"]):
            os.system("mkdir %s"%(args["o"]))
        if not os.path.exists(self.outdir):
            os.system("mkdir %s"%(self.outdir))

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
    optional.add_argument("-g", type=int, help="Whether use G group resolution annotation [0|1].", metavar="\b", default=0)
    optional.add_argument("-m", type=int, help="1 represents typing, 0 means only read assignment", metavar="\b", default=1)
    optional.add_argument("-k", type=int, help="The mean depth in a window lower than this value will be masked by N, set 0 to avoid masking", metavar="\b", default=5)
    optional.add_argument("-a", type=str, help="Prefix of filtered fastq file.", metavar="\b", default="long_read")
    optional.add_argument("-y", type=str, help="Read type, [nanopore|pacbio].", metavar="\b", default="pacbio")
    optional.add_argument("--minimap_index", type=int, help="Whether build Minimap2 index for the reference [0|1]. Using index can reduce memory usage.", metavar="\b", default=1)
    optional.add_argument("--db", type=str, help="db dir.", metavar="\b", default=sys.path[0] + "/../db/")
    optional.add_argument("--strand_bias_pvalue_cutoff", type=float, help="Remove a variant if the allele observations are biased toward one strand (forward or reverse). Recommand setting 0 to high-depth data.", metavar="\b", default=0.01)
    # optional.add_argument("-u", type=str, help="Choose full-length or exon typing. 0 indicates full-length, 1 means exon.", metavar="\b", default="0")
    optional.add_argument("--seed", type=int, help="seed to generate random numbers", metavar="\b", default=8)
    optional.add_argument("--max_depth", type=int, help="maximum depth for each HLA locus. Downsample if exceed this value.", metavar="\b", default=2000)
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    parameter = Parameters()
    # Min_score = 0.1  #the read is too long, so the score can be very low.
    Min_score = 0  #the read is too long, so the score can be very low.
    Min_diff = args["d"]  #0.001

    gene_list, interval_dict =  get_focus_gene(args)
    my_db = My_db(args)

    minimap_para = ''
    if args["y"] == "pacbio":
        minimap_para = " -x map-pb "
    elif args["y"] == "nanopore":
        minimap_para = " -x map-ont "


    ###assign reads
    if args["m"] == 10086:
        print ("skip assignment, just for testing")
    else:
        pbin = Pacbio_Binning()
        pbin.read_bam()        

    print ("Finished.")




