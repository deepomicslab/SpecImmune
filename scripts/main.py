
import sys
import os
import pysam
import gzip
import argparse
from db_objects import My_db


def main(args):
    if not os.path.exists(args["o"]):
        os.system("mkdir %s"%(args["o"]))
    outdir = args["o"] + "/" + args["n"]
    if not os.path.exists(outdir):
        os.system("mkdir %s"%(outdir))

    if args['i'] != "IG_TR":
        command = f"""
        ## first: read binning
        python3 {sys.path[0]}/long_read_typing.py -r {args["r"]} -n {args["n"]} -i {args["i"]} -o {args["o"]} -j {args["j"]} -k {args["k"]} -y {args["y"]} -m 0 --db {args["db"]}

        ## second: find a pair of alleles for each HLA locus
        # python3 {sys.path[0]}/select_best_reference_alleleV2.py -r {args["r"]} -n {args["n"]}  -i {args["i"]} -o {args["o"]} -j {args["j"]} -y {args["y"]} --db {args["db"]}
        python3 {sys.path[0]}/select_best_reference_alleleV2.py -r {args["r"]} -n {args["n"]}  -i {args["i"]} -o {args["o"]} -j {args["j"]} -y {args["y"]} --db {args["db"]}

        """
        os.system(command)

        # build individual ref when first run
        
        # my_db = My_db(args)
        # command = f"""
        # ## third: build individual reference for each HLA locus
        # python3 {sys.path[0]}/get_ref.py -n {args["n"]} -o {args["o"]} -j {args["j"]}    
        # python3 {sys.path[0]}/build_ref.py {args["o"]}/{args["n"]}/{args["n"]}.map.txt {my_db.full_db} 
        # """
        # if args["first_run"]:
        #     os.system(command)
        
        
        # command = f"""
        # ## forth: call & phasing variant & typing
        # python3 {sys.path[0]}/long_read_typing.py -r {args["r"]} -n {args["n"]} -o {args["o"]} -j {args["j"]} -k {args["k"]} -y {args["y"]} -m 2 --db {args["db"]} -i {args["i"]}

        # """
        # os.system(command)
    
    else:
        command = f"""
        mkdir {args["o"]}/{args["n"]}/tmp
        bash {sys.path[0]}/run.phase.IG.TR.sh {args["n"]} {args["r"]} {args["o"]}/{args["n"]} {args["db"]} {args["j"]} {args["k"]}
        """
        os.system(command)

        
    





if __name__ == "__main__":   

    parser = argparse.ArgumentParser(description="HLA Typing with only long-read data.", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    required.add_argument("-r", type=str, help="Long-read fastq file. PacBio or Nanopore.", metavar="\b")
    required.add_argument("-n", type=str, help="Sample ID", metavar="\b")
    required.add_argument("-o", type=str, help="The output folder to store the typing results.", metavar="\b", default="./output")
    required.add_argument("-i", type=str, help="HLA,KIR,CYP,IG_TR",metavar="\b", default="HLA")
    # optional.add_argument("-p", type=str, help="The population of the sample [Asian, Black, Caucasian, Unknown, nonuse] for annotation. Unknown means use mean allele frequency in all populations. nonuse indicates only adopting mapping score and considering zero-frequency alleles.", metavar="\b", default="Unknown")
    optional.add_argument("-j", type=int, help="Number of threads.", metavar="\b", default=5)
    # optional.add_argument("-d", type=float, help="Minimum score difference to assign a read to a gene.", metavar="\b", default=0.001)
    # optional.add_argument("-g", type=int, help="Whether use G group resolution annotation [0|1].", metavar="\b", default=0)
    # optional.add_argument("-m", type=int, help="1 represents typing, 0 means only read assignment", metavar="\b", default=1)
    optional.add_argument("-k", type=int, help="The mean depth in a window lower than this value will be masked by N, set 0 to avoid masking", metavar="\b", default=5)
    # optional.add_argument("-a", type=str, help="Prefix of filtered fastq file.", metavar="\b", default="long_read")
    optional.add_argument("-y", type=str, help="Read type, [nanopore|pacbio].", metavar="\b", default="pacbio")
    optional.add_argument("--db", type=str, help="db dir.", metavar="\b", default=sys.path[0] + "/../db/")
    # optional.add_argument("-dr", "--db_ref", type=str, help="database reference", metavar="\b", \
    #                        default=sys.path[0] + "/../db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta")
    optional.add_argument("-f", "--first_run", type=bool, help="set False for rerun", metavar="\b", default=True)
    # optional.add_argument("-u", type=str, help="Choose full-length or exon typing. 0 indicates full-length, 1 means exon.", metavar="\b", default="0")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)
    
    main(args)