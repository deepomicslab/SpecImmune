
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
        python3 {sys.path[0]}/read_binning.py -r {args["r"]} -n {args["n"]} -i {args["i"]} -o {args["o"]} -j {args["j"]} -k {args["k"]} -y {args["y"]} \
            --db {args["db"]} --min_identity {args["min_identity"]} --seq_tech {args["seq_tech"]} --RNA_type {args["RNA_type"]}
        ## second: find a pair of alleles for each HLA locus
        """
        if args["mode"] >= 4:
            os.system(command)

        command = f"""
        ## second: find a pair of alleles for each HLA locus
        python3 {sys.path[0]}/select_best_reference_alleleV2.py --max_read_num {args["max_read_num"]} --candidate_allele_num {args["candidate_allele_num"]} \
            --hete_p {args["hete_p"]} --align_method minimap2 -r {args["r"]} -n {args["n"]}  -i {args["i"]} -o {args["o"]} -j {args["j"]} -y {args["y"]} \
            --db {args["db"]} --seq_tech {args["seq_tech"]} --RNA_type {args["RNA_type"]}
        """
        if args["mode"] >= 3:
            os.system(command)

        # return

        # build individual ref when first run
        my_db = My_db(args)
        # print(f"""python3 {sys.path[0]}/get_ref.py -n {args["n"]} -o {args["o"]} -j {args["j"]}""")
        # command = f"""
        # ## third: build individual reference for each HLA locus
        # python3 {sys.path[0]}/get_ref.py -n {args["n"]} -o {args["o"]} -j {args["j"]}    
        # python3 {sys.path[0]}/build_ref.py {args["o"]}/{args["n"]}/{args["n"]}.map.txt {my_db.full_db} {my_db.individual_ref_dir}
        # """
        command = f"""
        ## third: build individual reference for each HLA locus, two ref version
        python3 {sys.path[0]}/get_2ref_align.py {args["n"]} {my_db.full_db} {my_db.individual_ref_dir} {args["o"]} {args["y"]} {args["j"]} {args["i"]} \
            {args["seq_tech"]} {args["RNA_type"]}
        """
        # if args["first_run"]:
        if args["mode"] >= 2:
            print (f"<<<<get_2ref_align.py\n {command}", flush=True)
            os.system(command)


        
        if args["analyze_method"] == "phase":
            command = f"""
            ## forth: call & phasing variant & typing
            python3 {sys.path[0]}/long_read_typing.py -r {args["r"]} -n {args["n"]} -o {args["o"]} -j {args["j"]} -k {args["k"]} -y {args["y"]} --db {args["db"]} \
                -i {args["i"]} --seq_tech {args["seq_tech"]} --RNA_type {args["RNA_type"]}

            """
            if args["mode"] >= 1:
                os.system(command)

        elif args["analyze_method"] == "assembly":
            command = f"""
            ## forth: assembly
            python3 {sys.path[0]}/assembly.py -i {args["i"]} -o {args["o"]} -n {args["n"]} -j {args["j"]} -y {args["y"]}
            """
            if args["mode"] >= 1:
                os.system(command)
        else:
            print("Please choose phase or assembly as analyze method.", flush=True)
            return
    
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
    optional.add_argument("-j", type=int, help="Number of threads.", metavar="\b", default=5)
    # optional.add_argument("-g", type=int, help="Whether use G group resolution annotation [0|1].", metavar="\b", default=0)
    optional.add_argument("--mode", type=int, help="4 represents all steps, 3 skip first, 2 skip two, 3, skipt three", metavar="\b", default=4)
    optional.add_argument("--analyze_method", type=str, help="phase/assembly", metavar="\b", default="phase")
    optional.add_argument("-k", type=int, help="The mean depth in a window lower than this value will be masked by N, set 0 to avoid masking", metavar="\b", default=5)
    # optional.add_argument("-a", type=str, help="Prefix of filtered fastq file.", metavar="\b", default="long_read")
    optional.add_argument("-y", type=str, help="Read type, [nanopore|pacbio|pacbio-hifi].", metavar="\b", default="pacbio")
    optional.add_argument("--db", type=str, help="db dir.", metavar="\b", default=sys.path[0] + "/../db/")
    optional.add_argument("-f", "--first_run", type=bool, help="set False for rerun", metavar="\b", default=True)
    optional.add_argument("--min_identity", type=float, help="Minimum identity to assign a read.", metavar="\b", default=0.85)
    optional.add_argument("--hete_p", type=float, help="Hete pvalue.", metavar="\b", default=0.3) 
    optional.add_argument("--candidate_allele_num", type=int, help="Maintain this number of alleles for ILP step.", metavar="\b", default=200)
    optional.add_argument("--min_read_num", type=int, help="min support read number for each locus.", metavar="\b", default=2)
    optional.add_argument("--max_read_num", type=int, help="max support read number for each locus.", metavar="\b", default=500)
    optional.add_argument("-rt", "--RNA_type", type=str, help="traditional,2D,Direct,SIRV",metavar="\b", default="traditional")
    optional.add_argument("--seq_tech", type=str, help="Amplicon sequencing or WGS sequencing [wgs|amplicon|rna].", metavar="\b", default="wgs")

    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)
    
    main(args)