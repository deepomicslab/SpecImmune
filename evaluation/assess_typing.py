
import sys
import argparse

sys.path.insert(0, sys.path[0]+'/../scripts/')

from four_field_compare import assess_sim_module
from determine_gene import get_focus_gene


# python3 ../evaluation/assess_typing.py -i HLA --true $outdir/$sample/$sample.HLA.hap.alleles.txt --infer $outdir/$sample/${sample}.HLA.type.result.txt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="HLA Typing with long-read data.", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    required.add_argument("-i", type=str, help="HLA,KIR,CYP",metavar="\b", default="HLA")
    required.add_argument("--true", type=str, help="true", metavar="\b")
    required.add_argument("--infer", type=str, help="infer.", metavar="\b")

    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    gene_list, interval_dict =  get_focus_gene(args)

    assess_sim_module(args["true"], args["infer"], gene_list)


