### count depth of CYP2D6 of each sample

import pandas as pd
import os
import sys
import numpy as np

cyp2d6_len = 4367  ## length in hg38


def count_depth(depth_file):
    ## read samtools depth file using pandas
    df = pd.read_csv(depth_file, sep="\t", header=None)
    ## get total depth
    total_depth = df[2].sum()
    ## get average depth
    avg_depth = total_depth/cyp2d6_len
    return round(avg_depth,2)

def for_all_samples_dp(folder):
    ## read all samples
    # folder = "/home/shuaiw/methylation/data/hla/CYP_rerun/out_CYP/"
    sample_list = []
    sample_depth_dict = {}
    for sample in os.listdir(folder):
        sample_list.append(sample)
        depth_file = f"{folder}/{sample}/{sample}/Genes_step2/CYP2D6.depth"
        dp1 = f"{folder}/{sample}/{sample}/Genes_step2/CYP2D6.0.depth"
        dp2 = f"{folder}/{sample}/{sample}/Genes_step2/CYP2D6.1.depth"
        if os.path.exists(depth_file):
            avg_depth = count_depth(depth_file)
            sample_depth_dict[sample] = avg_depth
            print (sample, avg_depth)
        elif os.path.exists(dp1) and os.path.exists(dp2):
            avg_depth1 = count_depth(dp1)
            avg_depth2 = count_depth(dp2)
            avg_depth = avg_depth1 + avg_depth2
            avg_depth = round(avg_depth, 2)
            print (sample, avg_depth)
            sample_depth_dict[sample] = avg_depth
    return sample_depth_dict


# for_all_samples_dp()
# depth_file = "/home/shuaiw/methylation/data/hla/CYP_rerun/out_CYP/HG01928/HG01928/Genes_step2/CYP2D6.depth"
# print (depth_file)
# avg_depth = count_depth(depth_file)
# print (avg_depth)