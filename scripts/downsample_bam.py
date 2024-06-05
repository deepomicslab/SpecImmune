import pysam
import os
import numpy as np
import random
import sys

from get_allele_depth import Get_depth


    

def downsample_bam(input_bam, output_bam, downsample_ratio, seed):
    random.seed(seed)
    # Open the input BAM file for reading
    input_samfile = pysam.AlignmentFile(input_bam, "rb")

    # Create a new BAM file for writing the downsampled records
    output_samfile = pysam.AlignmentFile(output_bam, "wb", header=input_samfile.header)

    # Iterate over alignment records and downsample based on the ratio
    for record in input_samfile:
        if downsample_ratio >= 1.0 or random.random() <= downsample_ratio:
            # Write the record to the output BAM file
            output_samfile.write(record)

    # Close the BAM files
    input_samfile.close()
    output_samfile.close()

def handle_bam(output_bam, output_depth):
    cmd = f"""
    samtools index {output_bam}
    samtools depth -d 1000000 -aa {output_bam} > {output_depth}
    """
    os.system(cmd)


def downsample_func(input_bam, output_bam, input_depth, output_depth, max_depth, seed):

    get_depth = Get_depth(input_depth)
    get_depth.record_depth()

    depth_list = []
    for gene in get_depth.depth_dict:
        for allele in get_depth.depth_dict[gene]:
            depth_list = get_depth.depth_dict[gene][allele]

    mean_depth = np.mean(depth_list)
    # print ("mean_depth", mean_depth, len(depth_list), depth_list, input_depth)

    if mean_depth <= max_depth:
        return 1
    else:
        downsample_ratio = float(max_depth)/mean_depth
        downsample_bam(input_bam, output_bam, downsample_ratio, seed)
        handle_bam(output_bam, output_depth)
        return downsample_ratio


# # Example usage:
# input_bam = "input.bam"  # Path to the input BAM file
# output_bam = "output.bam"  # Path for the output BAM file
# downsample_ratio = 0.5  # Downsample ratio (e.g., 0.5 for 50% downsample)

# downsample_bam(input_bam, output_bam, downsample_ratio)