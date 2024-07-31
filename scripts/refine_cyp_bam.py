"""
in cyp, a read can be mapped to 2d6 and 2d7, for such reads, check which gene the read should map, 
output in a new bam
"""

import pysam
import os
import sys

from read_objects import My_read

def refine_cyp_bam(bam_file, new_bam):
    bam = pysam.AlignmentFile(bam_file, 'rb')
    
    read_record = {}
    tolerance = 3000
    for read in bam:

        my_read = My_read()
        my_read.load_bam(read)
        if my_read.identity < min_identity:
            continue

        if read.query_name not in read_record:
            read_record[read.query_name] = [read]
        else:
            read_record[read.query_name].append(read)
    new_record = {}
    for read_name in read_record:
        # if read_name == 'SRR15476211.2802':
        #     print (read_name, len(read_record[read_name]))
        if len(read_record[read_name]) == 2:
            record_1 = read_record[read_name][0]
            record_2 = read_record[read_name][1]
            # if read_name == 'SRR15476211.2802':
            #     print (record_1.reference_start, record_2.reference_start, record_1.is_secondary, record_2.is_secondary, record_1.mapping_quality, record_2.mapping_quality)
            ## if one record is primary and the other is 0 mapping quality
            if not record_1.is_secondary and record_2.mapping_quality == 0:
                pass
            elif not record_2.is_secondary and record_1.mapping_quality == 0:
                pass
            else:
                continue

            # if read_name == 'SRR15476211.2802':
            #     print (read_name, record_1.reference_start, record_2.reference_start)
            ### if the two reads donot have clipped reads, and they are mapped to 2d6 and 2d7,
            ### check which gene the read should map
            if record_1.reference_name == 'chr22' and record_2.reference_name == 'chr22':

                # if read_name == 'SRR15476211.2802':
                #     print (read_name, record_1.reference_start, record_2.reference_start)
                if 42126499-tolerance <= record_1.reference_start <= 42130865+tolerance and 42139576-tolerance <= record_2.reference_start <= 42149500+tolerance:
                    my_read1 = My_read()
                    my_read1.load_bam(record_1)
                    my_read2 = My_read()
                    my_read2.load_bam(record_2)
                    # if read_name == 'SRR15476211.2802':
                    #     print (my_read1.reference_start,my_read1.identity, my_read2.reference_start, my_read2.identity)
                    # if my_read1.identity > my_read2.identity:
                    #     read_record[read_name] = [record_1]
                    # else:
                    #     read_record[read_name] = [record_2]
                    # if my_read1.reference_start < my_read2.reference_start:
                    if my_read1.identity > my_read2.identity:
                        # if read_name == 'SRR15476211.2802':
                        #     print ('1')
                        read_record[read_name] = [record_1]
                    else:
                        # if read_name == 'SRR15476211.2802':
                        #     print ('2')
                        read_record[read_name] = [record_2]
                    ## set the record as primary
                    read_record[read_name][0].is_secondary = False

        else:
            new_record[read_name] = read_record[read_name]

    
    new_bam = pysam.AlignmentFile(new_bam, 'wb', template=bam)
    for read in new_record:
        for i in range(len(new_record[read])):
            new_bam.write(new_record[read][i])
    new_bam.close()
    bam.close()
    ## sort and index the new bam using samtools
    # cmd = f"""samtools sort -o {new_bam}.sorted.bam {new_bam}\nsamtools index {new_bam}.sorted.bam"""
    # os.system(cmd)


# bam_file = '/mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/amplicon/NA17300/NA17300.bam'
# new_bam = '/mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/amplicon/NA17300/NA17300_refine.bam'
bam_file = sys.argv[1]
new_bam = sys.argv[2]
min_identity = 0.85
print ("refine bam...")
refine_cyp_bam(bam_file, new_bam)