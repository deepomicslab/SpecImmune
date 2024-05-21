"""
Choose best-mapped allele as reference

wangshuai, wshuai294@gmail.com
"""


import os
import numpy as np
import pickle
import sys
import argparse
import random

# import pulp
# from pulp import LpProblem, LpMinimize, LpMaximize, LpVariable, lpSum, PULP_CBC_CMD, value
from collections import defaultdict
import pysam

gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
# gene_list = ['A', 'C', 'DRB1']

def run_depth():
    cmd = f"""
    samtools view -bS -F 0x800  fredhutch-hla-1408-1012.db.sam | samtools sort - >fredhutch-hla-1408-1012.db.bam
    samtools depth -aa fredhutch-hla-1408-1012.db.bam>fredhutch-hla-1408-1012.db.depth
    """
    os.system(cmd)

def mapping_p(MAPQ, base = 10):
    ## MAPQ: MAPping Quality. It equals −10 log10 Pr{mapping position is wrong}, rounded to the nearest integer. A value 255 indicates that the mapping quality is not available.
    ## if MAPQ=0, base=1, 0.25. base=2,0.36  base=3, 0.49
    return 1 - 10**(-(base + MAPQ)/10)   

class Get_depth():

    def __init__(self, depth_file):
        self.depth_file = depth_file
        self.depth_dict = {}

    def record_depth(self):
        f = open(self.depth_file)
        for line in f:
            array = line.strip().split()

            allele = array[0]

            gene = allele.split("*")[0]

            depth = int(array[2])
            if gene not in self.depth_dict:
                self.depth_dict[gene] = {}
            if allele not in self.depth_dict[gene]:
                self.depth_dict[gene][allele] = []

            self.depth_dict[gene][allele].append(depth)
    
    def select(self, sort_depth_file):
        f = open(sort_depth_file, 'w')
        print ("Gene\tAllele\tDepth\tAllele_length", file = f)

        record_candidate_alleles = defaultdict(set)
        record_allele_length = {}
        for gene in self.depth_dict:
            if gene not in gene_list:
                continue
            record_allele_depth = {}
            
            record_allele_info = {}
            for allele in self.depth_dict[gene]:
                mean_depth = np.mean(self.depth_dict[gene][allele])
                median_depth = np.median(self.depth_dict[gene][allele])

                over_0 = 0
                for e in self.depth_dict[gene][allele]:
                    if e > 0:
                       over_0 += 1
                coverage =  float(over_0)/len(self.depth_dict[gene][allele])
                record_allele_length[allele] = len(self.depth_dict[gene][allele])

                record_allele_depth[allele] = mean_depth
                record_allele_info[allele] = [mean_depth, median_depth, coverage]

                # print (allele, mean_depth, median_depth, coverage)
            sorted_dict = sorted(record_allele_depth.items(), key=lambda x: x[1], reverse=True)
            for i in range(min([args["m"], len(sorted_dict)])):
                print (sorted_dict[i][0], round(sorted_dict[i][1],2), record_allele_length[sorted_dict[i][0]] , sep = "\t")
                record_candidate_alleles[gene].add(sorted_dict[i][0])

            ## output
            for i in range(len(sorted_dict)):
                print (gene, sorted_dict[i][0], round(sorted_dict[i][1],2), record_allele_length[sorted_dict[i][0]], sep = "\t", file = f)
        f.close()

        return record_candidate_alleles, record_allele_length


def find_primary_alignments(sam):

    primary_match_len_dict = {}
    bamfile = pysam.AlignmentFile(sam, 'r')  
    for read in bamfile:
        if read.is_unmapped:
            continue

        if read.is_secondary:
            continue
        read_name = read.query_name
        match_num, mismatch, alignment_score, identity = get_match_length(read)
        primary_match_len_dict[read_name] = [match_num, mismatch, alignment_score]
    
    bamfile.close()
    return primary_match_len_dict

def find_bad_reads(sam):

    bad_reads_dict = {} 
    bamfile = pysam.AlignmentFile(sam, 'r')  
    for read in bamfile:
        # if read.is_unmapped:
        #     continue

        if read.is_secondary:
            continue

        if read.is_unmapped:
            bad_reads_dict[read_name] = 1
            continue

        read_name = read.query_name
        match_num, mismatch, alignment_score, identity = get_match_length(read)

        # if identity < 0.92:
        #     bad_reads_dict[read_name] = 1

        # if match_num < 5000:
        #     bad_reads_dict[read_name] = 1
    
    bamfile.close()
    return bad_reads_dict

def construct_matrix(args, gene, sam, record_candidate_alleles, record_allele_length, mapping_quality_cutoff = 0 ):
    # sam = "/mnt/d/HLAPro_backup/Nanopore_optimize/output/fredhutch-hla-1408-1012/fredhutch-hla-1408-1012.db.bam"
    # primary_match_len_dict = find_primary_alignments(sam)
    # bad_reads_dict = find_bad_reads(sam)
    bad_reads_dict = {}
    bamfile = pysam.AlignmentFile(sam, 'r')  

    record_read_index = defaultdict(dict)
    record_allele_index = defaultdict(dict)

    record_read_index_base = defaultdict(dict)
    record_allele_index_base = defaultdict(dict)

    for read in bamfile:
        if read.is_unmapped:
            continue
        # print (read.mapping_quality)
        # if read.mapping_quality < mapping_quality_cutoff:
        #     continue

        # if read.is_secondary:
        #     continue
        # if random.randint(0, 10) < 2:
        #     continue

        read_name = read.query_name
        allele_name = read.reference_name
        gene = allele_name.split("*")[0]
        # print (read_name)

        if read_name in bad_reads_dict:
            continue

        if gene not in record_read_index:
            record_read_index_base[gene] = set()
            record_allele_index_base[gene] = set()
        
        if allele_name not in record_candidate_alleles[gene]:
            continue 

        if read_name not in record_read_index_base[gene]:
            record_read_index[gene][read_name] = len(record_read_index_base[gene])
            record_read_index_base[gene].add(read_name)
        if allele_name not in record_allele_index_base[gene]:
            record_allele_index[gene][allele_name] = len(record_allele_index_base[gene])
            record_allele_index_base[gene].add(allele_name)


    bamfile.close()
    # print (record_allele_index)

    read_matrix_dict = {}
    read_identity_matrix_dict = {}
    read_map_len_matrix_dict = {}
    read_mismatch_matrix_dict = {}

    for gene in record_read_index:
        x = np.empty(shape=(len(record_allele_index[gene]),len(record_read_index[gene])))
        x.fill(0)
        y = np.empty(shape=(len(record_allele_index[gene]),len(record_read_index[gene])))
        y.fill(0)
        z = np.empty(shape=(len(record_allele_index[gene]),len(record_read_index[gene])))
        z.fill(0)
        e = np.empty(shape=(len(record_allele_index[gene]),len(record_read_index[gene])))
        e.fill(0)

        read_matrix_dict[gene] = x
        read_map_len_matrix_dict[gene] = y
        read_identity_matrix_dict[gene] = z
        read_mismatch_matrix_dict[gene] = e

    bamfile = pysam.AlignmentFile(sam, 'r')  
    for read in bamfile:
        if read.is_unmapped:
            continue
        # if read.mapping_quality < mapping_quality_cutoff:
        #     continue


        # if read.mapping_quality < 1:
        #     delta = 0.1
        # else:
        #     delta = 1

        read_name = read.query_name
        allele_name = read.reference_name
        gene = allele_name.split("*")[0]
        match_num, mismatch, alignment_score, identity = get_match_length(read)

        # if  primary_match_len_dict[read_name][0] - match_num > 20 :   # 10
        #     continue
        # print ("alignment_score", alignment_score)

        # if gene == "C":
        #     if  primary_match_len_dict[read_name][2] - alignment_score > 50:   # 10
        #         continue

        # match_num = float(match_num)/record_allele_length[allele_name]
        # print (match_num)

        # if match_num < 2000:
        #     continue

        # if read_name == "9265b58c-82d8-4c1b-9c03-4c8a1bc05a15":
        #     print (read_name, allele_name, match_num)

        if read_name in record_read_index[gene] and allele_name in record_allele_index[gene]:
            a = record_allele_index[gene][allele_name]
            r = record_read_index[gene][read_name]
            read_matrix_dict[gene][a][r] = 1

            ## select the longest
            if read_map_len_matrix_dict[gene][a][r] < match_num:  # sometimes a read can have two alignments on a same allele, just choose the longest alignment length
                # 3af99479-6bb5-434c-89a6-45405b3716fd in fredhutch-hla-KT17.C.db.sam
                read_map_len_matrix_dict[gene][a][r] = match_num  #   #  round(match_num * mapping_p(read.mapping_quality))
                read_identity_matrix_dict[gene][a][r] = identity
                read_mismatch_matrix_dict[gene][a][r] = mismatch
            
            ## select the highest identity
            # if read_map_len_matrix_dict[gene][a][r] != 0:
            #     if read_map_len_matrix_dict[gene][a][r]/(read_map_len_matrix_dict[gene][a][r]+read_identity_matrix_dict[gene][a][r]) < match_num/(match_num + identity):
            #         read_map_len_matrix_dict[gene][a][r] = match_num  
            #         read_identity_matrix_dict[gene][a][r] = identity
            # else:
            #     read_map_len_matrix_dict[gene][a][r] = match_num  
            #     read_identity_matrix_dict[gene][a][r] = identity





    bamfile.close()
    # print (read_matrix_dict)
    print_read_matrix(args, gene, record_allele_index, record_read_index, read_map_len_matrix_dict, read_identity_matrix_dict)
    return read_matrix_dict, read_map_len_matrix_dict, record_allele_index, record_read_index, read_identity_matrix_dict,read_mismatch_matrix_dict

def print_read_matrix(args, gene, record_allele_index, record_read_index, read_map_len_matrix_dict,read_identity_matrix_dict):
    outdir = args["o"] + "/" + args["n"]
    out = open(f"""{outdir}/{args["n"]}.{gene}.read.matrix.csv""", 'w')
    first_line = "allele,"
    for allele in record_allele_index[gene]:
        first_line += allele + ","
    # first_line += "\n"
    print (first_line, file = out)
    for read_name in record_read_index[gene]:
        line = read_name + ","
        for allele_name in record_allele_index[gene]:
            a = record_allele_index[gene][allele_name]
            r = record_read_index[gene][read_name]
            line += str(read_map_len_matrix_dict[gene][a][r]) + "/" + str(round(read_identity_matrix_dict[gene][a][r], 4)) + ","
        # line += "\n"
        print (line, file = out)
    out.close()

def get_match_length(read, identity_cutoff = 0.85):
    match_num = 0        
    indel_num = 0
    unmap_part_num = 0
    for ci in read.cigar:
        if ci[0] == 0:
        # if ci[0] == 0 or ci[0] == 1 or ci[0] == 2:
            match_num += ci[1]
        elif ci[0] == 1 or ci[0] == 2:
            indel_num += ci[1]
        else:
            unmap_part_num += 1

    mis_NM = 0
    alignment_score = 0
    for ta in read.get_tags():
        if ta[0] == 'NM':
            mis_NM = ta[1]  
        if ta[0] == 'AS':
            alignment_score = ta[1]  
    
    if match_num == 0:
        print (read.query_name, read.cigar)
        sys.exit(0)
    
    identity = (match_num - mis_NM)/match_num

    # if identity < identity_cutoff:
    #     return 0 
    # if unmap_part_num > 2:
    #     print (read.cigar)
    #     return 0, 0

    return match_num - mis_NM,  mis_NM, alignment_score, identity
    # return alignment_score,  mis_NM, alignment_score 
    # return match_num, identity

def save_pkl(pkl_name, my_dict):
    with open(pkl_name, 'wb') as handle:
        pickle.dump(my_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

def load_pkl(pkl_name):
    with open(pkl_name, 'rb') as handle:
        my_dict = pickle.load(handle)
        return my_dict


def main(args):

    if not os.path.exists(args["o"]):
        os.system("mkdir %s"%(args["o"]))
    outdir = args["o"] + "/" + args["n"]
    if not os.path.exists(outdir):
        os.system("mkdir %s"%(outdir))
    # outdir = args["o"]
    
    
    result_dict = {}
    allele_match_dict = {}
    

    # for gene in read_matrix_dict:
    for gene in gene_list:
        # if gene not in gene_list:
        #     continue

        sam, depth_file, sort_depth_file = map2db(args, gene)

        get_depth = Get_depth(depth_file)
        get_depth.record_depth()
        record_candidate_alleles, record_allele_length = get_depth.select(sort_depth_file)

        print (sam)

        read_matrix_dict, read_map_len_matrix_dict, record_allele_index, record_read_index, read_identity_matrix_dict,read_mismatch_matrix_dict = \
            construct_matrix(args, gene, sam, record_candidate_alleles, record_allele_length)

        # read_matrix_dict_file = outdir + "/" + args["n"] + "." + gene + ".read_matrix_dict.pkl"
        # read_map_len_matrix_dict_file = outdir + "/" + args["n"] + "." + gene + ".read_map_len_matrix_dict.pkl"
        # record_allele_index_file = outdir + "/" + args["n"] + "." + gene + ".record_allele_index.pkl"
        # read_identity_matrix_file = outdir + "/" + args["n"] + "." + gene + ".read_identity_matrix.pkl"

        # save_pkl(read_matrix_dict_file, read_matrix_dict)
        # save_pkl(read_map_len_matrix_dict_file, read_map_len_matrix_dict)
        # save_pkl(record_allele_index_file, record_allele_index)
        # save_pkl(read_identity_matrix_file, read_identity_matrix_dict)

        # read_matrix_dict = load_pkl(read_matrix_dict_file)
        # read_map_len_matrix_dict = load_pkl(read_map_len_matrix_dict_file)
        # record_allele_index = load_pkl(record_allele_index_file)
        # read_identity_matrix_dict = load_pkl(read_identity_matrix_file)

        # print (read_matrix_dict_file)
        print ("finish matrix construction")


        record_index_allele = {}  #{index: allele}
        for allele in record_allele_index[gene]:
            record_index_allele[record_allele_index[gene][allele]] = allele
 
        read_matrix = read_matrix_dict[gene]
        read_map_len_matrix = read_map_len_matrix_dict[gene]
        read_identity_matrix = read_identity_matrix_dict[gene]
        read_mismatch_matrix = read_mismatch_matrix_dict[gene]
        
        # print (read_matrix)
        # print (record_allele_index[gene])
        for i in range(len(read_matrix)):
            allele = list(record_allele_index[gene].keys())[i]
            match_len = np.sum(read_map_len_matrix[i])
            allele_match_dict[allele] = match_len
            print (i, list(record_allele_index[gene].keys())[i], np.sum(read_matrix[i]), np.sum(read_map_len_matrix[i]), \
                round(np.sum(read_map_len_matrix[i])/record_allele_length[allele]), sep = "\t")

        # type_result, objective_value = model(read_matrix, read_map_len_matrix, allele_match_dict)
        type_result, objective_value, type_allele_result = model2(read_matrix, read_map_len_matrix, allele_match_dict, \
            record_index_allele, record_read_index, gene, read_identity_matrix, read_mismatch_matrix,record_allele_length)
        # print (gene, len(type_result), len(record_allele_index[gene]), record_allele_index[gene])
        # type_allele_result = []
        # for allele in record_allele_index[gene]:
        #     # print (allele, record_allele_index[gene][allele])
        #     if type_result[record_allele_index[gene][allele]] == 1:
        #         type_allele_result.append(allele)
        print (gene, type_allele_result, "\n\n")
        
        # homo_hete_ratio = (2 * allele_match_dict[type_allele_result[0]] - objective_value )/ allele_match_dict[type_allele_result[0]]
        # print (2 * allele_match_dict[type_allele_result[0]], objective_value, homo_hete_ratio)
        # if homo_hete_ratio  >  args["b"]:
        #     print (2 * allele_match_dict[type_allele_result[0]], objective_value, homo_hete_ratio)
        #     type_allele_result = [type_allele_result[0]]

        homo_hete_ratio = (objective_value - allele_match_dict[type_allele_result[0].split(",")[0]])/ allele_match_dict[type_allele_result[0].split(",")[0]]
        
        homo_hete_ratio_cutoff = args["b"]

        if gene == "DPA1":
            homo_hete_ratio_cutoff = 0.001
        # if gene == "DQA1":
        #     homo_hete_ratio_cutoff = 0.01

        print (homo_hete_ratio, homo_hete_ratio_cutoff)

        if gene == "C":
            if cal_sim_of_alleles(type_allele_result[0].split(",")[0], type_allele_result[1].split(",")[0]) != 6:
                if homo_hete_ratio  <  homo_hete_ratio_cutoff:
                    type_allele_result = [type_allele_result[0]]
        
        else:
            if homo_hete_ratio  <  homo_hete_ratio_cutoff:
                type_allele_result = [type_allele_result[0]]         


        result_dict[gene] = type_allele_result
        print (gene, type_allele_result, "\n\n")
    
    out = open(f"{outdir}/hla.new.result.txt", 'w')
    print ("#", file = out)
    print ("sample", end = "\t", file = out)
    for gene in gene_list:
        for i in range(2):
            print (f"HLA_{gene}_{i+1}", end = "\t", file = out)

    print (f"\n{args['n']}", end = "\t", file = out)
    for gene in gene_list:
        if len(result_dict[gene]) == 1:
            result_dict[gene].append(result_dict[gene][0])
 
        for i in range(2):
            print (result_dict[gene][i], end = "\t", file = out)

    # print (f"\nMatch_length", end = "\t", file = out)
    # for gene in gene_list:
    #     if len(result_dict[gene]) == 1:
    #         result_dict[gene].append(result_dict[gene][0])
    #     for i in range(2):
    #         print (allele_match_dict[result_dict[gene][i]], end = "\t", file = out)

    # print (f"\nObjective value:\t", end = "\t", file = out)
    out.close()


def model_bk(C, M):
    ## C: traditional read hit matrix   M: read map len matrix
    # allele_num = 3
    # read_num = 4
    allele_num = len(C)
    read_num = len(C[0])
    t_max = 2
    t_min = 1

    # C = [[1,1,1,1], [0,0,0,0], [1,1,1,1]]
    prob = LpProblem('HLA typing', LpMaximize)
    my_x =LpVariable.dicts("X", range(allele_num), cat = "Binary", lowBound = 0, upBound = 1)
    my_y =LpVariable.dicts("Y", range(read_num), cat = "Binary", lowBound = 0, upBound = 1)
    my_z =LpVariable.dicts("Z", range(allele_num*read_num), cat = "Binary", lowBound = 0, upBound = 1)  ## get x*y
    my_A =LpVariable.dicts("A", range(allele_num*read_num), cat = "Binary", lowBound = 0, upBound = 1)  ## get x*y
    my_L =LpVariable.dicts("L", range(read_num), lowBound = 0, cat='Integer')  # total match len from all reads

    # print (my_z)
    sum_x = 0
    for a in range(allele_num):
        sum_x += my_x[a]

    prob += sum_x <= t_max
    prob += sum_x >= t_min

    for r in range(read_num):
        sum_xc = 0
        sum_xm = 0
        for a in range(allele_num):
            sum_xc += my_x[a] * C[a][r]
            # print (C[a][r], M[a][r] )
            prob += my_z[a+allele_num*r] <= my_x[a]
            prob += my_z[a+allele_num*r] <= my_y[r]
            prob += my_z[a+allele_num*r] >= my_x[a] + my_y[r] - 1

            sum_xm += my_z[a+allele_num*r] * M[a][r] 

        prob += sum_xc >= my_y[r]
        prob += sum_xm >= my_L[r]

    # obj = 0
    # for r in range(read_num):
    #     # obj += my_y[r]
    #     obj +=  my_L[r]


    obj = lpSum([my_L[r] for r in range(read_num)])
    prob += obj
    # print(prob)

    solver = PULP_CBC_CMD()
    prob.solve(solver)

    type_result = []
    for i in prob.variables():
        if 'X' in i.name:
            type_result.append(i.varValue)
        #     print (i)
        # print (i, i.varValue)
    # return alpha
    print (type_result)
    return type_result

    # cost=LpVariable.dicts("cost",range(len(geno_set)*4),lowBound = 0)
    # beta_cost=LpVariable.dicts("beta_cost",range(len(geno_set)),lowBound = 0)

def downsample_matrix(M, max_num=500):
    # downsample reads
    # if len(M[0]) < max_num:
    #     return M
    # else:
    #     index_array = range(len(M[0]))
    #     random.shuffle(index_array)
    #     new_M = []
    #     for i in index_array[:max_num]:
    #         new_M.append(M[i])
    #     return new_M
    return M

def cal_zero_map(M):
    zero_read_dict = {}
    for z in range(len(M[0])):
        for a in range(len(M)):
            if M[a][z] == 0:
                zero_read_dict[z] = 1
    return zero_read_dict

def model2(C, M, allele_match_dict, record_index_allele, record_read_index, gene, read_identity_matrix,read_mismatch_matrix,record_allele_length):
    ## enumeration
    read_name_list = [0] * len(record_read_index[gene])
    for read_name in record_read_index[gene]:
        read_name_list[record_read_index[gene][read_name]] = read_name
    # M = downsample_matrix(M, 1000)
    print (len(M[0]))
    zero_read_dict = cal_zero_map(M)
    # C_sum, M_sum = 0, 0
    # for a in range(len(C)):
    #     # print (a, sum(C[a]), sum(M[a]))
    #     C_sum += sum(C[a])
    #     M_sum += sum(M[a])
    # print (read_num, "C_sum, M_sum", C_sum, M_sum, len(M))
    record_allele_pair_match_len = {}
    record_allele_pair_mismatch = defaultdict(int)
    record_allele_pair_identity = defaultdict(float)
    record_allele_pair_sep_match = {}
    record_allele_assign_reads_count = {}

    max_map_sum = 0
    allele_num = len(C)
    read_num = len(C[0])

    for i in range(allele_num):
        for j in range(i+1, allele_num):
            map_sum = 0
            tag = record_index_allele[i] + "&" + record_index_allele[j]
            # record_allele_pair_sep_match[tag] = {"a1_match":0, "a1_mismatch":0, "a1_identity":0, "a1_depth":0, "a2_match":0, "a2_mismatch":0, "a2_identity":0, "a2_depth":0}
            record_allele_pair_sep_match[tag] = {}
            record_allele_pair_sep_match[tag][record_index_allele[i]] = {"match":0, "mismatch":0, "identity":0, "depth":0}
            record_allele_pair_sep_match[tag][record_index_allele[j]] = {"match":0, "mismatch":0, "identity":0, "depth":0}
            
            allele_assigned_reads = [0, 0]

            for z in range(len(M[0])):

                ### adjustment, how to assign reads, not only by match length, also consider identity
                larger_index = determine_largest(read_identity_matrix[i][z], read_identity_matrix[j][z])
                if larger_index == 0:
                # if read_identity_matrix[i][z] > read_identity_matrix[j][z]:
                # if M[i][z] >= M[j][z]:
                    map_sum += M[i][z]
                    record_allele_pair_mismatch[tag]  += read_mismatch_matrix[i][z]
                    allele_assigned_reads[0] += 1
                    record_allele_pair_sep_match[tag][record_index_allele[i]]["match"] += M[i][z]
                    record_allele_pair_sep_match[tag][record_index_allele[i]]["mismatch"] += read_mismatch_matrix[i][z]
                else:
                    map_sum += M[j][z]
                    record_allele_pair_mismatch[tag]  += read_mismatch_matrix[j][z]
                    allele_assigned_reads[1] += 1
                    record_allele_pair_sep_match[tag][record_index_allele[j]]["match"] += M[j][z]
                    record_allele_pair_sep_match[tag][record_index_allele[j]]["mismatch"] += read_mismatch_matrix[j][z]


            # print (record_index_allele[i], record_index_allele[j], allele_assigned_reads, map_sum, max_map_sum, sep="\t")
            allele_names = [record_index_allele[i], record_index_allele[j]]
            for k in range(2):
                if allele_names[k] not in record_allele_assign_reads_count:
                    record_allele_assign_reads_count[allele_names[k]] = allele_assigned_reads[k]
                elif allele_assigned_reads[k] < record_allele_assign_reads_count[allele_names[k]]:
                    record_allele_assign_reads_count[allele_names[k]] = allele_assigned_reads[k]


            record_allele_pair_match_len[tag] = map_sum
            record_allele_pair_identity[tag] = map_sum/(map_sum + record_allele_pair_mismatch[tag])

            # print ("###", record_allele_pair_sep_match[tag])
            if (record_allele_pair_sep_match[tag][record_index_allele[i]]["match"] + record_allele_pair_sep_match[tag][record_index_allele[i]]["mismatch"]) > 0:
                record_allele_pair_sep_match[tag][record_index_allele[i]]["identity"] = record_allele_pair_sep_match[tag][record_index_allele[i]]["match"]/(record_allele_pair_sep_match[tag][record_index_allele[i]]["match"] + record_allele_pair_sep_match[tag][record_index_allele[i]]["mismatch"])
            if (record_allele_pair_sep_match[tag][record_index_allele[j]]["match"] + record_allele_pair_sep_match[tag][record_index_allele[j]]["mismatch"]) > 0:
                record_allele_pair_sep_match[tag][record_index_allele[j]]["identity"] = record_allele_pair_sep_match[tag][record_index_allele[j]]["match"]/(record_allele_pair_sep_match[tag][record_index_allele[j]]["match"] + record_allele_pair_sep_match[tag][record_index_allele[j]]["mismatch"])
            record_allele_pair_sep_match[tag][record_index_allele[i]]["depth"] = (record_allele_pair_sep_match[tag][record_index_allele[i]]["match"] + record_allele_pair_sep_match[tag][record_index_allele[i]]["mismatch"])/record_allele_length[record_index_allele[i]]
            record_allele_pair_sep_match[tag][record_index_allele[j]]["depth"] = (record_allele_pair_sep_match[tag][record_index_allele[j]]["match"] + record_allele_pair_sep_match[tag][record_index_allele[j]]["mismatch"])/record_allele_length[record_index_allele[j]]
            
            if map_sum >= max_map_sum:
                # print (">>>>>>>>")
                max_map_sum = map_sum
                type_result = [0] * allele_num
                type_result[i] = 1
                type_result[j] = 1
    print (type_result, max_map_sum)


    # if gene == "DQA1":
    #     alleles_after_read_assignment = select_alleles_by_read_assignment(record_allele_assign_reads_count, 0.5)
    # # else:
    # #     alleles_after_read_assignment = select_alleles_by_read_assignment(record_allele_assign_reads_count)
    #     record_allele_pair_match_len, record_allele_pair_identity = discard_alleles(alleles_after_read_assignment, record_allele_pair_match_len, record_allele_pair_identity)

    if gene in ["DPA1", 'DRB1'] :  # not DPB1
        allele_pair_after_identity_diff = select_alleles_by_identity_diff(record_allele_pair_sep_match, 0.1, 0.1)
        record_allele_pair_match_len, record_allele_pair_identity = discard_alleles(allele_pair_after_identity_diff, record_allele_pair_match_len, record_allele_pair_identity, "pair")
    

        

    tag_list, highest_score = choose_best_alleles(gene, record_allele_pair_match_len, record_allele_pair_identity,record_allele_pair_sep_match)
    tag_list = order_result_pair(tag_list, record_allele_pair_sep_match)
    type_allele_result =  generate_output(tag_list)          
    # for i in range(len(sorted_record_allele_pair_match_len)):
    #     if sorted_record_allele_pair_match_len[i][0] == 'C*04:01:01:01&C*03:03:01:01':
    #         print (sorted_record_allele_pair_match_len[i][0].split("&"), sorted_record_allele_pair_match_len[i][1], record_allele_pair_mismatch[sorted_record_allele_pair_match_len[i][0]], sorted_record_allele_pair_match_len[i][1]/(record_allele_pair_mismatch[sorted_record_allele_pair_match_len[i][0]]+sorted_record_allele_pair_match_len[i][1]), sep="\t")
    
    return type_result, highest_score, type_allele_result

def determine_largest(a, b):
    if a > b:
        return 0
    elif a < b:
        return 1
    else:
        return random.randint(0, 1)

def select_alleles_by_read_assignment(record_allele_assign_reads_count, cutoff=0.5):
    sorted_record_allele_assign_reads_count = sorted(record_allele_assign_reads_count.items(), key=lambda x: x[1], reverse=True)
    highest_reads_num = sorted_record_allele_assign_reads_count[0][1]
    alleles_after_read_assignment = {}

    for i in range(min([100, len(sorted_record_allele_assign_reads_count)])):
        print (sorted_record_allele_assign_reads_count[i])
    
    for i in range(len(sorted_record_allele_assign_reads_count)):
        alleles_after_read_assignment[sorted_record_allele_assign_reads_count[i][0]] = sorted_record_allele_assign_reads_count[i][1]
        if sorted_record_allele_assign_reads_count[i][1]/highest_reads_num < cutoff:
            break
    return alleles_after_read_assignment

def select_alleles_by_identity_diff(record_allele_pair_sep_match, identity_cutoff=0.07, depth_cutoff= 0.25):
    allele_pair_after_identity_diff = {}
    for tag in record_allele_pair_sep_match:
        allele_list = tag.split("&")
        if abs(record_allele_pair_sep_match[tag][allele_list[1]]["identity"]  - record_allele_pair_sep_match[tag][allele_list[0]]["identity"]) > identity_cutoff:
            continue
        depth_l = [record_allele_pair_sep_match[tag][allele_list[0]]["depth"], record_allele_pair_sep_match[tag][allele_list[1]]["depth"]]
        if min(depth_l)/max(depth_l) < depth_cutoff:
            continue
        allele_pair_after_identity_diff[tag] = 1
    print ("alleles_after_identity_diff pair", len(allele_pair_after_identity_diff))
    return allele_pair_after_identity_diff


def discard_alleles(alleles_after_read_assignment, record_allele_pair_match_len, record_allele_pair_identity, flag="allele"):
    new_record_allele_pair_match_len = {}
    new_record_allele_pair_identity = {}
    for tag in record_allele_pair_match_len:
        if flag == "allele":
            allele_list = tag.split("&")
            if allele_list[0] in alleles_after_read_assignment and  allele_list[1] in alleles_after_read_assignment :
                new_record_allele_pair_match_len[tag] = record_allele_pair_match_len[tag]
                new_record_allele_pair_identity[tag] = record_allele_pair_identity[tag]
        elif flag == "pair":
            if tag in alleles_after_read_assignment:
                new_record_allele_pair_match_len[tag] = record_allele_pair_match_len[tag]
                new_record_allele_pair_identity[tag] = record_allele_pair_identity[tag]

    print (len(record_allele_pair_match_len), "remain allele pair num is", len(new_record_allele_pair_identity))
    return new_record_allele_pair_match_len, new_record_allele_pair_identity


def print_match_results(sorted_record_allele_pair_match_len, record_allele_pair_sep_match, gene, record_allele_pair_identity):
    outdir = args["o"] + "/" + args["n"]
    out = open(f"""{outdir}/{args["n"]}.{gene}.allele.match.csv""", 'w')

    data = []
    for i in range(len(sorted_record_allele_pair_match_len)):
        tag = sorted_record_allele_pair_match_len[i][0]
        allele_list = tag.split("&")
        print (i, int(sorted_record_allele_pair_match_len[i][1]), round(record_allele_pair_identity[tag],3), allele_list[0], \
            round(record_allele_pair_sep_match[tag][allele_list[0]]["identity"],3),round(record_allele_pair_sep_match[tag][allele_list[0]]["depth"]),\
             allele_list[1], round(record_allele_pair_sep_match[tag][allele_list[1]]["identity"],3), round(record_allele_pair_sep_match[tag][allele_list[1]]["depth"]), sep = ",", file = out)
    out.close()



def choose_best_alleles(gene, record_allele_pair_match_len, record_allele_pair_identity,record_allele_pair_sep_match):
    sorted_record_allele_pair_match_len = sorted(record_allele_pair_match_len.items(), key=lambda x: x[1], reverse=True)
    
    print_match_results(sorted_record_allele_pair_match_len, record_allele_pair_sep_match, gene, record_allele_pair_identity)

    highest_match_score = sorted_record_allele_pair_match_len[0][1]
    

    len_diff_cutoff = 1e-4 
    ide_diff_cutoff = 0.001
    # if gene  in ["DQA1", "DRB1", "DPA1"]:
    #     len_diff_cutoff = 1e-3
    if gene  in ["DPA1"]:
        len_diff_cutoff = 1e-2
    if gene  in ["DRB1"]:
        len_diff_cutoff = 5e-2
    # if gene  in ["DPB1"]:
    #     len_diff_cutoff = 1e-2

    good_length_dict = {}
    for i in range(len(sorted_record_allele_pair_match_len)):
        tag = sorted_record_allele_pair_match_len[i][0]
        allele_list = tag.split("&")
        print ("#pass len", int(sorted_record_allele_pair_match_len[i][1]), round(record_allele_pair_identity[tag],3), allele_list[0], \
            round(record_allele_pair_sep_match[tag][allele_list[0]]["identity"],3),round(record_allele_pair_sep_match[tag][allele_list[0]]["depth"]),\
             allele_list[1], round(record_allele_pair_sep_match[tag][allele_list[1]]["identity"],3), round(record_allele_pair_sep_match[tag][allele_list[1]]["depth"]), sep="\t")
        if (highest_match_score - sorted_record_allele_pair_match_len[i][1])/highest_match_score <= len_diff_cutoff:
            good_length_dict[tag] = record_allele_pair_identity[tag]
        else:
            break

    for i in range(len(sorted_record_allele_pair_match_len)):
        tag = sorted_record_allele_pair_match_len[i][0]
        allele_list = tag.split("&")
        print ("#", int(sorted_record_allele_pair_match_len[i][1]), round(record_allele_pair_identity[tag],3), allele_list[0], \
            round(record_allele_pair_sep_match[tag][allele_list[0]]["identity"],3),round(record_allele_pair_sep_match[tag][allele_list[0]]["depth"]),\
             allele_list[1], round(record_allele_pair_sep_match[tag][allele_list[1]]["identity"],3), round(record_allele_pair_sep_match[tag][allele_list[1]]["depth"]), sep="\t")

        if i > 100:
            break
    identity_sorted_list = sorted(good_length_dict.items(), key=lambda x: x[1], reverse=True)
    match_len_with_max_identity = identity_sorted_list[0][1]
    full_result_list = []
    for i in range(len(identity_sorted_list)):
        if (match_len_with_max_identity - identity_sorted_list[i][1])/match_len_with_max_identity <= ide_diff_cutoff:
            print ("identity", identity_sorted_list[i])
            full_result_list.append(identity_sorted_list[i][0])
            
            if len(full_result_list) >= 30:
                break


    # first consider identity, then consider match length
    # if gene == "DQA1":
    #     len_diff_cutoff = 1e-2
    #     ide_diff_cutoff = 4e-3
    #     identity_sorted_list = sorted(record_allele_pair_identity.items(), key=lambda x: x[1], reverse=True)
    #     match_len_with_max_identity = identity_sorted_list[0][1]
    #     good_identity_dict = {}
    #     full_result_list = []
    #     for i in range(len(identity_sorted_list)):
    #         if (match_len_with_max_identity - identity_sorted_list[i][1])/match_len_with_max_identity <= ide_diff_cutoff:
    #             good_identity_dict[identity_sorted_list[i][0]] = record_allele_pair_match_len[identity_sorted_list[i][0]]
    #             print (identity_sorted_list[i], record_allele_pair_match_len[identity_sorted_list[i][0]])
    #         else:
    #             break
    #     match_len_sorted_list = sorted(good_identity_dict.items(), key=lambda x: x[1], reverse=True)
    #     match_len_with_max_identity = match_len_sorted_list[0][1]
    #     full_result_list = []
    #     for i in range(len(match_len_sorted_list)):
    #         if (match_len_with_max_identity - match_len_sorted_list[i][1])/match_len_with_max_identity <= len_diff_cutoff:
    #             print ("match length", match_len_sorted_list[i])
    #             full_result_list.append(match_len_sorted_list[i][0])
                
    #         if len(full_result_list) >= 30:
    #             break



    # print ("\nhomo\n")
    # for i in range(len(sorted_record_allele_pair_match_len)):
    #     tag = sorted_record_allele_pair_match_len[i][0]
    #     if tag.split("&")[0] == tag.split("&")[1]:
    #         print (tag.split("&"), int(sorted_record_allele_pair_match_len[i][1]), record_allele_pair_identity[tag], sep="\t")

    # print ("\ntruth\n")
    # tag = 'DQA1*01:03:01:01'  + "&" + "DQA1*01:04:02"
    # print (tag, record_allele_pair_match_len[tag], record_allele_pair_identity[tag], sep="\t")

    

    return full_result_list, highest_match_score


def choose_best_alleles_highest(gene, record_allele_pair_match_len, record_allele_pair_identity,record_allele_pair_sep_match):
    ## only report alleles with longest match length
    sorted_record_allele_pair_match_len = sorted(record_allele_pair_match_len.items(), key=lambda x: x[1], reverse=True)
    
    print_match_results(sorted_record_allele_pair_match_len, record_allele_pair_sep_match, gene, record_allele_pair_identity)

    highest_match_score = sorted_record_allele_pair_match_len[0][1]
    tag_list = []
    for i in range(len(sorted_record_allele_pair_match_len)):
        if i > 100:
            break

        tag = sorted_record_allele_pair_match_len[i][0]
        print (tag.split("&"), int(sorted_record_allele_pair_match_len[i][1]), record_allele_pair_identity[tag], sep="\t")
        
        if sorted_record_allele_pair_match_len[i][1] == highest_match_score:
            tag_list.append(tag)
        if gene == "C" and i < 2:
            tag_list.append(tag)
    return tag_list, highest_match_score

def generate_output(tag_list):
    type_allele_result = ['', '']
    for tag in tag_list:
        pair = tag.split("&")
        if len(type_allele_result[0]) == 0:
            type_allele_result = pair
            if type_allele_result[0] == 'C*04:01:01:11':
                type_allele_result[0] = 'C*04:01:01:01,C*04:01:01:11'
            if type_allele_result[1] == 'C*04:01:01:11':
                type_allele_result[1] = 'C*04:01:01:01,C*04:01:01:11'
        else:
            
            for p in pair:
                max_digit_num = [0, 0]
                keep = True
                for i in range(2):
                    for allele in type_allele_result[i].split(","):
                        if p == allele: 
                            keep = False
                        same_digit_num = cal_sim_of_alleles(p, allele)
                        # print (p, allele, same_digit_num)
                        if same_digit_num > max_digit_num[i]:
                            max_digit_num[i] = same_digit_num
                if keep:
                    if p == 'C*04:01:01:11':
                        p = 'C*04:01:01:01,C*04:01:01:11'
                    if max_digit_num[0] > max_digit_num[1]:
                        type_allele_result[0] += "," + p
                    elif max_digit_num[0] < max_digit_num[1]:
                        type_allele_result[1] += "," + p
                    else:
                        type_allele_result[0] += "," + p
                        type_allele_result[1] += "," + p

    return type_allele_result

def order_result_pair(type_allele_result, record_allele_pair_sep_match):
    sort_type_allele_result = []
    for tag in type_allele_result:
        allele_list = tag.split("&")

        if record_allele_pair_sep_match[tag][allele_list[1]]["depth"] > record_allele_pair_sep_match[tag][allele_list[0]]["depth"]:
            tag = allele_list[1] + "&" + allele_list[0]

        # if record_allele_pair_sep_match[tag][allele_list[1]]["identity"] > record_allele_pair_sep_match[tag][allele_list[0]]["identity"]:
        #     tag = allele_list[1] + "&" + allele_list[0]

        sort_type_allele_result.append(tag)
    return sort_type_allele_result



def cal_sim_of_alleles(allele1, allele2):
    clean_allele1 = allele1.split("*")[1].split(":")      
    clean_allele2 = allele2.split("*")[1].split(":")   
    same_digit_num = 0
    for i in range(min([len(clean_allele1), len(clean_allele2)])): 
        if clean_allele1[i] == clean_allele2[i]:
            same_digit_num += 2
        else:
            break
    return same_digit_num

def model(C, M, allele_match_dict):

    ## C: traditional read hit matrix   M: read map len matrix
    # allele_num = 3
    # read_num = 4

    M = downsample_matrix(M)

    allele_num = len(C)
    read_num = len(C[0])
    t_max = 2
    t_min = 1

    C_sum, M_sum = 0, 0
    for a in range(len(C)):
        # print (a, sum(C[a]), sum(M[a]))
        C_sum += sum(C[a])
        M_sum += sum(M[a])
    print (read_num, "C_sum, M_sum", C_sum, M_sum, len(M))

    # C = [[1,1,1,1], [0,0,0,0], [1,1,1,1]]
    prob = LpProblem('HLA typing', LpMaximize)
    my_x =LpVariable.dicts("X", range(allele_num), cat = "Binary", lowBound = 0, upBound = 1)
    my_z =LpVariable.dicts("Z", range(allele_num*read_num), cat = "Binary", lowBound = 0, upBound = 1)  ## get x*A
    # my_e =LpVariable.dicts("E", range(allele_num*read_num), cat = "Binary", lowBound = 0, upBound = 1)  ## get x*y
    my_A =LpVariable.dicts("A", range(allele_num*read_num), cat = "Binary", lowBound = 0, upBound = 1)  ## get x*y
    my_L =LpVariable.dicts("L", range(read_num), lowBound = 0, cat='Integer')  # total match len from all reads
    # my_g =LpVariable.dicts("G", range(allele_num*read_num), lowBound = 0, cat='Integer')  # whether two alleles

    # print (my_z)
    sum_x = 0   ## allele number
    for a in range(allele_num):
        sum_x += my_x[a]

    prob += sum_x <= t_max
    prob += sum_x >= t_min

    for r in range(read_num):
        sum_a = 0  # a read should be only assigned to a allele
        sum_xm = 0
        for a in range(allele_num):

            prob += my_z[a+allele_num*r] <= my_x[a]
            prob += my_z[a+allele_num*r] <= my_A[a+allele_num*r]
            prob += my_z[a+allele_num*r] >= my_x[a] + my_A[a+allele_num*r] - 1

            sum_xm += my_z[a+allele_num*r] * M[a][r] 
            # sum_a += my_z[a+allele_num*r] * C[a][r]
            sum_a += my_A[a+allele_num*r]

            # prob += my_g[a+allele_num*r] <= my_z[a+allele_num*r]
            # prob += my_g[a+allele_num*r] <= sum_x - 1
            # prob += my_g[a+allele_num*r] >= sum_x - 1 - (1-my_A[a+allele_num*r])
            # prob += sum_xm - my_g[a+allele_num*r]*beta*M[a][r]

        prob += sum_a <= 1 # constraint the number of alleles a read can be assigned to
        prob += sum_xm >= my_L[r]

    obj = 0
    for r in range(read_num):
        # obj += my_y[r]
        obj +=  my_L[r]
    # obj -= sum_x * gamma
    prob += obj
    


    # obj = lpSum([my_L[r] for r in range(read_num)])
    # prob += obj
    # print(prob)

    solver = PULP_CBC_CMD()
    prob.solve(solver)

    final_A_sum = 0
    objective_value = 0
    type_result = np.empty(allele_num)
    for i in prob.variables():
        if 'X' in i.name:
            index = int(i.name.split("_")[1])
            # print (i.name, i.varValue)
            type_result[index] = i.varValue
        elif 'L' in i.name:
            objective_value += i.varValue
        elif 'A' in i.name:
            final_A_sum += i.varValue
    print ("read num", read_num, "final_A_sum", final_A_sum)
    print (type_result)
    return type_result, objective_value

    # cost=LpVariable.dicts("cost",range(len(geno_set)*4),lowBound = 0)
    # beta_cost=LpVariable.dicts("beta_cost",range(len(geno_set)),lowBound = 0)

def map2db(args, gene):

    minimap_para = ''
    if args["y"] == "pacbio":
        minimap_para = " -x map-pb "
    elif args["y"] == "nanopore":
        minimap_para = " -x map-ont "

    outdir = args["o"] + "/" + args["n"]
    sam = outdir + "/" + args["n"] + "." + gene + ".db.sam"
    bam = outdir + "/" + args["n"] + "." + gene + ".db.bam"
    depth_file = outdir + "/" + args["n"] + "." + gene + ".db.depth"
    sort_depth_file = outdir + "/" + args["n"] + "." + gene + ".db.sort.depth.txt"
    # ref={args["f"] }
    # map raw reads to database
    alignDB_order = f"""
    fq={args["r"] }
    
    outdir={args["o"]}/{args["n"] }
    sample={args["n"] }

    fq=$outdir/{gene}.long_read.fq.gz
    ref={args["db"]}/HLA/whole/HLA_{gene}.fasta

    minimap2 -t {args["j"] } {minimap_para} -p 0.1 -N 100000 -a $ref $fq > {sam}
    samtools view -bS -F 0x800  {sam} | samtools sort - >{bam}
    samtools depth -aa {bam}>{depth_file}
    echo alignment done.
    """
    os.system(alignDB_order)
    return sam, depth_file, sort_depth_file


if __name__ == "__main__":
    # depth_file = "/mnt/d/HLAPro_backup/Nanopore_optimize/output/fredhutch-hla-1408-1012/fredhutch-hla-1408-1012.db.depth"
    # get_depth = Get_depth(depth_file)
    # get_depth.record_depth()
    # record_candidate_alleles = get_depth.select()

    # model()
    # construct_matrix()
    # sam = "/mnt/d/HLAPro_backup/Nanopore_optimize/output/fredhutch-hla-1408-1012/fredhutch-hla-1408-1012.db.bam"
    # main(sam)

    parser = argparse.ArgumentParser(description="HLA Typing with long-read data.", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    required.add_argument("-f", type=str, help="IMGT reference.", metavar="\b")
    required.add_argument("-r", type=str, help="Long-read fastq file. PacBio or Nanopore.", metavar="\b")
    required.add_argument("-n", type=str, help="Sample ID", metavar="\b")
    required.add_argument("-o", type=str, help="The output folder to store the typing results.", metavar="\b", default="./output")

    optional.add_argument("-j", type=int, help="Number of threads.", metavar="\b", default=5)
    optional.add_argument("-m", type=int, help="Maintain this number of alleles for ILP step.", metavar="\b", default=100)
    optional.add_argument("-b", type=float, help="The match length increase ratio lower than this value is homo [0-1].", metavar="\b", default=0.0007)
    optional.add_argument("--db", type=str, help="db dir.", metavar="\b", default=sys.path[0] + "/../db/")
    # optional.add_argument("-g", type=int, help="Whether use G group resolution annotation [0|1].", metavar="\b", default=0)
    # optional.add_argument("-m", type=int, help="1 represents typing, 0 means only read assignment", metavar="\b", default=1)
    optional.add_argument("-y", type=str, help="Read type, [nanopore|pacbio].", metavar="\b", default="pacbio")
    # optional.add_argument("-u", type=str, help="Choose full-length or exon typing. 0 indicates full-length, 1 means exon.", metavar="\b", default="0")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)




    
    # ref = "/mnt/d/HLAPro_backup/Nanopore_optimize/SpecHLA/db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta"

    main(args)





                # if max([M[i][z], M[j][z]]) < 500:
                # identity_i = 1 - read_identity_matrix[i][z]/(M[i][z] + read_identity_matrix[i][z])
                # identity_j = 1 - read_identity_matrix[j][z]/(M[j][z] + read_identity_matrix[j][z])
                # if max([identity_i, identity_j]) < 0.85:
                # #     continue
                # # if min([M[i][z], M[j][z]]) == 0:
                #     if record_index_allele[i] == 'C*16:01:01:01' and record_index_allele[j] == 'C*06:02:01:01' :#:'C*07:02:01:15':
                #         print (z, read_name_list[z], M[i][z], M[j][z], identity_i, identity_j)
                #     # continue
                # if z in zero_read_dict:
                #     continue
                
                # if record_index_allele[i] == 'C*04:01:01:01' and record_index_allele[j] == 'C*03:03:01:01' :#:'C*07:02:01:15':
                #     if M[i][z] > M[j][z] :#and min([M[i][z], M[j][z]]) ==0:
                #         print (read_name_list[z], record_index_allele[i], record_index_allele[j], int(M[i][z]), int(M[j][z]), int(max([M[i][z], M[j][z]])), read_identity_matrix[i][z],  read_identity_matrix[j][z],sep="\t")

                # if record_index_allele[i] == 'C*04:01:01:11' and record_index_allele[j] == 'C*03:03:01:01' :#:'C*07:02:01:15':
                #     if M[i][z] > M[j][z] :#and min([M[i][z], M[j][z]]) ==0:
                #         print (read_name_list[z], record_index_allele[i], record_index_allele[j], int(M[i][z]), int(M[j][z]), int(max([M[i][z], M[j][z]])), read_identity_matrix[i][z],  read_identity_matrix[j][z],sep="\t")
                # if M[i][z] >= M[j][z]:
                #     map_sum += M[i][z]
                #     record_allele_pair_mismatch[tag]  += read_identity_matrix[i][z]
                # else:
                #     map_sum += M[j][z]
                #     record_allele_pair_mismatch[tag]  += read_identity_matrix[j][z]
                # if record_index_allele[i] == 'DPA1*02:02:02:01' and record_index_allele[j] == 'DPA1*01:03:01:04' :
                #     print (read_identity_matrix[i][z], M[i][z], read_identity_matrix[j][z], M[j][z], sep = "\t")