import random
import numpy as np

def determine_largest(a, b, seed = 8):
    if a > b:
        return 0
    elif a < b:
        return 1
    else:
        return random.randint(0, 1)


class Align_Obj(object):

    def __init__(self, name):
        self.name = name
        self.alignment_length = 0
        self.match_num = 0
        self.mismatch_num = 0
        self.depth = 0
        self.identity = 0
        self.assigned_reads = []
        self.identity_list = []

        self.reference_start = float("inf")
        self.reference_end = 0
        self.coverage = 0
        self.median_identity = 0

    def add_read(self, match_num, mismatch_num):
        self.match_num += match_num
        self.mismatch_num  += mismatch_num
    
    def add_read_identity(self, identity):
        self.identity_list.append(identity)
    
    def add_read_ref(self, start, end):
        if start < self.reference_start:
            self.reference_start = start
        if end > self.reference_end:
            self.reference_end = end

    def get_identity(self):
        if self.match_num + self.mismatch_num == 0:
            self.identity = 0
            # print("Warning: Identity is 0", self.name)
        else:
            self.identity = self.match_num / (self.match_num + self.mismatch_num)
            # if self.name == "HLA-A*02:01:01:35&HLA-A*03:08:01:02":
            #     print (self.name, self.match_num, self.mismatch_num, self.identity)
    
    def get_median_identity(self):
        self.median_identity = np.median(self.identity_list)
    
    def get_depth(self, allele_length):
        self.depth = (self.match_num + self.mismatch_num)/allele_length

    def get_coverage(self, allele_length):
        self.coverage = (self.reference_end - self.reference_start)/allele_length


class My_allele_pair():

    def __init__(self, allele_1, allele_2):
        self.tag =  allele_1 + "&" + allele_2
        self.allele_1 = allele_1
        self.allele_2 = allele_2

        self.allele_1_obj = Align_Obj(allele_1)
        self.allele_2_obj = Align_Obj(allele_2)

        self.pair_obj = Align_Obj(self.tag)

        self.match_len_list = []
        self.identity_list = []
    
    def assign_reads(self, record_read_allele_dict):
        read_assign_dict = {}
        # test_tag = "HLA-S*01:01:01:01&HLA-S*01:02:01:03"
        # test_tag = "HLA-A*30:01:01:01&HLA-A*68:278"
        for read_name in record_read_allele_dict:

            ### read not support both allele
            if self.allele_1 not in record_read_allele_dict[read_name] and self.allele_2 not in record_read_allele_dict[read_name]:
                continue
            if self.allele_1 not in record_read_allele_dict[read_name]:
                larger_index = 1
            elif self.allele_2 not in record_read_allele_dict[read_name]:
                larger_index = 0
            else: 
                
                a1_identity = record_read_allele_dict[read_name][self.allele_1].identity
                a2_identity = record_read_allele_dict[read_name][self.allele_2].identity

                ## adjust identity based on match length
                a1_identity = a1_identity * 1.0001 ** (record_read_allele_dict[read_name][self.allele_1].match_num- record_read_allele_dict[read_name][self.allele_2].match_num)
                larger_index = determine_largest(a1_identity, a2_identity)

            if larger_index == 0:
                read_assign_dict[read_name] = [self.allele_1]
                self.allele_1_obj.add_read(record_read_allele_dict[read_name][self.allele_1].match_num, record_read_allele_dict[read_name][self.allele_1].mismatch_num)
                self.pair_obj.add_read(record_read_allele_dict[read_name][self.allele_1].match_num, record_read_allele_dict[read_name][self.allele_1].mismatch_num)
                self.pair_obj.add_read_identity(record_read_allele_dict[read_name][self.allele_1].identity)
                self.allele_1_obj.add_read_ref(record_read_allele_dict[read_name][self.allele_1].reference_start, record_read_allele_dict[read_name][self.allele_1].reference_end)
            else:
                read_assign_dict[read_name] = [self.allele_2]
                self.allele_2_obj.add_read(record_read_allele_dict[read_name][self.allele_2].match_num, record_read_allele_dict[read_name][self.allele_2].mismatch_num)
                self.pair_obj.add_read(record_read_allele_dict[read_name][self.allele_2].match_num, record_read_allele_dict[read_name][self.allele_2].mismatch_num)
                self.pair_obj.add_read_identity(record_read_allele_dict[read_name][self.allele_2].identity)
                self.allele_2_obj.add_read_ref(record_read_allele_dict[read_name][self.allele_2].reference_start, record_read_allele_dict[read_name][self.allele_2].reference_end)
            
            # if self.tag == test_tag:
            #     # print (read_name, larger_index, record_read_allele_dict[read_name][self.allele_1].identity, record_read_allele_dict[read_name][self.allele_2].identity)
            #     # print (read_name, larger_index, self.pair_obj.match_num)
            #     if record_read_allele_dict[read_name][self.allele_2].match_num > record_read_allele_dict[read_name][self.allele_1].match_num:
            #         if larger_index == 0:
            #             print (read_name, self.allele_1, record_read_allele_dict[read_name][self.allele_1].match_num,  self.allele_2, record_read_allele_dict[read_name][self.allele_2].match_num)

        # if self.tag == test_tag:
        #     print (self.pair_obj.match_num, len(read_assign_dict))

        self.allele_1_obj.get_identity()
        self.allele_2_obj.get_identity()
        self.pair_obj.get_identity()
        self.pair_obj.get_median_identity()

        return read_assign_dict