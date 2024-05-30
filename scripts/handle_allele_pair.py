import random

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

    def add_read(self, match_num, mismatch_num):
        self.match_num += match_num
        self.mismatch_num  += mismatch_num

    def get_identity(self):
        if self.match_num + self.mismatch_num == 0:
            self.identity = 0
        else:
            self.identity = self.match_num/ (self.match_num + self.mismatch_num)
    
    def get_depth(self, allele_length):
        self.depth = (self.match_num + self.mismatch_num)/allele_length


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
        for read_name in record_read_allele_dict:

            ### read not support both allele
            if self.allele_1 not in record_read_allele_dict[read_name] and self.allele_2 not in record_read_allele_dict[read_name]:
                continue
            if self.allele_1 not in record_read_allele_dict[read_name]:
                larger_index = 1
            elif self.allele_2 not in record_read_allele_dict[read_name]:
                larger_index = 0
            else: 
                larger_index = determine_largest(record_read_allele_dict[read_name][self.allele_1].identity, record_read_allele_dict[read_name][self.allele_2].identity)

            if larger_index == 0:
                self.allele_1_obj.add_read(record_read_allele_dict[read_name][self.allele_1].match_num, record_read_allele_dict[read_name][self.allele_1].mismatch_num)
                self.pair_obj.add_read(record_read_allele_dict[read_name][self.allele_1].match_num, record_read_allele_dict[read_name][self.allele_1].mismatch_num)
            else:
                self.allele_2_obj.add_read(record_read_allele_dict[read_name][self.allele_2].match_num, record_read_allele_dict[read_name][self.allele_2].mismatch_num)
                self.pair_obj.add_read(record_read_allele_dict[read_name][self.allele_2].match_num, record_read_allele_dict[read_name][self.allele_2].mismatch_num)

            # if self.allele_1 == "HLA-A*02:01:01:35" and self.allele_2 == "HLA-A*03:08:01:02":
            #     print (read_name, larger_index)

        self.allele_1_obj.get_identity()
        self.allele_2_obj.get_identity()
        self.pair_obj.get_identity()