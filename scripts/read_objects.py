




class My_read():

    def __init__(self, read):
        self.alignment_len = 0
        self.match_num = 0
        self.mismatch_num = 0
        self.indel_num = 0
        self.identity = 0
        self.alignment_score = 0
        self.read_name = None
        self.allele_name = None
        self.match_start_pos = 0
        self.match_end_pos = 0
        self.loci_name = None
        self.read_length = 0
        self.primary = False

        self.get_match_length(read)
        self.match_rate = self.identity
        self.mismatch_rate = 1 - self.match_rate
        


    def get_match_length(self, read):
   
        unmap_part_num = 0
        pos = 0
        for ci in read.cigar:
            self.read_length += ci[1]
            if ci[0] == 0:  # M
                self.alignment_len += ci[1]
                if self.match_start_pos == 0:
                    self.match_start_pos = pos
                pos += ci[1]
                self.match_end_pos = pos
            elif ci[0] == 1 or ci[0] == 2:  # I and D
                self.indel_num += ci[1]
                pos += ci[1]
                self.match_end_pos = pos
            else:
                unmap_part_num += 1
                pos += ci[1]
            
        if not read.is_secondary:
            self.primary = True

        for ta in read.get_tags():
            if ta[0] == 'NM':
                self.mismatch_num = ta[1]  
            if ta[0] == 'AS':
                self.alignment_score = ta[1]  
        
        if self.alignment_len == 0:
            print (read.query_name, read.cigar)
            sys.exit(0)
        
        self.match_num = self.alignment_len - self.mismatch_num
        self.identity = self.match_num/self.alignment_len

        self.read_name = read.query_name
        self.allele_name = read.reference_name
        self.loci_name = self.allele_name.split("*")[0]

        if self.loci_name == "KIR2DL5A" or self.loci_name == "KIR2DL5B":
            self.loci_name = "KIR2DL5"
        if self.read_length < 400 and self.allele_name[0:3] == "KIR":
            self.match_rate = 0


class My_locus():  # the match for a single read in all the alleles of a locus

    def __init__(self):
        self.read_name = None
        self.start = float('inf')
        self.end = 0
        self.represent_identity = 0
        self.represent_match_num = 0
        self.loci_name = None

    def add_record(self, read_obj):

        if self.loci_name == None:
            self.loci_name = read_obj.loci_name
        else:
            if self.loci_name != read_obj.loci_name:
                print (self.loci_name, read_obj.loci_name, "read name error")
                sys.exit(0)

        self.read_name = read_obj.read_name
        # if read_obj.read_name == "m54329U_200715_194535/18942351/ccs":
        #     print (read_obj.allele_name, read_obj.identity, read_obj.match_num, read_obj.match_start_pos, read_obj.match_end_pos)
        if read_obj.identity > self.represent_identity:
            self.represent_identity = read_obj.identity
            self.represent_match_num = read_obj.match_num

            if read_obj.match_start_pos < self.start:
                self.start = read_obj.match_start_pos
            if read_obj.match_end_pos > self.end:
                self.end = read_obj.match_end_pos


class Read_bin():  # the match for a single read in all the alleles of a locus

    def __init__(self, loci_object_dict):
        self.loci_object_dict = loci_object_dict
    
    def assign_multiple(self, identity_diff = 0.01):

        record_identity = {}
        for loci_name in self.loci_object_dict:
            record_identity[loci_name] = self.loci_object_dict[loci_name].represent_identity
            # print (loci_name, self.loci_object_dict[loci_name].represent_identity, self.loci_object_dict[loci_name].represent_match_num, self.loci_object_dict[loci_name].start,self.loci_object_dict[loci_name].end )
        assigned_locus = []
        record_interval = []
        gene_score = sorted(record_identity.items(), key=lambda item: item[1], reverse = True)
        highest_identity = gene_score[0][1]
        for i in range(len(gene_score)):
            loci_name = gene_score[i][0]
            if highest_identity - self.loci_object_dict[loci_name].represent_identity > identity_diff:
                break
            my_interval = [self.loci_object_dict[loci_name].start,self.loci_object_dict[loci_name].end]  

            if len(record_interval) >= 2:
                store_intervalue = [min(record_interval), max(record_interval)]
                # print (store_intervalue)
                if self.intervals_overlap(store_intervalue, my_interval):
                    print ("overlap", store_intervalue, my_interval, assigned_locus, loci_name, highest_identity, self.loci_object_dict[loci_name].represent_identity)
                    continue
                # else:
                #     print ("non overlap", store_intervalue, my_interval)
            # print (loci_name, self.loci_object_dict[loci_name].represent_identity, my_interval)
            assigned_locus.append(loci_name)
            record_interval += my_interval
        # print (assigned_locus)
        # print ("#####\n")
        return assigned_locus

    # def intervals_overlap(self, interval1, interval2):
    #     # Check if intervals overlap
    #     if interval1[1] >= interval2[0] and interval2[1] >= interval1[0]:
    #         return True
    #     else:
    #         return False

    def intervals_overlap(self, interval1, interval2):
        """
        Check if two intervals overlap.
        interval1 and interval2 are tuples of the form (start, end)
        Returns True if the intervals overlap, False otherwise.
        """
        start1, end1 = interval1
        start2, end2 = interval2
        
        # Check if one interval is completely to the left of the other
        if end1 < start2 or end2 < start1:
            return False
        
        # Otherwise, the intervals overlap
        return True





"""
M

BAM_CMATCH

0

I

BAM_CINS

1

D

BAM_CDEL

2

N

BAM_CREF_SKIP

3

S

BAM_CSOFT_CLIP

4

H

BAM_CHARD_CLIP

5

P

BAM_CPAD

6

=

BAM_CEQUAL

7

X

BAM_CDIFF

8

B

BAM_CBACK

9
"""
