import sys

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
        self.reference_start = 0
        self.reference_end = 0
        self.loci_name = None
        self.read_length = 0
        self.read_match_ratio = 0
        self.primary = False

        self.get_match_length(read)
        self.match_rate = self.identity
        self.mismatch_rate = 1 - self.match_rate
        
    def get_match_length(self, read):
        ### the NM tag consists all insertion, deletion and mismatches in the alignment

        # Get the alignment length in the read
        self.alignment_len = read.query_alignment_length

        self.read_length = read.query_length
        # self.read_match_ratio = self.alignment_len/self.read_length

        # Get the alignment start position in the read
        self.match_start_pos = read.query_alignment_start

        # Get the alignment end position in the read
        self.match_end_pos = read.query_alignment_end

        self.reference_start = read.reference_start
        self.reference_end = read.reference_end
            
        if not read.is_secondary:
            self.primary = True

        for ta in read.get_tags():
            if ta[0] == 'NM':
                self.mismatch_num = ta[1]  
        
        if self.alignment_len == 0:
            print ("unmapped read", read.query_name, read.cigar)
            sys.exit(0)
        
        self.match_num = self.alignment_len - self.mismatch_num
        self.identity = self.match_num/self.alignment_len

        self.read_name = read.query_name
        self.allele_name = read.reference_name
        self.loci_name = self.allele_name.split("*")[0]

        if self.loci_name == "KIR2DL5A" or self.loci_name == "KIR2DL5B":
            self.loci_name = "KIR2DL5"
        if self.alignment_len < 400 and self.allele_name[0:3] == "KIR":
            self.identity = 0

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

## given two intervals, get the distance between them
def get_interval_distance(interval1, interval2):
    return max(interval1[0], interval2[0]) - min(interval1[1], interval2[1])

class Read_bin():  # the match for a single read in all the alleles of a locus

    def __init__(self, loci_object_dict):
        self.loci_object_dict = loci_object_dict
    
    def assign_multiple(self, distance_matrix, read_name, identity_cutoff, identity_diff, dist_cutoff):

        record_identity = {}
        for loci_name in self.loci_object_dict:
            record_identity[loci_name] = self.loci_object_dict[loci_name].represent_identity
            # print (loci_name, self.loci_object_dict[loci_name].represent_identity, self.loci_object_dict[loci_name].represent_match_num, self.loci_object_dict[loci_name].start,self.loci_object_dict[loci_name].end )
        assigned_locus = []
        record_interval = []
        gene_score = sorted(record_identity.items(), key=lambda item: item[1], reverse = True)

        # if read_name == "e55f5884-50a1-43bc-a3eb-1003f364caa5":
        #     print ("###########", gene_score)

        highest_identity = gene_score[0][1]
        if highest_identity < identity_cutoff:   ## if the identity is too low, then skip
            return assigned_locus
        
        for i in range(len(gene_score)):
            loci_name = gene_score[i][0]
            if highest_identity - self.loci_object_dict[loci_name].represent_identity > identity_diff:
                break
            my_interval = [self.loci_object_dict[loci_name].start,self.loci_object_dict[loci_name].end]  

            confict_flag = False
            if len(record_interval) >= 2: # if there is at least a start and a end, then check overlap
                store_intervalue = [min(record_interval), max(record_interval)]
                # print (store_intervalue)
                best_locus = assigned_locus[0]

                ### get the distance between store_intervalue and my_interval
                local_distance = get_interval_distance(store_intervalue, my_interval)

                ### get the distance between best_locus and loci_name, check if the key is in the matrix
                if (best_locus, loci_name) in distance_matrix:
                    distance = distance_matrix[(best_locus, loci_name)]

                    # skip if the distance in read is shorter than real distance
                    if distance - local_distance > dist_cutoff:
                        confict_flag = True
                        # continue
                if local_distance <= 0:
                    # print ("#overlap", store_intervalue, my_interval, assigned_locus, loci_name, highest_identity, self.loci_object_dict[loci_name].represent_identity)
                    confict_flag = True
                    # continue
                # else:
                #     print ("non overlap", store_intervalue, my_interval)
            if confict_flag:
                ## if the match len of this locus is far long than the first, assign the read to this locus
                if self.loci_object_dict[loci_name].represent_match_num / self.loci_object_dict[assigned_locus[0]].represent_match_num > 2:
                    assigned_locus = [loci_name]
                    record_interval = my_interval
                    print (read_name, 
                           loci_name, 
                           self.loci_object_dict[loci_name].represent_match_num, 
                           self.loci_object_dict[loci_name].represent_identity, 
                           best_locus, 
                           self.loci_object_dict[best_locus].represent_match_num,
                           self.loci_object_dict[best_locus].represent_identity)
                else:
                    continue
            # print (loci_name, self.loci_object_dict[loci_name].represent_identity, my_interval)
            assigned_locus.append(loci_name)
            record_interval += my_interval
        # print (assigned_locus)
        # print ("#####\n")
        return assigned_locus

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
