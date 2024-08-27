import sys

MATCH_LEN_THRESHOLD = 2

class My_read():

    def __init__(self):
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

        self.match_rate = None
        self.mismatch_rate = None
        self.gap_ends_flag = False
        self.SA_flag = False
        self.XA_flag = False
        self.first_SA = '' 

        # self.get_match_length(read)

        # if self.read_name == "5d767760-9f0d-43a8-a115-0d4535e88218":
        #     if self.allele_name == "HLA-A*68:278" or self.allele_name == "HLA-A*30:01:01:01":
        #         print ("***", self.allele_name, self.match_num, self.mismatch_num, self.alignment_len, self.identity, self.match_start_pos, self.match_end_pos, self.reference_start, self.reference_end, self.loci_name, self.read_length, self.read_match_ratio, self.primary)
        
    def load_bam(self, read):
        ### the NM tag consists all insertion, deletion and mismatches in the alignment

        # Get the alignment length in the read
        # self.alignment_len = read.query_alignment_length
        gap_num, map_len, long_gap = self.count_cigar(read)
        self.alignment_len = map_len

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
            if ta[0] == 'AS':
                self.alignment_score = ta[1]
            if ta[0] == 'SA':
                # print (ta)
                self.SA_flag = True
                self.first_SA = ta[1].split(",")[0]
            if ta[0] == 'XA':
                # print (ta)
                self.XA_flag = True
        
        if self.alignment_len == 0:
            # print ("exit as the read is unmapped", read.query_name, read.cigar)
            # sys.exit(0)
            raise ValueError("unmapped read")
            
        
        self.match_num = self.alignment_len - self.mismatch_num 
        self.identity = self.match_num/(self.alignment_len - long_gap)

        self.read_name = read.query_name
        self.allele_name = read.reference_name
        self.loci_name = self.allele_name.split("*")[0]

        # if self.loci_name == "KIR2DL5A" or self.loci_name == "KIR2DL5B":
        #     self.loci_name = "KIR2DL5"
        if self.alignment_len < 400 and self.allele_name[0:3] == "KIR":
            self.identity = 0
        # if self.SA_flag and self.loci_name == "KIR3DP1" and self.first_SA.split("*")[0] == "KIR2DL4":
        #     self.identity = 0
        # if self.XA_flag and self.loci_name == "KIR3DP1":
        #     self.identity = 0

        self.match_rate = self.identity
        self.mismatch_rate = 1 - self.match_rate

        if self.match_num < 0:
            print ("negative match num", self.match_num, self.mismatch_num, self.alignment_len, self.identity, self.read_name, self.allele_name, self.loci_name)
            print (gap_num, map_len, long_gap)
            sys.exit(0)
        

    def count_cigar(self, read, long_cutoff = 50):
        gap_num = 0
        map_len = 0
        long_gap = 0
        

        for ci in read.cigar:
            if ci[0] == 0:
                map_len += ci[1]
            elif ci[0] == 1 or ci[0] == 2:
                gap_num += ci[1]
                map_len += ci[1]
                if ci[1] > long_cutoff:
                    long_gap += ci[1]
        self.gap_ends(read)
        return gap_num, map_len, long_gap
    
    ## check if there is long gap in the two ends
    def gap_ends(self, read, super_long_cutoff=300):
        # if read.query_name == "m64076_200603_055852/5440181/ccs":
        #     print (read.cigar, read.cigar[0][0], read.cigar[0][1], read.cigar[-1][0], read.cigar[-1][1])
        ## check if there is more than three cigars
        if len(read.cigar) >= 3 and read.cigar[0][1] > super_long_cutoff and  read.cigar[-1][1] > super_long_cutoff\
              and read.cigar[0][0] in [4,5] and read.cigar[-1][0] in [4,5]: # 4:S 5:H
            # print ("gap ends", read.query_name, read.cigar)
            self.gap_ends_flag = True



    def load_blast(self, line): # format 7
        field = line.strip().split("\t")

        # Get the alignment length in the read
        self.alignment_len = int(field[3])

        # Get the alignment start position in the read
        self.match_start_pos = int(field[6])

        # Get the alignment end position in the read
        self.match_end_pos = int(field[7])

        ## reverse match_start_pos and match_end_pos if match_end_pos < match_start_pos
        if self.match_end_pos < self.match_start_pos:
            self.match_start_pos, self.match_end_pos = self.match_end_pos, self.match_start_pos


        self.reference_start = int(field[8])
        self.reference_end = int(field[9])
        ## reverse reference_start and reference_end if reference_end < reference_start
        if self.reference_end < self.reference_start:
            self.reference_start, self.reference_end = self.reference_end, self.reference_start

        self.primary = False

        
        if self.alignment_len == 0:
            print ("unmapped read")
            sys.exit(0)

        # self.identity = float(field[2])/100
        # self.match_num = round(self.identity * self.alignment_len)
        # self.mismatch_num = self.alignment_len - self.match_num

        self.mismatch_num = int(field[4]) + int(field[5])
        self.match_num = self.alignment_len - self.mismatch_num
        self.identity = self.match_num/self.alignment_len

        self.read_name = field[0]
        self.allele_name = field[1]
        self.loci_name = self.allele_name.split("*")[0]

        # if self.loci_name == "KIR2DL5A" or self.loci_name == "KIR2DL5B":
        #     self.loci_name = "KIR2DL5"
        if self.alignment_len < 400 and self.allele_name[0:3] == "KIR":
            self.identity = 0

        self.match_rate = self.identity
        self.mismatch_rate = 1 - self.match_rate

    def load_second_blast(self, line): # format 7
        field = line.strip().split("\t")

        # Get the alignment length in the read
        self.alignment_len += int(field[3])

        # Get the alignment start position in the read
        match_start_pos = int(field[6])

        # Get the alignment end position in the read
        match_end_pos = int(field[7])

        ## reverse match_start_pos and match_end_pos if match_end_pos < match_start_pos
        if match_end_pos < match_start_pos:
            match_start_pos, match_end_pos = match_end_pos, match_start_pos
        if match_start_pos < self.match_start_pos:
            self.match_start_pos = match_start_pos
        if match_end_pos > self.match_end_pos:
            self.match_end_pos = match_end_pos


        reference_start = int(field[8])
        reference_end = int(field[9])
        ## reverse reference_start and reference_end if reference_end < reference_start
        if reference_end < reference_start:
            reference_start, reference_end = reference_end, reference_start
        if reference_start < self.reference_start:
            self.reference_start = reference_start
        if reference_end > self.reference_end:
            self.reference_end = reference_end

        self.primary = False

        
        if self.alignment_len == 0:
            print ("unmapped read")
            sys.exit(0)
        
        self.mismatch_num += (int(field[4]) + int(field[5]))
        self.match_num = self.alignment_len - self.mismatch_num
        self.identity = self.match_num/self.alignment_len

        self.read_name = field[0]
        self.allele_name = field[1]
        self.loci_name = self.allele_name.split("*")[0]

        # if self.loci_name == "KIR2DL5A" or self.loci_name == "KIR2DL5B":
        #     self.loci_name = "KIR2DL5"
        if self.alignment_len < 400 and self.allele_name[0:3] == "KIR":
            self.identity = 0

        self.match_rate = self.identity
        self.mismatch_rate = 1 - self.match_rate

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
        ### a new record should have a higher identity and match length not too short (>0.5*old record match length)
        if read_obj.identity > self.represent_identity and self.represent_match_num/read_obj.match_num < MATCH_LEN_THRESHOLD:
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

        highest_identity = gene_score[0][1]
        if highest_identity < identity_cutoff:   ## if the identity is too low, then skip
            return assigned_locus
        
        #### avoid the reads with short match length to be assigned to HLA-DRB1
        # if gene_score[0][0] == "HLA-DRB1" and self.loci_object_dict[gene_score[0][0]].represent_match_num < 300:
        #     print ("skip", read_name, gene_score, self.loci_object_dict[gene_score[0][0]].represent_match_num)
        #     return assigned_locus
        
        for i in range(len(gene_score)):
            loci_name = gene_score[i][0]

            if highest_identity - self.loci_object_dict[loci_name].represent_identity > identity_diff:
                break

            #### avoid the reads with short match length to be assigned to HLA-DRB1
            if loci_name == "HLA-DRB1" and self.loci_object_dict[loci_name].represent_match_num < 300:
                # print ("skip", read_name, loci_name, self.loci_object_dict[loci_name].represent_match_num, gene_score)
                continue
            # if loci_name in ['CYP2A13', 'CYP2A6', 'CYP2B6', 'CYP2C19', 'CYP2C8', 'CYP2C9', 'CYP3A4', 'CYP3A5', 'CYP4F2', 'NAT2', 'NUDT15', 'SLCO1B1']:
            #     if self.loci_object_dict[loci_name].represent_match_num < 1000 \
            #         or self.loci_object_dict[loci_name].represent_identity < 0.92:
            #         continue

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
                else:  # cannot find the distance in the matrix
                    confict_flag = True
                # if best_locus in ['CYP2D6', 'CYP2D7'] and loci_name in ['CYP2D6', 'CYP2D7']:
                #     confict_flag = False

                if local_distance <= 0:
                    # print ("#overlap", store_intervalue, my_interval, assigned_locus, loci_name, highest_identity, self.loci_object_dict[loci_name].represent_identity)
                    confict_flag = True
                    # continue
                # else:
                #     print ("non overlap", store_intervalue, my_interval)
            if confict_flag:
                ## if the match len of this locus is far long than the first, assign the read to this locus
                if self.loci_object_dict[loci_name].represent_match_num / self.loci_object_dict[assigned_locus[0]].represent_match_num > MATCH_LEN_THRESHOLD:
                    assigned_locus = [loci_name]
                    record_interval = my_interval
                    # print (read_name, 
                    #        loci_name, 
                    #        self.loci_object_dict[loci_name].represent_match_num, 
                    #        self.loci_object_dict[loci_name].represent_identity, 
                    #        best_locus, 
                    #        self.loci_object_dict[best_locus].represent_match_num,
                    #        self.loci_object_dict[best_locus].represent_identity)
                else:
                    continue
            else:
                assigned_locus.append(loci_name)
            record_interval += my_interval
        # print (assigned_locus)
        # print ("#####\n")
        # if len(assigned_locus) > 1 and assigned_locus[0] == assigned_locus[1]:
        #     print ("###########", gene_score, 
        #             assigned_locus, 
        #             read_name,
        #             self.loci_object_dict[gene_score[0][0]].represent_match_num, 
        #             self.loci_object_dict[gene_score[1][0]].represent_match_num)
        # assigned_locus = self.remove_reads(assigned_locus)
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

    def remove_reads(self, assigned_locus):
        ## independently remove the reads assigned to some genes
        cyp_others = ['CYP2A13', 'CYP2A6', 'CYP2B6', 'CYP2C19', 'CYP2C8', 'CYP2C9', 'CYP3A4', 'CYP3A5', 'CYP4F2', 'NAT2', 'NUDT15', 'SLCO1B1']
        new = []
        for loci_name in assigned_locus:
            if loci_name in cyp_others:
                if self.loci_object_dict[loci_name].represent_match_num > 3000 \
                    and self.loci_object_dict[loci_name].represent_identity > 0.95:
                    new.append(loci_name)
            else:
                new.append(loci_name)
        return new





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
