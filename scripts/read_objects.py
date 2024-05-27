




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

        self.get_match_length(read)
        


    def get_match_length(self, read):
   
        unmap_part_num = 0
        pos = 0
        for ci in read.cigar:
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
