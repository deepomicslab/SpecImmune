#!/usr/bin/env python3
from my_imports import *

class Reads():
    def __init__(self,vcffile,bamfile,conn_file):
        self.vcffile=vcffile
        self.bamfile=bamfile
        self.conn_file=conn_file
    def extract_snp(self):
        #need to filter the vcf file, maybe in the calling step.
        snp_list,freq_set=[],[]
        called_num,dp_sum=0,0
        f=open(self.vcffile,'r')
        for line in f:
            if line[0] != "#":
                called_num+=1
                line=line.strip()
                array=line.split()
                # snp_tag=str(array[0])+"_"+str(array[1])
                sub_array=array[-1].split(":")
                # print (sub_array,line)
                mini_array=sub_array[1].split(",")
                alt_dp=mini_array[1]
                dp=sub_array[2]
                dp_sum+=float(dp)
                if float(dp) == 0:
                    alt_freq=0
                else:
                    alt_freq=float(alt_dp)/float(dp)            
                freq_set.append(round(alt_freq,6))
                single_snp=[]
                single_snp.append(array[0])
                single_snp.append(array[1])
                single_snp.append(array[3])
                single_snp.append(array[4])
                snp_list.append(single_snp)
        print ("mean depth is",int(float(dp_sum)/len(snp_list)))
        return freq_set,snp_list
    def extract_delta(self):
        freq_set,snp_list=self.extract_snp()
        delta_set=init_delta(len(freq_set))
        f=open(self.conn_file,'r')
        for line in f:
            line=line.strip()
            array=line.split()
            if array[0] == '1': 
                delta_index=int(array[2])-1
                geno_type=array[3]
                for i in range(len(geno_type)-1):
                # delta_set[delta_index][int(geno_type[0])][int(geno_type[1])]+=1
                    delta_set[delta_index+i][int(geno_type[i])][int(geno_type[i+1])]+=1
                #print (geno_type[0],geno_type[1])
        # print (delta_set[:10])
        for delta in delta_set:
            sum_dp=sum(delta[0])+sum(delta[1])
            # print (sum_dp)
            if sum_dp>4:
                delta[0][0]=float(delta[0][0])/sum_dp
                delta[0][1]=float(delta[0][1])/sum_dp
                delta[1][0]=float(delta[1][0])/sum_dp
                delta[1][1]=float(delta[1][1])/sum_dp
            else:
                delta[:]=[[0,0],[0,0]]
        return freq_set,delta_set

def reads_info(bamfile,snp_set,cutoff):
    delta_set,share_set=[],[]
    # snp_set=[['AHZD01000001', '376', 'C', 'G'],['AHZD01000001', '376', 'G', 'C'],['AHZD01000001', '377', 'A', 'G']]
    samfile = pysam.AlignmentFile(bamfile, "rb")
    no_conn=0
    for i in range(len(snp_set)-1):
        first=snp_set[i]
        second=snp_set[i+1]
        # print (int(second[1])-int(first[1]))
        first_reads=reads_support(samfile,first)
        second_reads=reads_support(samfile,second)
        connect_dp=[[0,0],[0,0]]
        for i in range(2):
            for j in range(2):
                connect_dp[i][j]=list_overlap(first_reads[i],second_reads[j])
                # if connect_dp[i][j] <3:
                #     connect_dp[i][j]=0
        dp_sum=sum(connect_dp[0])+sum(connect_dp[1])
        share_set.append(dp_sum)
        # print (connect_dp)
        if dp_sum>cutoff:
            connect_ratio=[[0,0],[0,0]]
            for m in range(2):
                for n in range(2):
                    connect_ratio[m][n]=round(float(connect_dp[m][n])/dp_sum,6)
        else:
            connect_ratio=[[0.25,0.25],[0.25,0.25]]
            no_conn+=1
        #### gap too long can not trust
        # if abs(int(first[1]) - int(second[1])) >80:
        #     connect_ratio=[[0,0],[0,0]]
        delta_set.append(connect_ratio)
    # print ("%s snps in total, %s truncated points."%(len(snp_set),no_conn),no_conn/len(snp_set))
    # print ("edge sum",np.mean(four_edge),len(delta_set))
    return delta_set,share_set
def list_overlap(a,b):
    # print (list(set(a).intersection(set(b))))
    return len(set(a).intersection(set(b)))    
def reads_support(samfile,first):
    num=0
    save_reads=[[],[]]
    # for read in samfile.fetch('AHZD01000001', first, first+1):  #6,7 means locus 7 in bam file.
    for read in samfile.fetch(first[0],int(first[1])-1,int(first[1])):
        num+=1
        # print (read)
        # if int(first[1])-1 not in read.get_reference_positions(full_length=True):
        if int(first[1])-1 in read.get_reference_positions(full_length=False):
            # print (len(read.get_reference_positions(full_length=True)))
            # print (int(first[1])-1 )
            # print ('_______________')
            reads_index=read.get_reference_positions(full_length=True).index(int(first[1])-1)
            if read.query_sequence[reads_index] == first[2]:
                save_reads[0].append(read.query_name)
            elif read.query_sequence[reads_index] == first[3]:
                save_reads[1].append(read.query_name)
        # print (read.query_name)
        # print (read.query_sequence[reads_index])
        # print (read.query_sequence)
        # print (read.get_reference_positions(full_length=True))
    # print (save_reads,num)
    return save_reads