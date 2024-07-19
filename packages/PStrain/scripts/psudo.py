#!/usr/bin/env python3

from argparse import ArgumentParser
from my_imports import *
import pipeline

Usage = \
"""%(prog)s [options] -b/--bam <bam file> -v/--vcf <vcf file> -o/--outdir <output directory>

Before using this script, you should map the reads to your reference, and call SNPs. Then you can provide the bam and vcf file. \
Once you have the prior knowledge of the strains number, you can provie it the to script by -k/--strainsNum parameter. Otherwise,\
it will choose the strains number by elbow method. 
Help information can be found by %(prog)s -h/--help, additional information can be found in README.MD or www.xxx.com.
"""
scripts_dir=sys.path[0]+'/'
parser = ArgumentParser(description="Needle: identify strains in samples with single species.",prog='python3 single_species.py',usage=Usage)
optional=parser._action_groups.pop()
required=parser.add_argument_group('required arguments')
#necessary parameter
required.add_argument("-b", "--bam",help="The bam file of the input samples.",dest='bamfile',metavar='')
required.add_argument("-v", "--vcf",help="The vcf file of the input samples.",dest='vcffile',metavar='')
required.add_argument("-o", "--outdir",help="The output directory.",dest='outdir',metavar='')
#alternative parameter
optional.add_argument("-k", "--strainsNum",help="The number of strains in the sample (default is chosen by elbow method)."\
    ,dest='k',metavar='', type=int,default='0')
optional.add_argument("--snp_dp",help="The minimum depth of SNPs to be considered in strain profiling step (default is 5).",\
    dest='snp_dp',metavar='',default=6, type=int)
optional.add_argument("--extractHAIRS",help="Path to extractHAIRS binary (default version is in the packages path).",dest='extractHAIRS',metavar='',\
    default=scripts_dir+'../packages/HAPCUT2/extractHAIRS')
#added parameters
optional.add_argument("--prior",help="The file of prior knowledge of genotype frequencies in the population.",dest='prior',metavar='',\
    default=scripts_dir+'/prior_beta.pickle')
optional.add_argument("--elbow",help="The cutoff of elbow method while identifying strains number. \
If the loss reduction ratio is less than the cutoff, then the strains number is determined.",\
    dest='elbow',metavar='',default=0.24, type=float)
optional.add_argument("--qual",help="The minimum quality score of SNVs to be considered in strain profiling step (default is 20).",\
    dest='qual',metavar='',default=20, type=int)  
optional.add_argument("-w", "--weight",help="The weight of genotype frequencies while computing loss, then the weight of\
 linked read type frequencies is 1-w. The value is between 0~1. (default is 0.0)",dest='weight',metavar='',default=0.3, type=float)
optional.add_argument( "--lambda1",help="The weight of prior knowledge while rectifying genotype frequencies. \
The value is between 0~1. (default is 0.0)",dest='lambda1',metavar='',default=0.0, type=float)
optional.add_argument( "--lambda2",help="The weight of prior estimation while rectifying second order genotype frequencies. \
The value is between 0~1. (default is 0.0)",dest='lambda2',metavar='',default=0.1, type=float)    

parser._action_groups.append(optional)
args = parser.parse_args()



def read_vcf(vcffile,outdir,snp_dp):
    if not os.path.exists(outdir):
        os.system('mkdir '+outdir)
    snp_list,beta_set=[],[]
    locus_list=[]
    alt_homo=[]
    myvcf=VariantFile(vcffile)
    vcf_out=VariantFile('%s/filter.vcf'%(outdir),'w',header=myvcf.header)
    sample=list(myvcf.header.samples)[0]
    to_sort=[]
    for record in myvcf.fetch():
        if record.info['DP'] <1:
            continue
        geno=record.samples[sample]['GT']        
        depth=record.samples[sample]['AD']
        dp=depth[0]+depth[1]   
        if len(record.ref)==1 and len(record.alts)==1 and len(record.alts[0])==1 and dp>=snp_dp  and record.qual > 20:
            snp=[record.chrom,record.pos,record.ref,record.alts[0]]
            if geno==(0,1):
                to_sort.append(snp)
                vcf_out.write(record)
                beta=round(float(depth[1])/dp,6)
                snp_list.append(snp)
                beta_set.append(beta)
                locus_list.append(int(record.pos))
            elif geno==(1,1):
                to_sort.append(snp)
                alt_homo.append(snp)
    return snp_list,beta_set,to_sort
def delta(outdir,extractHAIRS,bamfile,beta_set):
    hapcut_order='%s --bam %s --VCF %s/filter.vcf --out %s/vcf.conn'%(extractHAIRS,bamfile,outdir,outdir)
    os.system(hapcut_order)
    snp_num=len(beta_set)
    delta_set=[]
    for i in range(snp_num-1):
        delta_set.append([[0,0],[0,0]])
    for line in open('%s/vcf.conn'%(outdir)):
        line=line.strip()
        array=line.split()
        if array[0] == '1': 
            delta_index=int(array[2])-1
            geno_type=array[3]
            for i in range(len(geno_type)-1):
                delta_set[delta_index+i][int(geno_type[i])][int(geno_type[i+1])]+=1
    for delta in delta_set:
        sum_dp=sum(delta[0])+sum(delta[1])
        if sum_dp>0:
            delta[0][0]=float(delta[0][0])/sum_dp
            delta[0][1]=float(delta[0][1])/sum_dp
            delta[1][0]=float(delta[1][0])/sum_dp
            delta[1][1]=float(delta[1][1])/sum_dp
        else:
            delta[:]=[[0.25,0.25],[0.25,0.25]]
    return delta_set
def output(outdir,final_alpha,seq_list,to_sort,snp_list):
    ra_file=open(outdir+'/strain_RA.txt','w')
    print ('# Strain_ID\tStrain_RA',file=ra_file)
    seq=seq_list
    snp=snp_list  
    alpha=final_alpha    
    for j in range(len(alpha)):
        print ('str-'+str(j+1),alpha[j],file=ra_file)
    sequence=open(outdir+'/strain_seq.txt','w')
    print ('# Gene\tLocus\tRef\tAlt\t',end='',file=sequence)
    for i in range(len(alpha)):
        print ('str-'+str(i+1),end='\t',file=sequence)
    print ('',file=sequence)
    he_num=len(snp)
    he=0
    ho=0       
    for point in to_sort:
        for info in point:
            print (info,end='\t',file=sequence)
        if he<he_num and point[0] == snp[he][0] and point[1] == snp[he][1]:
            for j in range(len(alpha)):
                print (seq[j][he],end='\t',file=sequence)
            he+=1
        else:
            for j in range(len(alpha)):
                print (1,end='\t',file=sequence)   
            ho+=1
        print ('',file=sequence)  
    sequence.close()
    print ('done')

if __name__ == "__main__":
    if len(sys.argv)==1:
        print (Usage%{'prog':sys.argv[0]})
    else:
        bamfile,vcffile,outdir,strainsNum,extractHAIRS,snp_dp=args.bamfile,args.vcffile,args.outdir,args.k,args.extractHAIRS,args.snp_dp
        qual,weight,lambda1,lambda2,prior,elbow=args.qual,args.weight,args.lambda1,args.lambda2,args.prior,args.elbow
        with open(prior, 'rb') as f:
            popu = pickle.load(f)
        f.close()
        print ('Database is loaded.')
        data,alt_homo,to_sort,hete_species=pipeline.read_vcf(outdir,[],snp_dp,qual,vcffile)
        data=pipeline.delta(outdir,extractHAIRS,data,bamfile)
        sp=data[0]
        species=sp[0]
        snp_list=sp[1]
        dps_set=[]
        for dps in sp[2]:
            hat=dps[1]/sum(dps)
            dps_set.append(hat)
        beta_set,delta_set,share_set=pipeline.rectify(sp,lambda1,lambda2,popu)
        for i in range(len(beta_set)):
            print (dps_set[i],beta_set[i])
'''
        if strainsNum==0:
            wo=Workflow(beta_set,delta_set,share_set,weight,elbow)
            final_alpha,seq_list=wo.choose_k()
        else:
            wo=Workflow(beta_set,delta_set,share_set,weight,elbow)
            final_alpha,seq_list=wo.given_k(strainsNum)   
        output(outdir,final_alpha,seq_list,to_sort,snp_list)    
'''
