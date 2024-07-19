#!/usr/bin/env python3

from my_imports import *
import multiprocessing
import logging
from reads_info import reads_info


def run_metaphlan2(metaphlan2dir,metaphlan2,fq1,fq2,nproc,bowtie2):
    print ("run metaphlan...")
    metaphlan2_file=metaphlan2dir+'/metaphlan2_output.txt'
    if not os.path.isfile(metaphlan2_file) or os.stat(metaphlan2_file).st_size == 0:
        if fq2=='single_end':
            metaphlan2_order='python3 %s %s --input_type fastq --bowtie2_exe %s --nproc %s --bowtie2out %s/bowtie.out.bz2 >%s/metaphlan2_output.txt'%(metaphlan2,fq1,bowtie2,nproc,metaphlan2dir,metaphlan2dir)
        else:
            metaphlan2_order='python3 %s %s,%s --input_type fastq --bowtie2_exe %s --nproc %s --bowtie2out %s/bowtie.out.bz2 >%s/metaphlan2_output.txt'%(metaphlan2,fq1,fq2,bowtie2,nproc,metaphlan2dir,metaphlan2dir)
        # metaphlan2_order='%s %s,%s --input_type fastq --bowtie2_exe %s --nproc %s --bowtie2out %s/bowtie.out.bz2 >%s/metaphlan2_output.txt'%(metaphlan2,fq1,fq2,bowtie2,nproc,metaphlan2dir,metaphlan2dir)
        print (metaphlan2_order)
        os.system(metaphlan2_order)
        os.system('rm %s/bowtie.out.bz2'%(metaphlan2dir))
        logging.info('Metaphlan2 is done.') 
    else:
        logging.info('Metaphlan2 result exists already.')
def read_metaphlan2(metaphlan2dir):
    metaphlan2_file=metaphlan2dir+'/metaphlan2_output.txt'
    species_set=[]
    sp_ra={}
    for line in open(metaphlan2_file,'r'):
        line=line.strip()
        name=re.search('s__(.+)$',line)
        if name and not re.search('t__',line) and not re.search('unclassified',line):
            array=name.group(1).split()
            if float(array[1]) > 0.05:
                species_set.append(array[0])
                #sp_ra.append(float(array[1]))
                sp_ra[array[0]]=float(array[1]) 
    return species_set,sp_ra
def extract_ref(species_set,refdir,dbdir):
    # dbdir='../db_v20/'
    gene_dict={}
    # species_set=['s__Escherichia_coli']
    for line in gzip.open(dbdir+'/species_markers.txt.gz','rt'):
        line=line.strip()
        array=line.split()
        for sp in species_set:
            if array[1]=='s__'+sp:
                gene_dict['>'+array[0]]=sp
    flag=False
    merged_fh=open(refdir+'/merged_ref.fa','w')
    for line in gzip.open(dbdir+'/marker_gene.fna.gz','rt'):
        line=line.strip()       
        if line[0] == '>':
            if line in gene_dict.keys() and not re.search(',',line):
                line=line+'|'+gene_dict[line]
                flag=True
            else:
                flag=False
        if flag == True:
            print (line,file=merged_fh)
def index_ref(refdir,picard,samtools,bowtie2_build):
    # picard_index='java -jar %s CreateSequenceDictionary REFERENCE=%s/merged_ref.fa OUTPUT=%s/merged_ref.dict'%(picard,refdir,refdir)
    dict_index='%s dict %s/merged_ref.fa > %s/merged_ref.dict'%(samtools,refdir,refdir)
    samtools_index='%s faidx %s/merged_ref.fa'%(samtools,refdir)
    bowtie2_index='%s -f %s/merged_ref.fa %s/merged_ref'%(bowtie2_build,refdir,refdir)
    index_order=dict_index+'\n'+samtools_index+'\n'+bowtie2_index+'\n'
    os.system(index_order)
    logging.info('Sample specific reference is finished.')
def bowtie2_map(bowtie2,samtools,refdir,mapdir,fq1,fq2,nproc,picard):
    if fq2=='single_end':
        bowtie2_order='%s -x %s/merged_ref -p %s -U %s |grep -v "XS:i:"|%s view -bS -F 4 |%s sort -o %s/mapped.sort.bam'%(bowtie2,refdir,nproc,fq1,samtools,samtools,mapdir)
    else:
        bowtie2_order='%s -x %s/merged_ref -p %s -1 %s -2 %s |grep -v "XS:i:"|%s view -bS -F 4 |%s sort -o %s/mapped.sort.bam'%(bowtie2,refdir,nproc,fq1,fq2,samtools,samtools,mapdir)
    addheader='java -jar %s AddOrReplaceReadGroups I=%s/mapped.sort.bam O=%s/mapped.bam LB=whatever PL=illumina PU=whatever SM=whatever'%(picard,mapdir,mapdir)
    os.system(bowtie2_order)
    os.system(addheader)
    os.system('rm %s/mapped.sort.bam'%(mapdir))
    os.system('%s index %s/mapped.bam'%(samtools,mapdir))
    os.system('%s depth %s/mapped.bam >%s/mapped.depth'%(samtools,mapdir,mapdir))
    logging.info('Sample reads mapped to its specific reference.')
def call_snp(gatk,mapdir,refdir,bamfile):
    gatk_order='java -Xmx5g -jar %s \
    -T HaplotypeCaller -R %s/merged_ref.fa -allowPotentiallyMisencodedQuals \
    -I %s -o %s/mapped.vcf.gz'%(gatk,refdir,bamfile,mapdir)
    os.system(gatk_order)
    logging.info('SNPs calling is done.')

    # freebayes_order = 'freebayes -f -R %s/merged_ref.fa -b %s -v %s/mapped.vcf\n \
    # gzip %s/mapped.vcf'%(refdir,bamfile,mapdir,mapdir)
    # os.system(freebayes_order)
def read_vcf(mapdir,removed_gene,snp_dp,qual,vcffile,dict_species_depth):
    # vcffile='%s/mapped.vcf.gz'%(mapdir)
    if not os.path.exists(mapdir):
        os.system('mkdir '+mapdir)
    hete_species=[]
    snp_list,beta_set=[],[]
    locus_list=[]
    alt_homo=[]
    myvcf=VariantFile(vcffile)
    vcf_out=VariantFile('%s/mapped.filter.vcf'%(mapdir),'w',header=myvcf.header)
    sample=list(myvcf.header.samples)[0]
    data=[]
    species=''
    to_sort=[]
    flag=False
    for record in myvcf.fetch():
        if record.info['DP'] <1:
            continue
        geno=record.samples[sample]['GT']        
        depth=record.samples[sample]['AD']
        dp=depth[0]+depth[1]   
        my_snp_dp = snp_dp * dict_species_depth[record.chrom.split('|')[-1]]
        # print (my_snp_dp)
        if len(record.ref)==1 and len(record.alts)==1 and len(record.alts[0])==1 and dp>=my_snp_dp  and record.qual >= qual and record.chrom not in removed_gene:
        # if len(record.ref)==1 and len(record.alts)==1 and len(record.alts[0])==1 and dp>=snp_dp  and record.qual >= qual and record.chrom not in removed_gene:
            snp=[record.chrom,str(record.pos),record.ref,record.alts[0]]
            if geno==(0,1):
                flag=True
                to_sort.append(snp)
                vcf_out.write(record)
                array=record.chrom.split('|')
                sp=array[-1]
                if sp != species:
                    hete_species.append(species)
                    data.append([species,snp_list,beta_set])
                    species=sp
                    snp_list=[]
                    beta_set=[]
                beta=[depth[0],depth[1]]  # depth of 0 and 1
                # beta=np.array(depth)
                snp_list.append(snp)
                beta_set.append(beta)
                locus_list.append(int(record.pos))
            elif geno==(1,1):
                to_sort.append(snp)
                alt_homo.append(snp)
    hete_species.append(species)
    data.append([species,snp_list,beta_set])
    data=data[1:]
    hete_species=hete_species[1:]
    logging.info('SNPs is filtered.')
    if flag == False:
        logging.warning('There is no heterozygous loci at all.')
    return data,alt_homo,to_sort,hete_species
def depth_of_each_species(mapdir):
    dict_species_depth = {}
    for line in open('%s/mapped.depth'%(mapdir)):
        line=line.strip()
        array=line.split()   
        species = array[0].split('|')[-1]
        if species not in dict_species_depth.keys():
            dict_species_depth[species] = [float(array[2])]
        else:
            dict_species_depth[species].append(float(array[2]))
    for species in dict_species_depth.keys():
        dict_species_depth[species] = np.mean(np.array(dict_species_depth[species]))
    return (dict_species_depth)
def copy_number(mapdir,species_dp,samtools):
    gene_depth=[]
    chrom=''
    dp_set=[0]
    depth_file = '%s/mapped.depth'%(mapdir)
    if not os.path.isfile(depth_file):
        os.system('%s depth %s/mapped.bam >%s/mapped.depth'%(samtools,mapdir,mapdir))
    for line in open(depth_file, 'r'):
        line=line.strip()
        array=line.split()
        if array[0] != chrom:
            gene_depth.append([chrom,np.mean(dp_set),dp_set])
            chrom=array[0]
            dp_set=[]
        dp_set.append(int(array[2]))
    gene_depth.append([chrom,np.mean(dp_set),dp_set])
    gene_depth=gene_depth[1:]
    dict={}
    for gene in gene_depth:
        array=gene[0].split('|')
        sp=array[-1]
        if sp in dict.keys():
            dict[sp].append(gene[1])
        else:
            dict[sp]=[gene[1]]
    sp_mean={}
    low_dp_sp=[]
    # print (dict.keys())
    for sp in dict.keys():
        sp_mean[sp]=[round(np.mean(dict[sp]),6),round(np.std(dict[sp]),6)]
        if sp_mean[sp][0] < species_dp:
            low_dp_sp.append(sp)

    removed_gene=[]
    num=3
    for gene in gene_depth:
        array=gene[0].split('|')
        sp=array[-1]  
        mean_std= sp_mean[sp]
        # if mean_std[0] < species_dp:
        #     removed_gene.append(gene[0])
        if gene[1]<mean_std[0]-num*mean_std[1] or gene[1] > mean_std[0]+num*mean_std[1]:#or len(gene[2])<100 or max(gene[2])>200 :
            removed_gene.append(gene[0])
    return removed_gene,low_dp_sp
def delta(mapdir,data,bamfile):
    
    for i in range(len(data)):
        # share_reads=[]
        sp = data[i]
        beta_set = sp[2]
        snp_list = sp[1]
        delta_set,share_set=reads_info(bamfile,snp_list,0)
        data[i].append(delta_set)  
        data[i].append(share_set)  
        # print (delta_set)
    return data
    '''
    hapcut_order='%s --bam %s --VCF %s/mapped.filter.vcf --out %s/mapped.conn'%(extractHAIRS,bamfile,mapdir,mapdir)
    os.system(hapcut_order)
    snp_num=0
    for sp in data:
        snp_num+=len(sp[-1])
    delta_set=[]
    for i in range(snp_num-1):
        delta_set.append([[0,0],[0,0]])
    for line in open('%s/mapped.conn'%(mapdir)):
        line=line.strip()
        array=line.split()
        if array[0] == '1': 
            delta_index=int(array[2])-1
            geno_type=array[3]
            for i in range(len(geno_type)-1):
                delta_set[delta_index+i][int(geno_type[i])][int(geno_type[i+1])]+=1
    base=0
    for i in range(len(data)):
        share_reads=[]
        sp=data[i]
        beta_set=sp[2]
        sp_snp_num=len(sp[-1])
        start=base
        end=base+sp_snp_num-1
        base=end+1
        for j in range(start,end):
            delta=np.array(delta_set[j])
            sum_dp=sum(delta[0])+sum(delta[1])
            share_reads.append(sum_dp)
            if sum_dp >0:
                delta=delta/sum_dp
            else:
                delta=np.array([[0.25,0.25],[0.25,0.25]])
            delta_set[j]=delta.tolist()
        # print (delta_set[start:end])
        data[i].append(delta_set[start:end])  # hat beta^2
        data[i].append(share_reads)  #reads number that cover the two loci
    logging.info('The input for core algorithm is done.')
    return data
    '''  
def prior_locus(popu,nucleotide,point):
    name=str(point[0])+'_'+str(point[1])
    if name in popu.keys():
        allele_freq=popu[name]
        ref_freq=float(allele_freq[nucleotide[point[2]]])
        alt_freq=float(allele_freq[nucleotide[point[3]]])
        if ref_freq+alt_freq > 0:
            prior_beta=round(alt_freq/(ref_freq+alt_freq),6)
        else:
            prior_beta=0.5
    else:
        prior_beta=0.5
    return prior_beta
def rectify(sp,lambda1,lambda2,popu):
    reads_set=sp[2]  # reads number of 0 and 1 in locus j  
    delta_hat=sp[3]  # init beta^2
    share_set=sp[4]  # share reads number in loci j and j+1 
    snp_list=sp[1]   #snp keys

    nucleotide={'A':0,'T':1,'C':2,'G':3}

    # rectify beta
    beta_set=[]
    for i in range(len(reads_set)):
        prior_f=prior_locus(popu,nucleotide,snp_list[i])
        dp=int(reads_set[i][0]) + int(reads_set[i][1])
        hat_beta=round(float(reads_set[i][1])/dp,6)
        if lambda1 > 1+dp:
            lambda1 = 1+dp
        beta=hat_beta*(1-lambda1/(1+dp))+prior_f*(lambda1/(1+dp))
        # print ('psudocount',hat_beta,prior_f,beta)
        beta_set.append(beta)
    # rectify beta^2
    delta_set=[]
    for i in range(len(share_set)):       
        c=int(share_set[i])
        hat_delta=delta_hat[i]
        inde_delta=[[(1-beta_set[i])*(1-beta_set[i+1]),(1-beta_set[i])*(beta_set[i+1])],[(beta_set[i])*(1-beta_set[i+1]),(beta_set[i])*(beta_set[i+1])]]
        hat_delta=np.array(hat_delta)
        inde_delta=np.array(inde_delta)
        if lambda2 > 1+c:
            lambda2 = 1+c
        # if c < 5:
        #     delta = np.array([[0, 0], [0, 0]])
        # else:
        #     delta = inde_delta*(lambda2/(1+c)) + hat_delta*(1-lambda2/(1+c))  #retify
        delta=inde_delta*(lambda2/(1+c)) + hat_delta*(1-lambda2/(1+c))  #retify
        delta_set.append(delta.tolist())
    return beta_set,delta_set,share_set
def profiling(data,weight,lambda1,lambda2,popu,elbow,low_dp_sp):
    species_alpha,species_seq,species_snp=[],[],[]
    for sp in data:
        species=sp[0]
        snp_list=sp[1]
        beta_set,delta_set,share_set=rectify(sp,lambda1,lambda2,popu)    
        # print (beta_set[:10],delta_set[:10])
        # if len(beta_set) == 0:
        #     #for the species that consists no hete locus.
        #     species_alpha.append([1.0])
        #     species_seq.append([])
        #     species_snp.append([])
        if species not in low_dp_sp:# and len(beta_set) > 0:
            wo=Workflow(beta_set,delta_set,share_set,weight,elbow)
            final_alpha,seq_list,final_loss=wo.choose_k()
            species_alpha.append(final_alpha)
            species_seq.append(seq_list)
            species_snp.append(snp_list)
        else:
            #for species with low abundance, we profile the consensus sequence.
            wo=Workflow(beta_set,delta_set,share_set,weight,elbow)
            final_alpha,seq_list,final_loss=wo.given_k(1)
            species_alpha.append(final_alpha)
            species_seq.append(seq_list)
            species_snp.append(snp_list)        
    print ('The core algorithm is done.')
    logging.info('The core algorithm is done.')
    return species_alpha,species_seq,species_snp
def single_run(sample,outdir,fq1,fq2,arg_list,popu):
    nproc=arg_list[0] 
    species_dp=arg_list[1]
    snp_dp=arg_list[2]
    dbdir=arg_list[3]
    picard=arg_list[4]
    samtools=arg_list[5]
    bowtie2_build=arg_list[6]
    bowtie2=arg_list[7]
    # extractHAIRS=arg_list[8]
    gatk=arg_list[9]
    metaphlan2=arg_list[10]
    weight=arg_list[12]
    lambda1=arg_list[13]
    lambda2=arg_list[14]
    prior=arg_list[15]
    elbow=arg_list[16]
    qual=arg_list[17]
    outdir=outdir+'/%s/'%(sample)
    refdir=outdir+'/ref/'
    metaphlan2dir=outdir+'/metaphlan2/'
    mapdir=outdir+'/map'
    resultdir=outdir+'/result/'
    if not os.path.exists(outdir):
        os.system('mkdir '+outdir)
    if not os.path.exists(refdir):
        os.system('mkdir '+refdir)
    if not os.path.exists(metaphlan2dir):
        os.system('mkdir '+metaphlan2dir)
    if not os.path.exists(mapdir):    
        os.system('mkdir '+mapdir)
    if not os.path.exists(resultdir):
        os.system('mkdir '+resultdir)
    if not os.path.exists(resultdir+'/seq/'):
        os.system('mkdir '+resultdir+'/seq/')
    print ('##############')
    logging.basicConfig(level=logging.DEBUG,  
                        format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',  
                        datefmt='%a, %d %b %Y %H:%M:%S',  
                        filename=outdir+'/run.log',  
                        filemode='w')  
    vcffile,bamfile='%s/mapped.vcf.gz'%(mapdir),'%s/mapped.bam'%(mapdir)
    if os.path.isfile(vcffile) and os.path.isfile(bamfile):
        species_set,sp_ra=read_metaphlan2(metaphlan2dir)
        removed_gene,low_dp_sp=copy_number(mapdir,species_dp,samtools)
        dict_species_depth = depth_of_each_species(mapdir)
    else:
        run_metaphlan2(metaphlan2dir,metaphlan2,fq1,fq2,nproc,bowtie2)
        species_set,sp_ra=read_metaphlan2(metaphlan2dir)
        extract_ref(species_set,refdir,dbdir)
        index_ref(refdir,picard,samtools,bowtie2_build)
        bowtie2_map(bowtie2,samtools,refdir,mapdir,fq1,fq2,nproc,picard)    
        removed_gene,low_dp_sp=copy_number(mapdir,species_dp,samtools)
        dict_species_depth = depth_of_each_species(mapdir)
        call_snp(gatk,mapdir,refdir,bamfile) 
    print ('......................................')
    data,alt_homo,to_sort,hete_species=read_vcf(mapdir,removed_gene,snp_dp,qual,vcffile, dict_species_depth)
    print ('after reads vcf, start delta')
    data=delta(mapdir,data,bamfile)  #add delta_set
    print ('start profiling.')
    species_alpha,species_seq,species_snp=profiling(data,weight,lambda1,lambda2,popu,elbow,low_dp_sp)
    output(species_alpha,species_seq,species_snp,species_set,sp_ra,alt_homo,to_sort,resultdir,hete_species,low_dp_sp)
    print ('Sample %s is done.'%(sample))
    logging.info('Sample %s is done.'%(sample))
def output(species_alpha,species_seq,species_snp,species_set,sp_ra,alt_homo,to_sort,resultdir,hete_species,low_dp_sp):
    ra_file=open(resultdir+'/strain_RA.txt','w')
    print ('# Species\tSpecies_RA\tStrain_ID\tStrain_Freq\tStrain_RA',file=ra_file)
    for i in range(len(hete_species)):
        species=hete_species[i]
        #ra=sp_ra[i]
        ra=sp_ra[species]
        alpha=species_alpha[i]
        seq=species_seq[i]
        snp=species_snp[i]  
        if species not in low_dp_sp:  
            for j in range(len(alpha)):
                print (species,ra,'str-'+str(j+1),round(alpha[j],6),round(ra*alpha[j],5),file=ra_file)
        else:
            print (species,ra,'Consensus_Seq','1.0',ra,file=ra_file)
        sequence=open(resultdir+'/seq/'+species+'_seq.txt','w')
        print ('# Gene\tLocus\tRef\tAlt\t',end='',file=sequence)
        if species not in low_dp_sp: 
            for i in range(len(alpha)):
                print ('str-'+str(i+1),end='\t',file=sequence)
        else:
            print ('Consensus_Seq',end='\t',file=sequence)
        print ('',end='\n',file=sequence)
        he_num=len(snp)
        he=0
        ho=0       
        for point in to_sort:
            match=re.search('(.*)(\|%s)'%(species),point[0])
            if match:
                point[0]=match.group(1)
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
                point[0] = point[0] + '|' +  species
        sequence.close()
    print ('hete species is done')
    finish_species=hete_species[:]
    #the species that only have 1/1 snvs
    ####################################
    pre_species=''
    # print ('relative abundance', sp_ra, species_set, hete_species)
    for point in to_sort:
        # match=re.search('(.*)(\|%s)'%(species),point[0])
        match=re.search('(.*\|)(.*?)$', point[0])
        if match:
            species=match.group(2)   
            new_match=re.search('(.*)(\|%s)'%(species),point[0])  
            point[0]=new_match.group(1)
            if species in hete_species:
                continue
            if species != pre_species:
                finish_species.append(species)
                if species not in species_set:
                    print ('There is no species called %s!!'%(species))
                index=species_set.index(species)
                #ra=sp_ra[index]
                ra=sp_ra[species]
                # print (species, index, ra)
                pre_species=species
                # print (resultdir+'/seq/'+species+'_seq.txt')
                sequence=open(resultdir+'/seq/'+species+'_seq.txt','w')
                print ('# Gene\tLocus\tRef\tAlt\tstr-1',end='\n',file=sequence)
                print (species,ra,'str-1','1.0',ra,file=ra_file)
            for info in point:
                print (info,end='\t',file=sequence)
            print ('1',file=sequence)  
    print ('homo species is done')
    #the species that have no snv
    for index in range(len(species_set)):
        species=species_set[index]
        ra=sp_ra[species]
        if species not in finish_species:
            # print (species,'no valid SNV')
            sequence=open(resultdir+'/seq/'+species+'_seq.txt','w')
            print ('# Gene\tLocus\tRef\tAlt\tNo_Valid_SNV',end='',file=sequence)
            print (species,ra,'No_Valid_SNV','1.0',ra,file=ra_file)
            sequence.close()
    print ('output is done')
def multi_samples(outdir,cfgfile,arg_list,popu):

    if not os.path.exists(outdir):
        os.system('mkdir '+outdir)
    cfg_list=[]
    for line in open(cfgfile,'r'):
        line=line.strip()
        if line != '':
            cfg_list.append(line)
    sample_num=int(len(cfg_list)/4)
    for i in range(sample_num):
        sample_name=cfg_list[4*i+1].split(':')[1].strip()
        fq1=cfg_list[4*i+2].split(':')[1].strip()
        fq2_array=cfg_list[4*i+3].split(':')
        if len(fq2_array) == 2 and len(fq2_array[1].strip())>0:
            fq2=fq2_array[1].strip()
        else:
            fq2='single_end'
        single_run(sample_name,outdir,fq1,fq2,arg_list,popu)
        print ('Sample %s is done.'%(sample_name))
def multiproc(outdir,cfgfile,arg_list, args):
    prior=arg_list[15]

    if os.path.isfile(prior) and args.skip_prior:
        with open(prior, 'rb') as f:
            popu = pickle.load(f)
        f.close()
        # popu = {}
        print ('Database is loaded.')
        logging.info('Database is loaded.')
    else:
        popu = {}

    if not os.path.exists(outdir):
        os.system('mkdir '+outdir)
    cfg_list=[]
    for line in open(cfgfile,'r'):
        line=line.strip()
        if line != '':
            cfg_list.append(line)
    sample_num=int(len(cfg_list)/4)
    pool=multiprocessing.Pool(processes=arg_list[11])
    pool_list=[]    
    for i in range(sample_num):
        sample_name=cfg_list[4*i+1].split(':')[1].strip()
        # print (f"config {sample_name}")
        fq1=cfg_list[4*i+2].split(':')[1].strip()
        fq2_array=cfg_list[4*i+3].split(':')
        if len(fq2_array) == 2 and len(fq2_array[1].strip())>0:
            fq2=fq2_array[1].strip()
        else:
            fq2='single_end'
        # if not os.path.isfile('/mnt/disk2_workspace/wangshuai/00.strain/01.real_strains/real_data/CRC_new/%s/result/strain_RA.txt'%(sample_name)):
        pool_list.append(pool.apply_async(single_run,(sample_name,outdir,fq1,fq2,arg_list,popu,)))
    pool.close()
    pool.join()
        
# if __name__ == "__main__":
#     multi_samples(sys.argv[1],sys.argv[2])
