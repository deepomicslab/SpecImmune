#!/usr/bin/env python3

from my_imports import *

def read_cfg(cfgfile):
    cfg_list=[]
    for line in open(cfgfile,'r'):
        line=line.strip()
        if line != '':
            cfg_list.append(line)
    sample_num=int(len(cfg_list)/4)
    return cfg_list, sample_num

def makedir(outdir, similarity):
    mergedir=outdir+'/merge_%s/'%(similarity)  
    if not os.path.exists(mergedir):
        os.system('mkdir '+ mergedir)
    if not os.path.exists(mergedir+'/seq/'):
        os.system('mkdir '+ mergedir+'/seq/')
    return mergedir

def read_rafile(sample_num, outdir, cfg_list):
    species_dict={}
    population_species=[]
    for i in range(sample_num):
        sample_species={}
        sample_name=cfg_list[4*i+1].split(':')[1].strip()
        RAdir=outdir+'/%s/'%(sample_name)+'/result/strain_RA.txt'
        if os.path.isfile(RAdir):  
            for line in open(RAdir,'r'):
                line=line.strip()
                if line[0] != '#':
                    array=line.split()   
                    if array[0] not in sample_species.keys():
                        sample_species[array[0]] = 1
                    else:  
                        sample_species[array[0]] += 1  
                    if array[0] not in population_species:
                        population_species.append(array[0])  
        else:
            print ('cannot find %s'%(RAdir))  
        species_dict[sample_name] = sample_species
    return species_dict, population_species

def output_strain_number_population(mergedir, population_species, species_dict):
    strain_number=open(mergedir+'strain_number.txt','w')

    print ('#samples',end = '\t',file=strain_number)
    for sp in population_species:
        print (sp,end='\t',file=strain_number)
    print ('',file=strain_number)
    for sample in species_dict.keys():
        print (sample,end='\t',file=strain_number)
        single_sample=species_dict[sample]
        for sp in population_species:            
            if sp in single_sample.keys():
                print (single_sample[sp],end='\t',file=strain_number)
            else:
                print (0,end='\t',file=strain_number)
        print ('',file=strain_number)
    strain_number.close()

def read_species_seq_population(species_dict, outdir, sp):
    locus_dict = {}
    allele_dict={}
    num_dict = {}
    name_dict={}
    for sample in species_dict.keys():
        sample_dict={}
        seq_file='%s/%s/result/seq/%s_seq.txt'%(outdir,sample,sp)
        if not os.path.isfile(seq_file): 
            # print ('sample: %s doesnot have %s!'%(sample, sp))
            continue
        for line in open(seq_file,'r'):
            line=line.strip()
            if line[0] != '#':
                array=line.split()
                tag='%s,%s'%(array[0],array[1])
                sample_dict[tag]=array[4:]
                if tag in locus_dict.keys():
                    if array[2] not in locus_dict[tag]:
                        locus_dict[tag].append(array[2])
                    elif array[3] not in locus_dict[tag]:
                        locus_dict[tag].append(array[3])
                else:
                    allele=[array[2],array[3]]
                    locus_dict[tag]=allele
            elif line[0] == '#':
                array=line.split()
                num_dict[sample]=len(array)-5
                strain_name=array[5:]
                name_dict[sample]=strain_name
        allele_dict[sample]=sample_dict
    # print (num_dict)
    return locus_dict, allele_dict, num_dict, name_dict

def generate_all_points(locus_dict, allele_dict, num_dict):
    all_points=[]
    for lo in locus_dict.keys():
        if len(locus_dict[lo]) != 2:
            continue
        point_list=[]
        for sample in allele_dict.keys():
            str_num=0
            if sample in num_dict.keys():
                str_num=num_dict[sample]

            if lo in allele_dict[sample].keys():
                if len(allele_dict[sample][lo]) != int(str_num):
                    print (sample, int(str_num), len(allele_dict[sample][lo]), 'warning')
                point_list+=allele_dict[sample][lo]
            else:
                point_list+=int(str_num)*[0]
        point_list=list(map(int,point_list))
        all_points.append(point_list)
    return all_points

def species_samples(outdir,cfgfile,similarity):
    cfg_list, sample_num = read_cfg(cfgfile)
    if sample_num < 2:
        print ('No need to merge in case of one sample.')
        return 0
    mergedir = makedir(outdir, similarity)
    species_dict, population_species = read_rafile(sample_num, outdir, cfg_list)
    output_strain_number_population(mergedir, population_species, species_dict)
    #merge points
    final_dict={}
    for sample in species_dict.keys():
        final_dict[sample]={}
    for sp in population_species:
        locus_dict, allele_dict, num_dict, name_dict = read_species_seq_population(species_dict, outdir, sp)
        all_points = generate_all_points(locus_dict, allele_dict, num_dict)
 
        # print ('merge samples', sp)
        # if sp == 'Streptococcus_salivarius':
        #     print (len(locus_dict), len(all_points))
        cluster_index,maj_seq=cluster(all_points,similarity)
        #output all seq 
        sample_seq=open(mergedir+'/seq/'+'%s_cluster_%s_seq.txt'%(sp, similarity),'w')
        clu_seq=open(mergedir+'/seq/'+'%s_clu.txt'%(sp),'w')
        print ('# Gene\tLocus\tRef\tAlt\t',end='',file=sample_seq)    
        print ('# Gene\tLocus\tRef\tAlt\t',end='',file=clu_seq)        
        for sample in allele_dict.keys():
            if sample in num_dict.keys():
                str_num=num_dict[sample]
                str_name=name_dict[sample]
                for i in range(str_num):
                    column='%s/%s'%(sample,str_name[i])
                    print (column,end='\t',file=sample_seq)
        print ('',file=sample_seq)
        # for the species without snvs
        if len(maj_seq) == 0:
            maj_seq=[[0]]
            sample_numer=len(list(num_dict.keys()))
            cluster_index=[1]*sample_numer
        for i in range(len(maj_seq[0])):
            clu_name='clu-%s'%(i+1)
            print (clu_name,end='\t',file=clu_seq)
        print ('',file=clu_seq)

        i=0
        for lo in locus_dict.keys():
            if len(locus_dict[lo]) == 2:
                array=lo.split(',')
                print ('%s\t%s\t%s\t%s'%(array[0],array[1],locus_dict[lo][0],locus_dict[lo][1]),end='\t',file=sample_seq)
                print ('%s\t%s\t%s\t%s'%(array[0],array[1],locus_dict[lo][0],locus_dict[lo][1]),end='\t',file=clu_seq)
                for p in all_points[i]:
                    print (p,end='\t',file=sample_seq)
                # print (lo,all_points[i])
                
                print ('',file=sample_seq)
                #cluster output
                for maj in maj_seq[i]:
                    print (maj,end='\t',file=clu_seq)
                print ('',file=clu_seq)

                i+=1
        sample_seq.close()
        clu_seq.close()
        
        start=0
        for sample in species_dict.keys():
            if sample in num_dict.keys():
                str_num=num_dict[sample]
                sample_index=cluster_index[start:start+str_num]
                start=start+str_num
                
                final_dict[sample][sp]=sample_index
        # break
    for sample in final_dict.keys():
        RAdir=outdir+'/%s/'%(sample)+'/result/strain_RA.txt'
        if os.path.isfile(RAdir):
            clusterdir=outdir+'/%s/'%(sample)+'/result/strain_merged_RA_%s.txt'%(similarity)
            cout=open(clusterdir,'w')
            for line in open(RAdir,'r'):
                line=line.strip()
                if line[0] == '#':
                    print ('# Species\tSpecies_RA\tStrain_ID\tCluster_ID\tStrain_Freq\tStrain_RA',file=cout)
                else:
                    array=line.split()
                    # print (array)
                    if array[0] in final_dict[sample].keys():# and :
                        # print (final_dict[sample][array[0]])
                        index=final_dict[sample][array[0]][0]
                        cluster_ID='clu-%s'%(index)
                        print (array[0],array[1],array[2],cluster_ID,array[3],array[4],file=cout)
                        final_dict[sample][array[0]]=np.delete(final_dict[sample][array[0]],0)
                    else:
                        print (array[0],array[1],array[2],'NA',array[3],array[4],file=cout)
        else:
            print ('no abundance file!')
            cout.close()
    print ('merge done')

def cluster(all_points,similarity):
    # similarity=0.95
    cutoff=1-similarity
    points_number=len(all_points)
    all_points=np.array(all_points)
    seq_strs=all_points.transpose()
    str_num=len(seq_strs)
    if str_num>1:
        y=pdist(seq_strs,'cityblock')
        Z = average(y)
        cluster_index=fcluster(Z, int(cutoff*points_number), criterion='distance')        
    else:
        cluster_index=np.array([1])
    clu_num=max(cluster_index)
    # print (clu_num,cluster_index)
    maj_seq=[]
    for points in all_points:
        major_list=np.zeros((clu_num,2))
        for i in range(str_num):
            str_allele=points[i]
            major_list[cluster_index[i]-1][str_allele]+=1
        major_allele=[]
        for clu_list in major_list:
            if clu_list[0] >= clu_list[1]:
                major_allele.append(0)
            else:
                major_allele.append(1)
        maj_seq.append(major_allele)
        # print (major_allele,clu_num)
    return cluster_index,maj_seq

Usage = \
"""%(prog)s [options] -c/--conf <config file> -o/--outdir <output directory>

The config file should follow the format:

//
sample: [sample1_ID]
fq1: [forward reads fastq]
fq2: [reverse/mate reads fastq]
//
sample: [sample2_ID]
fq1: [forward reads fastq]
fq2: [reverse/mate reads fastq]
...

Help information can be found by %(prog)s -h/--help, config file format for single end reads and additional information can be found in README.MD or www.xxx.com.
"""
def input(Usage):

    scripts_dir=sys.path[0]+'/'
    parser = ArgumentParser(description="Merge the results from all of the samples.",prog='python3 merge.py',usage=Usage)
    optional=parser._action_groups.pop()
    required=parser.add_argument_group('required arguments')    
    #necessary parameter
    required.add_argument("-c", "--conf",help="The configuration file of the input samples.",dest='cfgfile',metavar='')
    required.add_argument("-o", "--outdir",help="The output directory.",dest='outdir',metavar='')
    optional.add_argument("--similarity",help="The similarity cutoff of hierachical clustering in merge step (default is 0.95).",\
        dest='similarity',metavar='',default=0.95, type=float)
    parser._action_groups.append(optional)
    args = parser.parse_args()
    return args.cfgfile, args.outdir, args.similarity


if __name__ == "__main__":
    if len(sys.argv)==1:
        print (Usage%{'prog':sys.argv[0]})
    else:
        # species_samples(sys.argv[1],sys.argv[2],0.95)
        cfgfile, outdir, similarity = input(Usage)
        species_samples(outdir,cfgfile,similarity)
