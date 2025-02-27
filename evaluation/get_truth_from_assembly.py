import os, re

def get_phased_assemblies():
    record_truth_file_dict = {}
    
    for file in os.listdir(inpath):
        if file[-5:] != "fasta":
            continue
        sample = file.split("_")[1]
        if sample not in record_truth_file_dict:
            record_truth_file_dict[sample] = ['', '']
        full_file = inpath + file
        if re.search(".h1-", full_file):
            record_truth_file_dict[sample][0] = full_file
        else:
            record_truth_file_dict[sample][1] = full_file
    # print (record_truth_file_dict)
    return  record_truth_file_dict

def get_phased_assemblies_hprc():
    record_truth_file_dict = {}
    
    for file in os.listdir(inpath):
        if file[-3:] != ".fa":
            continue
        # print (file)
        sample = file.split("_")[0]
        # print (sample)
        if sample not in record_truth_file_dict:
            record_truth_file_dict[sample] = ['', '']
        full_file = inpath + file
        if re.search("hap1", full_file):
            record_truth_file_dict[sample][0] = full_file
        else:
            record_truth_file_dict[sample][1] = full_file
    # print (record_truth_file_dict)
    return  record_truth_file_dict

# inpath = "/mnt/d/my_HLA/assembly/"
# outdir = "/mnt/d/HLAPro_backup/Nanopore_optimize/hgscv2_truth_bwa/"
# record_truth_file_dict = get_phased_assemblies()


inpath = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/HPRC_assembly/"
outdir = "/mnt/d/HLAPro_backup/Nanopore_optimize/hprc_truth_bwa/"
record_truth_file_dict = get_phased_assemblies_hprc()

for sample in record_truth_file_dict:
    # if sample != "HG00096":
    #     continue
    # print ("xx", record_truth_file_dict[sample])
    ## check if the bwa index exsits for record_truth_file_dict[sample][0], if not. index
    # for i in range(2):
    #     prefix = record_truth_file_dict[sample][i]
    #     bwt = prefix + ".bwt"
    #     if not os.path.exists(bwt):
    #         print (f"indexing {record_truth_file_dict[sample][i]}, {bwt} does not exist")
    #         cmd = f"""bwa index {record_truth_file_dict[sample][i]}"""
    #         os.system(cmd)
    #     else:
    #         print (f"{record_truth_file_dict[sample][i]} already indexed")
    # for gene_class in ["HLA", "CYP", "KIR", "IG_TR"]:
    for gene_class in ["HLA"]:
        cmd = f"""
        python3 ../scripts/typing_from_assembly.py -j 15 -1 {record_truth_file_dict[sample][0]} -2 {record_truth_file_dict[sample][1]} -n {sample} -i {gene_class} -o {outdir}  
        """
        print (cmd)
        os.system(cmd)
    # break

    # cmd = f"""
    # python3 ../scripts/anno.IG.TR.py {sample} {record_truth_file_dict[sample][0]} {record_truth_file_dict[sample][1]}  {outdir} ../db/IG_TR/IG_TR.fasta 15
    # """
    # print (cmd)
    # os.system(cmd)
    # break
