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

inpath = "/mnt/d/my_HLA/assembly/"
outdir = "/mnt/d/HLAPro_backup/Nanopore_optimize/pacbio_truth/"

record_truth_file_dict = get_phased_assemblies()
for sample in record_truth_file_dict:
    if sample == "HG00514":
        continue
    print (record_truth_file_dict[sample])
    for gene_class in ["HLA", "CYP", "KIR"]:
        cmd = f"""
        python3 ../scripts/typing_from_assembly.py -1 {record_truth_file_dict[sample][0]} -2 {record_truth_file_dict[sample][1]} -n {sample} -i {gene_class} -o {outdir} -j 15 --db /mnt/d/HLAPro_backup/Nanopore_optimize/SpecComplex/db/
        """
        print (cmd)
        os.system(cmd)
    # break
