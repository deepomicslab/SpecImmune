import sys
sys.path.insert(0, sys.path[0]+'/../scripts/')

from four_field_compare import main_cyp_hprc, load_GeT_RM4, main_all_cyp
from four_field_compare import assess_sim_module
from determine_gene import get_focus_gene

def evaluate_sim():
    ## simumation for all genes
    gene_list, interval_dict =  get_focus_gene({'i':'CYP'})
    outdir="/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/reads/"
    resultdir="/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/results/"
    data = []
    dict = {'90':'90', '98':'95'}
    for prefix in ['90', '98']:
        total_accuracy_dict = {}
        total = [0, 0]
        for i in range(1, 51):
            sample = f"CYP_dp50_acc{prefix}_{i}"

            
            truth = f"{outdir}/{sample}/{sample}.CYP.hap.alleles.txt"
            # infer = f"{resultdir}/{sample}/{sample}.CYP.type.result.txt"
            infer = f"{resultdir}/{sample}/{sample}.CYP.final.type.1.result.txt"
            spec_gene_accuracy_dict = assess_sim_module(truth, infer, gene_list, 'CYP', mode="one_guess")
            print(spec_gene_accuracy_dict)
            for gene in spec_gene_accuracy_dict:
                if gene not in total_accuracy_dict:
                    total_accuracy_dict[gene] = [0, 0]
                total_accuracy_dict[gene][0] += spec_gene_accuracy_dict[gene][1]
                total_accuracy_dict[gene][1] += spec_gene_accuracy_dict[gene][2]
                total[0] += spec_gene_accuracy_dict[gene][1]
                total[1] += spec_gene_accuracy_dict[gene][2]
        
        for gene in total_accuracy_dict:
            accuracy = total_accuracy_dict[gene][0]/total_accuracy_dict[gene][1]
            print(f"{gene}: {accuracy}")
            data.append([gene, accuracy, total_accuracy_dict[gene][0], total_accuracy_dict[gene][1], dict[prefix]+"%"])
        print(f"Total: {total[0]/total[1]}", total)
        break
    ## convert to df, and save it in a csv file
    import pandas as pd
    df = pd.DataFrame(data, columns = ['Gene', 'Accuracy', 'Correct', 'Total', 'Seq_acc'])
    df.to_csv("cyp_results/sim_all_gene.csv", index=False)

def evaluate_real():
    real_truth = "cyp_results/cyp_all_gene.csv"

if __name__ == "__main__":
    #### CYP

    pangu_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/pangu_hprc/"
    spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_hprc/"
    result_file = "cyp_results/hprc_hifi_cyp_depth_cutoff.csv"
    main_cyp_hprc(pangu_dir, spec_dir, result_file, 'hprc')

    # pangu_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/pangu_hprc_ont/"
    # spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_hprc_ont2/"
    # result_file = "cyp_results/hprc_ont_cyp_depth_cutoff.csv"
    # main_cyp_hprc(pangu_dir, spec_dir, result_file, 'hprc')

    # pangu_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/pangu_1k/"
    # spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_1k3/"
    # result_file = "cyp_results/cyp_depth_cutoff.csv"
    # main_cyp_hprc(pangu_dir, spec_dir, result_file)

    # spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_1k_all2/"
    # result_file = "cyp_results/cyp_all_gene2.csv"
    # main_all_cyp(spec_dir, result_file)

    # load_GeT_RM4()

    evaluate_sim()


