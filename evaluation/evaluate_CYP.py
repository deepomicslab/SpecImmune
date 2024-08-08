from four_field_compare import main_cyp_hprc

if __name__ == "__main__":
    #### CYP
    # pangu_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/pangu_hprc/"
    # spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_hprc/"

    # pangu_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/pangu_hprc_ont/"
    # spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_hprc_ont/"
    
    # read_pangu_result("/mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/amplicon2/NA17246-SRR15476220/NA17246-SRR15476220_report.json")
    # read_spec_result("/mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/amplicon2/NA17246-SRR15476220/NA17246-SRR15476220.CYP.merge.type.result.txt")
    # load_HPRC_CYP_truth()
    # print (validate_star_allele('*41x2', '*41x3'))
    pangu_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/pangu_1k/"
    spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_1k3/"
    main_cyp_hprc(pangu_dir, spec_dir)