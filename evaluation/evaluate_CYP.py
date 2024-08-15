from four_field_compare import main_cyp_hprc, load_GeT_RM4, main_all_cyp



if __name__ == "__main__":
    #### CYP

    #pangu_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/pangu_hprc/"
    #spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_hprc/"
    #result_file = "cyp_results/hprc_hifi_cyp_depth_cutoff.csv"
    #main_cyp_hprc(pangu_dir, spec_dir, result_file, 'hprc')

    # pangu_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/pangu_hprc_ont/"
    # spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_hprc_ont2/"
    # result_file = "cyp_results/hprc_ont_cyp_depth_cutoff.csv"
    # main_cyp_hprc(pangu_dir, spec_dir, result_file, 'hprc')

    # pangu_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/pangu_1k/"
    # spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_1k3/"
    # result_file = "cyp_results/cyp_depth_cutoff.csv"
    # main_cyp_hprc(pangu_dir, spec_dir, result_file)

    spec_dir = "/home/wangshuai/00.hla/long/experiments/cyp/cyp_results/spec_1k3/"
    result_file = "cyp_results/cyp_all_gene.csv"
    main_all_cyp(spec_dir, result_file)

    #load_GeT_RM4()
