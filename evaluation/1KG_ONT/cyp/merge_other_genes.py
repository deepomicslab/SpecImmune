import os
import pandas as pd

def main(spec_dir, dataset, result):
    merged_df = pd.DataFrame()
    data = []
    for folder in os.listdir(spec_dir):
        ## check if the folder is a folder
        if not os.path.isdir(os.path.join(spec_dir, folder)):
            continue
        sample = folder
        spec_result = os.path.join(spec_dir, folder, f"{folder}.CYP.merge.type.result.txt")
        if not os.path.exists(spec_result):
            spec_result = os.path.join(spec_dir, folder, folder, f"{folder}.CYP.merge.type.result.txt")
            print (spec_result)
            # continue
        if not os.path.exists(spec_result):
            print (f"{spec_result} does not exist")
            continue
        ## read the spec result using pandas, skipping lines starting with '#'
        spec_result_df = pd.read_csv(spec_result, sep="\t", comment='#')
        ## add the sample name to the dataframe and dataset
        spec_result_df['Sample'] = sample
        spec_result_df['Dataset'] = dataset
        ## add spec_result_df to the merged dataframe
        merged_df = pd.concat([merged_df, spec_result_df], ignore_index=True)

    merged_df.to_csv(result, index=False)
    print (f"Saved the merged result to {result}")

if __name__ == "__main__":

    spec_dir = "/home/shuaiw/methylation/data/hla/CYP_rerun/out_1kg_CYP_ont/"
    result = "all_results/test"
    dataset = "1KGP"
    main(spec_dir, dataset, result)

