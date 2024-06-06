import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Define output directory and error list
outdir = "/scratch/project/cs_shuaicli/wxd/app/SpecLong/simulation/batch_sim_pacbio_CLR"
err_list = [round(x * 0.01, 2) for x in range(78, 80)]

# Function to parse accuracy from log files
def parse_acc(file): 
    res = {}
    try:
        with open(file, 'r') as f:
            for line in f:
                gene, acc = line.strip().split(""")
                res[gene] = float(acc)
    except FileNotFoundError:
        print(f"File not found: {file}")
    return res

# Initialize dictionaries to store accuracies
hlala_dict = {}
step1_dict = {}
step2_dict = {}

# Populate dictionaries with accuracy data
for err in err_list:
    formatted_err = f"{err:.2f}"
    for depth in range(5, 55, 5):
        sample = f"pacbio_dp{depth}_acc{err}"
        step1_log = os.path.join(outdir, sample, "acc.step1.log")
        step2_log = os.path.join(outdir, sample, "acc.step2.log")
        hlala_log = os.path.join(outdir, sample, "acc.hlala.log")
        print(formatted_err, depth)
        
        step1_dict[(formatted_err, depth)] = parse_acc(step1_log)
        step2_dict[(formatted_err, depth)] = parse_acc(step2_log)
        hlala_dict[(formatted_err, depth)] = parse_acc(hlala_log)

# Convert dictionaries to dataframes for plotting
def dict_to_df(data_dict, step_name):
    records = []
    for (err, depth), accs in data_dict.items():
        for gene, acc in accs.items():
            records.append((err, depth, gene, acc))
    return pd.DataFrame(records, columns=["Error Rate", "Depth", "Gene", f"{step_name} Accuracy"])

step1_df = dict_to_df(step1_dict, "Step 1")
step2_df = dict_to_df(step2_dict, "Step 2")
hlala_df = dict_to_df(hlala_dict, "Hlala")

# Melt dataframes for seaborn
def melt_df(df, step_name):
    return df.melt(id_vars=["Error Rate", "Depth", "Gene"], value_name=f"{step_name} Accuracy", var_name="Metric")

step1_melted = melt_df(step1_df, "Step 1")
step2_melted = melt_df(step2_df, "Step 2")
hlala_melted = melt_df(hlala_df, "Hlala")

# Combine all melted dataframes
combined_df = pd.concat([step1_melted, step2_melted, hlala_melted])

# Plotting
plt.figure(figsize=(12, 8))
sns.set(style="whitegrid")

# Create bar plot using seaborn
sns.barplot(x="Depth", y="value", hue="Metric", data=combined_df)
plt.title("Accuracies of Different Steps at Various Depths and Error Rates")
plt.xlabel("Sequencing Depth")
plt.ylabel("Accuracy")
plt.legend(title="Step")

# Save the plot to a file
plt.savefig("accuracies_bar_plot.png")
plt.close()