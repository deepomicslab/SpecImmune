import csv
import os

# Path to your CSV file
csv_file_path = 'SRP272207.csv'
data_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/VDJ/sra/"
results_dir = "/mnt/d/HLAPro_backup/Nanopore_optimize/vdj_results_tcr3/"

run_dict = {}
# Open the CSV file and read the top 10 lines (excluding the header)
with open(csv_file_path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    count = 0
    for row in reader:
        if count < 11:  # Only process the top 10 lines
            run_name = row['Run']
            sample_name = row['Sample Name']
            print(f"{run_name}\t{sample_name}")
            run_dict[run_name] = sample_name
            count += 1
        else:
            break

for run_name, sample_name in run_dict.items():
    if sample_name == "NA18506":
        continue
    print(f"Run: {run_name}, Sample: {sample_name}")
    ## download the sra and convert to fastq
    # cmd = f"""
    # prefetch -O {data_dir} {run_name}
    # fastq-dump --outdir {data_dir} --split-3 {data_dir}/{run_name}/{run_name}.sra --gzip
    # mv {data_dir}/{run_name}.fastq.gz {data_dir}/{sample_name}.fastq.gz
    # """
    # os.system(cmd)
    # break

    ## run typing
    cmd = f"""
        # echo python3 ../scripts/main.py -n {sample_name} -o {results_dir} -j 15 -y nanopore -i IG_TR  -r {data_dir}/{sample_name}.fastq.gz
        python3 ../scripts/main.py --hg38 /mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/hg38_no_alt.fa -n {sample_name} -o {results_dir} -j 15 -y nanopore -i IG_TR  -r {data_dir}/{sample_name}.fastq.gz
    """
    os.system(cmd)
