import os
import subprocess
import matplotlib.pyplot as plt
import random
from tqdm import tqdm

def generate_depth_file(bam_file, output_folder):
    depth_file = os.path.join(output_folder, os.path.basename(bam_file) + ".depth")
    with open(depth_file, "w") as out:
        subprocess.run(["samtools", "depth", "-aa", bam_file], stdout=out)
    return depth_file

def plot_depth_curves(depth_files, sample_names, output_folder):
    plt.figure(figsize=(12, 6))
    
    for depth_file, sample_name in zip(depth_files, sample_names):
        positions = []
        depths = []

        with open(depth_file, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                positions.append(int(fields[1]))
                depths.append(int(fields[2]))

        color = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])  # 随机颜色
        plt.plot(positions, depths, label=sample_name, color=color)

    plt.xlabel("Position")
    plt.ylabel("Depth")
    plt.title("Depth Curves for Samples")
    plt.legend()
    plt.grid(True)

    plot_file = os.path.join(output_folder, "combined_depth_plot.png")
    plt.savefig(plot_file)
    plt.close()

def process_bam_files(sample_file, bam_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    depth_files = []
    sample_names = []

    with open(sample_file, "r") as f:
        for line in tqdm(f):
            sample_name = line.strip()
            bam_file = os.path.join(bam_folder,f"{sample_name}", f"{sample_name}.HLA.sorted.bam")
            if os.path.exists(bam_file):
                print(f"Processing {bam_file}...")

                depth_file = generate_depth_file(bam_file, output_folder)
                depth_files.append(depth_file)
                sample_names.append(sample_name)
            else:
                print(f"BAM file for sample {bam_file} does not exist!")

    if depth_files:
        plot_depth_curves(depth_files, sample_names, output_folder)
        print(f"Combined depth curve saved in {output_folder}.")
    else:
        print("No valid BAM files were processed.")

if __name__ == "__main__":
    sample_file = "/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/1KGP_ONT/all.samples"  # 修改为你的样本名文件路径
    bam_folder = "/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/1KGP_ONT/readsout/HLA/"  # 修改为你存放BAM文件的目录
    output_folder = "/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/1KGP_ONT/HLA_depth"  # 修改为你希望保存输出文件的目录

    process_bam_files(sample_file, bam_folder, output_folder)