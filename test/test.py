import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

def plot_depth(file_path):
    # Check if the file exists
    if not os.path.isfile(file_path):
        print(f"Error: The file {file_path} does not exist.")
        return

    # Load depth data
    depth_data = pd.read_csv(file_path, sep='\t', header=None, names=['chromosome', 'position', 'depth'])

    # Plot using seaborn
    plt.figure(figsize=(12, 6))
    sns.lineplot(data=depth_data, x='position', y='depth', hue='chromosome', linewidth=1)
    # plot a horizontal line at y=30
    plt.axhline(y=24.5, color='r', linestyle='--', label='Threshold (30)')


    plt.title('Depth across multiple chromosomes')
    plt.xlabel('Position')
    plt.ylabel('Depth')
    plt.legend(title='Chromosome')
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_depth.py <path_to_depth_file>")
    else:
        file_path = sys.argv[1]
        plot_depth(file_path)