import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

# Data
data = {
    "HG002-IsoSeq": [2, 2, 2, 2, 2, 2],
    "HG002-MAS-Seq": [2, 2, 2, 2, 2, 1],
    "HG002-ONT-Direct": [2, 2, 2, 2, 2, 2],
    "HG004-IsoSeq": [2, 2, 2, 2, 2, 2],
    "HG004-MAS-Seq": [2, 2, 2, 2, 2, 2],
    "HG004-ONT-Direct": [2, 2, 2, 2, 2, 2],
    "HG005-IsoSeq": [2, 2, 2, 2, 2, 2],
    "HG005-MAS-Seq": [2, 2, 2, 2, 2, 1],
    "HG005-ONT-Direct": [2, 2, 1, 2, 2, 2],
}

# Row labels
index = ['A', 'B', 'C', 'DQB1', 'DPB1', 'DRB1']

# Create DataFrame
df = pd.DataFrame(data, index=index)

# Custom colormap for values 0, 1, and 2
colors = ['#A8C9E0', '#468EBF', '#14C0CC']
cmap = ListedColormap(colors)

# Normalization to ensure correct mapping of colors to 0, 1, and 2
norm = BoundaryNorm([0, 1, 2, 3], cmap.N)

# Plot heatmap
plt.figure(figsize=(12, 6))
heatmap = sns.heatmap(df, cmap=cmap, cbar=False, linewidths=1, vmin=0, vmax=2, norm=norm)
# set font size for labels
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Create custom legend using three colored cells
legend_labels = [Patch(facecolor='#A8C9E0', edgecolor='black', label='Unresolved'),
                 Patch(facecolor='#468EBF', edgecolor='black', label='Half resolved'),
                 Patch(facecolor='#14C0CC', edgecolor='black', label='Resolved')]

# Add the custom legend to the plot
plt.legend(handles=legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=14)

# Set labels, title, and show plot
# plt.title("Heatmap of HG002, HG004, and HG005 Data with Custom Legend")
plt.tight_layout()  # Adjust layout to make room for the legend
# save to svg
plt.savefig("heatmap_custom_legend.svg", format="svg", dpi=600, transparent=True)
plt.show()