import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.transforms as mtransforms

# Define the data for the heatmap
data = {
    'Locus': ['A', 'B', 'C', 'E', 'F', 'G', 'H', 'J', 'K', 'L', 'N', 'P', 'S', 'T', 'U', 'V', 'W', 'Y',
              'DRA', 'DQA1', 'DQA2', 'DQB1', 'DQB2', 'DPA1', 'DPA2', 'DPB1', 'DPB2', 'DMA', 'DMB', 'DOA', 'DOB',
              'DRB1', 'DRB3', 'DRB4', 'DRB5', 'MICA', 'MICB', 'TAP1', 'TAP2'],
    'SpecLong': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                 1, 1, 1, 1, 1, 1, 1, 1],
    'SpecHLA':  [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 
                 1, 0, 0, 0, 0, 0, 0, 0],
    'HLA*LA':    [1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 
                 1, 1, 1, 0, 0, 0, 0, 0]
}

# Convert to a pandas DataFrame
df = pd.DataFrame(data)

# Transpose the DataFrame so that systems are rows and genes are columns
df_transposed = df.T
print(df_transposed)

# Define a custom discrete colormap (e.g., light blue for 0 and green for 1)
cmap = ListedColormap(['#DFF3F8', '#519D78'])  # '#DFF3F8' for 0, '#519D78' for 1

# Create the heatmap
plt.figure(figsize=(12, 16))  # 增加高度，确保显示所有基因
ax = sns.heatmap(df_transposed, cmap=cmap, cbar=False, linewidths=0.5, linecolor='gray')

# Loop through the DataFrame to add diagonal lines on cells with a value of 0
for i in range(df_transposed.shape[0]):  # rows (systems)
    for j in range(df_transposed.shape[1]):  # columns (genes)
        if df_transposed.iloc[i, j] == 0:
            # Draw a diagonal line from top left to bottom right
            ax.plot([j, j + 1], [i + 1, i], color='gray', lw=1.5)

# **确保纵坐标显示Locus基因名称**
# 使用df.index确保显示Locus基因
ax.set_yticklabels(df.index, rotation=0, fontsize=10)

# **手动设置ylim以确保显示所有基因**
ax.set_ylim(len(df_transposed), 0)  # 设置ylim范围，确保完整显示

# 调整y轴标签的可见性，确保不自动裁剪
ax.tick_params(axis='y', which='both', length=0)  # 移除刻度线

# **手动调整子图边距，以确保标签不会被剪裁**
plt.subplots_adjust(left=0.3, right=0.9, top=0.95, bottom=0.1)  # left 增大使得Y轴标签有足够空间

# Create custom patches for the legend
handled_patch = mpatches.Patch(facecolor='#519D78', edgecolor='gray', label='Handled')  # Solid green

# Custom patch for 'Not Handled' including diagonal slash
class DiagonalPatch(mpatches.Patch):
    def __init__(self, *args, **kwargs):
        super(DiagonalPatch, self).__init__(*args, **kwargs)

    def draw(self, renderer):
        super(DiagonalPatch, self).draw(renderer)
        # Add diagonal line on the patch
        x0, y0, w, h = self.get_bbox().bounds
        diagonal_line = mtransforms.Affine2D().rotate_deg_around(x0 + w / 2, y0 + h / 2, -45)
        line = Line2D([x0, x0 + w], [y0, y0 + h], transform=diagonal_line + self.get_transform(), color='gray', lw=1.5)
        line.draw(renderer)

# Create a custom patch for "Not Handled" with diagonal slash
not_handled_patch = DiagonalPatch(facecolor='#DFF3F8', edgecolor='gray', label='Not Handled')

# Add the custom legend
plt.legend(handles=[handled_patch, not_handled_patch], loc='upper right', bbox_to_anchor=(1.15, 1), fontsize=8, title_fontsize=12, frameon=False)

# 调整布局，确保所有元素不重叠
plt.tight_layout()

# Save the plot to a PDF file
plt.savefig(f"handle_hla.svg", format='svg', bbox_inches='tight', dpi=600)

# Show plot
plt.show()