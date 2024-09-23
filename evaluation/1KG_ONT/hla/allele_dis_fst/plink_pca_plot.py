import pandas as pd
import matplotlib.pyplot as plt

# 读取 PCA 结果
pca_data = pd.read_csv('pca_out.eigenvec', delim_whitespace=True, header=None)

# 绘制前两主成分的散点图
plt.scatter(pca_data[2], pca_data[3])
plt.title('PCA: First two principal components')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.show()