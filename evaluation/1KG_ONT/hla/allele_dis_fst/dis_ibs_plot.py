import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# 读取个体间的距离矩阵
distance_matrix = np.loadtxt('dis_matrix_A.mibs')

# 读取个体的 ID 文件
ids = pd.read_csv('dis_matrix_A.mibs.id', delim_whitespace=True, header=None, names=['FID', 'IID'])

# 读取种群信息文件
pop_info = pd.read_csv('pop.txt', delim_whitespace=True)

# 将种群信息合并到 IDs 数据框中
ids = ids.merge(pop_info, on=['FID', 'IID'], how='left')

# 获取唯一种群
unique_pops = ids['Population'].unique()

# 初始化种群间平均距离矩阵
pop_mean_dist = pd.DataFrame(index=unique_pops, columns=unique_pops)

# 计算种群间的平均距离
for pop1 in unique_pops:
    for pop2 in unique_pops:
        # 找到属于这两个种群的个体的索引
        indices_pop1 = ids[ids['Population'] == pop1].index
        indices_pop2 = ids[ids['Population'] == pop2].index
        
        # 提取这两个种群之间的距离
        distances = []
        for id1 in indices_pop1:
            for id2 in indices_pop2:
                distances.append(distance_matrix[id1, id2])
        
        # 计算平均距离
        mean_distance = np.mean(distances)
        pop_mean_dist.loc[pop1, pop2] = mean_distance

# 输出种群间的平均距离矩阵
print(pop_mean_dist)

# 保存为 CSV 文件
pop_mean_dist.to_csv('population_distance_matrix.csv', index=True)

# 将距离矩阵转为浮点数格式
pop_mean_dist = pop_mean_dist.astype(float)

# 绘制热图
plt.figure(figsize=(10, 8))
sns.heatmap(pop_mean_dist, annot=True, cmap="YlGnBu", fmt=".5f")
plt.title('Population Distance Matrix')
plt.xlabel('Population')
plt.ylabel('Population')
plt.show()