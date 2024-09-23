from ete3 import Tree, TreeStyle, NodeStyle
import pandas as pd

# 读取树文件
tree = Tree("p_dis.nwk")

# 读取样本-人口注释文件
annotations = pd.read_csv("pop.txt", sep="\t")

# 为每个 population 定义一种颜色
population_colors = {
    "EUR": "red",
    "AFR": "green",
    "EAS": "blue",
    "SAS": "purple",
    "AMR": "orange",
}

def layout(node):
    if node.is_leaf():
        name = node.name
        if name in sample_to_population:
            population = sample_to_population[name]
            color = population_colors.get(population, "black")
            
            # 设置叶节点样式
            nstyle = NodeStyle()
            nstyle["fgcolor"] = color  # 设置颜色
            nstyle["size"] = 10  # 设置节点大小
            node.set_style(nstyle)
# 样本名到 population 的映射
sample_to_population = dict(zip(annotations['Individual'], annotations['Population']))

# 设置节点样式，包括叶子节点和分支颜色
def set_leaf_branch_style(node, color):
    nstyle = NodeStyle()
    nstyle["fgcolor"] = color  # 设置节点颜色
    nstyle["size"] = 10        # 设置节点大小 (叶节点)
    
    # 设置连接到叶节点的分支的颜色和宽度
    nstyle["vt_line_color"] = color  # 设置垂直线颜色
    nstyle["hz_line_color"] = color  # 设置水平线颜色
    nstyle["vt_line_width"] = 4      # 设置垂直线宽度
    nstyle["hz_line_width"] = 4      # 设置水平线宽度
    
    node.set_style(nstyle)

# 递归地为叶节点及其分支设置颜色
def color_leaf_branches(node):
    if node.is_leaf():
        name = node.name
        if name in sample_to_population:
            population = sample_to_population[name]
            color = population_colors.get(population, "black")
            set_leaf_branch_style(node, color)

# 定义树的布局
ts = TreeStyle()
ts.show_leaf_name = True  # 显示叶节点的名字
ts.mode = "c"  # 设置为环状布局

# 为树中的每个叶节点及其相连的分支设置颜色
for leaf in tree.iter_leaves():
    color_leaf_branches(leaf)

# 绘制树
tree.show(tree_style=ts)