import pandas as pd

# 读取CIF文件，提取_atom_site.label_asym_id列
cif_file = "HLA-A-11-01-01-01-all/fold_2024_12_08_20_38_model_0.cif"
cif_file = "KIR2DL1-0030204-all/fold_2024_12_08_2_model_0.cif"
from Bio.PDB import MMCIFParser

# 加载 mmCIF 文件
# cif_file = "example.cif"  # 替换为你的文件名
parser = MMCIFParser(QUIET=True)

# 解析结构
try:
    structure = parser.get_structure("model", cif_file)
except Exception as e:
    print(f"Error parsing CIF file: {e}")
    exit()

# 遍历结构并提取链 ID
print("Chain IDs in the structure:")
for model in structure:  # 遍历所有模型
    for chain in model:  # 遍历所有链
        print(f"Chain ID: {chain.id}")