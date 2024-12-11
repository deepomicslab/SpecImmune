from Bio.PDB import PDBParser, PDBIO
from string import ascii_uppercase

def reassign_chain_ids(structure, start_index):
    """
    为结构中的链重新分配链 ID，从指定的起始索引开始。
    """
    chain_ids = list(ascii_uppercase)  # 链 ID 可用范围：A, B, C, ..., Z
    i = start_index
    for model in structure:
        for chain in model:
            chain.id = chain_ids[i]
            i += 1
    return i  # 返回下一个可用的链 ID 索引

# 定义解析器和输出器
parser = PDBParser(QUIET=True)
io = PDBIO()

# 读取两个 PDB 文件
structure1 = parser.get_structure("structure1", "cyp2d6.1-all/cyp2d6.1.pdb")
structure2 = parser.get_structure("structure2", "HLA-A-01010101-all/esm.pdb")

# 给两个结构重新分配链 ID
next_chain_index = reassign_chain_ids(structure1, 0)  # 从 A 开始
reassign_chain_ids(structure2, next_chain_index)      # 接着分配链 ID

# 获取第一个 Model
model1 = next(iter(structure1))  # 获取 structure1 的第一个 Model
model2 = next(iter(structure2))  # 获取 structure2 的第一个 Model

# 将 structure2 的链添加到 structure1 的 Model 中
for chain in model2:
    model1.add(chain)  # 将第二个结构的链添加到第一个结构中

# 保存为新的 PDB 文件
io.set_structure(structure1)
io.save("cyp2d6.1_hla-a-01010101-merge/ESM_merged_with_new_chain_ids.pdb")

print("两个文件已成功合并为 cyp2d6.1_hla-a-01010101-merge/ESM_merged_with_new_chain_ids.pdb")