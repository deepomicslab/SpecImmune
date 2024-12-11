from Bio.PDB import MMCIFParser, PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model

def merge_cif_to_pdb(cif_file1, cif_file2, output_file):
    """
    合并两个 mmCIF 文件，并按链 ID (A, B, C, ...) 重新编号，保存为 PDB 格式
    """
    parser = MMCIFParser(QUIET=True)
    
    # 加载两个 CIF 文件
    structure1 = parser.get_structure("structure1", cif_file1)
    structure2 = parser.get_structure("structure2", cif_file2)
    
    # 合并链信息
    merged_chains = []
    for model in structure1:
        for chain in model:
            merged_chains.append(chain)
    for model in structure2:
        for chain in model:
            merged_chains.append(chain)
    
    # 按顺序重新分配链 ID (A, B, C, ...)
    chain_id_list = [chr(65 + i) for i in range(len(merged_chains))]  # A, B, C, ...
    for new_id, chain in zip(chain_id_list, merged_chains):
        chain.id = new_id  # 修改链 ID
    
    # 创建一个新的 Structure 对象
    merged_structure = Structure("merged_structure")
    new_model = Model(0)
    merged_structure.add(new_model)
    
    for chain in merged_chains:
        new_model.add(chain)

    # 保存到 PDB 文件
    io = PDBIO()
    io.set_structure(merged_structure)
    io.save(output_file)
    print(f"Merged structure saved to PDB file: {output_file}")

# # 示例：合并两个 mmCIF 文件并保存为 PDB 格式
# cif_file1 = "file1.cif"  # 替换为你的第一个文件路径
# cif_file2 = "file2.cif"  # 替换为你的第二个文件路径
# output_file = "merged_structure.pdb"  # 输出文件名

# merge_cif_to_pdb(cif_file1, cif_file2, output_file)


# 示例：合并两个文件并保存
cif_file1 = "HLA-A-11-01-01-01-all/fold_2024_12_08_20_38_model_0.cif"
cif_file2 = "KIR2DL1-0030204-all/fold_2024_12_08_2_model_0.cif"
output_file = "hla-a-kir2dl1-merge/merged.pdb"  # 输出文件名

merge_cif_to_pdb(cif_file1, cif_file2, output_file)