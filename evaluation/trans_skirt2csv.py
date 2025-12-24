#!/usr/bin/env python
import os
import re
import argparse
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="解析 merged.hap.csv 文件，将每个基因的两个 allele 调用解析后输出到 CSV 文件中"
    )
    parser.add_argument(
        "-i", "--input_file",
        type=str,
        default="out_dir/NA21309/NA21309/out/merged.hap.csv",
        help="输入 CSV 文件路径，默认: out_dir/NA21309/NA21309/out/merged.hap.csv"
    )
    parser.add_argument(
        "-o", "--output_file",
        type=str,
        default="kir_alleles_parsed.csv",
        help="输出 CSV 文件路径，默认: kir_alleles_parsed.csv"
    )
    return parser.parse_args()

def process_cell(cell_text, gene):
    """
    对单个 cell 的文本进行解析处理：
      1. 优先查找括号中是否有形如 "of <数字>" 的信息；
         如果有，则返回： gene + "*" + <数字> 。
      2. 否则，先去除所有括号及括号中的内容，
         合并空白、去除尾部特殊字符（如 #, $, +），
         并检查结果中是否缺少星号 "*" —— 如果缺少，则在 gene 与末尾数字中间插入星号。
    """
    # ① 尝试查找包含 "of <数字>" 的模式
    m = re.search(r"of\s*(\d+)", cell_text)
    if m:
        allele_num = m.group(1)
        return gene + "*" + allele_num

    # ② 否则处理：去除括号及括号内的内容
    cleaned = re.sub(r"\([^)]*\)", "", cell_text)
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    # 如果有 "~" 分隔，则只处理第一个token（可根据需要调整，此处优先取第一个）
    if "~" in cleaned:
        tokens = cleaned.split("~")
        token = tokens[0].strip()
    else:
        token = cleaned.strip()
    # 去除尾部特殊符号
    token = token.rstrip("#$+ ").strip()
    
    # 如果 token 中没有 "*"，尝试在 gene 后面插入 "*" —— 前提是 token以 gene 开头后跟数字
    if "*" not in token:
        # 针对类似 "KIR3DL3005" 这样的格式，提取 gene部分与数字部分
        m2 = re.match(re.escape(gene) + r'(\d+)', token)
        if m2:
            return gene + "*" + m2.group(1)
    return token

def process_row(row, gene):
    """
    对于同一行（对应一个 allele 调用），依次处理从第二列开始的每个 non-empty cell：
      - 如果某个 cell 的原始文本中含有 "of" 信息，则优先采用该 cell 处理后的结果；
      - 否则，按顺序使用第一个非空 cell 的处理结果。
    返回该行处理后的 allele 字符串。
    """
    candidate = None
    candidate_with_of = None
    # 从第二列开始处理（索引1之后）
    for cell in row.iloc[1:]:
        if pd.notnull(cell) and str(cell).strip() != "":
            cell_text = str(cell).strip()
            proc = process_cell(cell_text, gene)
            # 如果原始 cell 包含 "of" 信息则记为优先候选
            if re.search(r"of\s*\d+", cell_text):
                candidate_with_of = proc
            elif candidate is None:
                candidate = proc
    # 如果优先候选存在，则返回优先候选，否则返回第一个候选
    if candidate_with_of is not None:
        return candidate_with_of
    else:
        return candidate if candidate is not None else ""

def parse_kir_csv(input_file):
    """
    解析 CSV 文件：
      - 第一列为 contig，如 "KIR3DL3_0" 或 "KIR3DL3_1" 用于确定 gene 及 allele 编号。
      - 从第二列开始的所有非空 cell，使用 process_row() 进行分行处理，得到具体 allele 调用。
    返回字典格式： { gene: { 0: allele_call, 1: allele_call } }
    """
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"未找到文件: {input_file}")

    # 此处假设 CSV 无表头，故 header=None
    df = pd.read_csv(input_file, header=None)
    results = {}
    for _, row in df.iterrows():
        contig_value = row.iloc[0]
        if not isinstance(contig_value, str) or "_" not in contig_value:
            continue
        try:
            gene, index_str = contig_value.rsplit("_", 1)
            allele_index = int(index_str)
        except ValueError:
            continue

        allele_call = process_row(row, gene)
        if gene not in results:
            results[gene] = {}
        results[gene][allele_index] = allele_call
    return results

def save_results(results, output_file):
    """
    将 results 构建为 DataFrame，包含列：Gene、Allele1、Allele2，并输出到 CSV 文件中。
    """
    output_list = []
    for gene, alleles in results.items():
        allele1 = alleles.get(0, "")
        allele2 = alleles.get(1, "")
        output_list.append({
            "Gene": gene,
            "Allele1": allele1,
            "Allele2": allele2
        })
    df_out = pd.DataFrame(output_list)
    df_out.sort_values("Gene", inplace=True)
    df_out.to_csv(output_file, index=False)
    print(f"解析结果已保存到: {output_file}")

def main():
    args = parse_arguments()
    results = parse_kir_csv(args.input_file)
    save_results(results, args.output_file)

if __name__ == "__main__":
    main()