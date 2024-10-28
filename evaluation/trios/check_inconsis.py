import pysam

# 定义孟德尔遗传定律的检查函数
def mendelian_inheritance(parent1_alleles, parent2_alleles, child_alleles):
    """
    检查是否符合孟德尔遗传定律
    parent1_alleles, parent2_alleles, child_alleles：分别是父母和孩子的等位基因列表，例如 ['A', 'T'], ['A', 'A'], ['A', 'T']
    返回值：True 表示符合，False 表示不符合
    """
    # 获取父母的所有可能的等位基因组合，并进行排序
    possible_alleles = set([''.join(sorted([a, b])) for a in parent1_alleles for b in parent2_alleles])
    
    # 检查孩子的等位基因是否在父母的可能组合中，并对孩子的基因型排序
    return ''.join(sorted(child_alleles)) in possible_alleles

def get_alleles(record, gt):
    """
    根据给定的record和基因型（GT），返回该基因型对应的等位基因
    record: VCF中的变异记录
    gt: 基因型（如 '0/1'）
    返回值：基因型对应的等位基因列表
    """
    alleles = []
    for allele_index in gt:
        if allele_index == 0:
            alleles.append(record.ref)  # 参考等位基因
        elif allele_index > 0:
            alleles.append(record.alts[allele_index - 1])  # 替代等位基因
        else:
            alleles.append('.')
    return alleles

def count_mendelian_errors(vcf_file_parent1, vcf_file_parent2, vcf_file_child, regions):
    """
    统计在给定区域内不符合孟德尔遗传定律的突变数量
    vcf_file_parent1, vcf_file_parent2, vcf_file_child: 父母和孩子的VCF文件路径
    regions: 指定的基因组区域列表，格式为 ['chr:start-end', 'chr2:start2-end2', ...]
    """
    # 打开VCF文件
    parent1_vcf = pysam.VariantFile(vcf_file_parent1)
    parent2_vcf = pysam.VariantFile(vcf_file_parent2)
    child_vcf = pysam.VariantFile(vcf_file_child)
    
    # 初始化错误计数
    mendelian_errors = 0
    total_variants = 0

    # 遍历每个区域
    for gene, region in regions.items():
        # 解析区域
        region_chr, region_range = region.split(':')
        region_start, region_end = map(int, region_range.split('-'))

        # 在每个区域内遍历变异
        print(f"正在检查区域 {region}")
        # for parent1_record in parent1_vcf.fetch(region_chr, region_start, region_end):
        for child_record in child_vcf.fetch(region_chr, region_start, region_end):

            total_variants += 1
            parent1_record = parent1_vcf.fetch(region_chr, child_record.pos - 1, child_record.pos)
            parent2_record = parent2_vcf.fetch(region_chr, child_record.pos - 1, child_record.pos)
            # child_record = child_vcf.fetch(region_chr, child_record.pos - 1, child_record.pos)
            
            parent2_record = next(parent2_record, None)
            parent1_record = next(parent1_record, None)

           
            
            if parent2_record is None or parent1_record is None:
                continue



            

            # 获取父母和孩子的基因型
            parent1_gt = parent1_record.samples[0]['GT']
            parent2_gt = parent2_record.samples[0]['GT']
            child_gt = child_record.samples[0]['GT']

            # print(parent1_record.ref, parent1_record.alts, parent1_gt)
            # print(parent2_record.ref, parent2_record.alts, parent2_gt)
            # print(child_record.ref, child_record.alts, child_gt)


            

            # 确保三者的变异位置一致
            if parent1_record.pos != parent2_record.pos or parent1_record.pos != child_record.pos:
                continue

             # if record is not snvs continue
            if len(parent1_record.ref) > 1 or len(parent1_record.alts[0]) > 1:
                continue
            if len(parent2_record.ref) > 1 or len(parent2_record.alts[0]) > 1:
                continue
            if len(child_record.ref) > 1 or len(child_record.alts[0]) > 1:
                continue


            # 如果基因型中有缺失值（例如 './.'），直接跳过
            if None in parent1_gt or None in parent2_gt or None in child_gt:
                continue

            # 获取实际的等位基因
            parent1_alleles = get_alleles(parent1_record, parent1_gt)
            parent2_alleles = get_alleles(parent2_record, parent2_gt)
            child_alleles = get_alleles(child_record, child_gt)

            # 检查是否符合孟德尔遗传定律
            if not mendelian_inheritance(parent1_alleles, parent2_alleles, child_alleles):
                mendelian_errors += 1
                # print(f"不符合孟德尔遗传定律的突变：{region_chr}:{parent1_record.pos} 父母({parent1_alleles}, {parent2_alleles}) 孩子({child_alleles})")

    return mendelian_errors, total_variants
# 示例使用 c_HG00420

# HLA-T: chr6:29896443-29898947
# HLA-U: chr6:29933764-29934880
# HLA-DPA1: chr6:33064569-33080775
# HLA-G: chr6:29826474-29831125
# MICA: chr6:31399784-31415315
# HLA-DQB1: chr6:32659467-32668383
# HLA-E: chr6:30489509-30494194
# HLA-W: chr6:29955834-29959058
# HLA-DRB5: chr6:32517353-32530287
region_name_dict={
    "HLA-T": "chr6:29896443-29898947",
    "HLA-U": "chr6:29933764-29934880",
    "HLA-DPA1": "chr6:33064569-33080775",
    "HLA-G": "chr6:29826474-29831125",
    "MICA": "chr6:31399784-31415315",
    "HLA-DQB1": "chr6:32659467-32668383",
    "HLA-E": "chr6:30489509-30494194",
    "HLA-W": "chr6:29955834-29959058",
    "HLA-DRB5": "chr6:32517353-32530287"
}

vcf_file_parent1 = "c_HG00420/vcf/merged_HG00418_bcftools.sorted.vcf.gz"
vcf_file_parent2 = "c_HG00420/vcf/merged_HG00419_bcftools.sorted.vcf.gz"
vcf_file_child = "c_HG00420/vcf/merged_HG00420_bcftools.sorted.vcf.gz"
regions = ["chr6:29896443-29898947", "chr6:29933764-29934880", "chr6:33064569-33080775", "chr6:29826474-29831125", "chr6:31399784-31415315", "chr6:32659467-32668383", "chr6:30489509-30494194", "chr6:29955834-29959058", "chr6:32517353-32530287"]

errors, total_variants = count_mendelian_errors(vcf_file_parent1, vcf_file_parent2, vcf_file_child, region_name_dict)
print(f"不符合孟德尔遗传定律的突变数量: {errors}")
print(f"总共检查的突变数量: {total_variants}")
print(f"错误率: {errors / total_variants:.2%}")

# c_HG01258

# HLA-T: chr6:29896443-29898947
# TAP2: chr6:32821831-32838739
# HLA-DQA2: chr6:32741391-32747198
# HLA-DRB1: chr6:32577902-32589848
# HLA-DRB5: chr6:32517353-32530287
# HLA-L: chr6:30259562-30266951
# HLA-DPB2: chr6:33112516-33129113

region_name_dict={
    "HLA-T": "chr6:29896443-29898947",
    "TAP2": "chr6:32821831-32838739",
    "HLA-DQA2": "chr6:32741391-32747198",
    "HLA-DRB1": "chr6:32577902-32589848",
    "HLA-DRB5": "chr6:32517353-32530287",
    "HLA-L": "chr6:30259562-30266951",
    "HLA-DPB2": "chr6:33112516-33129113"
}

vcf_file_parent1 = "c_HG01258/vcf/merged_HG01256_bcftools.sorted.vcf.gz"
vcf_file_parent2 = "c_HG01258/vcf/merged_HG01257_bcftools.sorted.vcf.gz"
vcf_file_child = "c_HG01258/vcf/merged_HG01258_bcftools.sorted.vcf.gz"
regions = ["chr6:29896443-29898947", "chr6:32821831-32838739", "chr6:32741391-32747198", "chr6:32577902-32589848", "chr6:32517353-32530287", "chr6:30259562-30266951", "chr6:33112516-33129113"]

errors, total_variants = count_mendelian_errors(vcf_file_parent1, vcf_file_parent2, vcf_file_child, region_name_dict)
print(f"不符合孟德尔遗传定律的突变数量: {errors}")
print(f"总共检查的突变数量: {total_variants}")
print(f"错误率: {errors / total_variants:.2%}")

# c_NA12877

# HLA-DQA1: chr6:32628179-32655272
# HLA-DPA2: chr6:33091482-33097295
# MICA: chr6:31399784-31415315
# HLA-DQB1: chr6:32659467-32668383
# HLA-DQA2: chr6:32741391-32747198
# DRA: chr6:32439878-32445046
# HLA-L: chr6:30259562-30266951
# HLA-DRB5: chr6:32517353-32530287

region_name_dict={
    "HLA-DQA1": "chr6:32628179-32655272",
    "HLA-DPA2": "chr6:33091482-33097295",
    "MICA": "chr6:31399784-31415315",
    "HLA-DQB1": "chr6:32659467-32668383",
    "HLA-DQA2": "chr6:32741391-32747198",
    "DRA": "chr6:32439878-32445046",
    "HLA-L": "chr6:30259562-30266951",
    "HLA-DRB5": "chr6:32517353-32530287"
}


vcf_file_parent1 = "c_NA12877/vcf/merged_NA12889_bcftools.sorted.vcf.gz"
vcf_file_parent2 = "c_NA12877/vcf/merged_NA12890_bcftools.sorted.vcf.gz"
vcf_file_child = "c_NA12877/vcf/merged_NA12877_bcftools.sorted.vcf.gz"
regions = ["chr6:32628179-32655272", "chr6:33091482-33097295", "chr6:31399784-31415315", "chr6:32659467-32668383", "chr6:32741391-32747198", "chr6:32439878-32445046", "chr6:30259562-30266951", "chr6:32517353-32530287"]

errors, total_variants = count_mendelian_errors(vcf_file_parent1, vcf_file_parent2, vcf_file_child, region_name_dict)
print(f"不符合孟德尔遗传定律的突变数量: {errors}")
print(f"总共检查的突变数量: {total_variants}")
print(f"错误率: {errors / total_variants:.2%}")


# c_NA12878

# HLA-T: chr6:29896443-29898947
# MICA: chr6:31399784-31415315
# MICB: chr6:31494881-31511124
# HLA-DOA: chr6:33004182-33009591
# HLA-DRB5: chr6:32517353-32530287

region_name_dict={
    "HLA-T": "chr6:29896443-29898947",
    "MICA": "chr6:31399784-31415315",
    "MICB": "chr6:31494881-31511124",
    "HLA-DOA": "chr6:33004182-33009591",
    "HLA-DRB5": "chr6:32517353-32530287"
}


vcf_file_parent1 = "c_NA12878/vcf/merged_NA12891_bcftools.sorted.vcf.gz"
vcf_file_parent2 = "c_NA12878/vcf/merged_NA12892_bcftools.sorted.vcf.gz"
vcf_file_child = "c_NA12878/vcf/merged_NA12878_bcftools.sorted.vcf.gz"
regions = ["chr6:29896443-29898947", "chr6:31399784-31415315", "chr6:31494881-31511124", "chr6:33004182-33009591", "chr6:32517353-32530287"]

errors, total_variants = count_mendelian_errors(vcf_file_parent1, vcf_file_parent2, vcf_file_child, region_name_dict)
print(f"不符合孟德尔遗传定律的突变数量: {errors}")
print(f"总共检查的突变数量: {total_variants}")
print(f"错误率: {errors / total_variants:.2%}")


# c_NA19129

# HLA-DMA: chr6:32948613-32969094
# HLA-DQA1: chr6:32628179-32655272
# MICA: chr6:31399784-31415315
# HLA-DOA: chr6:33004182-33009591
# HLA-DQA2: chr6:32741391-32747198
# DRA: chr6:32439878-32445046
# HLA-DRB1: chr6:32577902-32589848
# HLA-DRB4: chr6:32530288-32543197
# HLA-DRB5: chr6:32517353-32530287

region_name_dict={
    "HLA-DMA": "chr6:32948613-32969094",
    "HLA-DQA1": "chr6:32628179-32655272",
    "MICA": "chr6:31399784-31415315",
    "HLA-DOA": "chr6:33004182-33009591",
    "HLA-DQA2": "chr6:32741391-32747198",
    "DRA": "chr6:32439878-32445046",
    "HLA-DRB1": "chr6:32577902-32589848",
    "HLA-DRB4": "chr6:32530288-32543197",
    "HLA-DRB5": "chr6:32517353-32530287"
}

vcf_file_parent1 = "c_NA19129/vcf/merged_NA19127_bcftools.sorted.vcf.gz"
vcf_file_parent2 = "c_NA19129/vcf/merged_NA19128_bcftools.sorted.vcf.gz"
vcf_file_child = "c_NA19129/vcf/merged_NA19129_bcftools.sorted.vcf.gz"
regions = ["chr6:32948613-32969094", "chr6:32628179-32655272", "chr6:31399784-31415315", "chr6:33004182-33009591", "chr6:32741391-32747198", "chr6:32439878-32445046", "chr6:32577902-32589848", "chr6:32530288-32543197", "chr6:32517353-32530287"]

errors, total_variants = count_mendelian_errors(vcf_file_parent1, vcf_file_parent2, vcf_file_child, region_name_dict)
print(f"不符合孟德尔遗传定律的突变数量: {errors}")
print(f"总共检查的突变数量: {total_variants}")
print(f"错误率: {errors / total_variants:.2%}")

# c_NA19828

# HLA-A: chr6:29941260-29949572
# HLA-DQA1: chr6:32628179-32655272
# HLA-E: chr6:30489509-30494194
# TAP1: chr6:32838740-32855648
# HLA-DRB1: chr6:32577902-32589848
# HLA-DRB4: chr6:32530288-32543197
region_name_dict={
    "HLA-A": "chr6:29941260-29949572",
    "HLA-DQA1": "chr6:32628179-32655272",
    "HLA-E": "chr6:30489509-30494194",
    "TAP1": "chr6:32838740-32855648",
    "HLA-DRB1": "chr6:32577902-32589848",
    "HLA-DRB4": "chr6:32530288-32543197"
}
vcf_file_parent1= "c_NA19828/vcf/merged_NA19819_bcftools.sorted.vcf.gz"
vcf_file_parent2 = "c_NA19828/vcf/merged_NA19828_bcftools.sorted.vcf.gz"
vcf_file_child = "c_NA19828/vcf/merged_NA19828_bcftools.sorted.vcf.gz"
regions= region_name_dict.values()
regions = ["chr6:29941260-29949572", "chr6:32628179-32655272", "chr6:30489509-30494194", "chr6:32838740-32855648", "chr6:32577902-32589848", "chr6:32530288-32543197"]

errors, total_variants = count_mendelian_errors(vcf_file_parent1, vcf_file_parent2, vcf_file_child, region_name_dict)
print(f"不符合孟德尔遗传定律的突变数量: {errors}")
print(f"总共检查的突变数量: {total_variants}")
print(f"错误率: {errors / total_variants:.2%}")


