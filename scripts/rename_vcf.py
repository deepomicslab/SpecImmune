import sys
import pysam
import tempfile

def rename_sample(input_vcf, new_sample_name, output_vcf):
    # 1. 读取原VCF头部，文本替换样本名
    with pysam.VariantFile(input_vcf, "r") as invcf:
        old_samples = list(invcf.header.samples)
        if len(old_samples) != 1:
            print("警告：输入VCF中的样本数为{}，将全部重命名为 {}。".format(len(old_samples), new_sample_name))
        header_lines = []
        for line in str(invcf.header).splitlines():
            if line.startswith("#CHROM"):
                parts = line.strip().split("\t")
                header_lines.append("\t".join(parts[:9] + [new_sample_name]))
            else:
                header_lines.append(line)
        new_header_str = "\n".join(header_lines) + "\n"

    # 2. 用修改后的header创建新的VariantFile header对象
    with tempfile.NamedTemporaryFile("wt", delete=False) as tmp_header:
        tmp_header.write(new_header_str)
        tmp_header.flush()
        header_path = tmp_header.name
    with pysam.VariantFile(header_path) as header_in:
        new_header = header_in.header.copy()

    # 3. 写新VCF，record直接写出即可
    with pysam.VariantFile(input_vcf, "r") as invcf, \
         pysam.VariantFile(output_vcf, "w", header=new_header) as outvcf:
        for record in invcf:
            outvcf.write(record)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("用法: python rename_vcf_sample.py <input.vcf[.gz]> <new_sample_name> <output.vcf>")
        sys.exit(1)
    rename_sample(sys.argv[1], sys.argv[2], sys.argv[3])