#!/bin/bash
# 检查是否提供了参考基因组
REFERENCE="/scratch/project/cs_shuaicli/hg38_ref/BWA_GATK_index/hg38.fa"
if [ ! -f "$REFERENCE" ]; then
  echo "参考基因组 $REFERENCE 不存在，请确保已下载并放置在正确路径。"
  exit 1
fi

# 家系信息
six_trios=(
  "2418 NA19828 NA19818 NA19819"
  "CLM16 HG01258 HG01256 HG01257"
  "1463-Paternal NA12877 NA12889 NA12890"
  "1463-Maternal NA12878 NA12891 NA12892"
  "SH006 HG00420 HG00418 HG00419"
  "Y077 NA19129 NA19128 NA19127"
)

# 创建输出目录
mkdir -p ./output

# 遍历每个家系
for trio in "${six_trios[@]}"; do
  set -- $trio   # 解析家系信息
  trio_id=$1
  child=$2
  father=$3
  mother=$4

  echo "处理家系: $trio_id - 儿童: $child, 父亲: $father, 母亲: $mother"

  # 定义样本数组
  samples=($child $father $mother)

  # 遍历每个样本并进行比对和变异检测
  for sample in "${samples[@]}"; do
    echo "处理样本: $sample"

    # 定义FASTQ文件路径
    fastq1="./fastq_files/${sample}_R1.fastq"
    fastq2="./fastq_files/${sample}_R2.fastq"

    if [[ ! -f "$fastq1" || ! -f "$fastq2" ]]; then
      echo "FASTQ文件缺失: $fastq1 或 $fastq2"
      continue
    fi

    # 比对FASTQ文件到参考基因组，生成SAM文件
    sam_output="./output/${sample}.sam"
    minimap2 -ax sr "$REFERENCE" "$fastq1" "$fastq2" > "$sam_output"

    # 将SAM文件转换为BAM文件，进行排序并生成索引
    bam_output="./output/${sample}.sorted.bam"
    samtools view -Sb "$sam_output" | samtools sort -o "$bam_output"
    samtools index "$bam_output"

    # 删除中间的SAM文件以节省空间
    rm "$sam_output"

    # 使用longshot进行变异检测
    vcf_output="./output/${sample}.vcf"
    longshot -P -F -r "$REFERENCE" -o "$vcf_output" -bam "$bam_output"

    echo "样本 $sample 的比对和变异检测已完成，结果保存在 $vcf_output"
  done
done

echo "所有家系的样本处理已完成。"
