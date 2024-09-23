import re
import pandas as pd
import sys

# GTF 文件路径
gtf_file = sys.argv[1]
# 输出 BED 文件路径
bed_file = sys.argv[2]

# 定义读取GTF文件的函数
def parse_gtf(gtf_file):
    transcripts = {}

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            
            fields = line.strip().split('\t')
            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            # 解析属性字段，提取 transcript_id
            transcript_match = re.search('transcript_id "([^"]+)"', attributes)
            if not transcript_match:
                continue
            transcript_id = transcript_match.group(1)

            # 初始化字典
            if transcript_id not in transcripts:
                transcripts[transcript_id] = {
                    'chrom': chrom,
                    'strand': strand,
                    'thick_start': None,
                    'thick_end': None,
                    'exons': []
                }

            # 处理 CDS 特征
            if feature == 'CDS':
                start = int(start)
                end = int(end)
                if transcripts[transcript_id]['thick_start'] is None or start < transcripts[transcript_id]['thick_start']:
                    transcripts[transcript_id]['thick_start'] = start
                if transcripts[transcript_id]['thick_end'] is None or end > transcripts[transcript_id]['thick_end']:
                    transcripts[transcript_id]['thick_end'] = end

            # 处理 exon 特征
            if feature == 'exon':
                start = int(start) - 1  # BED 是 0-based，所以 GTF 的起始位置要减 1
                end = int(end)
                transcripts[transcript_id]['exons'].append((start, end))

    return transcripts

# 转换为 BED-like 格式
def gtf_to_bed(transcripts, bed_file):
    with open(bed_file, 'w') as f:
        for transcript_id, data in transcripts.items():
            chrom = data['chrom']
            strand = data['strand']
            exons = sorted(data['exons'])

            # 获取转录本的开始和结束位置
            tx_start = exons[0][0]
            tx_end = exons[-1][1]

            # 处理 thickStart 和 thickEnd (CDS 区域)
            thick_start = data['thick_start'] if data['thick_start'] is not None else tx_start
            thick_end = data['thick_end'] if data['thick_end'] is not None else tx_end

            # 计算 blockCount, blockSizes, blockStarts
            block_count = len(exons)
            block_sizes = ",".join([str(end - start) for start, end in exons]) + ','
            block_starts = ",".join([str(start - tx_start) for start, _ in exons]) + ','

            # 输出格式: chrom, tx_start, tx_end, transcript_id, score (0), strand, thick_start, thick_end, reserved (0), block_count, block_sizes, block_starts
            f.write(f"{chrom}\t{tx_start}\t{tx_end}\t{transcript_id}\t0\t{strand}\t{thick_start}\t{thick_end}\t0\t{block_count}\t{block_sizes}\t{block_starts}\n")

# 解析 GTF 文件并转换为 BED 格式
transcripts = parse_gtf(gtf_file)
gtf_to_bed(transcripts, bed_file)

print(f"Conversion complete. BED file saved to {bed_file}")