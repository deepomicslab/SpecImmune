import argparse
import pysam
import subprocess, argparse, os, re
import numpy as np

def bam_average_depth(bam_file):
    # 打开BAM文件
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        total_depth = 0
        total_positions = 0

        # 遍历BAM文件中的每一个染色体（reference）
        for ref in bam.references:
            # 遍历pileup中的每个列
            for pileupcolumn in bam.pileup(ref):
                total_depth += pileupcolumn.n
                total_positions += 1

    # 计算并返回平均深度
    average_depth = total_depth / total_positions
    return round(average_depth)

def check_svtype(sv_type_dict, chrom, start, end):
    
    if (chrom, start, end) in sv_type_dict:
        return sv_type_dict[(chrom, start, end)]
    else:
        regions = sv_type_dict.keys()
        for region in regions:
            region_chr, region_start, region_end = region
            if region_chr == chrom:
                if start >= int(region_start) and end <= int(region_end):
                    return sv_type_dict[region]
        else:
            return 'normal'

def parse_BND_end(text):
    pattern = r'chr6:(\d+)'
    matches = re.findall(pattern, text)
    print(matches)
    return [int(match) for match in matches]  

def vcf2seg():
    sv_vcf = pysam.VariantFile(args.sv_vcf)
    sv_file =  open('{}_seg.txt'.format(args.prefix), 'w')
    ori_sv_file = open('{}_ori_seg.txt'.format(args.prefix), 'w')
    # sv_file.write('chr\tpos1\tpos2\tsvtype\n')
    all_segs = []

    sv_dict = {}
    all_poses = []
    sv_juncs = []
    seg_direct = {}
    dup_tags = {}
    sv_pairs = []
    for rec in sv_vcf:
        print(rec.samples[0]['GT'])
        if rec.samples[0]['GT'][args.hap_idx] == 1:
            # chr1, pos1 = rec.chrom, str(rec.pos)
            # chr2, pos2 = rec.chrom, str(rec.stop)
            chrom = rec.chrom
            svtype = rec.info['SVTYPE']
            if svtype == 'INS':
                continue
            if ']' in rec.alts[0] or '[' in rec.alts[0]:
                pos1, pos2 = min(rec.pos, parse_BND_end(rec.alts[0])[0]), max(rec.pos, parse_BND_end(rec.alts[0])[0])
                # pos1, pos2 = rec.pos, parse_BND_end(rec.alts[0])[0]
            else:
                pos1, pos2 = min(rec.pos, rec.stop), max(rec.pos, rec.stop)
                # pos1, pos2 = rec.pos, rec.stop

            sv_dict[(chrom, pos1, pos2)] = svtype
            all_poses.append(pos1)
            all_poses.append(pos2)
            str1, str2 = '+', '+'
            if 'STRANDS' in rec.info:
                if rec.info['STRANDS'] == '++':
                    str1, str2 = '+', '-'
                elif rec.info['STRANDS'] == '--':
                    str1, str2 = '-', '+'
                elif rec.info['STRANDS'] == '+-':
                    str1, str2 = '+', '+'
                elif rec.info['STRANDS'] == '-+':
                    str1, str2 = '-', '-'
            sv_juncs.append([chrom, pos1, str1, pos2, str2, svtype])
            sv_pairs.append((pos1, pos2))
            # depth = str(rec.samples[0]['DV'])
            # svtype = rec.info['SVTYPE']
            # sv.append([chr1, pos1, str1, chr2, pos2, str2, depth, svtype])
            if svtype == 'DUP':
                dup_tags[(chrom, pos1, pos2)] = rec.info['STRANDS'][1]
    all_poses.append(mhc_start)
    all_poses.append(mhc_end)
    all_poses = sorted(list(set(all_poses)))
    print("sv_dict:", sv_dict)
    if len(sv_dict) == 0:
        all_segs.append([chrom_name, mhc_start, mhc_end, 'normal', mhc_end-mhc_start])
        return all_segs, sv_juncs, {}, {}, {}, []
    start_seg_dict = {}
    end_seg_dict = {}
    for start, end in zip(all_poses, all_poses[1:]):
        svtype = check_svtype(sv_dict, chrom, start, end)
        start_seg_dict[start] = [chrom, start, end, svtype, end-start]
        end_seg_dict[end] = [chrom, start, end, svtype, end-start]
        sv_file.write("{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, svtype, end-start))
        if svtype != 'DEL':
            all_segs.append([chrom, start, end, svtype, end-start])
    sv_file.close()

    print("all_pos:",all_poses)
    print("all_segs:", all_segs)
    print("start_seg_dict:", start_seg_dict)
    print("end_seg_dict:", end_seg_dict)
    print("dup_tags:", dup_tags)
    print("sv_pairs:", sv_pairs)


    return all_segs, sv_juncs, start_seg_dict, end_seg_dict, dup_tags, sv_pairs


def generate_sv_juncs_new2(all_segs):
    # assign sequence direction by sv type
    JUNC_TMP="JUNC {} {} {} {} {}"
    all_junc_strs = []
    for idx, seg in enumerate(all_segs):
        if seg[3] == "INV":
            all_segs[idx].append("-")
        else:
            all_segs[idx].append("+")
    CN_segs = []
    for idx, seg in enumerate(all_segs):
        if seg[3] == "DUP":
            CN_segs.append(seg)
            CN_segs.append(seg)
        else:
            CN_segs.append(seg)
    #generate junctions
    for idx, seg in enumerate(CN_segs):
        if idx == len(CN_segs)-1:
            break
        seg1_str = "{}-{}".format(seg[1], seg[2])
        seg2_str = "{}-{}".format(CN_segs[idx+1][1], CN_segs[idx+1][2])
        print(idx, seg1_str, seg2_str)
        sv_junc = JUNC_TMP.format(seg1_str, seg[5], seg2_str, CN_segs[idx+1][5], round(dp/2))
        all_junc_strs.append(sv_junc)
    return all_junc_strs

def generate_seg_fasta(segs_list):
    out_fa=open(args.prefix+"_seg.fa", 'w')
    ref_fasta = pysam.FastaFile(args.ref_fasta)
    for seg in segs_list:
        chrom, start, end, svtype, length = seg
        seq = ref_fasta.fetch(chrom, start, end)
        region="{}:{}-{}".format(chrom, start, end)
        out_fa.write(">{}\n{}\n".format(region, seq))


def get_depth(bam_file, region_chromosome, region_start, region_end):
    # for key, arr in pos.items():
    # for n in range(1, len(arr)):
    cnt = bam_file.count_coverage(region_chromosome, region_start, region_end, quality_threshold = 0)
    posDepth = []
    for i in range(len(cnt[0])):
        temp = 0
        for j in range(4):
            temp += cnt[j][i]
        posDepth.append(temp)
    # name = key+':'+str(arr[n-1])+'-'+str(arr[n])
    avg_dp = sum(posDepth)/len(posDepth)
    return avg_dp

def write_graph(final_juncs, final_segs):
    out_file = open(args.prefix+"_graph.txt", 'w')
    for seg in final_segs:
        out_file.write(seg+'\n')
    for junc in final_juncs:
        out_file.write(junc+'\n')
    out_file.close()

def generate_seg_cn(all_segs, start_seg_dict, end_seg_dict, dup_tags):
    ori_bam_f = pysam.AlignmentFile(args.ori_bam)
    hap0_bam_f = pysam.AlignmentFile(args.hap0_bam)
    hap1_bam_f = pysam.AlignmentFile(args.hap1_bam)
    final_segs = []
    SEG_TMP = "SEG {}-{} {} {}"
    for seg in all_segs:
        print("seg:", seg)
        if seg[3] == "DEL":
            seg_str = SEG_TMP.format(seg[1], seg[2], 0, 0)
            final_segs.append(seg_str)
        elif seg[3] == "DUP":
            seg_str = SEG_TMP.format(seg[1], seg[2], dp, 2)
            final_segs.append(seg_str)
        elif seg[3] == "INV":
            seg_str = SEG_TMP.format(seg[1], seg[2], dp/2, 1)
            final_segs.append(seg_str)
        elif seg[3] == "INS":
            seg_str = SEG_TMP.format(seg[1], seg[2], dp/2, 1)
            final_segs.append(seg_str)
        elif seg[3] == "TRA":
            seg_str = SEG_TMP.format(seg[1], seg[2], dp/2, 1)
            final_segs.append(seg_str)
        elif seg[3] == "normal":
            seg_str = SEG_TMP.format(seg[1], seg[2], dp/2, 1)
            final_segs.append(seg_str)
        else:
            ori_dp = get_depth(ori_bam_f, seg[0], seg[1], seg[2])
            ori_cn = round(ori_dp/(0.5*dp))
            print(seg[0], seg[1], seg[2], ori_dp, ori_cn)
            h0_dp = get_depth(hap0_bam_f, seg[0], seg[1], seg[2])
            h1_dp = get_depth(hap1_bam_f, seg[0], seg[1], seg[2])
            if h0_dp + h1_dp == 0:
                h0_cn = round(ori_cn/2) if round(ori_cn/2) > 1 else 1
                h1_cn = round(ori_cn/2) if round(ori_cn/2) > 1 else 1
            else:
                h0_h1_ratio = h0_dp/(h0_dp+h1_dp)
                if ori_cn ==1 and round(h0_h1_ratio) == 1:
                    if h0_dp <= h1_dp:
                        h0_cn = 0
                        h1_cn = 1
                    else:
                        h0_cn = 1
                        h1_cn = 0
                else:
                    h0_cn = round(ori_cn*h0_h1_ratio)
                    h1_cn = ori_cn - h0_cn
            if h0_cn == 0:
                h0_cn = 1
            if h1_cn == 0:
                h1_cn = 1
            seg_str = SEG_TMP.format(seg[1], seg[2], [h0_dp,h1_dp][args.hap_idx], [h0_cn,h1_cn][args.hap_idx])
            final_segs.append(seg_str)
    return final_segs


def main():
    all_segs, sv_juncs, start_seg_dict, end_seg_dict, dup_tags, sv_pair = vcf2seg()
    print(sv_juncs)
    generate_seg_fasta(all_segs)
    print("sv juncs:", sv_juncs)
    if len(sv_juncs) == 0:
        with open(args.prefix+"_nosv.flag", 'w') as out_file:
            out_file.write("no sv\n")
        exit()
    final_juncs =  generate_sv_juncs_new2(all_segs)
    print(final_juncs)
    final_segs = generate_seg_cn(all_segs, start_seg_dict, end_seg_dict, dup_tags)
    write_graph(final_juncs, final_segs)
    print("dp:", dp)
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Extract VCF records within a specified chromosome coordinate range.')
    parser.add_argument('-sv', '--sv_vcf', required=False, help='Output bed1 file')
    parser.add_argument('-p', '--prefix', required=False, help='Output bed1 file')
    parser.add_argument('-r', '--ref_fasta', required=False, help='Output bed1 file')
    parser.add_argument('-hid', '--hap_idx', required=False, type=int, help='Output bed1 file')
    parser.add_argument('-ob', '--ori_bam', required=False, help='Output bed1 file')
    parser.add_argument('-h0b', '--hap0_bam', required=False, help='Output bed1 file')
    parser.add_argument('-h1b', '--hap1_bam', required=False, help='Output bed1 file')
    parser.set_defaults(float=False)
    args = parser.parse_args()


        
    # get length of ref_fasta use pysam.FastaFile
    ref_fasta = pysam.FastaFile(args.ref_fasta)
    mhc_start = 0
    # get end
    chrom_name=ref_fasta.references[0]
    # TypeError: get_reference_length() takes exactly 1 positional argument (0 given)
    mhc_end = ref_fasta.get_reference_length(chrom_name)

    dp=bam_average_depth(args.ori_bam)


    main()