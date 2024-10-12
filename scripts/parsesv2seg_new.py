import argparse
import pysam
import subprocess, argparse, os, re
import numpy as np

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
        all_segs.append(['chr6', mhc_start, mhc_end, 'normal', mhc_end-mhc_start])
        return all_segs, sv_juncs, {}, {}, {}
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

    return all_segs, sv_juncs, start_seg_dict, end_seg_dict, dup_tags, sv_pairs

def run_command(command):
    # python 3.7 >
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        shell=True
    )
    return result


def check_dup_tag(dup_tags, chrom, start, end):
    all_keys = dup_tags.keys()
    for key in all_keys:
        if chrom == key[0]:
            if start >= key[1] and end <= key[2]:
                return dup_tags[key]

def in_dup_region(start, end, dup_regions):
    for region in dup_regions:
        if start >= region[1] and end <= region[2]:
            return True
    return False
    

def generate_sv_juncs_new(all_segs, sv_juncs, start_seg_dict, end_seg_dict, dup_tags, sv_pairs):
    dup_regions = list(dup_tags.keys())
    JUNC_TMP="JUNC {} {} {} {} {}"
    all_junc_strs = []
    visited_junc_pairs = set()
    for sv_jun in sv_juncs:
        chrom, pos1, sign1, pos2, sign2, svtype = sv_jun[0], sv_jun[1], sv_jun[2], sv_jun[3], sv_jun[4], sv_jun[5]
        seg1, seg2 = end_seg_dict[pos1], start_seg_dict[pos2]
        if seg1 not in all_segs or seg2 not in all_segs:
            continue
        print("svtype:", svtype)
        print(pos1, pos2)
        if (chrom, pos1, pos2) in dup_tags:
            seg1, seg2 = start_seg_dict[pos1], end_seg_dict[pos2]
            print("raw seg1:", seg1)
            print("raw seg2:", seg2)
            if seg1 != seg2:
                junc_str = JUNC_TMP.format("{}-{}".format(seg1[1], seg1[2]), "+", "{}-{}".format(seg1[1], seg1[2]), "+", round(dp/2))
                all_junc_strs.append(junc_str)
                visited_junc_pairs.add(("{}-{}".format(seg1[1], seg1[2]), "{}-{}".format(seg1[1], seg1[2])))
                junction = JUNC_TMP.format("{}-{}".format(seg2[1], seg2[2]), "+", "{}-{}".format(seg2[1], seg2[2]), "+", round(dp/2))
                all_junc_strs.append(junction)
                visited_junc_pairs.add(("{}-{}".format(seg2[1], seg2[2]), "{}-{}".format(seg2[1], seg2[2])))
            else:
                junc_str = JUNC_TMP.format("{}-{}".format(seg1[1], seg1[2]), "+", "{}-{}".format(seg1[1], seg1[2]), "+", round(dp/2))
                all_junc_strs.append(junc_str)
                visited_junc_pairs.add(("{}-{}".format(seg1[1], seg1[2]), "{}-{}".format(seg1[1], seg1[2])))

            # junc_str = JUNC_TMP.format("{}-{}".format(seg1[1], seg1[2]), '+', "{}-{}".format(seg2[1], seg2[2]), '+', round(dp/2))
            # all_junc_strs.append(junc_str)
            # visited_junc_pairs.add(("{}-{}".format(seg1[1], seg1[2]), "{}-{}".format(seg2[1], seg2[2])))
        elif svtype == 'INV':
            print("inv1 seg1:", seg1)
            print("inv1 seg2:", seg2)
            seg1, seg2 = end_seg_dict[pos1], end_seg_dict[pos2]
            junc_str = JUNC_TMP.format("{}-{}".format(seg1[1], seg1[2]), "+", "{}-{}".format(seg2[1], seg2[2]), "-", round(dp/2))
            all_junc_strs.append(junc_str)
            visited_junc_pairs.add(("{}-{}".format(seg1[1], seg1[2]), "{}-{}".format(seg2[1], seg2[2])))
            seg1, seg2 = start_seg_dict[pos1], start_seg_dict[pos2]
            print("inv2 seg1:", seg1)
            print("inv2 seg2:", seg2)
            junc_str = JUNC_TMP.format("{}-{}".format(seg1[1], seg1[2]), "-", "{}-{}".format(seg2[1], seg2[2]), "+", round(dp/2))
            all_junc_strs.append(junc_str)
            visited_junc_pairs.add(("{}-{}".format(seg1[1], seg1[2]), "{}-{}".format(seg2[1], seg2[2])))
        elif svtype == 'DUP':
            seg1, seg2 = end_seg_dict[pos2], start_seg_dict[pos1]
            print("dup1 seg1:", seg1)
            print("dup1 seg2:", seg2)
            junc_str = JUNC_TMP.format("{}-{}".format(seg1[1], seg1[2]), "+", "{}-{}".format(seg2[1], seg2[2]), "+", round(dp/2))
            all_junc_strs.append(junc_str)
            visited_junc_pairs.add(("{}-{}".format(seg1[1], seg1[2]), "{}-{}".format(seg2[1], seg2[2])))
            seg1, seg2 = start_seg_dict[pos1], end_seg_dict[pos2]
            print("dup2 seg1:", seg1)
            print("dup2 seg2:", seg2)
            junc_str = JUNC_TMP.format("{}-{}".format(seg1[1], seg1[2]), "+", "{}-{}".format(seg2[1], seg2[2]), "+", round(dp/2))
            all_junc_strs.append(junc_str)
            visited_junc_pairs.add(("{}-{}".format(seg1[1], seg1[2]), "{}-{}".format(seg2[1], seg2[2])))
        else:
            print("else seg1:", seg1)
            print("else seg2:", seg2)
            junc_str = JUNC_TMP.format("{}-{}".format(seg1[1], seg1[2]), sign1, "{}-{}".format(seg2[1], seg2[2]), sign2, round(dp/2))
            all_junc_strs.append(junc_str)
            visited_junc_pairs.add(("{}-{}".format(seg1[1], seg1[2]), "{}-{}".format(seg2[1], seg2[2])))

    # remove DEL segs from all_segs, store res in nodel_segs
    nodel_segs = []
    for seg in all_segs:
        if seg[3] != 'DEL':
            nodel_segs.append(seg)
    # check missed normal junctions
    visited_loss_dup_juncs = set() # for those normal junctions in dup regions
    for idx, seg in enumerate(nodel_segs):
        if idx == len(nodel_segs)-1:
            break
        pos1, pos2 = seg[2], nodel_segs[idx+1][1]
        # if (pos1, pos2) in visited_junc_pairs:
        #     continue
        print(seg, nodel_segs[idx+1])
        seg1_str = "{}-{}".format(seg[1], seg[2])
        seg2_str = "{}-{}".format(nodel_segs[idx+1][1], nodel_segs[idx+1][2])
        print("init seg1:", seg1_str)
        print("init seg2:", seg2_str)
        if (seg1_str, seg2_str) in visited_junc_pairs:
            continue
        print("unvisited seg1:", seg1_str)
        print("unvisited seg2:", seg2_str)
        if in_dup_region(pos1, pos2, dup_regions):
            print("in_dup")
            print("seg1:", seg1_str)
            print("seg2:", seg2_str)
        if in_dup_region(pos1, pos2, dup_regions) and (pos1, pos2) not in sv_pairs:
            print("in_dup and not in sv_pairs")
            print("seg1:", seg1_str)
            print("seg2:", seg2_str)
            junc_str = JUNC_TMP.format(seg1_str, '+', seg2_str, '+', round(dp))
            all_junc_strs.append(junc_str)
        else:
            # print("add normal")
            junc_str = JUNC_TMP.format(seg1_str, '+', seg2_str, '+', round(dp/2))
        all_junc_strs.append(junc_str)
        print(idx, seg1_str, seg2_str)
    #for self loop
    for idx, seg in enumerate(nodel_segs):
        if (seg[1], seg[2]) not in sv_pairs and seg[3]=="DUP":
            seg1_str = "{}-{}".format(seg[1], seg[2])
            seg2_str = "{}-{}".format(seg[1], seg[2])
            junc_str = JUNC_TMP.format(seg1_str, '+', seg2_str, '+', round(dp/2))
            all_junc_strs.append(junc_str)

    return all_junc_strs


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


def generate_sv_juncs(all_segs, sv_juncs, start_seg_dict, end_seg_dict, dup_tags):
    # generate sv juncs first
    final_sv_juncs = []
    # remove DEL segs from all_segs, store res in nodel_segs
    nodel_segs = []
    for seg in all_segs:
        if seg[3] != 'DEL':
            nodel_segs.append(seg)
    print("nodel seg:", nodel_segs)
    JUNC_TMP="JUNC {} {} {} {} {}"
    for idx, seg in enumerate(nodel_segs):
        if idx == len(nodel_segs)-1:
            break
        seg1_str = "{}-{}".format(seg[1], seg[2])
        seg2_str = "{}-{}".format(nodel_segs[idx+1][1], nodel_segs[idx+1][2])
        print(idx, seg1_str, seg2_str)
        if seg[3] == 'normal' and nodel_segs[idx+1][3] == 'normal':
            final_junc = JUNC_TMP.format(seg1_str, '+', seg2_str, '+', round(dp/2))
            final_sv_juncs.append(final_junc)
        elif seg[3] == 'normal' and nodel_segs[idx+1][3] == 'INV':
            final_junc = JUNC_TMP.format(seg1_str, '+', seg2_str, '-', round(dp/2))
            final_sv_juncs.append(final_junc)
        elif seg[3] == 'normal' and nodel_segs[idx+1][3] == 'DUP':
            dup_start = nodel_segs[idx+1][1]
            dup_end = nodel_segs[idx+1][2]
            dup_tag = check_dup_tag(dup_tags, seg[0], dup_start, dup_end)
            if dup_tag == "+":
                final_junc = JUNC_TMP.format(seg1_str, '+', seg2_str, '+', round(dp/2))
                final_sv_juncs.append(final_junc)
                final_junc = JUNC_TMP.format(seg2_str, '+', seg2_str, '+', round(dp/2))
                final_sv_juncs.append(final_junc)                
            else:
                final_junc = JUNC_TMP.format(seg1_str, '+', seg2_str, '+', round(dp/2))
                final_sv_juncs.append(final_junc)
                final_junc = JUNC_TMP.format(seg2_str, '+', seg2_str, '-', round(dp/2))
                final_sv_juncs.append(final_junc)
        elif seg[3] == 'INV' and nodel_segs[idx+1][3] == 'normal':
            final_junc = JUNC_TMP.format(seg1_str, '-', seg2_str, '+', round(dp/2))
            final_sv_juncs.append(final_junc)
        elif seg[3] == 'INV' and nodel_segs[idx+1][3] == 'INV':
            final_junc = JUNC_TMP.format(seg1_str, '-', seg2_str, '-', round(dp/2))
            final_sv_juncs.append(final_junc)
        elif seg[3] == 'INV' and nodel_segs[idx+1][3] == 'DUP':
            dup_start = nodel_segs[idx+1][1]
            dup_end = nodel_segs[idx+1][2]
            dup_tag = check_dup_tag(dup_tags, seg[0], dup_start, dup_end)
            if dup_tag == "+":
                final_junc = JUNC_TMP.format(seg1_str, '-', seg2_str, '+', round(dp/2))
                final_sv_juncs.append(final_junc)
                final_junc = JUNC_TMP.format(seg2_str, '+', seg2_str, '+', round(dp/2))
                final_sv_juncs.append(final_junc)
            else:
                final_junc = JUNC_TMP.format(seg1_str, '-', seg2_str, '+', round(dp/2))
                final_sv_juncs.append(final_junc)
                final_junc = JUNC_TMP.format(seg2_str, '+', seg2_str, '-', round(dp/2))
                final_sv_juncs.append(final_junc)
        elif seg[3] == 'DUP' and nodel_segs[idx+1][3] == 'normal':
            dup_start = seg[1]
            dup_end = seg[2]
            dup_tag = check_dup_tag(dup_tags, seg[0], dup_start, dup_end)
            if dup_tag == "+":
                final_junc = JUNC_TMP.format(seg1_str, '+', seg2_str, '+', round(dp/2))
                final_sv_juncs.append(final_junc)
            else:
                final_junc = JUNC_TMP.format(seg1_str, '-', seg2_str, '+', round(dp/2))
                final_sv_juncs.append(final_junc)
        elif seg[3] == 'DUP' and nodel_segs[idx+1][3] == 'INV':
            dup_start = seg[1]
            dup_end = seg[2]
            dup_tag = check_dup_tag(dup_tags, seg[0], dup_start, dup_end)
            if dup_tag == "+":
                final_junc = JUNC_TMP.format(seg1_str, '+', seg2_str, '-', round(dp/2))
                final_sv_juncs.append(final_junc)
            else:
                final_junc = JUNC_TMP.format(seg1_str, '-', seg2_str, '-', round(dp/2))
                final_sv_juncs.append(final_junc)
        elif seg[3] == 'DUP' and nodel_segs[idx+1][3] == 'DUP':
            seg1_dup_start = seg[1]
            seg1_dup_end = seg[2]
            seg1_dup_tag = check_dup_tag(dup_tags, seg[0], seg1_dup_start, seg1_dup_end)
            seg2_dup_start = nodel_segs[idx+1][1]
            seg2_dup_end = nodel_segs[idx+1][2]
            seg2_dup_tag = check_dup_tag(dup_tags, seg[0], seg2_dup_start, seg2_dup_end)
            if seg1_dup_tag == "+" and seg2_dup_tag == "+":
                final_junc = JUNC_TMP.format(seg1_str, '+', seg2_str, '+', round(dp/2))
                final_sv_juncs.append(final_junc)
                final_junc = JUNC_TMP.format(seg2_str, '+', seg2_str, '+', round(dp/2))
                final_sv_juncs.append(final_junc)
            elif seg1_dup_tag == "+" and seg2_dup_tag == "-":
                final_junc = JUNC_TMP.format(seg1_str, '+', seg2_str, '+', round(dp/2))
                final_sv_juncs.append(final_junc)
                final_junc = JUNC_TMP.format(seg2_str, '+', seg2_str, '-', round(dp/2))
                final_sv_juncs.append(final_junc)
            elif seg1_dup_tag == "-" and seg2_dup_tag == "+":
                final_junc = JUNC_TMP.format(seg1_str, '-', seg2_str, '+', round(dp/2))
                final_sv_juncs.append(final_junc)
                final_junc = JUNC_TMP.format(seg2_str, '+', seg2_str, '+', round(dp/2))
                final_sv_juncs.append(final_junc)
            elif seg1_dup_tag == "-" and seg2_dup_tag == "-":
                final_junc = JUNC_TMP.format(seg1_str, '-', seg2_str, '+', round(dp/2))
                final_sv_juncs.append(final_junc)
                final_junc = JUNC_TMP.format(seg2_str, '+', seg2_str, '-', round(dp/2))
                final_sv_juncs.append(final_junc)

    return final_sv_juncs


            



def generate_seg_fasta(segs_list):
    out_fa=open(args.prefix+"_seg.fa", 'w')
    for seg in segs_list:
        region="{}:{}-{}".format(seg[0], seg[1], seg[2])
        command = "bcftools view -r {}  {} -Oz -o {}".format(region, args.snp_vcf, args.prefix+'_tmp.vcf.gz')
        print(command)
        os.system(command)
        command = "tabix {} -f".format(args.prefix+'_tmp.vcf.gz')
        print(command)
        os.system(command)
        # samtools faidx ref.fa 8:11870-11890 | bcftools consensus in.vcf.gz > out.fa
        command = "samtools faidx {} {} | bcftools consensus {}_tmp.vcf.gz -H {}".format(args.ref_fasta, region, args.prefix, args.hap_idx+1)
        print(command)
        result = run_command(command)
        # print(result.stdout)
        # seq = result.stdout.split('\n')
        # print(seq)
        print(result.stdout)
        out_fa.write('{}'.format(result.stdout))
        # print(result.stdout)
        # print(result.stderr)

def cal_copy_number():
    f = open(args.depth_file, 'r')
    depth_dict = {}
    all_depth = []
    for line in f:
        array = line.strip().split()
        chrom = array[0]
        pos = int(array[1])
        depth = int(array[2])
        if chrom not in depth_dict:
            depth_dict[chrom] = []
        depth_dict[chrom].append(depth)
        all_depth.append(depth)
    f.close()
    median_depth = np.median(all_depth)
    chrom_copy = {} 
    chrom_depth = {}
    for chrom in depth_dict:
        chrom_depth[chrom] = np.median(depth_dict[chrom])
        chrom_copy[chrom] = round(float(np.median(depth_dict[chrom]))/median_depth)
    return chrom_copy, median_depth, chrom_depth

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

def write_graph(final_juncs, final_segs):
    out_file = open(args.prefix+"_graph.txt", 'w')
    for seg in final_segs:
        out_file.write(seg+'\n')
    for junc in final_juncs:
        out_file.write(junc+'\n')
    out_file.close()

def main():
    all_segs, sv_juncs, start_seg_dict, end_seg_dict, dup_tags, sv_pair = vcf2seg()
    print(sv_juncs)
        
    generate_seg_fasta(all_segs)
    print("sv juncs:", sv_juncs)
    if len(sv_juncs) == 0:
        with open(args.prefix+"_nosv.flag", 'w') as out_file:
            out_file.write("no sv\n")
        exit()

    # generate_sv_juncs()
    # chrom_copy, median_depth, chrom_depth = cal_copy_number()
    # print("chrom_copy:", chrom_copy)
    # print("median_depth:", median_depth)
    # print("chrom_depth:", chrom_depth)
    # print("sv_juncs:", sv_juncs)
    # maybe integrate fasta here
    
    final_juncs =  generate_sv_juncs_new2(all_segs)
    print(final_juncs)
    final_segs = generate_seg_cn(all_segs, start_seg_dict, end_seg_dict, dup_tags)
    write_graph(final_juncs, final_segs)
    print("dp:", dp)
    

if __name__ == '__main__':
    mhc_start = 28510120
    mhc_end = 33480577
    parser = argparse.ArgumentParser(description='Extract VCF records within a specified chromosome coordinate range.')
    parser.add_argument('-sv', '--sv_vcf', required=False, help='Output bed1 file')
    parser.add_argument('-snp', '--snp_vcf', required=False, help='Output bed1 file')
    parser.add_argument('-p', '--prefix', required=False, help='Output bed1 file')
    parser.add_argument('-r', '--ref_fasta', required=False, help='Output bed1 file')
    parser.add_argument('-hid', '--hap_idx', required=False, type=int, help='Output bed1 file')
    parser.add_argument('-ob', '--ori_bam', required=False, help='Output bed1 file')
    parser.add_argument('-h0b', '--hap0_bam', required=False, help='Output bed1 file')
    parser.add_argument('-h1b', '--hap1_bam', required=False, help='Output bed1 file')

    
    # parser.add_argument('-sv', '--sv_vcf', required=False, help='Output bed1 file')
    # parser.add_argument('-o', '--out_vcf', required=False, help='Output bed1 file')

    parser.set_defaults(float=False)
    args = parser.parse_args()

    dp=bam_average_depth(args.ori_bam)


    main()