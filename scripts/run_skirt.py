import subprocess
import shutil
import os


def run_skirt_pre_directly(fasta_path, output_path, python_bin):
    # 确保输出目录存在
    os.makedirs(output_path, exist_ok=True)

    # SKIRT 主目录
    script_dir = os.path.dirname(os.path.realpath(__file__))
    skirt_wd = os.path.realpath(os.path.join(script_dir, "SKIRT"))

    # 引用文件路径
    nuc_query = f"{skirt_wd}/IPDKIR/kir_nuc.fasta"
    gen_query = f"{skirt_wd}/IPDKIR/kir_gen.fasta"
    eds1_query = f"{skirt_wd}/IPDKIR/fasta/KIR3DS1_nuc.fasta"
    zds4_fusion = f"{skirt_wd}/IPDKIR/fasta/KIR2DS4-00101e124567_3DL1-03501e89_nuc.fasta"
    zdl3_fusion = f"{skirt_wd}/IPDKIR/fasta/KIR2DL3-00101e1245_2DP1-00201e6789_nuc.fasta"
    novel_allele = f"{skirt_wd}/novel.fa"

    # 输出前缀
    fastafile = os.path.basename(fasta_path)
    outputname = os.path.splitext(fastafile)[0]
    output_prefix = os.path.join(output_path, outputname)
    paf_path = f"{output_prefix}.paf"

    # 创建 novel.fa
    with open(novel_allele, "wb") as out_f:
        for path in [zds4_fusion, zdl3_fusion, nuc_query]:
            with open(path, "rb") as src:
                shutil.copyfileobj(src, out_f)

    # minimap2 运行
    print("[INFO] Running minimap2 alignments...")
    with open(paf_path, "w") as paf_out, open(f"{paf_path}.err", "w") as paf_err:
        subprocess.run([
            "minimap2", "-cx", "splice:hq", "-G16k", "-y", "--cs", "-t32", "-2", fasta_path, novel_allele
        ], stdout=paf_out, stderr=paf_err)

        subprocess.run([
            "minimap2", "-cx", "splice:hq", "-G16k", "-k8", "--end-seed-pen", "5", "-y", "--cs", "-t32", "-2", fasta_path, eds1_query
        ], stdout=paf_out, stderr=paf_err)

    gen_paf_path = f"{output_prefix}.gen.paf"
    with open(gen_paf_path, "w") as gen_out, open(f"{gen_paf_path}.err", "w") as gen_err:
        subprocess.run([
            "minimap2", "-cx", "splice:hq", "-G16k", "-y", "--cs", "-t32", "-2", fasta_path, gen_query
        ], stdout=gen_out, stderr=gen_err)
    return paf_path, gen_paf_path

    