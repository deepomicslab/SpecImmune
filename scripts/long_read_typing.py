"""
Function 2 : typing with only long reads
"""

import sys
import os
import pysam
import gzip
import argparse
from collections import defaultdict
import subprocess
from Bio import SeqIO


from downsample_bam import downsample_func
from read_objects import My_read, My_locus, Read_bin
from determine_gene import get_focus_gene, get_folder_list
from db_objects import My_db
from alignment_modules import Read_Type
from folder_objects import My_folder


class Pacbio_Binning():

    def __init__(self):
        self.db = my_db.full_db
        if args["seq_tech"] == "rna":
            self.cds_db=my_db.full_cds_db

        self.sam = f"""{parameter.outdir}/{parameter.sample}.db.bam"""

        self.bamfile = pysam.AlignmentFile(self.sam, 'rb')   
        self.assign_file = f"{parameter.outdir}/{parameter.sample}.assign.txt"


class Parameters():

    def __init__(self):

        self.sample = args["n"]
        self.raw_fq = args["r"]
        outdir = args["o"]
        self.population = args["p"]
        self.threads = args["j"]
        self.bin = "%s/../bin/"%(sys.path[0])      
        # self.outdir = "%s/%s/"%(outdir, self.sample)
        self.outdir = my_folder.step2_genes_dir
        # self.whole_dir = "%s/whole/"%(sys.path[0])



class Fasta():

    def alignment(self, gene):
        ### call and phase snps
        cmd = """
        sample=%s
        bin=%s
        outdir=%s
        hla=%s
        hla_ref=%s
        minimap2 -t %s %s -a $hla_ref $outdir/$hla.%s.fq.gz | samtools view -bS -F 0x804 -| samtools sort - >$outdir/$hla.bam
        samtools index $outdir/$hla.bam
        samtools depth -d 1000000 -aa $outdir/$hla.bam >$outdir/$hla.depth
        """%(parameter.sample, parameter.bin, parameter.outdir, gene, my_db.get_gene_alleles(gene), parameter.threads, minimap_para, args["a"])
        os.system(cmd)

        max_depth = args["max_depth"]
        seed = args["seed"]
        input_bam = f"{parameter.outdir}/{gene}.bam"
        input_depth = f"{parameter.outdir}/{gene}.depth"
        output_bam = f"{parameter.outdir}/{gene}.downsample.bam"
        output_depth = f"{parameter.outdir}/{gene}.downsample.depth"
        downsample_ratio = downsample_func(input_bam, output_bam, input_depth, output_depth, max_depth, seed)
        print (f"downsample ratio is {downsample_ratio} for {gene}", flush=True)
        if downsample_ratio < 1:
            os.system(f"rm {input_bam}")
            os.system(f"rm {input_depth}")
            os.system(f"mv {output_bam} {input_bam}")
            os.system(f"mv {output_depth} {input_depth}")

    def vcf2fasta(self, gene):
        ### map and downsample alignment
        self.alignment(gene)
        ### call and phase snps
        awk_script = '{{sum+=$3}} END {{ if (NR>0) print sum/NR; else print 0; }}'
        hla_ref=my_db.get_gene_alleles(gene)
        bam=f"{parameter.outdir}/{gene}.bam"
        # print("bam for vars:", bam, parameter.outdir)
        depth_file=f"{parameter.outdir}/{gene}.depth"
        # print ("xxx", bam)
        # print ("xxx", hla_ref)
        # call snp
        # mask_bed=f"{parameter.outdir}/low_depth.bed"
        ovcf=f"{parameter.outdir}/{parameter.sample}.{gene}.dv.vcf"
        ogvcf=f"{parameter.outdir}/{parameter.sample}.{gene}.dv.g.vcf"
        call_phase_cmd_longshot= f"""
            set_dp={args["k"]}
            avg_depth=$(samtools depth -a {bam} | awk '{awk_script}')
            if (( $(echo "$avg_depth < 5" | bc -l) )) && [ {args['y']} = "pacbio" ]; then
                set_dp=0
            else
                set_dp=5
            fi
            python3 {sys.path[0]}/mask_low_depth_region.py -f False -c {depth_file} -b {my_folder.step2_genes_dir}/{gene}.low_depth.bed -w 20 -d 5
            
            longshot -F -c 2 -C 100000 -P {args["strand_bias_pvalue_cutoff"]} -r {interval_dict[gene]} --bam {bam} --ref {hla_ref} --out {my_folder.step2_genes_dir}/{parameter.sample}.{gene}.longshot.vcf 
            bgzip -f {my_folder.step2_genes_dir}/{parameter.sample}.{gene}.longshot.vcf
            tabix -f {my_folder.step2_genes_dir}/{parameter.sample}.{gene}.longshot.vcf.gz
            zcat {my_folder.step2_genes_dir}/{parameter.sample}.{gene}.longshot.vcf.gz >{my_folder.step2_genes_dir}/{parameter.sample}.{gene}.phased.vcf          
            bgzip -f {my_folder.step2_genes_dir}/{parameter.sample}.{gene}.phased.vcf
            tabix -f {my_folder.step2_genes_dir}/{parameter.sample}.{gene}.phased.vcf.gz

            # bash {sys.path[0]}/run_dv.sh {hla_ref} {bam} {ovcf} {ogvcf} {parameter.threads} {interval_dict[gene]} 
        """
        # os.system(call_phase_cmd_longshot)
        call_phase_cmd_dv= f"""
            set_dp={args["k"]}
            avg_depth=$(samtools depth -a {bam} | awk '{awk_script}')
            if (( $(echo "$avg_depth < 5" | bc -l) )) && [ {args['y']} = "pacbio" ]; then
                set_dp=0
            else
                set_dp=5
            fi
            python3 {sys.path[0]}/mask_low_depth_region.py -f False -c {depth_file} -b {my_folder.step2_genes_dir}/{gene}.low_depth.bed -w 20 -d 5


            bash {sys.path[0]}/run_dv.sh {hla_ref} {bam} {ovcf} {ogvcf} {parameter.threads} {interval_dict[gene]} 
        """
        # os.system(call_phase_cmd_dv)
        if args["snv_tool"] == "longshot":
            os.system(call_phase_cmd_longshot)
        elif args["snv_tool"] == "deepvariant":
            if args["dv_sif"]:
                os.system(call_phase_cmd_dv)
            else:
                print ("Please provide the deepvariant sif file !")
                sys.exit(1)

        print ("call snp done")
        # call sv & phase snv-sv
        gene_work_dir=f"{parameter.outdir}/{gene}_work"
        if not os.path.exists(gene_work_dir):
            os.makedirs(gene_work_dir)
        sv_cmd_longshot = f"""
            bash {sys.path[0]}/refine_haplotype_pipe.sh {bam} {hla_ref} {gene} {interval_dict[gene]} {mask_bed} {gene_work_dir} {parameter.threads} {parameter.sample} 
        """
        sv_cmd_dv = f"""
            bash {sys.path[0]}/refine_haplotype_dv_pipe.sh {bam} {hla_ref} {gene} {interval_dict[gene]} {mask_bed} {gene_work_dir} {parameter.threads} {parameter.sample} 
        """
        # os.system(sv_cmd)
        if args["snv_tool"] == "longshot":
            os.system(sv_cmd_longshot)
        elif args["snv_tool"] == "deepvariant":
            os.system(sv_cmd_dv)
        print ("call sv done")
        ## reconstruct HLA sequence based on the phased snps & sv
        for index in range(2):
            order = f"""
            echo ">{gene}_{index}" >{parameter.outdir}/hla.allele.{index+1}.{gene}.fasta
            cat {gene_work_dir}/{gene}.{index+1}.raw.fa|grep -v ">" >>{parameter.outdir}/hla.allele.{index+1}.{gene}.fasta    
            samtools faidx {parameter.outdir}/hla.allele.{index+1}.{gene}.fasta    
            """
            os.system(order)

    def is_hom(self, gene):
        # if elf.individual_ref_dir + f"/{gene}/{gene}.{ref_idx+1}.fasta" and elf.individual_ref_dir + f"/{gene}/{gene}.{ref_idx+2}.fasta" exists
        if os.path.exists(my_db.get_gene_alleles_2ref(gene, 0)) and os.path.exists(my_db.get_gene_alleles_2ref(gene, 1)):
            return False
        return True

    def calculate_avg_depth(self, bam, awk_script):
        avg_depth_cmd = f"samtools depth  {bam} | awk '{awk_script}'"
        return float(subprocess.check_output(avg_depth_cmd, shell=True).strip())

    def run_command(self, cmd, description):
        # print(f"Executing {description}:\n{cmd}", flush=True)
        try:
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
        except subprocess.CalledProcessError as e:
            print(f"Error executing {description}: {e}", flush=True)
            raise

    def downsample_and_process(self, bam, output_bam, depth_file, output_depth, max_depth, seed, downsample_func, awk_script):
        downsample_ratio = downsample_func(bam, output_bam, depth_file, output_depth, max_depth, seed)
        if downsample_ratio < 1:
            os.remove(bam)
            os.remove(depth_file)
            os.rename(output_bam, bam)
            os.rename(output_depth, depth_file)
        avg_depth = self.calculate_avg_depth(bam, awk_script)
        return avg_depth, downsample_ratio



    def generate_sequence(self, gene, index, gene_work_dir, output_dir):
        if index is None:
            for idx in range(2):
                order = f"""
                echo ">{gene}_{idx}" > {output_dir}/{args['i']}.allele.{idx+1}.{gene}.fasta
                cat {gene_work_dir}/{gene}.{idx+1}.raw.fa | grep -v ">" >> {output_dir}/{args['i']}.allele.{idx+1}.{gene}.fasta
                samtools faidx {output_dir}/{args['i']}.allele.{idx+1}.{gene}.fasta
                """
                self.run_command(order, f"sequence generation for {gene}_{idx}")
        else:
            order = f"""
            echo ">{gene}_{index}" > {output_dir}/{args['i']}.allele.{index+1}.{gene}.fasta
            cat {gene_work_dir}/{gene}.{index+1}.raw.fa | grep -v ">" >> {output_dir}/{args['i']}.allele.{index+1}.{gene}.fasta
            samtools faidx {output_dir}/{args['i']}.allele.{index+1}.{gene}.fasta
            """
            self.run_command(order, f"sequence generation for {gene}_{index}")


    def remove_N_characters(self, gene):
        for idx in range(2):
            input_file = f"{my_folder.sequence_dir}/{args['i']}.allele.{idx+1}.{gene}.fasta"
            output_file = f"{my_folder.sequence_dir}/{args['i']}.allele.{idx+1}.{gene}.noN.fasta"
            with open(input_file, "r") as infile, open(output_file, "w") as outfile:
                sequences = SeqIO.parse(infile, "fasta")
                modified_sequences = []

                for seq_record in sequences:
                    seq_record.seq = seq_record.seq.replace("N", "")
                    modified_sequences.append(seq_record)
                
                SeqIO.write(modified_sequences, outfile, "fasta")

    def process_phase(self, gene, bam, depth_file, mask_bed, hla_ref, interval, set_dp, min_cov, args, parameter, index=None):
        sample = parameter.sample
        # outdir = parameter.outdir
        set_window = 10
        if args['seq_tech'] == "rna":
            set_window = 1


        suffix = f".{index}" if index is not None else ""
        call_phase_cmd = f"""
        samtools index {bam}
        python3 {sys.path[0]}/mask_low_depth_region.py -f False -c {depth_file} -b {my_folder.step2_genes_dir}/{gene}{suffix}.low_depth.bed -w {set_window} -d {int(set_dp)}
        
        longshot -F -c {min_cov} -C 100000 -P {args["strand_bias_pvalue_cutoff"]} -r {interval} --bam {bam} --ref {hla_ref} --out {my_folder.step2_genes_dir}/{sample}.{gene}{suffix}.longshot.vcf
        bgzip -f {my_folder.step2_genes_dir}/{sample}.{gene}{suffix}.longshot.vcf
        tabix -f {my_folder.step2_genes_dir}/{sample}.{gene}{suffix}.longshot.vcf.gz

        # samtools dict {hla_ref} > {hla_ref[:-6]}.dict
        # java -Xmx5g -jar {sys.path[0]}/../packages/GenomeAnalysisTK.jar -T HaplotypeCaller -R {hla_ref} -I {bam} -o {my_folder.step2_genes_dir}/{sample}.{gene}{suffix}.longshot.raw.vcf.gz  -L {interval}
        # whatshap phase -o {my_folder.step2_genes_dir}/{sample}.{gene}{suffix}.longshot.vcf.gz -r {hla_ref} --indels {my_folder.step2_genes_dir}/{sample}.{gene}{suffix}.longshot.raw.vcf.gz {bam}
        # tabix -f {my_folder.step2_genes_dir}/{sample}.{gene}{suffix}.longshot.vcf.gz
        

        zcat {my_folder.step2_genes_dir}/{sample}.{gene}{suffix}.longshot.vcf.gz > {my_folder.step2_genes_dir}/{sample}.{gene}{suffix}.phased.vcf
        bgzip -f {my_folder.step2_genes_dir}/{sample}.{gene}{suffix}.phased.vcf
        tabix -f {my_folder.step2_genes_dir}/{sample}.{gene}{suffix}.phased.vcf.gz
        """
        self.run_command(call_phase_cmd, f"call phase for {gene}{suffix}")
        mask_bed = f"{my_folder.step2_genes_dir}/{gene}{suffix}.low_depth.bed"
        return mask_bed

    def process_sv(self, gene, bam, hla_ref, interval, mask_bed, gene_work_dir, parameter, index=None):
        sample = parameter.sample
        threads = parameter.threads

        suffix = f"_{index}" if index is not None else ""
        refine_script = f"{sys.path[0]}/refine_haplotype_pipe.sh" if index is None else f"{sys.path[0]}/refine_haplotype_2ref_pipe.sh"
        sv_cmd = f"""bash {refine_script} {bam} {hla_ref} {gene} {interval} {mask_bed} {gene_work_dir} {threads} {sample}"""
        seq_tag=""
        if args['seq_tech'] == "rna":
            seq_tag=args["RNA_type"]
        else:
            seq_tag=args["y"]
        if index is not None:
            sv_cmd = f"{sv_cmd} {index}"
        else:
            sv_cmd = f"{sv_cmd} {seq_tag}"
        self.run_command(sv_cmd, f"sv for {gene}{suffix}")

    def process_gene(self, gene, hla_ref, interval, bam, depth_file, mask_bed, set_dp, min_cov, args, parameter, awk_script, index=None):
        mask_bed = self.process_phase(gene, bam, depth_file, mask_bed, hla_ref, interval, set_dp, min_cov, args, parameter, index)
        gene_work_dir = f"{my_folder.step2_genes_dir}/{gene}_work" if index is None else f"{my_folder.step2_genes_dir}/{gene}_{index}_work"
        if not os.path.exists(gene_work_dir):
            os.makedirs(gene_work_dir)
        self.process_sv(gene, bam, hla_ref, interval, mask_bed, gene_work_dir, parameter, index)
        self.generate_sequence(gene, index, gene_work_dir, my_folder.sequence_dir)
        if args['seq_tech'] == "rna":
            self.remove_N_characters(gene)

    def vcf2fasta2ref(self, gene):
        awk_script = '{sum+=$3} END { if (NR>0) print sum/NR; else print 0; }'
        max_depth = args["max_depth"]
        seed = args["seed"]
        set_dp = args['k']

        if self.is_hom(gene) or args['seq_tech'] == "rna":
            print(f"Processing RNAseq {gene}", flush=True) if args['seq_tech'] == "rna" else \
            print(f"Processing homo {gene}", flush=True)
            bam = f"{my_folder.step2_genes_dir}/{gene}.bam"
            if not os.path.exists(bam):
                print(f"{bam} does not exist", flush=True)
                return
            hla_ref = my_db.get_gene_alleles_ref(gene)
            depth_file = f"{my_folder.step2_genes_dir}/{gene}.depth"
            mask_bed = f"{my_folder.step2_genes_dir}/low_depth.bed"
            output_bam = f"{my_folder.step2_genes_dir}/{gene}.downsample.bam"
            output_depth = f"{my_folder.step2_genes_dir}/{gene}.downsample.depth"
            avg_depth, downsample_ratio = self.downsample_and_process(bam, output_bam, depth_file, output_depth, max_depth, seed, downsample_func, awk_script)
            if avg_depth < 2 and args['y']=="nanopore":
                print(f"Reads coverage is too low for {bam}, skip it")
                return
            elif avg_depth<1:
                print(f"Reads coverage is too low for {bam}, skip it")
                return
            print(f"downsample ratio is {downsample_ratio} for {gene}", flush=True)
            if args['y'] == "pacbio-hifi":
                set_dp == 1
            elif args['y'] != "pacbio-hifi":
                set_dp = int(max(0.3 * avg_depth, 0)) if avg_depth >= 5 else 0
            elif args['seq_tech'] == "rna" and args['RNA_type'] == "traditional":
                set_dp = 1



            # set_dp = int(max(0.3 * avg_depth, 0)) if args['y'] != 'pacbio-hifi' and avg_depth >= 5 else 0
            if args['y'] == "nanopore" and set_dp==0:
                set_dp=1
            min_cov = 2 if args['y'] == 'pacbio-hifi' else set_dp
            interval = f"{interval_dict[gene]}"
            self.process_gene(gene, hla_ref, interval, bam, depth_file, mask_bed, set_dp, min_cov, args, parameter, awk_script)
        else:
            print(f"Processing hete {gene}", flush=True)
            for index in range(2):
                bam = f"{my_folder.step2_genes_dir}/{gene}.{index}.bam"
                if not os.path.exists(bam):
                    print(f"{bam} does not exist", flush=True)
                    continue

                hla_ref = my_db.get_gene_alleles_2ref(gene, index)
                depth_file = f"{my_folder.step2_genes_dir}/{gene}.{index}.depth"
                mask_bed = f"{my_folder.step2_genes_dir}/low_depth.bed"
                output_bam = f"{my_folder.step2_genes_dir}/{gene}.{index}.downsample.bam"
                output_depth = f"{my_folder.step2_genes_dir}/{gene}.{index}.downsample.depth"
                avg_depth, downsample_ratio = self.downsample_and_process(bam, output_bam, depth_file, output_depth, max_depth, seed, downsample_func, awk_script)
                if avg_depth < 2 and args['y']=="nanopore":
                    print(f"Reads coverage is too low for {bam}, skip it")
                    return
                elif avg_depth<1:
                    print(f"Reads coverage is too low for {bam}, skip it")
                    return
                print(f"downsample ratio is {downsample_ratio} for {gene}", flush=True)
                if args['y'] == "pacbio-hifi":
                    set_dp == 1
                elif args['y'] != "pacbio-hifi":
                    set_dp = int(max(0.3 * avg_depth, 0)) if avg_depth >= 5 else 0
                elif args['seq_tech'] == "rna" and args['RNA_type'] == "traditional":
                    set_dp = 1
                # set_dp = int(max(0.3 * avg_depth, 0)) if args['y'] != 'pacbio' or avg_depth >= 5 else 0
                    
                    
                if args['y'] == "nanopore" and set_dp==0:
                    set_dp=1
                min_cov = 2 if args['y'] == 'pacbio' else set_dp


                interval = f"{interval_dict[gene]}_ref{index+1}"
                if ":" in interval_dict[gene]:
                    interval_split = interval_dict[gene].split(":")
                    interval = f"{interval_split[0]}_ref{index+1}:{interval_split[1]}"

                self.process_gene(gene, hla_ref, interval, bam, depth_file, mask_bed, set_dp, min_cov, args, parameter, awk_script, index)



    def get_fasta(self):
        for gene in gene_list:
            # self.vcf2fasta(gene)
            # print (gene)
            self.vcf2fasta2ref(gene)
        # self.annotation()
            
    def splice_align_db(self):
        for gene in gene_list:
            hap1=f"{my_folder.sequence_dir}/{args['i']}.allele.1.{gene}.noN.fasta"
            hap2=f"{my_folder.sequence_dir}/{args['i']}.allele.2.{gene}.noN.fasta"
            if os.path.exists(hap1) and os.path.exists(hap2):
                # map to gene ref 
                cmd = f"""
                minimap2 -t {parameter.threads} -ax splice {my_db.full_db} {hap1} | samtools view -bS -F 0x804 -| samtools sort - >{my_folder.step2_genes_dir}/{gene}.allele1.bam
                samtools index {my_folder.step2_genes_dir}/{gene}.allele1.bam
                minimap2 -t {parameter.threads} -ax splice {my_db.full_db} {hap2} | samtools view -bS -F 0x804 -| samtools sort - >{my_folder.step2_genes_dir}/{gene}.allele2.bam
                samtools index {my_folder.step2_genes_dir}/{gene}.allele2.bam
                """
                os.system(cmd)
            else:
                print (f"WARNING: {hap1} or {hap2} does not exist")
                # sys.exit(1)

    def get_g_dict(self, g_anno_file):
        g_anno_dict = {}
        ## check if g_anno_file exists
        if not os.path.exists(g_anno_file):
            print(f"WARNING: {g_anno_file} does not exist")
            return g_anno_dict
        with open(g_anno_file, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if "/" in line:
                    items = line.strip().split(";")
                    alleles = items[1].split("/")
                    allele_g = items[-1]
                    gene_name=items[0]
                    for allele in alleles:
                        g_anno_dict[f"{gene_name}{allele}"] = allele_g
        return g_anno_dict


    def get_type_from_bam(self, bam):
        bam=pysam.AlignmentFile(bam, "rb")
        types = []
        for read in bam:
            reference_name = bam.get_reference_name(read.reference_id)
            types.append(reference_name)
        if len(types) == 0:
            types = ["NA"]
        return types
    
    def get_g_group(self, types, g_anno_dict):
        g_groups = []
        for t in types:
            t_fmt=t.replace("HLA-", "")
            if t_fmt in g_anno_dict and g_anno_dict[t_fmt] not in g_groups:
                g_groups.append(g_anno_dict[t_fmt])
        return g_groups

    def annoRNA(self):
        g_anno_file=my_db.g_group_annotation
        g_anno_dict=self.get_g_dict(g_anno_file)
        allele_dict={}
        for gene in gene_list:
            allele1_bam=f"{my_folder.step2_genes_dir}/{gene}.allele1.bam"
            allele2_bam=f"{my_folder.step2_genes_dir}/{gene}.allele2.bam"
            if os.path.exists(allele1_bam) and os.path.exists(allele2_bam):
                allele1_types=self.get_type_from_bam(allele1_bam)
                allele2_types=self.get_type_from_bam(allele2_bam)
                allele_dict[gene] = [allele1_types, allele2_types]
        # write to file
        result_file=f"{my_folder.sample_prefix}.{args['i']}.final.rna.type.result.txt"
        g_result_file=f"{my_folder.sample_prefix}.{args['i']}.final.rna.type.result.g.txt"
        with open(result_file, 'w') as f:
            f.write(f"# {my_db.version_info}\n")
            f.write("Locus\tChromosome\tAllele\n")
            for gene in gene_list:
                if gene in allele_dict:
                    allele1_types, allele2_types = allele_dict[gene]
                    f.write(f"{gene}\t{1}\t{';'.join(allele1_types)}\n")
                    f.write(f"{gene}\t{2}\t{';'.join(allele2_types)}\n")
        # write g group to file
        with open(g_result_file, 'w') as f:
            f.write(f"# {my_db.version_info}\n")
            f.write("Locus\tChromosome\tAllele\n")

            for gene in gene_list:
                if gene in allele_dict:
                    allele1_types, allele2_types = allele_dict[gene]
                    g_groups_1=self.get_g_group(allele1_types, g_anno_dict)
                    g_groups_2=self.get_g_group(allele2_types, g_anno_dict)
                    f.write(f"{gene}\t{1}\t{';'.join(g_groups_1)}\n")
                    f.write(f"{gene}\t{2}\t{';'.join(g_groups_2)}\n")   

        print ("please check the result in %s"%(result_file), flush=True)


    def annotation(self):
        if args['seq_tech'] == "rna":
            rna_tag='1'
        else:
            rna_tag='0'
        if args['i'] == "HLA" and args['seq_tech'] != 'rna':
            anno = f"""
            perl {sys.path[0]}/annoHLA.pl -s {parameter.sample} -i {my_folder.sequence_dir} -p {parameter.population} -r tgs -g {args["g"]} -t {rna_tag} -d {args["db"]}
            """
            # print (anno)
            os.system(anno)
        elif args['i'] == "HLA" and args['seq_tech'] == 'rna':
            # remove N in seq and align haplotype to the reference by splicing
            self.splice_align_db()
            self.annoRNA()

        elif args['i'] == "KIR":
            anno = f"""
            perl {sys.path[0]}/annoKIR.pl -s {parameter.sample} -i {my_folder.sequence_dir} -p {parameter.population} -r tgs -d {args["db"]} -t {rna_tag} -g {args["g"]}

            """
            # print (anno)
            os.system(anno)
        elif args['i'] == "CYP":
            anno = f"""
            perl {sys.path[0]}/annoCYP.pl -s {parameter.sample} -i {my_folder.sequence_dir} -p {parameter.population} -r tgs -d {args["db"]} -t {rna_tag} -g {args["g"]}
            """
            # print (anno)
            os.system(anno)

        if args['seq_tech'] != 'rna':
            refine = f"""
            python3 {sys.path[0]}/refine_typing.py -n {args["n"]} -o {args["o"]}  --db {args["db"]}  -i {args["i"]} --seq_tech {args["seq_tech"]} --RNA_type {args["RNA_type"]} --iteration {args["iteration"]}
            """
            os.system(refine)





if __name__ == "__main__":   

    parser = argparse.ArgumentParser(description="HLA Typing with only long-read data.", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    required.add_argument("-r", type=str, help="Long-read fastq file. PacBio or Nanopore.", metavar="\b")
    required.add_argument("-n", type=str, help="Sample ID", metavar="\b")
    required.add_argument("-o", type=str, help="The output folder to store the typing results.", metavar="\b", default="./output")
    required.add_argument("-i", type=str, help="HLA,KIR,CYP",metavar="\b", default="HLA")
    optional.add_argument("-p", type=str, help="The population of the sample [Asian, Black, Caucasian, Unknown, nonuse] for annotation. Unknown means use mean allele frequency in all populations. nonuse indicates only adopting mapping score and considering zero-frequency alleles.", metavar="\b", default="Unknown")
    optional.add_argument("-j", type=int, help="Number of threads.", metavar="\b", default=5)
    optional.add_argument("-d", type=float, help="Minimum score difference to assign a read to a gene.", metavar="\b", default=0.001)
    optional.add_argument("-g", type=int, help="Whether use G group resolution annotation [0|1].", metavar="\b", default=0)
    optional.add_argument("-m", type=int, help="1 represents typing, 0 means only read assignment", metavar="\b", default=1)
    optional.add_argument("-k", type=int, help="The mean depth in a window lower than this value will be masked by N, set 0 to avoid masking", metavar="\b", default=5)
    optional.add_argument("-a", type=str, help="Prefix of filtered fastq file.", metavar="\b", default="long_read")
    optional.add_argument("-y", type=str, help="Read type, [nanopore|pacbio|pacbio-hifi].", metavar="\b", default="pacbio")
    optional.add_argument("--minimap_index", type=int, help="Whether build Minimap2 index for the reference [0|1]. Using index can reduce memory usage.", metavar="\b", default=1)
    optional.add_argument("--db", type=str, help="db dir.", metavar="\b", default=sys.path[0] + "/../db/")
    optional.add_argument("--strand_bias_pvalue_cutoff", type=float, help="Remove a variant if the allele observations are biased toward one strand (forward or reverse). Recommand setting 0 to high-depth data.", metavar="\b", default=0.01)
    # optional.add_argument("-u", type=str, help="Choose full-length or exon typing. 0 indicates full-length, 1 means exon.", metavar="\b", default="0")
    optional.add_argument("--seed", type=int, help="seed to generate random numbers", metavar="\b", default=8)
    optional.add_argument("--max_depth", type=int, help="maximum depth for each HLA locus. Downsample if exceed this value.", metavar="\b", default=2000)
    optional.add_argument("-rt", "--RNA_type", type=str, help="traditional,2D,Direct,SIRV",metavar="\b", default="traditional")
    optional.add_argument("--seq_tech", type=str, help="Amplicon sequencing or WGS sequencing [wgs|amplicon|rna].", metavar="\b", default="wgs")
    optional.add_argument("--iteration", type=int, help="Iteration count.", metavar="\b", default=1)
    optional.add_argument("--dv_sif", type=str, help="DeepVariant sif file", metavar="\b")
    optional.add_argument("--snv_tool", type=str, help="longshot or deepvariant", metavar="\b", default="longshot")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    
    # Min_score = 0.1  #the read is too long, so the score can be very low.
    Min_score = 0  #the read is too long, so the score can be very low.
    Min_diff = args["d"]  #0.001

    # gene_list, interval_dict =  get_focus_gene(args)
    my_db = My_db(args)
    my_folder = My_folder(args)
    parameter = Parameters()

    # db_folder=os.path.dirname(my_db.full_cds_db) if args["seq_tech"] == "rna" else os.path.dirname(my_db.full_db)
    db_folder = os.path.dirname(my_db.full_db)
    gene_list = get_folder_list(db_folder)
    interval_dict={}
    for gene in gene_list:
        interval_dict[gene]=gene

    read_type = Read_Type(args["seq_tech"], args["y"], args["RNA_type"])
    minimap_para = read_type.get_minimap2_param()



    print ("start variant pipeline ...")
    # pbin = Pacbio_Binning()
        
    fa = Fasta()
    fa.get_fasta()
    print ("Sequence is reconstructed, start annotation...")
    fa.annotation()

    print ("Finished.")




