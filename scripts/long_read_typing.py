"""
Function 2 : typing with only long reads
"""

import sys
import os
import pysam
import gzip
import argparse
from collections import defaultdict


from downsample_bam import main
from read_objects import My_read, My_locus, Read_bin
from determine_gene import get_focus_gene
from db_objects import My_db
from alignment_modules import Read_Type


class Pacbio_Binning():

    def __init__(self):
        self.db = my_db.full_db

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
        self.outdir = "%s/%s/"%(outdir, self.sample)
        # self.whole_dir = "%s/whole/"%(sys.path[0])

        if not os.path.exists(args["o"]):
            os.system("mkdir %s"%(args["o"]))
        if not os.path.exists(self.outdir):
            os.system("mkdir %s"%(self.outdir))


class Fasta():

    def alignment(self, gene):
        ### call and phase snps
        cmd = """
        sample=%s
        bin=%s
        outdir=%s
        hla=%s
        hla_ref=%s
        minimap2 -t %s %s -a $hla_ref $outdir/$hla.%s.fq.gz | samtools view -bS -F 0x800 -| samtools sort - >$outdir/$hla.bam
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
        downsample_ratio = main(input_bam, output_bam, input_depth, output_depth, max_depth, seed)
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
        print("bam for vars:", bam, parameter.outdir)
        depth_file=f"{parameter.outdir}/{gene}.depth"
        print ("xxx", bam)
        print ("xxx", hla_ref)
        # call snp
        mask_bed=f"{parameter.outdir}/low_depth.bed"
        ovcf=f"{parameter.outdir}/{parameter.sample}.{gene}.dv.vcf"
        ogvcf=f"{parameter.outdir}/{parameter.sample}.{gene}.dv.g.vcf"
        call_phase_cmd= f"""
            set_dp={args["k"]}
            avg_depth=$(samtools depth -a {bam} | awk '{awk_script}')
            if (( $(echo "$avg_depth < 5" | bc -l) )) && [ {args['y']} = "pacbio" ]; then
                set_dp=0
            else
                set_dp=5
            fi
            python3 {sys.path[0]}/mask_low_depth_region.py -f False -c {depth_file} -o {parameter.outdir} -w 20 -d 5
            
            cp {mask_bed} {parameter.outdir}/{gene}.low_depth.bed
            longshot -F -c 2 -C 100000 -P {args["strand_bias_pvalue_cutoff"]} -r {interval_dict[gene]} --bam {bam} --ref {hla_ref} --out {parameter.outdir}/{parameter.sample}.{gene}.longshot.vcf 
            bgzip -f {parameter.outdir}/{parameter.sample}.{gene}.longshot.vcf
            tabix -f {parameter.outdir}/{parameter.sample}.{gene}.longshot.vcf.gz
            zcat {parameter.outdir}/{parameter.sample}.{gene}.longshot.vcf.gz >{parameter.outdir}/{parameter.sample}.{gene}.phased.vcf          
            bgzip -f {parameter.outdir}/{parameter.sample}.{gene}.phased.vcf
            tabix -f {parameter.outdir}/{parameter.sample}.{gene}.phased.vcf.gz

            # bash {sys.path[0]}/run_dv.sh {hla_ref} {bam} {ovcf} {ogvcf} {parameter.threads} {interval_dict[gene]}
        """
        os.system(call_phase_cmd)

        print ("call snp done")
        # call sv & phase snv-sv
        gene_work_dir=f"{parameter.outdir}/{gene}_work"
        if not os.path.exists(gene_work_dir):
            os.makedirs(gene_work_dir)
        sv_cmd = f"""
            bash {sys.path[0]}/refine_haplotype_pipe.sh {bam} {hla_ref} {gene} {interval_dict[gene]} {mask_bed} {gene_work_dir} {parameter.threads} {parameter.sample} 
        """
        # sv_cmd = f"""
        #     bash {sys.path[0]}/refine_haplotype_dv_pipe.sh {bam} {hla_ref} {gene} {interval_dict[gene]} {mask_bed} {gene_work_dir} {parameter.threads} {parameter.sample} 
        # """
        os.system(sv_cmd)
        ## reconstruct HLA sequence based on the phased snps & sv
        for index in range(2):
            order = f"""
            echo ">{gene}_{index}" >{parameter.outdir}/hla.allele.{index+1}.{gene}.fasta
            cat {gene_work_dir}/{gene}.{index+1}.raw.fa|grep -v ">" >>{parameter.outdir}/hla.allele.{index+1}.{gene}.fasta    
            samtools faidx {parameter.outdir}/hla.allele.{index+1}.{gene}.fasta    
            """
            os.system(order)


    def vcf2fasta2ref(self, gene):
        awk_script = '{{sum+=$3}} END {{ if (NR>0) print sum/NR; else print 0; }}'
        for index in range(2):
            bam=f"{parameter.outdir}/{gene}.{index}.bam"
            hla_ref=my_db.get_gene_alleles_2ref(gene, index)
            depth_file=f"{parameter.outdir}/{gene}.{index}.depth"
            mask_bed=f"{parameter.outdir}/low_depth.bed"
            max_depth = args["max_depth"]
            seed = args["seed"]
            output_bam = f"{parameter.outdir}/{gene}.{index}.downsample.bam"
            output_depth = f"{parameter.outdir}/{gene}.{index}.downsample.depth"
            downsample_ratio = main(bam, output_bam, depth_file, output_depth, max_depth, seed)
            if ":" in interval_dict[gene]:
                interval_split=interval_dict[gene].split(":")
                interval=f"{interval_split[0]}_ref{index+1}:{interval_split[1]}"
            else:
                interval=f"{interval_dict[gene]}_ref{index+1}"
            print (f"downsample ratio is {downsample_ratio} for {gene}", flush=True)
            if downsample_ratio < 1:
                os.system(f"rm {bam}")
                os.system(f"rm {depth_file}")
                os.system(f"mv {output_bam} {bam}")
                os.system(f"mv {output_depth} {depth_file}")

            call_phase_cmd= f"""
            set_dp={args["k"]}
            avg_depth=$(samtools depth -a {bam} | awk '{awk_script}')
            if (( $(echo "$avg_depth < 5" | bc -l) )) && [ {args['y']} = "pacbio" ]; then
                set_dp=0
            else
                set_dp=5
            fi
            python3 {sys.path[0]}/mask_low_depth_region.py -f False -c {depth_file} -o {parameter.outdir} -w 20 -d $set_dp
            
            cp {mask_bed} {parameter.outdir}/{gene}.{index}.low_depth.bed
            longshot -F -c 2 -C 100000 -P {args["strand_bias_pvalue_cutoff"]} -r {interval} --bam {bam} --ref {hla_ref} --out {parameter.outdir}/{parameter.sample}.{gene}.{index}.longshot.vcf 
            bgzip -f {parameter.outdir}/{parameter.sample}.{gene}.{index}.longshot.vcf
            tabix -f {parameter.outdir}/{parameter.sample}.{gene}.{index}.longshot.vcf.gz
            zcat {parameter.outdir}/{parameter.sample}.{gene}.{index}.longshot.vcf.gz >{parameter.outdir}/{parameter.sample}.{gene}.{index}.phased.vcf          
            bgzip -f {parameter.outdir}/{parameter.sample}.{gene}.{index}.phased.vcf
            tabix -f {parameter.outdir}/{parameter.sample}.{gene}.{index}.phased.vcf.gz
            """
            os.system(call_phase_cmd)
            print ("call snp done", flush=True)
            # call sv & phase snv-sv
            gene_work_dir=f"{parameter.outdir}/{gene}_{index}_work"
            if not os.path.exists(gene_work_dir):
                os.makedirs(gene_work_dir)
            sv_cmd = f"""
                bash {sys.path[0]}/refine_haplotype_2ref_pipe.sh {bam} {hla_ref} {gene} {interval} {mask_bed} {gene_work_dir} {parameter.threads} {parameter.sample} {index}
            """
            # sv_cmd = f"""
            #     bash {sys.path[0]}/refine_haplotype_dv_pipe.sh {bam} {hla_ref} {gene} {interval_dict[gene]} {mask_bed} {gene_work_dir} {parameter.threads} {parameter.sample} 
            # """
            os.system(sv_cmd)
                
            # generate sequence
            order = f"""
            echo ">{gene}_{index}" >{parameter.outdir}/hla.allele.{index+1}.{gene}.fasta
            cat {gene_work_dir}/{gene}.{index+1}.raw.fa|grep -v ">" >>{parameter.outdir}/hla.allele.{index+1}.{gene}.fasta    
            samtools faidx {parameter.outdir}/hla.allele.{index+1}.{gene}.fasta    
            """
            os.system(order)




    def get_fasta(self):
        for gene in gene_list:
            # self.vcf2fasta(gene)
            self.vcf2fasta2ref(gene)
        # self.annotation()

    def annotation(self):
        print(f"""perl {sys.path[0]}/annoHLA.pl -s {parameter.sample} -i {parameter.outdir} -p {parameter.population} -r tgs -g {args["g"]} -d {args["db"]}""")
        print(f"""python3 {sys.path[0]}/refine_typing.py -n {parameter.sample} -o {parameter.outdir}  --db {args["db"]}""")
        anno = f"""
        perl {sys.path[0]}/annoHLA.pl -s {parameter.sample} -i {parameter.outdir} -p {parameter.population} -r tgs -g {args["g"]} -d {args["db"]}
        cat {parameter.outdir}/hla.result.txt
        python3 {sys.path[0]}/refine_typing.py -n {parameter.sample} -o {parameter.outdir}  --db {args["db"]}  -i {args["i"]}
        """
        # print (anno)
        os.system(anno)


    #python3 {parameter.whole_dir}/g_group_annotation.py -s {parameter.sample} -i {parameter.outdir} -p {parameter.population} -j {args["j"]} 
    # def phase(self):
    #     hairs = f"$bin/ExtractHAIRs --triallelic 1 --pacbio 1 --indels 1 --ref $ref\
    #      --bam $outdir/$sample.tgs.sort.bam --VCF %s --out $outdir/fragment.tgs.file"
    # order = '%s/../bin/SpecHap -P --window_size 15000 --vcf %s --frag %s/fragment.sorted.file\
    #  --out %s/%s.specHap.phased.vcf'%(sys.path[0],my_new_vcf, outdir, outdir,gene)

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
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    parameter = Parameters()
    # Min_score = 0.1  #the read is too long, so the score can be very low.
    Min_score = 0  #the read is too long, so the score can be very low.
    Min_diff = args["d"]  #0.001

    gene_list, interval_dict =  get_focus_gene(args)
    my_db = My_db(args)

    read_type = Read_Type(args["y"])
    minimap_para = read_type.get_minimap2_param()



    print ("start variant pipeline ...")
    pbin = Pacbio_Binning()
        
    fa = Fasta()
    fa.get_fasta()
    print ("Sequence is reconstructed, start annotation...")
    fa.annotation()

    print ("Finished.")




