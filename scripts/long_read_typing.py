"""
Function 1 : assign long reads to the gene
Function 2 : typing with only long reads
"""

import sys
import os
import pysam
import gzip
import argparse
from downsample_bam import main
from read_objects import My_read

# interval_dict = {"A":"HLA_A:1000-4503", "B":"HLA_B:1000-5081","C": "HLA_C:1000-5304","DPA1":"HLA_DPA1:1000-10775",\
#     "DPB1":"HLA_DPB1:1000-12468","DQA1":"HLA_DQA1:1000-7492","DQB1":"HLA_DQB1:1000-8480","DRB1":"HLA_DRB1:1000-12229" }
# gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
# gene_list = ['A']


# class Read_Obj():
#     # for each alignment record, extract information (identity, gene name)
#     def __init__(self, read):
#         mis_NM = 0
#         for ta in read.get_tags():
#             if ta[0] == 'NM':
#                 mis_NM = ta[1]  

#         self.match_num = 0  
#         self.read_length = 0      
#         for ci in read.cigar:
#             self.read_length += ci[1]
#             if ci[0] == 0:
#                 self.match_num += ci[1]
#             # if ci[0] == 0 or ci[0] == 1 or ci[0] == 2:
#             #     self.match_num += ci[1]
#             # if ci[0] == 1 or ci[0] == 2:
#             #     mis_NM += ci[1]
#         # print (read.query_name, read.reference_name, self.read_length, mis_NM)   
#         # self.read_length = len(read.query_sequence) 

#         self.read_name = read.query_name
#         self.allele_name = read.reference_name
#         self.mismatch_rate = round(float(mis_NM)/self.match_num, 6)
#         self.match_rate = 1 - self.mismatch_rate

#         # if self.match_num < 500:
#         #     self.match_rate = 0

#         self.loci_name = self.allele_name.split("*")[0]
#         if self.loci_name == "KIR2DL5A" or self.loci_name == "KIR2DL5B":
#             self.loci_name = "KIR2DL5"
#         if self.read_length < 400 and self.allele_name[0:3] == "KIR":
#             self.match_rate = 0

class Score_Obj():
    # determine which gene to assign
    def __init__(self):
        self.loci_score = {}
        self.loci_mismatch_score = {}
        self.read_loci = {}
    
    def add_read(self, read_obj):
        # score = round(read_obj.match_rate * (1 - read_obj.mismatch_rate), 6)
        score = read_obj.match_rate
        if read_obj.read_name not in self.loci_score:
            self.loci_score[read_obj.read_name] = {}
            self.loci_score[read_obj.read_name][read_obj.loci_name] = [score, read_obj.match_num]
        elif read_obj.loci_name not in self.loci_score[read_obj.read_name]:
            self.loci_score[read_obj.read_name][read_obj.loci_name] = [score, read_obj.match_num]
        else:
            if score > self.loci_score[read_obj.read_name][read_obj.loci_name][0]:
                self.loci_score[read_obj.read_name][read_obj.loci_name] = [score, read_obj.match_num]
    
    def assign(self, assign_file):
        f = open(assign_file, 'w')
        # print (len(self.loci_score))
        for read_name in self.loci_score: # for each read
            assigned_locus = []
            gene_score = sorted(self.loci_score[read_name].items(), key=lambda item: item[1][0], reverse = True)
            gene_match_len = sorted(self.loci_score[read_name].items(), key=lambda  x: x[1][1], reverse = True)
            # if len(gene_score) > 1 and (gene_score[0][0] == "DQB1"):
            #     print (read_name, gene_score[:2])
            if gene_score[0][1][0] <= Min_score:
                continue
            if len(gene_score) == 1: # mapped to only one gene, directly assign to that gene
                assigned_locus = [gene_score[0][0]]
            else:
                # real-data based adjustment
                

                # if gene_score[0][0] == "DRB1" or gene_score[1][0] == "DRB1":
                #     print (read_name, gene_score[0][0], gene_score[:5], gene_match_len[:5])
                if gene_score[0][0] in ["HLA-U"] and gene_score[1][0] == "HLA-A" :
                    assigned_locus = ["HLA-A"]                
                elif gene_score[0][0] == "HLA-DRB1" and gene_score[0][1][0] - gene_score[1][1][0] < 0.05:  # 0.02 0.05
                    continue
                elif gene_score[0][0] == "HLA-DQB1" and gene_score[0][1][0] < 0.9:
                    continue
                elif gene_score[0][0] == 'HLA-DPB2' and gene_score[1][0] == "HLA-DPA1":
                    assigned_locus = ["HLA-DPA1"]
                elif gene_score[0][0] in ['HLA-DPB1', "HLA-DPA1"] and gene_score[1][0] in ['HLA-DPB1', "HLA-DPA1"]:
                    assigned_locus = ['HLA-DPB1', "HLA-DPA1"]
                # map to more than one gene, check the score difference
                elif gene_score[0][1][0] - gene_score[1][1][0] >= Min_diff:
                    assigned_locus = [gene_score[0][0]]
                # score diff too small, can not determine which gene to assign
                # discard this read
                else:
                    continue
            # print ("assigned locus", assigned_locus)
            print (read_name, assigned_locus, file = f)
            self.read_loci[read_name] = assigned_locus
        f.close()
        return self.read_loci

class Pacbio_Binning():

    def __init__(self):
        
         
        # self.db = f"{sys.path[0]}/../db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta"
        # self.db = f"{sys.path[0]}/../db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta"
        # self.db = f"""{args["db"]}/ref/hla_gen.format.filter.extend.DRB.no26789.fasta"""

        if gene_class == "HLA":
            self.db = f"""{args["db"]}/HLA/ref/HLA.extend.fasta"""
        elif gene_class == "KIR":
            self.db = f"""{args["db"]}/KIR/ref/KIR.extend.select.fasta"""
            #self.db = f"{sys.path[0]}/../db/KIR/ref/KIR.extend.fasta"
        elif gene_class == "CYP":
            self.db = f"""{args["db"]}/CYP/ref/CYP.merge.fasta"""
        else:
            print ("wrong gene_class")

        self.sam = f"""{parameter.outdir}/{parameter.sample}.db.bam"""
        if args["m"] != 2:
            self.map2db()
        self.bamfile = pysam.AlignmentFile(self.sam, 'rb')   
        self.assign_file = f"{parameter.outdir}/{parameter.sample}.assign.txt"

    def index_db(self):
        ref_index = self.db[:-5] + args["y"] + ".mmi"
        # print ("search the reference index:", ref_index)
        if not os.path.isfile(ref_index):
            print ("start build Minimap2 index for the reference...")
            os.system(f"minimap2 {minimap_para} -d {ref_index} {self.db} ")
        else:
            print (f"Detect Minimap2 index for the reference: {ref_index}")
        self.db = ref_index

    def map2db(self):
        if args["minimap_index"] == 1:
            self.index_db()
        # map raw reads to database
        alignDB_order = f"""
        fq={parameter.raw_fq}
        ref={self.db}
        outdir={parameter.outdir}
        bin={sys.path[0]}/../bin
        sample={parameter.sample}
        minimap2 -t {parameter.threads} {minimap_para} -p 0.1 -N 100000 -a $ref $fq |samtools view -bS -o {self.sam}
        echo alignment done.
        """
        os.system(alignDB_order)

    def read_bam(self):
        # observe each read, assign it to gene based on alignment records
        scor = Score_Obj()
        for read in self.bamfile:
            if read.is_unmapped:
                continue
            # print (read)
            # read_obj = Read_Obj(read)
            read_obj = My_read(read)
            scor.add_read(read_obj)
            # print (read_obj.read_name, read_obj.mismatch_rate, read_obj.allele_name )
        read_loci = scor.assign(self.assign_file)
        for gene in gene_list:
            self.filter_fq(gene, read_loci)
        print ("reads-binning done.")

    def filter_fq(self, gene, dict):
        # output the assigned reads to the fastq file of each gene
        i = 0
        #gene = 'A'
        outfile = parameter.outdir + '/%s.%s.fq'%(gene, args["a"])
        out = open(outfile, 'w')
        flag = False
        if parameter.raw_fq.split(".")[-1] == "gz":
            f = gzip.open(parameter.raw_fq,'rt')
        else:
            f = open(parameter.raw_fq)
        for line in f:
            line = line.strip()
            if i % 4 == 0:
                read_name = line.split()[0][1:]
                if read_name in dict.keys() and gene in dict[read_name]:
                    flag = True
                    num = 1
                    print (line, file = out)
            elif flag:
                print (line, file = out)
                num += 1
                if num == 4:
                    flag = False
            i += 1
        f.close()
        out.close()
        os.system('gzip -f %s'%(outfile))


class Parameters():

    def __init__(self):

        self.sample = args["n"]
        self.raw_fq = args["r"]
        outdir = args["o"]
        self.population = args["p"]
        self.threads = args["j"]
        self.db = "%s/"%(args["db"])
        self.bin = "%s/../bin/"%(sys.path[0])      
        self.outdir = "%s/%s/"%(outdir, self.sample)
        self.whole_dir = "%s/whole/"%(sys.path[0])

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
        db=%s
        outdir=%s
        hla=%s
        hla_ref=$db/split_ref/$hla.fasta
        minimap2 -t %s %s -a $hla_ref $outdir/$hla.%s.fq.gz | samtools view -bS -F 0x800 -| samtools sort - >$outdir/$hla.bam
        samtools index $outdir/$hla.bam
        samtools depth -d 1000000 -aa $outdir/$hla.bam >$outdir/$hla.depth
        """%(parameter.sample, parameter.bin, args["db"], parameter.outdir, gene, parameter.threads, minimap_para, args["a"])
        os.system(cmd)

        max_depth = args["max_depth"]
        seed = args["seed"]
        input_bam = f"{parameter.outdir}/{gene}.bam"
        input_depth = f"{parameter.outdir}/{gene}.depth"
        output_bam = f"{parameter.outdir}/{gene}.downsample.bam"
        output_depth = f"{parameter.outdir}/{gene}.downsample.depth"
        downsample_ratio = main(input_bam, output_bam, input_depth, output_depth, max_depth, seed)
        print (f"downsample ratio is {downsample_ratio} for {gene}")
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
        hla_ref=f"{parameter.db}/split_ref/{gene}.fasta"
        bam=f"{parameter.outdir}/{gene}.bam"
        depth_file=f"{parameter.outdir}/{gene}.depth"
        print ("xxx", bam)
        # call snp
        mask_bed=f"{parameter.outdir}/low_depth.bed"
        call_phase_cmd= f"""
            set_dp={args["k"]}
            avg_depth=$(samtools depth -a {bam} | awk '{awk_script}')
            if (( $(echo "$avg_depth < 5" | bc -l) )) && [ {args['y']} = "pacbio" ]; then
                set_dp=0
            else
                set_dp=5
            fi
            python3 {sys.path[0]}/mask_low_depth_region.py -f False -c {depth_file} -o {parameter.outdir} -w 20 -d $set_dp
            
            cp {mask_bed} {parameter.outdir}/{gene}.low_depth.bed
            longshot -F -c 2 -C 100000 -P {args["strand_bias_pvalue_cutoff"]} -r {interval_dict[gene]} --bam ${bam} --ref {hla_ref} --out {parameter.outdir}/{parameter.sample}.{gene}.longshot.vcf 
            bgzip -f {parameter.outdir}/{parameter.sample}.{gene}.longshot.vcf
            tabix -f {parameter.outdir}/{parameter.sample}.{gene}.longshot.vcf.gz
            zcat {parameter.outdir}/{parameter.sample}.{gene}.longshot.vcf.gz >{parameter.outdir}/$sample.{gene}.phased.vcf          
            bgzip -f {parameter.outdir}/{parameter.sample}.{gene}.phased.vcf
            tabix -f {parameter.outdir}/{parameter.sample}.{gene}.phased.vcf.gz
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
        os.system(sv_cmd)

        ## reconstruct HLA sequence based on the phased snps & sv
        for index in range(2):
            order = f"""
            echo ">HLA_{gene}_{index}" >$outdir/hla.allele.$i.HLA_$hla.fasta
            cat {gene_work_dir}/{gene}.{index+1}.raw.fa|grep -v ">" >>$outdir/hla.allele.$i.HLA_$hla.fasta    
            samtools faidx $outdir/hla.allele.{index+1}.HLA_{gene}.fasta    
            """
            os.system(order)



    def get_fasta(self):
        for gene in gene_list:
            self.vcf2fasta(gene)
        # self.annotation()

    def annotation(self):
        anno = f"""
        perl {parameter.whole_dir}/annoHLA.pl -s {parameter.sample} -i {parameter.outdir} -p {parameter.population} -r tgs -g {args["g"]} -d {args["db"]}/HLA 
        cat {parameter.outdir}/hla.result.txt
        python3 {parameter.whole_dir}/../refine_typing.py -n {parameter.sample} -o {parameter.outdir}  --db {args["db"]}
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
    optional.add_argument("-y", type=str, help="Read type, [nanopore|pacbio].", metavar="\b", default="pacbio")
    optional.add_argument("--minimap_index", type=int, help="Whether build Minimap2 index for the reference [0|1]. Using index can reduce memory usage.", metavar="\b", default=0)
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

    gene_class = args["i"]
    if gene_class == "HLA":
        gene_list = [ 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DPB2', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-E', 'HLA-F', 'HLA-G', 'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L', 'HLA-P', 'HLA-V', 'HLA-DQA2', 'HLA-DPA2', 'HLA-N', 'HLA-S', 'HLA-T', 'HLA-U', 'HLA-W', 'MICA', 'MICB', 'TAP1', 'TAP2', 'HFE' ]
        interval_dict = {"HFE": "HFE:301-8261","HLA-A": "HLA-A:301-3802","HLA-B": "HLA-B:301-4381","HLA-C": "HLA-C:301-4618","HLA-DMA": "HLA-DMA:301-5310","HLA-DMB": "HLA-DMB:301-7040","HLA-DOA": "HLA-DOA:301-3953","HLA-DOB": "HLA-DOB:301-5086","HLA-DPA1": "HLA-DPA1:301-10075","HLA-DPA2": "HLA-DPA2:301-7043","HLA-DPB1": "HLA-DPB1:301-11826","HLA-DPB2": "HLA-DPB2:301-18134","HLA-DQA1": "HLA-DQA1:301-6784","HLA-DQA2": "HLA-DQA2:301-6152","HLA-DQB1": "HLA-DQB1:301-7402","HLA-DRA": "HLA-DRA:301-6005","HLA-DRB1": "HLA-DRB1:301-11380","HLA-DRB3": "HLA-DRB3:301-13888","HLA-DRB4": "HLA-DRB4:301-15764","HLA-DRB5": "HLA-DRB5:301-13745","HLA-E": "HLA-E:301-4122","HLA-F": "HLA-F:301-3848","HLA-G": "HLA-G:301-3438","HLA-H": "HLA-H:301-3810","HLA-J": "HLA-J:301-3844","HLA-K": "HLA-K:301-3852","HLA-L": "HLA-L:301-4070","HLA-N": "HLA-N:301-935","HLA-P": "HLA-P:301-3231","HLA-S": "HLA-S:301-1174","HLA-T": "HLA-T:301-2787","HLA-U": "HLA-U:301-1030","HLA-V": "HLA-V:301-2203","HLA-W": "HLA-W:301-3272","MICA": "MICA:301-13027","MICB": "MICB:301-12616","TAP1": "TAP1:301-9570","TAP2": "TAP2:301-10907"}
    if gene_class == "CYP":
        gene_list = [ 'CYP19A1', 'CYP1A1', 'CYP1B1', 'CYP26A1', 'CYP2A13', 'CYP2A6', 'CYP2B6', 'CYP2C19', 'CYP2C8', 'CYP2C9', 'CYP2D6', 'CYP2F1', 'CYP2J2', 'CYP2R1', 'CYP2S1', 'CYP2W1', 'CYP3A4', 'CYP3A43', 'CYP4A22', 'CYP4B1', 'CYP4F2', 'CYP8A1', 'CYP3A5', 'CYP3A7' ]
        interval_dict = {"CYP4B1": "CYP4B1:1-2170","CYP3A4": "CYP3A4:1-34205","CYP3A7": "CYP3A7:1-37162","CYP2F1": "CYP2F1:1-20929","CYP2A13": "CYP2A13:1-14732","CYP2C9": "CYP2C9:1-58934","CYP26A1": "CYP26A1:1-11410","CYP3A5": "CYP3A5:1-38805","CYP1A1": "CYP1A1:1-7878","CYP3A43": "CYP3A43:1-45538","CYP2A6": "CYP2A6:1-13910","CYP19A1": "CYP19A1:1-47775","CYP2S1": "CYP2S1:1-21330","CYP1B1": "CYP1B1:1-12177","CYP2B6": "CYP2B6:1-34098","CYP8A1": "CYP8A1:1-5603","CYP2C19": "CYP2C19:1-99871","CYP4F2": "CYP4F2:1-27051","CYP2D6": "CYP2D6:1-11312","CYP4A22": "CYP4A22:1-12568","CYP2R1": "CYP2R1:1-21197","CYP2C8": "CYP2C8:1-39726","CYP2W1": "CYP2W1:1-13442","CYP2J2": "CYP2J2:1-40444"}
    if gene_class == "KIR":
        gene_list = ['KIR2DL1', 'KIR2DL2', 'KIR2DL3', 'KIR2DL4', 'KIR2DL5', 'KIR2DP1', 'KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS4', 'KIR2DS5', 'KIR3DL1', 'KIR3DL2', 'KIR3DL3', 'KIR3DP1', 'KIR3DS1']
        interval_dict = {"KIR2DL1": "KIR2DL1:301-15041","KIR2DL2": "KIR2DL2:301-15082","KIR2DL3": "KIR2DL3:301-15068","KIR2DL4": "KIR2DL4:301-11434","KIR2DL5": "KIR2DL5:301-10072","KIR2DP1": "KIR2DP1:301-13426","KIR2DS1": "KIR2DS1:301-15020","KIR2DS2": "KIR2DS2:301-14878","KIR2DS3": "KIR2DS3:301-15403","KIR2DS4": "KIR2DS4:301-16392","KIR2DS5": "KIR2DS5:301-15548","KIR3DL1": "KIR3DL1:301-14846","KIR3DL2": "KIR3DL2:301-17301","KIR3DL3": "KIR3DL3:301-12699","KIR3DP1": "KIR3DP1:301-4540","KIR3DS1": "KIR3DS1:301-15232"}


    minimap_para = ''
    if args["y"] == "pacbio":
        minimap_para = " -x map-pb "
    elif args["y"] == "nanopore":
        minimap_para = " -x map-ont "


    ###assign reads
    if args["m"] == 10086:
        print ("skip assignment, just for testing")
    elif args["m"] == 2:
        print ("start variant pipeline ...")
        pbin = Pacbio_Binning()
    else:
        pbin = Pacbio_Binning()
        pbin.read_bam()        

    if args["m"] != 0:
        fa = Fasta()
        fa.get_fasta()
        print ("Sequence is reconstructed, start annotation...")
        fa.annotation()
    print ("Finished.")




