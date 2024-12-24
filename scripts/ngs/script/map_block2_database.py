"""
to obtain the linkage between blocks,
get block seq,
map seq to the database
get highest score of two blocks mapped to a same allele

wangshuai July 8, 2022
"""
import sys
import os




def blast_map(fragment1):
    command = f"""
    bin={sys.path[0]}/../bin
    $bin/samtools faidx {hla_ref} {fragment1}  | bcftools consensus -H 1 {vcf} >{outdir}/{fragment1}_hap1.fasta
    $bin/samtools faidx {hla_ref} {fragment1}  | bcftools consensus -H 2 {vcf} >{outdir}/{fragment1}_hap2.fasta

    $bin/blastn -query {outdir}/{fragment1}_hap1.fasta -out {outdir}/{fragment1}_hap1.fasta.out -subject {gene_db} -outfmt 6 -max_target_seqs 10000 -strand plus
    $bin/blastn -query {outdir}/{fragment1}_hap2.fasta -out {outdir}/{fragment1}_hap2.fasta.out -subject {gene_db} -outfmt 6 -max_target_seqs 10000 -strand plus

    """
    os.system(command)

class Construct_Graph():

    def __init__(self):
        self.fragments = []
        self.noise = 150
        self.break_point_list = [300]
        self.link_type = [[[1,1],[2,2]],  [[2,1],[1,2]] ]  
        self.dup_start = 3250 #3898
        self.dup_end = 3600 #4400

    def split_fragments(self):
        f = open(break_point_file)
        record_breakpoint_num = 0
        for line in f:
            if line[0] == "#":
                continue
            record_breakpoint_num += 1
            array = line.strip().split()
            # break_point = int(array[1])
            break_point = round((int(array[1]) + int(array[2]))/2)
            if break_point > self.dup_start - self.noise and break_point < self.dup_end + self.noise:
                continue
            self.break_point_list.append(break_point)
        self.break_point_list.append(gene_length)
        if gene == "HLA_DRB1":
            self.break_point_list += [self.dup_start, self.dup_end]
        self.break_point_list = sorted(self.break_point_list)
        self.save_fragments()

        if record_breakpoint_num == 0: # no breakpoint, thus no fragment splitted
            self.fragments = []
    
    def save_fragments(self):
        for i in range(len(self.break_point_list) - 1):
            start = self.break_point_list[i]+1
            end = self.break_point_list[i+1]
            fragment = "%s:%s-%s"%(gene, str(start), str(end))
            if fragment == "HLA_DRB1:%s-%s"%(self.dup_start+1, self.dup_end):
                continue
            self.fragments.append(fragment)
        # print (self.fragments)
    
    def get_edge(self):
        score_out = open(score_file, "w")
        print ("#frag1 frag2 00_edge_score;01_edge_score  00_allele;frag1_map_score;frag1_map_len;frag2_map_score;frag2_map_len \
        01_allele;frag1_map_score;frag1_map_len;frag2_map_score;frag2_map_len", file = score_out)
        if len(self.fragments) <= 1:
            print ("No need to phase block.")
            return 0
        for i in range(len(self.fragments)):
            blast_map(self.fragments[i])

        for i in range(len(self.fragments)):
            for j in range(i+1, len(self.fragments)):
                fragment1 = self.fragments[i]
                fragment2 = self.fragments[j]
                egde_info = []
                for x in range(2):
                    edge_score = 0
                    edge_link = None
                    for y in range(2):
                        # print (self.link_type[x][y])
                        blast_file_1 = f"{outdir}/{fragment1}_hap{self.link_type[x][y][0]}.fasta.out"
                        blast_file_2 = f"{outdir}/{fragment2}_hap{self.link_type[x][y][1]}.fasta.out"
                        analyze = Analyze_map()
                        link = analyze.main(blast_file_1, blast_file_2)
                        if link.high_score >= edge_score:
                            edge_link = link
                            edge_score = link.high_score
                    # print (fragment1,fragment2, x, edge_score, edge_link.support_allele)
                    egde_info += [edge_score, edge_link.support_allele]
                
                support_1 = ""
                for allele in egde_info[1]:  
                    for ele in allele:
                        support_1 = support_1 + str(ele) + ";" 
                support_2 = ""
                for allele in egde_info[3]:  
                    for ele in allele:
                        support_2 = support_2 + str(ele) + ";" 
                print (fragment1, fragment2, "%s;%s"%(egde_info[0], egde_info[2]), support_1, support_2, file = score_out)
        score_out.close()
   
class Analyze_map():
    # def __init__(self, gene):
    #     self.gene = gene
    def read_blast(self, blast_file): # in outfmt 6
        score_dict = {}
        f = open(blast_file)
        for line in f:
            array = line.strip().split()
            allele = array[1]
            score = round(float(array[2]),2)
            map_len = int(array[3])
            score_dict[allele] = [score, map_len]
        return score_dict
    
    def merge_score(self, score_dict_1, score_dict_2):
        merged_score_dict = {}
        for allele in score_dict_1:
            if allele in score_dict_2:
                merged_score_dict[allele] = score_dict_1[allele] + score_dict_2[allele]
        return merged_score_dict
    
    def main(self, blast_file_1, blast_file_2):
        # blast_file_1 = "/mnt/d/HLAPro_backup/haplotype/sample20/frag1_1.out"
        # blast_file_2 = "/mnt/d/HLAPro_backup/haplotype/sample20/frag1_1.out"
        score_dict_1 = self.read_blast(blast_file_1)
        score_dict_2 = self.read_blast(blast_file_2)
        merged_score_dict = self.merge_score(score_dict_1, score_dict_2)
        link = Linkage(merged_score_dict)
        # print (link.high_score, link.support_allele)
        return link

class Linkage():
    def __init__(self, merged_score_dict):
        high_score = 0
        support_allele = []
        if len(merged_score_dict) > 0:
            i = 0
            for allele in merged_score_dict:
                score_list = merged_score_dict[allele]
                weight_score = score_list[0] * score_list[1] + score_list[2] * score_list[3]
                if i == 0:
                    high_score = weight_score
                    support_allele.append([allele] + score_list )
                else:
                    if weight_score == high_score:
                        support_allele.append([allele] + score_list )

                i += 1
        self.high_score = high_score
        self.support_allele = support_allele


# analyze = Analyze_map()
# analyze.main()
# map = Map_database()
# map.blast_map("HLA_DRB1:1001-3950", "HLA_DRB1:4300-6378")
if __name__ == "__main__":

    # gene = "HLA_DRB1"
    # outdir = "/mnt/d/HLAPro_backup/haplotype/sample20/"
    gene = sys.argv[1]
    outdir = sys.argv[2] + "/"

    hla_ref = '%s/../db/immune.complex.gene.ref.extend.fasta'%(sys.path[0])
    gene_db = '%s/../db/whole/%s.fasta'%(sys.path[0], gene)
    vcf = "%s/%s.spechap.vcf.gz"%(outdir, gene)

    break_point_file = outdir + "/%s_break_points_spechap.txt"%(gene)
    score_file = outdir + "/%s_break_points_score.txt"%(gene)

    #gene_length_dict = {'HLA_A':[1000,4503],'HLA_B':[1000,5081],'HLA_C':[1000,5304],'HLA_DPA1':[1000,10775],'HLA_DPB1':[1000,12468],'HLA_DQA1':[1000,7492],'HLA_DQB1':[1000,8480],'HLA_DRB1':[1000,12229]}
    gene_length_dict = { 'HFE':[301, 8261], 'HLA-A':[301, 3802], 'HLA-B':[301, 4381], 'HLA-C':[301, 4618], 'HLA-DMA':[301, 5310], 'HLA-DMB':[301, 7040], 'HLA-DOA':[301, 3953], 'HLA-DOB':[301, 5086], 'HLA-DPA1':[301, 10075], 'HLA-DPA2':[301, 7043], 'HLA-DPB1':[301, 11826], 'HLA-DPB2':[301, 18134], 'HLA-DQA1':[301, 6784], 'HLA-DQA2':[301, 6152], 'HLA-DQB1':[301, 7402], 'HLA-DRA':[301, 6005], 'HLA-DRB1':[301, 11380], 'HLA-DRB3':[301, 13888], 'HLA-DRB4':[301, 15764], 'HLA-DRB5':[301, 13745], 'HLA-E':[301, 4122], 'HLA-F':[301, 3848], 'HLA-G':[301, 3438], 'HLA-H':[301, 3810], 'HLA-J':[301, 3844], 'HLA-K':[301, 3852], 'HLA-L':[301, 4070], 'HLA-N':[301, 935], 'HLA-P':[301, 3231], 'HLA-S':[301, 1174], 'HLA-T':[301, 2787], 'HLA-U':[301, 1030], 'HLA-V':[301, 2203], 'HLA-W':[301, 3272], 'MICA':[301, 13027], 'MICB':[301, 12616], 'TAP1':[301, 9570], 'TAP2':[301, 10907],\
            'KIR2DL1':[301, 15041], 'KIR2DL2':[301, 15082], 'KIR2DL3':[301,15068], 'KIR2DL4':[301, 11434], 'KIR2DL5':[301, 10072], 'KIR2DP1':[301, 13426], 'KIR2DS1':[301, 15020], 'KIR2DS2':[301, 14878], 'KIR2DS3':[301, 15403], 'KIR2DS4':[301, 16370], 'KIR2DS5':[301, 15548], 'KIR3DL1':[301, 14846], 'KIR3DL2':[301, 17301], 'KIR3DL3':[301, 12699], 'KIR3DP1':[301, 4540], 'KIR3DS1':[301, 15232],\
            'CYP19A1':[1,47775], 'CYP1A1':[1,7878], 'CYP1B1':[1,12177], 'CYP26A1':[1,11410], 'CYP2A13':[1,14732], 'CYP2A6':[1,13910], 'CYP2B6':[1,34098], 'CYP2C19':[1,99871], 'CYP2C8':[1,39726], 'CYP2C9':[1,58934], 'CYP2D6':[1,11312], 'CYP2F1':[1,20929], 'CYP2J2':[1,40444], 'CYP2R1':[1,21197], 'CYP2S1':[1,21330], 'CYP2W1':[1,13442], 'CYP3A4':[1,34205], 'CYP3A43':[1,45538], 'CYP4A22':[1,12568], 'CYP4B1':[1,20353], 'CYP4F2':[1,27051], 'CYP8A1':[1,64264], 'CYP3A5':[1,38805], 'CYP3A7':[1,37162] }
    gene_length = gene_length_dict[gene][1]


    cons = Construct_Graph()
    cons.split_fragments()
    cons.get_edge()
