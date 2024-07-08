import random
import sys, os
import subprocess
from collections import defaultdict
from Bio import SeqIO  # Add this import statement

sys.path.insert(0, sys.path[0]+'/../scripts/')

from determine_gene import get_focus_gene
from Bio import SeqIO

"""
# given a file like: IGHA1   chr14_igh       1025372 1025764 16
IGHA2   chr14_igh       905141  905533  16
IGHD    chr14_igh       1159103 1159426 16
extract the gene name and the gene interval, and save them in a dict
"""
def get_gene_interval(gene_file):
	gene_list = []
	interval_dict = {}
	strand_dict = {}
	with open(gene_file, "r") as f:
		for line in f:
			line = line.strip().split()
			gene_list.append(line[0])
			# print (line)
			strand_dict[line[0]] = line[4]
			if line[1] not in interval_dict:
				interval_dict[line[1]] = {}
			interval_dict[line[1]][line[0]] = (int(line[2]), int(line[3]))
	# print ("intervali_dict",interval_dict)
	for genome in interval_dict:
		# sort the gene on it by their start position
		interval_dict[genome] = dict(sorted(interval_dict[genome].items(), key=lambda item: item[1][0]))

	return interval_dict, strand_dict

def DNA_complement2(sequence):
    sequence = sequence[::-1]
    trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    string = sequence.translate(trantab)
    return string


## use biopython load the allele sequence from db_fasta, save to gene:allele_list, and use a dict to save allele:sequence, use the Seq IO
def load_allele(db_fasta):
	gene_allele_dict = defaultdict(list)
	allele_seq_dict = {}
	for record in SeqIO.parse(db_fasta, "fasta"):
		gene = record.id.split("*")[0]
		allele = record.id
		sequence = str(record.seq).upper()
		gene_allele_dict[gene].append(allele)
		allele_seq_dict[allele] = sequence
	return gene_allele_dict, allele_seq_dict


# for each gene in gene_list, randomly select two alleles from gene_allele_dict, and save them in a new dict
def choose_allele():
	allele_dict = {}
	# for gene in gene_list:
	for gene in gene_allele_dict:
		alleles = gene_allele_dict[gene]
		# print (gene, alleles)
		if len(alleles) == 0:
			allele_dict[gene] = []
		elif len(alleles) == 1:
			allele_dict[gene] = alleles + alleles
		else:
			allele_dict[gene] = random.sample(alleles, 2)
	return allele_dict

# for each gene in allele_dict, write the gene and two alleles to out_fa, and write the gene and two alleles to out_allele
def output():
	with open(out_fa, "w") as f:
		for gene in allele_dict:
			if len(allele_dict[gene]) == 0:
				print (f"{gene} has no alleles")
				continue
			for allele in allele_dict[gene]:
				# print (allele)
				f.write(f">{allele}\n{allele_seq_dict[allele]}\n")
	with open(out_allele, "w") as f:
		f.write("Gene\tHap1\tHap2\n")
		for gene in allele_dict:
			if len(allele_dict[gene]) == 0:
				print (f"{gene} has no alleles")
				continue
			f.write(f"{gene}\t{allele_dict[gene][0]}\t{allele_dict[gene][1]}\n")

def insert_allele_to_ref(interval_dict, allele_dict, ref_hap,strand_dict):
	# for each genome in interval_dict, in ref_hap, replace the gene sequence using the allele sequence, use biopython to do this

	# Load the reference haplotype sequence for each chr, and store them in a dict
	ref_seq = {}
	with open(ref_hap, "r") as f:
		for record in SeqIO.parse(f, "fasta"):
			ref_seq[record.id] = str(record.seq).upper()
	# focus_chr_segments = {"chr14":[21621838, 106875071]}
	hap1_seq = {}
	hap2_seq = {}
		
	f = open(out_allele, "w")
	f.write("Gene\tHap1\tHap2\n")
	for genome in interval_dict:
		# if genome != "chr14_104363198_108375071":
		# 	continue
		gap_start = 0
		hap1_seq[genome] = ''
		hap2_seq[genome] = ''
		for gene in interval_dict[genome]:
			if gene not in allele_dict:
				print (f"{gene} has no alleles")
				# break
				continue
			# if gene != "TRDC":
			# 	continue
			# continue
			start, end = interval_dict[genome][gene]

			allele1, allele2 = allele_dict[gene]
			if len(allele_seq_dict[allele1]) == 0 or len(allele_seq_dict[allele2]) == 0:
				print (f"{allele1} or {allele2} has no base")
				# break
				continue
			f.write(f"{gene}\t{allele1}\t{allele2}\n")
			print (gene, genome, start, end, len(allele_seq_dict[allele1]), len(allele_seq_dict[allele2]), gap_start, start, allele1)
			if strand_dict[gene] == "+":
				hap1_seq[genome] += ref_seq[genome][gap_start:start] + allele_seq_dict[allele1]
				hap2_seq[genome] += ref_seq[genome][gap_start:start] + allele_seq_dict[allele2]
			else:
				hap1_seq[genome] += ref_seq[genome][gap_start:start] + DNA_complement2(allele_seq_dict[allele1])
				hap2_seq[genome] += ref_seq[genome][gap_start:start] + DNA_complement2(allele_seq_dict[allele2])

			gap_start = end

		hap1_seq[genome] += ref_seq[genome][gap_start:]
		hap2_seq[genome] += ref_seq[genome][gap_start:]
	f.close()
	# hap1_seq["chr14"] = hap1_seq["chr14"][22450000:22480000]
	# hap2_seq["chr14"] = hap2_seq["chr14"][22450000:22480000]
	# hap1_seq["chr14"] = hap1_seq["chr14"][21621838:22480000]
	# hap2_seq["chr14"] = hap2_seq["chr14"][21621838:22480000]
	# output hap1_seq and hap2_seq to hap1 and hap2
	with open(hap1, "w") as f:
		for genome in hap1_seq:
			f.write(f">{genome}\n{hap1_seq[genome]}\n")
	with open(hap2, "w") as f:
		for genome in hap2_seq:
			f.write(f">{genome}\n{hap2_seq[genome]}\n")
	
	os.system(f"cat {hap1} {hap2} > {hap}")

if __name__ == "__main__":   
	sample, outdir = sys.argv[1], sys.argv[2]
	## construct output dir
	if not os.path.exists(outdir ):
		print ("make dir")
		os.makedirs(outdir )

	print ("hi", outdir)
	gene_list, xx =  get_focus_gene("IG_TR")

	db_fasta = "../db/IG_TR/IG_TR.fasta"
	ref_hap = "../VDJ_ref/IG_TR.segment.fa"
	gene_file = "../gene_dist/IG_TR.gene.lite.bed"


	out_fa = f"{outdir}/{sample}.IG_TR.sep.fa"
	out_allele = f"{outdir}/{sample}.IG_TR.hap.alleles.txt"

	hap1 = f"{outdir}/{sample}.IG.TR.hap1.fa"
	hap2 = f"{outdir}/{sample}.IG.TR.hap2.fa"
	hap = f"{outdir}/{sample}.IG.TR.hap.fa"

	gene_allele_dict, allele_seq_dict = load_allele(db_fasta)
	allele_dict = choose_allele()
	# print (allele_dict)
	# output()
	# print (gene_list)
	
	interval_dict, strand_dict = get_gene_interval(gene_file)
	# print (interval_dict)
	insert_allele_to_ref(interval_dict, allele_dict, ref_hap, strand_dict)

	# f.write("Gene\tHap1\tHap2\n")






