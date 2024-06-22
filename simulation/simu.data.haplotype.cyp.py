import random
import sys, os
import subprocess
from collections import defaultdict
from Bio import SeqIO  # Add this import statement

sys.path.insert(0, sys.path[0]+'/../scripts/')

from determine_gene import get_focus_gene


## use biopython load the allele sequence from db_fasta, save to gene:allele_list, and use a dict to save allele:sequence, use the Seq IO
def load_allele(db_fasta):
	gene_allele_dict = defaultdict(list)
	allele_seq_dict = {}
	for record in SeqIO.parse(db_fasta, "fasta"):
		gene = record.id.split("*")[0]
		allele = record.id
		sequence = str(record.seq)
		gene_allele_dict[gene].append(allele)
		allele_seq_dict[allele] = sequence
	return gene_allele_dict, allele_seq_dict


# for each gene in gene_list, randomly select two alleles from gene_allele_dict, and save them in a new dict
def choose_allele():
	allele_dict = {}
	for gene in gene_list:
		alleles = gene_allele_dict[gene]
		print (gene, alleles)
		if len(alleles) == 0:
			allele_dict[gene] = []
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

sample, outdir = sys.argv[1], sys.argv[2]
## construct output dir
if not os.path.exists(outdir + "/" + sample):
	os.makedirs(outdir + "/" + sample)
gene_list, interval_dict =  get_focus_gene("CYP")

db_fasta = "../db/CYP/CYP.full.fasta"
out_fa = f"{outdir}/{sample}.CYP.sep.fa"
out_allele = f"{outdir}/{sample}.CYP.hap.alleles.txt"

gene_allele_dict, allele_seq_dict = load_allele(db_fasta)
allele_dict = choose_allele()
output()

# f.write("Gene\tHap1\tHap2\n")






