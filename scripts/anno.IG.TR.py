import sys
import subprocess
import os

from bed_objects import Bed_db
from determine_gene import get_focus_gene, get_folder_list
from db_objects import My_db

## get the chrs from the input hap1 and hap2
def get_chr(fasta):
	chrs = []
	with open(fasta) as in_file:
		for line in in_file:
			if line.startswith(">"):
				chrs.append(line.strip().split()[0][1:])
	return chrs

def anno_method1():

	sample, hap1, hap2, outdir, db, threads = sys.argv[1:]
	# chrs = ["chr1", "chr2", "chr7", "chr14", "chr22"]

	hap_list = [hap1, hap2]
	## check if the tmp folder exists, if not, create it
	if not os.path.exists(f"{outdir}/tmp"):
		os.makedirs(f"{outdir}/tmp")

	with open(f"{outdir}/{sample}.IG.TR.allele.txt", "w") as out_file:
		out_file.write("sample\tgene\tallele\tscore\tlen\tstart\tend\tchr\thap\n")
			
		for i in range(1, 3):
			chrs = get_chr(hap_list[i-1])
			for chr in chrs:
				
				os.system(f"samtools faidx {hap_list[i-1]} '{chr}' > {outdir}/tmp/{sample}.{chr}.hap{i}.fasta")
				subprocess.run(["blastn", "-query", f"{outdir}/tmp/{sample}.{chr}.hap{i}.fasta", "-out", f"{outdir}/tmp/{sample}.{chr}.hap{i}.blast.txt", "-db", db, "-outfmt", "7", "-max_target_seqs", "3000", "-num_threads", threads])
				# delete the fasta
				os.system(f"rm {outdir}/tmp/{sample}.{chr}.hap{i}.fasta")
				hashg = {}
				hashs = {}
				hashl = {}
				hashp = {}
				
				with open(f"{outdir}/tmp/{sample}.{chr}.hap{i}.blast.txt") as in_file:
					for line in in_file:
						if line.startswith("#"):
							continue
						
						# Split the line into a list of values
						values = line.split()

						# Assign values to variables, converting numerical values to integers
						allele = values[1]
						identity = float(values[2])
						length = int(values[3])
						mis = int(values[4])
						gap = int(values[5])
						start = int(values[6])
						end = int(values[7])
						score = float(values[11])

						gene = allele.split("*")[0]
						
						if mis >= 5:
							continue
						
						if identity < 97 and "TRB" in gene:
							continue
						
						if identity < 98 and "TRA" in gene:
							continue
						
						if identity < 97 and "IGK" in gene:
							continue
						
						if identity < 97 and "IGL" in gene:
							continue
						
						if identity < 98 and "IGH" in gene and gene != "IGHV1-2":
							continue
						
						if "IGKV1/ORY-1" in allele:
							continue
						
						if "IGKV1-NL1" in allele:
							continue
						
						if chr == "chr22_igk" and start < 1000:
							continue
						
						
						if "IGHV" in allele and length < 150:
							continue
						
						if "TRAV" in allele or "TRBV" in allele or "TRDV" in allele or "TRGV" in allele:
							if length < 100:
								continue
						
						if "IGHV" in allele and "OR" in chr:
							continue
						
						if chr == "chr15_igh" and not allele.startswith("IGHV") and not allele.startswith("IGHV/OR15"):
							continue
						
						if chr == "chr16_igh" and not allele.startswith("IGHV") and not allele.startswith("IGHV/OR16"):
							continue
						
						if gene not in hashg:
							hashg[gene] = allele
							hashs[gene] = identity
							hashl[gene] = length
							hashp[gene] = f"{start}\t{end}"
						else:
							if length >= hashl[gene] - 5 and identity > hashs[gene]:
								hashg[gene] = allele
								hashs[gene] = identity
								hashl[gene] = length
								hashp[gene] = f"{start}\t{end}"
				
				for g in sorted(hashg.keys()):
					out_file.write(f"{sample}\t{g}\t{hashg[g]}\t{hashs[g]}\t{hashl[g]}\t{hashp[g]}\t{chr}\thap{i}\n")

if __name__ == "__main__":
	anno_method1()