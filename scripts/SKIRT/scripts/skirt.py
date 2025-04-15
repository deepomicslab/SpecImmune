import argparse
import pandas as pd
import re
import gzip
import os

from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio import SeqIO
from Bio.Seq import Seq

class PafParser:
	columns=['query_name', 'query_length', 'query_start', 'query_end',
		'strand', 'target_name', 'target_length', 'target_start',
		'target_end', 'matching_length', 'mapping_length',
		'mapping_quality']

	def __init__(self, paf_file_path):
		self.paf_file_path = paf_file_path
		self.df = self.parse()

	def get_initial_columns(self):
		return self.columns

	def parse(self):
		optional_fields = set()
		data = []
		with open(self.paf_file_path) as f:
			for line in f:
				fields = line.strip().split('\t')
				row = {col: self.field_value(col, fields[i]) for i, col in enumerate(self.columns)}
				# Parse optional fields
				if len(fields) > len(self.columns):
					for opt_field in fields[len(self.columns):len(fields)-1]:
						key, dtype, val = opt_field.split(':', 2)
						# Drop 's2' column if it exists
						if (key!='s2'):
							if (dtype=='i'):
								row[key] = int(val)
							elif (dtype=='f'):
								row[key] = float(val)
							else:
								row[key] = val
							optional_fields.add(key)
					row['kir_allele'] = fields[len(fields) - 1].split(' ')[0]
					kir_nomenclature = row['kir_allele'].split('*')
					row['kir_gene'] = kir_nomenclature[0][3:]
					row['kir_cds_3'] = kir_nomenclature[1][:3]
					row['kir_cds_5'] = kir_nomenclature[1][:5]
					row['kir_cds_7'] = kir_nomenclature[1][:7]
					row['num_exons'] = row['cs'].count('~') + 1
				row['matching_rate'] = 1 - (float(row['query_length'] - row['query_end'] + row['query_start'] + row['NM']) / float(max(row['mapping_length'], row['query_length'])))
				row['resolution'] = 5 if (row['matching_rate']==1) else 0
				# Handle special case for known fusion alleles
				row['resolution'] = len(row['kir_allele']) - len(kir_nomenclature[0]) - 1 if (row['query_name'].count("KIR") > 2) else row['resolution']
				data.append(row)
		df = pd.DataFrame(data)
		# Set optional fields with missing values to empty string
		for field in optional_fields:
			if field not in df.columns:
				df[field] = ''

		# Ignore lengths bigger than a normal KIR gene
		df = df[(abs(df['target_end'] - df['target_start']) <= 18000)]# & (df['matching_rate'] >= 0.9)]

		return df

	def filter_perfect_match(self):
		df = self.df
		self.df = df[(df['NM'] == 0) & (df['query_start'] == 0) & (df['query_length'] == df['query_end'])]

	def field_value(self, col_name, value):
		if (col_name.endswith('_length') or col_name.endswith('_start') or col_name.endswith('_end')
			or col_name.endswith('_bases') or col_name.endswith('_quality')):
			return int(value)
		else:
			return value

class FormatComposer:
	gene_colors = {
		"3DL3" : "138,43,226",
		"2DS2" : "0,100,255",
		"2DL2" : "0,100,255",
		"2DL3" : "0,100,255",
		"2DL5B": "0,100,255",
		"2DL5" : "0,100,255",
		"2DS3" : "0,100,255",
		"2DS5" : "0,100,255",
		"2DP1" : "77,88,88",
		"2DL1" : "0,100,255",
		"3DP1" : "138,43,226",
		"2DL4" : "138,43,226",
		"3DL1" : "0,100,255",
		"3DS1" : "0,100,255",
		"2DL5A": "0,100,255",
		"2DS1" : "0,100,255",
		"2DS4" : "0,100,255",
		"3DL2" : "138,43,226"
	}

	gene_counts = {
		"3DL3" : 0,
		"2DS2" : 0,
		"2DL2" : 0,
		"2DL3" : 0,
		"2DL5B": 0,
		"2DL5" : 0,
		"2DS3" : 0,
		"2DS5" : 0,
		"2DP1" : 0,
		"2DL1" : 0,
		"3DP1" : 0,
		"2DL4" : 0,
		"3DL1" : 0,
		"3DS1" : 0,
		"2DL5A": 0,
		"2DS1" : 0,
		"2DS4" : 0,
		"3DL2" : 0
	}

	bed_columns = ["chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"]
	vcf_columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

	bed_dfs = {}
	vcf_dfs = {}

	def __init__(self, df, out_path):
		self.df = df
		self.out_path = out_path
		self.process_alleles()

	def reset_color(self):
		self.gene_counts = {key : 0 for key in self.gene_counts} # reset all counts

	def set_color(self, gene_name):
		dup_color = "255,0,0"
		self.gene_counts[gene_name] += 1
		if (self.gene_counts[gene_name]>1):
			return dup_color
		else:
			return self.gene_colors[gene_name]

	def process_alleles(self):
		# Group the input dataframe by 'target_name' and sort by 'target_start' and 'strand'
		df_sorted = self.df.groupby("target_name").apply(sort_based_on_strand).reset_index(drop=True)
		kir_contigs = df_sorted.groupby("target_name")

		# Process each group in the sorted input dataframe
		for target, kir_contig in kir_contigs:

			# Create an empty output dataframe with the desired columns
			beddf = pd.DataFrame(columns=self.bed_columns)
			vcfdf = pd.DataFrame(columns=self.vcf_columns)
			bedrows = []
			vcfrows = []
			self.reset_color()
			for pos, (idx, kir) in enumerate(kir_contig.iterrows()):
				kir_gene = kir["kir_gene"]

				bedrow = {self.bed_columns[0]: target}
				bedrow[self.bed_columns[1]] = kir["target_start"]
				bedrow[self.bed_columns[2]] = kir["target_end"]
				bedrow[self.bed_columns[3]] = format_allele_info(kir)
				bedrow[self.bed_columns[4]] = 0 #"." isn't compatible with Genome Browser
				bedrow[self.bed_columns[5]] = kir["strand"]
				bedrow[self.bed_columns[6]] = kir["target_start"]
				bedrow[self.bed_columns[7]] = kir["target_end"]
				bedrow[self.bed_columns[8]] = self.set_color(kir_gene)
				exon_infos = cstr_cutter(kir["cs"])
				bedrow[self.bed_columns[9]] = len(exon_infos[0])
				bedrow[self.bed_columns[10]] = ",".join([str(size) for size in exon_infos[0]])
				bedrow[self.bed_columns[11]] = ",".join([str(start) for start in exon_infos[1]])
				bedrows.append(bedrow)
				#bed_columns = ["chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd",
				#	"itemRgb", "blockCount", "blockSizes", "blockStarts"]

				for start_offset, variant in zip(exon_infos[2], exon_infos[3]):

					variant_start = kir["target_start"] + start_offset

					bedrow = {self.bed_columns[0]: target}
					bedrow[self.bed_columns[1]] = variant_start
					bedrow[self.bed_columns[2]] = variant_start + 1
					bedrow[self.bed_columns[3]] = variant
					bedrow[self.bed_columns[4]] = 0	#"." isn't compatible with Genome Browser
					bedrow[self.bed_columns[5]] = kir["strand"]
					bedrow[self.bed_columns[6]] = variant_start
					bedrow[self.bed_columns[7]] = variant_start + 1
					bedrow[self.bed_columns[8]] = "255,255,0"
					bedrow[self.bed_columns[9]] = 1
					bedrow[self.bed_columns[10]] = 1
					bedrow[self.bed_columns[11]] = 0
					bedrows.append(bedrow)

					vcfrow = {key : "." for key in self.vcf_columns} #initialize vcfrow with default value
					#vcf_columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
					vcfrow[self.vcf_columns[0]] = target
					vcfrow[self.vcf_columns[1]] = variant_start + 1
					vcfrow[self.vcf_columns[2]] = kir["kir_allele"]
					if (len(variant)>1):
						vcfrow[self.vcf_columns[3]] = variant[0] if variant[1]=='>' else (variant[1::] if variant[0]=='-' else ".")
						vcfrow[self.vcf_columns[4]] = variant[2] if variant[1]=='>' else (variant[1::] if variant[0]=='+' else ".")
					vcfrows.append(vcfrow)

			# Add the processed row to the output dataframe
			beddf = pd.concat([beddf, pd.DataFrame(bedrows)], ignore_index=True)
			vcfdf = pd.concat([vcfdf, pd.DataFrame(vcfrows)], ignore_index=True)

			self.bed_dfs[target] = beddf
			self.vcf_dfs[target] = vcfdf

	def output_bed_file(self):
		for target in self.bed_dfs:
			with open(self.out_path + "." + target.replace("#", "-") + ".bed", "w+") as f:
				f.write("track name=KIR description=\"killer-cell immunoglobulin-like receptors\" visibility=2 itemRgb=On\n")
				# Write the BED file
				self.bed_dfs[target].to_csv(f, sep="\t", index=False, header=False)

	def output_vcf_file(self):
		for target in self.vcf_dfs:
			with open(self.out_path + "." + target.replace("#", "-") + ".vcf", "w+") as f:
				f.write("##fileformat=VCFv4.2\n")
				f.write("#" + "\t".join(self.vcf_columns) + "\n")
				# Write the VCF file
				self.vcf_dfs[target].to_csv(f, sep="\t", index=False, header=False)

def cstr_cutter(cstr = ''):
	index = 0
	snp = 0
	len_list = []
	start_list = []
	variant_content = []
	variant_start = []
	for m in re.findall(r'([-+])([acgtn]+)|(:)(\d+)|(\*)([acgtn][acgtn])|(~[acgtn]{2})(\d+)[acgtn]{2}', cstr):
		for i in range(len(m)):
			if (not m[i]):
				continue
			if(m[i]==':'):
				if(snp==1):
					len_list.append(len_list.pop()+int(m[i+1]))
				else:
					start_list.append(index)
					len_list.append(int(m[i+1]))
				index += int(m[i+1])
				i+=1
				snp = 0
			elif(m[i]=='*'):
				variant_start.append(index)
				variant_content.append(m[i+1][0].upper() + ">" + m[i+1][1].upper())
				len_list.append(len_list.pop()+1)
				index += 1
				i += 1
				snp = 1
			elif(m[i][0]=='~'):
				index += int(m[i+1])
				i+=1
				snp = 0
			elif(m[i]=='+'):
				variant_start.append(index)
				variant_content.append(m[i] + m[i+1].upper())
				i+=1
			elif(m[i]=='-'):
				variant_start.append(index)
				variant_content.append(m[i] + m[i+1].upper())
				index += len(m[i+1])
				i+=1
	return (len_list, start_list, variant_start, variant_content)

def filter_overlapping_rows(df):
	groups = df.groupby('target_name')
	rows_to_drop = []
	for target, group in groups:
		if len(group) <= 1:
			continue
		group = group.sort_values(['target_start', 'NM', 'matching_length', 'target_end'], ascending=[True, True, False, False])
		for i in range(len(group)):
			row1 = group.iloc[i]
			if (row1.name in rows_to_drop):
				continue
			for j in range(i + 1, len(group)):
				row2 = group.iloc[j]
				if (row2.name in rows_to_drop or row1.target_end < row2.target_start or row2.target_end < row1.target_start):
					continue
				if (abs(row1.target_start - row2.target_start) <= 3000 or abs(row1.target_end - row2.target_end) <= 3000):
					if ((row1.matching_rate == 1 or row1.NM == 0) and row1.matching_length > row2.matching_length):
						rows_to_drop.append(row2.name)
					elif ((row2.matching_rate == 1 or row2.NM == 0) and row1.matching_length < row2.matching_length):
						rows_to_drop.append(row1.name)
					# some alleles lost several bases at the head of exon1 but match perfectly
					elif (row1.kir_gene == row2.kir_gene and row1.num_exons != row2.num_exons and
							(row1.NM == 0 or row2.NM == 0) and row1.matching_rate >= 0.9 and row2.matching_rate >= 0.9):
						df.at[group.index[i], 'target_start'] = df.at[group.index[j], 'target_start'] = min(row1.target_start, row2.target_start)
						df.at[group.index[i], 'target_end'] = df.at[group.index[j], 'target_end'] = max(row1.target_end, row2.target_end)
					elif (row1.matching_rate < row2.matching_rate):
						rows_to_drop.append(row1.name)
					elif (row1.matching_rate > row2.matching_rate):
						rows_to_drop.append(row2.name)
					elif (row1.NM > row2.NM):
						rows_to_drop.append(row1.name)
					elif (row1.NM < row2.NM):
						rows_to_drop.append(row2.name)
					elif (row1.matching_length < row2.matching_length):
						rows_to_drop.append(row1.name)
					elif (row1.matching_length > row2.matching_length):
						rows_to_drop.append(row2.name)
					elif (row1.mapping_length < row2.mapping_length):
						rows_to_drop.append(row1.name)
					elif (row1.mapping_length > row2.mapping_length):
						rows_to_drop.append(row2.name)
	df.drop(rows_to_drop, inplace=True)
	return df

def check_genomic_match(cds_df, genomic_df):
	for idx, row in cds_df[cds_df['matching_rate']==1].iterrows():
		genomic_rows = genomic_df[(abs(genomic_df['target_start']-row['target_start'])<=3000) & (genomic_df['query_name'] == row['query_name'])].sort_values('matching_rate', ascending=False)
		if not genomic_rows.empty:
			if (genomic_rows.iloc[0]['matching_rate']==1 and
				genomic_rows.iloc[0]['target_end']-genomic_rows.iloc[0]['target_start'] == genomic_rows.iloc[0]['query_end']-genomic_rows.iloc[0]['query_start']):
				cds_df.loc[idx, 'resolution'] = 7
			else:
				cds_df.loc[idx, 'resolution'] = 6
	return cds_df

def check_synonymous(subject_fasta_path, query_fasta_path):

	# Set the command line arguments for tblastn with the specified filters
	tblastn_cline = NcbitblastnCommandline(query=query_fasta_path, subject=subject_fasta_path, outfmt=5,
						evalue=1e-10, out=subject_fasta_path + "_results.xml")

	# Run the tblastn search
	stdout, stderr = tblastn_cline()
	# Parse the records from the XML output file as a list
	blast_results = list(NCBIXML.parse(open(subject_fasta_path + "_results.xml")))
	all_hsps = []
	# Loop through the blast results and extract the hits that meet the specified filters
	for result in blast_results:
		result.max_hsp = None
		for alignment in result.alignments:
			#for hsp in alignment.hsps:
			sorted_hsps = sorted((hsp for hsp in alignment.hsps if hsp.gaps == 0 and hsp.expect == 0),
					key=lambda hsp: (-hsp.identities, -hsp.align_length))
			if (len(sorted_hsps) > 0):
				max_hsp = sorted_hsps[0]
				if (result.max_hsp is None or result.max_hsp.identities < max_hsp.identities or
					(result.max_hsp.identities==max_hsp.identities and result.max_hsp.align_length < max_hsp.align_length)):
					result.max_hsp = max_hsp
	sorted_results = sorted((result for result in blast_results if result.max_hsp is not None),
				key=lambda result: (-result.max_hsp.identities, -(result.max_hsp.identities/result.query_length), -result.max_hsp.align_length))
	max_identities = sorted_results[0].max_hsp.identities if len(sorted_results) > 0 else 0
	syn_or_closest_matchings = sorted((result for result in sorted_results if result.max_hsp.identities == max_identities),
					key=lambda result: (-(result.max_hsp.identities/result.query_length), result.query.split(' ')[1]))
	if (len(syn_or_closest_matchings)):
		return syn_or_closest_matchings[0].query.split(' ')[1], syn_or_closest_matchings[0].max_hsp.identities/syn_or_closest_matchings[0].query_length
	else:
		return None, 0

def locus_grouping(x):
	# Group the rows based on the target_start in 3000 bp range as one gene locus
	groups = []
	for i, row in x.iterrows():
		found_group = False
		for group in groups:
			if (len(group) > 0 and (abs(group[0]['target_start'] - row['target_start']) <= 3000
					or abs(group[0]['target_end'] - row['target_end']) <= 3000)):
				group.append(row)
				found_group = True
				break
		if (not found_group):
			groups.append([row])

	# Sort each sub-group by resolution in descending order
	for i, group in enumerate(groups):
		groups[i] = pd.DataFrame(group).sort_values(by=['resolution', 'matching_rate', 'kir_allele'], ascending=[False, False, True])

	# Return the first row of each subgroup
	return [group.iloc[0] for group in groups]

def find_contig_by_id(contigs, contig_id):
	# Look for the sequence with the specified ID
	for contig in contigs:
		if (contig.id == contig_id):
			return contig
	return None

def handle_variant_alleles(df, contigs, output_path):
	current_contig = None
	wd = os.getenv('SKIRT_WD', os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
	for idx, group in df[(df['matching_rate'] > 0.99) & (df['matching_rate'] < 1)].groupby('target_name').apply(locus_grouping).items():
		for row in group:
			if (current_contig is None or current_contig.id != idx):
				current_contig = find_contig_by_id(contigs, idx)
			if (current_contig is not None and current_contig.id == idx):
				# Extract the desired substring
				genomic_seq = str(current_contig.seq[row.target_start:row.target_end])
				# get lists of coordinates of exons
				exon_lens, exon_starts, vars_start, vars_nuc = cstr_cutter(row.cs)
				exons = [genomic_seq[start:start+length] for start, length in zip(exon_starts, exon_lens)]
				coding_seq = "".join(exons)
				#print(coding_seq)
				extract_fasta_path = output_path + "." + current_contig.id.replace("#", "-") + "." + row.kir_allele.replace("*","-")  + ".fa"
				fp = open(extract_fasta_path, "w")
				fp.write(">" + row.kir_allele + "_CDS_" + row.query_name + "_" + current_contig.id + "\n")
				fp.write(coding_seq)
				fp.write("\n")
				fp.close()
				if (row.kir_gene not in ["2DP1", "3DP1"]):
					syn_or_closest, percent_identity = check_synonymous(extract_fasta_path, wd + "/IPDKIR/fasta/KIR" + row.kir_gene + "_prot.fasta")
					if (syn_or_closest is None):
						print("none matching", row.kir_allele)
						continue
					for i, local_row in df[(df['target_name']==idx) & (df['kir_gene']==row.kir_gene) &
							((abs(df['target_start'] - row.target_start) <= 3000) |
							(abs(df['target_end'] - row.target_end) <= 3000)) & (df['matching_rate'] < 1)].iterrows():
						len_coding = syn_or_closest.find('*')+4
						if (local_row.kir_allele[:len_coding] == syn_or_closest[:len_coding]):
							if (percent_identity == 1):
								df.at[i, 'resolution'] = 3
							elif (percent_identity > 0):
								df.at[i, 'resolution'] = 2
						else:
							if (local_row.kir_gene=="3DL2"):
								print("-" + local_row.kir_allele, "+" + syn_or_closest, percent_identity)
							df.drop(i, inplace=True)

	return df

def drop_redundant_alleles(df):
	for idx, group in df.groupby('target_name').apply(locus_grouping).items():
		for row in group:
			for i, local_row in df[(df['target_name']==idx) & (df['kir_gene']==row.kir_gene) &
					((abs(df['target_start'] - row.target_start) <= 3000) |
					(abs(df['target_end'] - row.target_end) <= 3000))].iterrows():
				if (local_row.kir_allele != row.kir_allele):
					df.drop(i, inplace=True)
	return df

def format_allele_info(kir):
	kir_gene = kir["kir_gene"]
	kir_allele = kir["kir_allele"]
	kir_resolution = kir["resolution"]
	kir_mat_rate = kir["matching_rate"]

	if (kir_resolution in [0,2] and kir_gene not in ["2DP1","3DP1"]):
		return kir_allele.split('*')[0] + "# (%.2f%% of " % (kir_mat_rate * 100) + kir["kir_cds_5"] + ")"
	elif (kir_resolution != 0):
		suffix = "=" if kir_mat_rate<1 and kir_resolution==3 else ("$" if kir_resolution==6 else ("+" if kir_resolution <=5 else ""))
		kir_resolution = (((kir_resolution-1)>>1)<<1)+1 # round down to 3,5,7-digit
		return kir_allele[:kir_allele.find('*')+kir_resolution+1] + suffix

	return kir_allele[:kir_allele.find('*')+6] + ("(%.2f%%" % (kir_mat_rate * 100) + ")" if kir_mat_rate<1 else "")

def sort_based_on_strand(group):
	group_pos = group[group['strand'] == '+'].sort_values('target_start')
	group_neg = group[group['strand'] == '-'].sort_values('target_start', ascending=False)
	return pd.concat([group_pos, group_neg])

def generate_haplotype_layout(df, start_asc = True):
	# Sort the input dataframe by 'target_name' and 'target_start'
	df = df.sort_values(by=["target_name", "target_start"], ascending=[True, start_asc])

	# Define the order of the KIR genes for the output dataframe
	ordered_kir_genes = [
		"3DL3", "2DS2", "2DL2", "2DL3", "2DL5B", "2DS3",
		"2DS5", "2DP1", "2DL1", "3DP1", "central region", "2DL4",
		"3DL1", "3DS1", "2DL5A", "2DS3", "2DS5", "2DS1", "2DS4", "3DL2"
	]

	# Create a new list of column names with suffix for duplicate names
	unique_kir_genes = []
	old_to_new_names = {}
	new_to_old_names = {}
	for gene in ordered_kir_genes:
		if gene in unique_kir_genes:
			new_name = gene + "t"
			unique_kir_genes.append(new_name)
			new_to_old_names[new_name] = gene
			old_to_new_names[gene] = new_name
		else:
			unique_kir_genes.append(gene)

	# Split the KIR genes into three parts
	centromeric_motif = ordered_kir_genes[:10]
	telomeric_motif = ordered_kir_genes[11:]

	centromeric_len = len(centromeric_motif)
	telomeric_len = len(telomeric_motif)

	# Group the input dataframe by 'target_name' and 'kir_gene'
	df_sorted = df.groupby("target_name").apply(sort_based_on_strand).reset_index(drop=True)
	kir_contigs = df_sorted.groupby("target_name")

	# Create an empty output dataframe with the desired columns
	result_df = pd.DataFrame(columns=["contig"] + unique_kir_genes)

	# Process each group in the sorted input dataframe
	for target_name, kir_contig in kir_contigs:
		row = {"contig": target_name}

		# Initialize the 'central region' column
		row["central region"] = ""
		central_region = []

		current_centro_idx = 0
		current_telo_idx = len(kir_contig) - 1

		# Fill the row with KIR alleles from centromeric
		checked_centromeric_genes = []
		for pos, (idx, kir) in enumerate(kir_contig.iterrows()):
			kir_gene = kir["kir_gene"]

			# Check if the KIR gene is in the first part of the ordered list
			if (kir_gene in centromeric_motif):
				# Check if this KIR gene can be added at this point
				if (kir_gene not in checked_centromeric_genes):
					for i in range(len(checked_centromeric_genes), centromeric_motif.index(kir_gene) + 1):
						checked_centromeric_genes.append(centromeric_motif[i])
					row[kir_gene] = format_allele_info(kir)
					current_centro_idx = pos + 1
				else:
					break
			else:
				break

		# Fill the row with KIR alleles from part3
		checked_telomeric_genes = []
		for pos, (idx, kir) in enumerate(kir_contig.iloc[::-1].iterrows()):
			kir_gene = kir["kir_gene"]

			# Check if the KIR gene is in the third part of the ordered list
			if (kir_gene in telomeric_motif):
				# Check if this KIR gene can be added at this point
				if (kir_gene not in checked_telomeric_genes):
					for i in range(telomeric_len - len(checked_telomeric_genes) - 1, telomeric_motif.index(kir_gene) - 1, -1):
						checked_telomeric_genes.insert(0, telomeric_motif[i])
					row[kir_gene if kir_gene not in list(new_to_old_names.values()) else old_to_new_names[kir_gene] ] = format_allele_info(kir)
					current_telo_idx = len(kir_contig) - pos - 2
				else:
					break
			else:
				break

			# Stop the loop at the end of centromeric
			if (current_telo_idx <= current_centro_idx-1):
				break

		# Fill the row with KIR alleles that don't belong to part1 or part3
		for _, kir in kir_contig.iloc[current_centro_idx:current_telo_idx+1].iterrows():
			# Add KIR alleles to the 'central region' separated by "~"
			central_region.append(format_allele_info(kir))
		row["central region"] = "~".join(central_region)

		# Add the processed row to the output dataframe
		result_df = pd.concat([result_df, pd.DataFrame([row])], ignore_index=True)

	# Rename the columns back to original names
	result_df = result_df.rename(columns=new_to_old_names)

	# Save the output dataframe to a CSV file with a header
	return result_df

def main():
	parser = argparse.ArgumentParser(description='Parse a PAF file and output a tab-separated table.')
	parser.add_argument("-b", help="Output BED file", action="store_true")
	parser.add_argument("-v", help="Output VCF file", action="store_true")
	parser.add_argument("-p", help="Output HAP file", action="store_false")
	parser.add_argument("-g", help="path to the input Genomic Alleles Mapping PAF file", nargs='?', const='', default=None)
	parser.add_argument("-asm", help="path to the input assembly fasta (gz) file", nargs='?', const='', default=None)
	parser.add_argument('paf_file', type=str, help='path to the input PAF file')
	parser.add_argument('output_file', type=str, help='path to the output file')
	args = parser.parse_args()

	cds_paf = PafParser(args.paf_file)
	cds_df = filter_overlapping_rows(cds_paf.df)
	genomic_paf = PafParser(args.g)
	result_df = check_genomic_match(cds_df, genomic_paf.df)
	genomic_paf.filter_perfect_match()
	if (args.asm is not None):
		# Use gzip to open the fasta file
		gzipped = 0
		with open(args.asm, 'rb') as f:
			gzipped = f.read(2) == b'\x1f\x8b'

		if gzipped:
			# File is gzipped, open using gzip
			with gzip.open(args.asm, 'rt') as f:
				# Use Biopython's SeqIO module to read the fasta file
				records = SeqIO.parse(f, 'fasta')
				result_df = handle_variant_alleles(result_df, records, args.output_file)
		else:
			# File is not gzipped, open normally
			with open(args.asm, 'r') as f:
				# Use Biopython's SeqIO module to read the fasta file
				records = SeqIO.parse(f, 'fasta')
				result_df = handle_variant_alleles(result_df, records, args.output_file)
	target_order = True if result_df.loc[result_df.first_valid_index(), 'strand']=="+" else False
	result_df = result_df.sort_values(['target_name','target_start','matching_rate','matching_length','target_end'], ascending=[True, target_order, False, False, target_order])
	result_df = drop_redundant_alleles(result_df)
	result_df = result_df.drop_duplicates(cds_paf.get_initial_columns())
	result_df.to_csv(args.output_file + ".allele.csv", columns=cds_paf.get_initial_columns()+['NM','cs','kir_allele','matching_rate','resolution'], sep=',')
	genomic_paf.df.to_csv(args.output_file + ".g.csv")
	df = result_df.reset_index(drop=True)
	hap_df = generate_haplotype_layout(df, True if df.loc[df.first_valid_index(), 'strand']=="+" else False)
	hap_df.to_csv(args.output_file + ".hap.csv", index=False)

	formater = FormatComposer(df, args.output_file)
	formater.output_bed_file()
	formater.output_vcf_file()

if __name__ == '__main__':
	main()

