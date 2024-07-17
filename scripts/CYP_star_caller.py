"""
python CYP_star_caller.py cyp2d6 /mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/dp_vcf/NA19239.dv.small.phased.vcf.gz chr22:42126499-42130865
"""

import copy
import datetime
import operator
import os
import statistics
import subprocess
import sys
import timeit
import types
import shutil
import logging
import gzip
import pkgutil
import time

import pandas as pd
import numpy as np

class VCF:
    def __init__(self):
        self.meta = []
        self.header = []
        self.data = []

def copy_vcf(original_vcf, items):
    copied_vcf = VCF()
    if 'meta' in items:
        copied_vcf.meta = copy.deepcopy(original_vcf.meta)
    if 'header' in items:
        copied_vcf.header = copy.deepcopy(original_vcf.header)
    if 'data' in items:
        copied_vcf.data = copy.deepcopy(original_vcf.data)
    return copied_vcf

class SNPAllele:
    def __init__(self):
        self.pos = '' # reference genome position
        self.wt = '' # wild type (*1) allele
        self.var = '' # variant allele
        self.rs = '' # rs ID
        self.het = False # heterozygous
        self.ad = 0 # allelic depth
        self.td = 0 # total depth
        self.n = '' # SNP table number
        self.hg = '' # reference genome allele
        self.so = '' # sequence ontology
        self.effect = '' # coding effect
        self.impact = '' # variant impact
        self.rev = False # reverting variant

    @property
    def key(self):
        return (self.pos, self.wt, self.var)

    @property
    def af(self):
        return 0 if self.td == 0 else self.ad / self.td

    def __eq__(self, other):
        return self.key == other.key

    def __hash__(self):
        return hash(self.key)

    def summary(self):
        return '<{}:{}>{}:{}/{}:{:.2f}:{}:{}:{}>'.format(self.pos, self.wt, self.var, self.ad, self.td, self.af, self.so, self.impact, self.effect)


class StarAllele:
    def __init__(self):
        self.name = ''
        self.score = -100.0
        self.core = []
        self.tag = []
        self.sv = ''

    @property
    def ranked_as(self):
        """
        Unknown function alleles should be broken ties with normal function
        alleles using attributes other than activity score. Increased
        function alleles should come before normal function alleles.
        """

        if self.score < 0:
            return 1.0
        elif self.score > 1:
            return 0.99
        else:
            return self.score

    @property
    def rank(self):
        return (self.ranked_as, -1 * int(bool(self.sv)), -1 * len(self.core))

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)


def get_snp_list(target_gene, snp_table, genome_build):
    # df = pd.read_table(f"{program_dir}/snp_table.tsv")
    df = pd.read_table(snp_table)
    df = df[df['gene'] == target_gene]
    def func(row):
        snp_allele = SNPAllele()
        snp_allele.n = row['sg_id']
        snp_allele.effect = row['functional_effect']
        snp_allele.pos = row[f'{genome_build}_pos']
        snp_allele.id = row['rs_id']
        snp_allele.hg = row[f'{genome_build}_allele']
        snp_allele.var = row['var_allele']
        snp_allele.wt = row['wt_allele']
        snp_allele.so = row['sequence_ontology']
        snp_allele.impact = row['variant_impact']
        snp_allele.rev = row[f'{genome_build}_revertant'] == 'revertant'
        return snp_allele
    snp_alleles = df.apply(func, axis=1).to_list()

    return snp_alleles


def get_star_dict(target_gene, snp_list, program_dir, genome_build, star_table):
    #df = pd.read_table(f"{program_dir}/star_table.tsv")
    df = pd.read_table(star_table)
    df = df[df["gene"] == target_gene]
    def func(row):
        star_allele = StarAllele()
        star_allele.name = row["name"]
        star_allele.score = row["score"]
        if row[f"{genome_build}_core"] in ['.', 'ref']:
            star_allele.core = []
        else:
            star_allele.core = copy.deepcopy([x for x in snp_list if '{}:{}>{}'.format(x.pos, x.wt, x.var) in row[f"{genome_build}_core"].split(',')])

        if row[f"{genome_build}_tag"] == '.':
            star_allele.tag = []
        else:
            star_allele.tag = copy.deepcopy([x for x in snp_list if '{}:{}>{}'.format(x.pos, x.wt, x.var) in row[f"{genome_build}_tag"].split(',')])

        if row["sv"] == '.':
            star_allele.sv = ''
        else:
            star_allele.sv = row["sv"]

        return star_allele

    star_alleles = df.apply(func, axis=1).to_list()

    star_dict = {}

    for star_allele in star_alleles:
        print (star_allele.name, star_allele.core, star_allele.tag, star_allele.sv)
        star_dict[star_allele.name] = star_allele


    return star_dict

class Haplotype:
    def __init__(self):
        self.cand = []
        self.obs = []

    @property
    def sv(self):
        sv = 'no_sv'
        sv_list = []
        for star_allele in self.cand:
            if star_allele.sv and star_allele.sv not in sv_list:
                sv_list.append(star_allele.sv)
        if len(sv_list) > 1:
            raise ValueError('haplotype contains multiple structural variant calls')
        if len(sv_list) == 1:
            sv = sv_list[0]
        return sv

    @property
    def af(self):
        return [0] if not self.obs else [x.af for x in self.obs]

    @property
    def af_mean_main(self):
        filtered = [x.af for x in self.obs if x.td > 10 and x.het and x in [y for y in self.cand[0].core]]
        return -1 if not filtered else statistics.mean(filtered)

    def af_mean_gene(self, start, end):
        filtered = [x.af for x in self.obs if x.td > 10 and x.het and start <= x.pos <= end]
        return -1 if not filtered else statistics.mean(filtered)

    def fit_data(self, total_cn, start, end):
        """Return the fit MAF and CN."""
        maf_choices = []
        for i in range(1, total_cn):
            maf_choices.append(i / total_cn)
        fit_maf = min(maf_choices, key = lambda x: abs(x - self.af_mean_gene(start, end)))
        fit_cn = maf_choices.index(fit_maf) + 1
        return fit_maf, fit_cn

    def remove_star(self, sX):
        """Remove the given star allele from the candidates list."""
        for i, star in enumerate(self.cand):
            if star.name == sX.name:
                del self.cand[i]
                break

    def add_dup(self, cn):
        """Duplicate the main star allele by the given CN."""
        if cn == 1:
            return
        if cn > 10:
            cn = 10
        sX = self.cand[0]
        score = sX.score * cn
        name = sX.name + 'x' + str(cn)
        sY = StarAllele()

        sY.name = name; sY.score = score; sY.core = copy.deepcopy(sX.core); sY.sv = 'cnv{}'.format(cn)

        self.cand.insert(0, sY)
        self.remove_star(sX)


class Person:
    def __init__(self):
        self.name = '' # sample ID
        self.gt = False # true if genotyped
        self.sv = ['', ''] # SV calls
        self.pt = '' # predicted phenotype
        self.ssr = '' # sum of squared residuals
        self.dip_cand = [] # candidate stars
        self.hap = [Haplotype(), Haplotype()]
        self.bad = False # true if QC failed
        self.af_df = None

        

    def get_status(self):
        return 'g' if self.gt else 'ng'

    def get_hap1_main(self):
        if self.gt:
            return self.hap[0].cand[0].name
        else:
            return '.'

    def get_hap2_main(self):
        if self.gt:
            return self.hap[1].cand[0].name
        else:
            return '.'

    def get_hap1_cand(self):
        return ','.join([x.name for x in self.hap[0].cand])

    def get_hap2_cand(self):
        return ','.join([x.name for x in self.hap[1].cand])

    def get_hap1_score(self):
        if self.gt:
            return str(self.hap[0].cand[0].score)
        else:
            return '.'

    def get_hap2_score(self):
        if self.gt:
            return str(self.hap[1].cand[0].score)
        else:
            return '.'

    def get_dip_score(self):
        if self.gt:
            return str(self.hap[0].cand[0].score + self.hap[1].cand[0].score)
        else:
            return '.'

    def get_phenotype(self):
        if self.gt:
            return self.pt
        else:
            return '.'

    def get_dip_sv(self):
        return ','.join(self.sv)

    def get_hap1_sv(self):
        if self.gt:
            return self.hap[0].sv
        else:
            return '.'

    def get_hap2_sv(self):
        if self.gt:
            return self.hap[1].sv
        else:
            return '.'

    def get_ssr(self):
        return self.ssr

    def get_dip_cand(self):
        return ','.join([x.name for x in self.dip_cand])

    def get_hap1_main_core(self):
        if self.gt:
            return [x.summary() for x in self.hap[0].obs if x in self.hap[0].cand[0].core]
        else:
            return '.'

    def get_hap2_main_core(self):
        if self.gt:
            return [x.summary() for x in self.hap[1].obs if x in self.hap[1].cand[0].core]
        else:
            return '.'

    def get_hap1_main_tag(self):
        if self.gt:
            return [x.summary() for x in self.hap[0].obs if x in self.hap[0].cand[0].tag]
        else:
            return '.'

    def get_hap2_main_tag(self):
        if self.gt:
            return [x.summary() for x in self.hap[1].obs if x in self.hap[1].cand[0].tag]
        else:
            return '.'


    def get_hap1_af_mean_gene(self, target_locus):
        if self.gt:
            return "{:.2f}".format(self.hap[0].af_mean_gene(target_locus.gene_start, target_locus.gene_end))
        else:
            return '.'

    def get_hap2_af_mean_gene(self, target_locus):
        if self.gt:
            return "{:.2f}".format(self.hap[1].af_mean_gene(target_locus.gene_start, target_locus.gene_end))
        else:
            return '.'

    def get_hap1_af_mean_main(self):
        if self.gt:
            return "{:.2f}".format(self.hap[0].af_mean_main)
        else:
            return '.'

    def get_hap2_af_mean_main(self):
        if self.gt:
            return "{:.2f}".format(self.hap[1].af_mean_main)
        else:
            return '.'



    def summarize(self, target_locus):
        list2str = lambda x: '.' if not x else ','.join([str(x) for x in x])

        results = [
            self.name,
            self.get_status(),
            self.get_hap1_main(),
            self.get_hap2_main(),
            self.get_hap1_cand(),
            self.get_hap2_cand(),
            self.get_hap1_score(),
            self.get_hap2_score(),
            self.get_dip_score(),
            self.get_phenotype(),
            self.get_dip_sv(),
            self.get_hap1_sv(),
            self.get_hap2_sv(),
            self.get_ssr(),
            self.get_dip_cand(),
            list2str(self.get_hap1_main_core()),
            list2str(self.get_hap2_main_core()),
            list2str(self.get_hap1_main_tag()),
            list2str(self.get_hap2_main_tag()),
            self.get_hap1_af_mean_gene(target_locus),
            self.get_hap2_af_mean_gene(target_locus),
            self.get_hap1_af_mean_main(),
            self.get_hap2_af_mean_main(),
        ]
        
        return results


def vcf2samples(vcf):
    samples = []
    for name in vcf.header[9:]:
        sample = Person()
        sample.name = name
        i = vcf.header.index(name)
        for fields in vcf.data:
            
            pos, rs, ref, alt, inf, fmt = int(fields[1]), fields[2], fields[3], fields[4].split(','), fields[7].split(';'), fields[8]
            print (pos, rs, ref, alt, inf, fmt)

            # if not any(['PS=D' in x for x in inf]):
            #     continue

            #gt = [int(x) for x in fields[i].split(":")[0].split("|")]
            gt = [x for x in fields[i].split(':')[0].split('|')]
            if gt[0] == '0/0':
                continue
            if gt[0] == '1/1':
                gt = [1, 1]
            else:
                gt = [int(y) if y != "." else y for y in gt]
            print (gt)

            al = [ref] + alt
            # vi_list = ['no_change'] + [x for x in inf if 'VI=' in x][0].replace('VI=', '').split(',')
            # so_list = ['no_change'] + [x for x in inf if 'SO=' in x][0].replace('SO=', '').split(',')
            # fe_list = ['no_change'] + [x for x in inf if 'FE=' in x][0].replace('FE=', '').split(',')

            vi_list = ['no_change', 'no_change']
            so_list = ['no_change', 'no_change']
            fe_list = ['no_change', 'no_change']

            for j in [0, 1]:

                snp = SNPAllele()
                if gt[j] != ".":
                    #print(pos,gt[j])
                    #print(gt)
                    snp.pos, snp.wt, snp.var, snp.rs, snp.het, snp.so, snp.impact, snp.effect = pos, ref, al[gt[j]], rs, gt[0] != gt[1], so_list[gt[j]], vi_list[gt[j]], fe_list[gt[j]]

                if 'AD' in fmt:
                    ad_ind=fmt.split(":").index("AD")
                    ad_list=[x for x in fields[i].split(':')[ad_ind].split(',')]
                    ad_list = [int(y) if y != "." else 0 for y in ad_list]
                    #ad_list = [int(x) for x in fields[i].split(':')[1].split(',')]
                    #print(ad_list)

                    if not ad_list:
                         snp.ad=0; snp.td=0
                    elif gt[j] == ".":
                         snp.ad=0; snp.td=0
                    else:
                     #    print(gt[j],ad_list)
                         snp.ad = ad_list[gt[j]]; snp.td = sum(ad_list)

                sample.hap[j].obs.append(snp)
        samples.append(sample)
    return samples



def sort_star_names(names):
    def f(x):
        cn = 1
        if '*' not in x or x == '*DEL':
            n = 999
        else:
            n = int(''.join([y for y in x.split('+')[0].split('x')[0] if y.isdigit()]))
            if 'x' in x.split('+')[0]:
                cn = int(x.split('+')[0].split('x')[1])
        return (n, cn, len(x))

    return sorted(names, key = f)


def parse_region(region):
    return {'chr': region.split(':')[0].replace('chr', ''), 'start': int(region.split(':')[1].split('-')[0]), 'end': int(region.split(':')[1].split('-')[1])}

def read_vcf_simple(file):
    f = gzip.open(file, 'rt') if '.gz' in file else open(file)
    vcf = VCF()
    for line in f:
        if '##' in line:
            vcf.meta.append(line)
            continue
        fields = line.strip().split('\t')
        if fields[0] == '#CHROM':
            vcf.header = fields
            continue
        chr = fields[0].replace('chr', '')
        vcf.data.append([chr] + fields[1:])
    f.close()
    return vcf

def read_vcf_region(file, region):
    vcf = VCF()
    region_dict = parse_region(region)
    f = gzip.open(file, 'rt') if '.gz' in file else open(file)
    for line in f:
        if '##' in line:
            vcf.meta.append(line)
            continue
        fields = line.strip().split('\t')
        if fields[0] == '#CHROM':
            vcf.header = fields
            continue
        chr, pos = fields[0].replace('chr', ''), int(fields[1])
        if chr != region_dict['chr'] or pos < region_dict['start']:
            continue
        if pos > region_dict['end']:
            break
        vcf.data.append([chr] + fields[1:])
    f.close()
    return vcf

def read_sv_table(fn):
    result = {}

    with open(fn) as f:
        header = next(f).strip().split()
        for line in f:
            #print(line)
            fields = line.strip().split()
            gene = fields[0]
            name = fields[2]

            if gene not in result:
                result[gene] = {}

            result[gene][name] = dict(zip(header, fields))

    return result

# Import the data files.
snp_table_file = sys.path[0] + "/../CYP_data/snp_table.tsv"  #"/mnt/d/HLAPro_backup/Nanopore_optimize/stargazer-grc38-v.2.0.2/stargazer/snp_table.tsv"
star_table_file = sys.path[0] + "/../CYP_data/star_table.tsv"
sv_table = read_sv_table(sys.path[0] + "/../CYP_data/sv_table.tsv")

# gene = "cyp2d6"
# phased_vcf = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38/dp_vcf/NA19239.dv.small.phased.vcf.gz"
# region = 'chr22:42126499-42130865'

gene = sys.argv[1]
phased_vcf = sys.argv[2]
region = sys.argv[3]


ref= 'grc38'

gene = gene.lower()


snp_list = get_snp_list(gene, snp_table_file, ref)
# print (snp_list)
star_dict = get_star_dict(gene, snp_list,
    '', ref, star_table_file)
# print (star_dict)

sv_dict = sv_table[gene]

input_vcf = read_vcf_region(phased_vcf, region)
# print (input_vcf)
persons = vcf2samples(input_vcf)
for sample in persons:
    sample.sv = ['no_sv', 'no_sv']
    sample.ssr = '.'

    # print (sample.name, sample.hap[0].obs, sample.hap[1].obs)
    # f = lambda x: sorted([v for k, v in star_dict.items() if set(v.core).issubset(x) and not (v.sv and v.sv not in sample.sv)], key = lambda x: x.rank)

    def f(x):
        filtered_stars = []
        for k, v in star_dict.items():
            if set(v.core).issubset(x) and (not v.sv or v.sv in sample.sv):
            # if len(set(v.core) & set(x)) >= round(len(v.core) * 1)  and (not v.sv or v.sv in sample.sv):
                # for y in x:
                #     print (len(x), y.pos, y.wt, y.var)
                print (v.name, v.rank)
                filtered_stars.append(v)
        return sorted(filtered_stars, key=lambda star: star.rank)

    hap1_snp = [x for x in sample.hap[0].obs if x.wt != x.var]
    hap2_snp = [x for x in sample.hap[1].obs if x.wt != x.var]
    sample.hap[0].cand = f(hap1_snp)
    sample.hap[1].cand = f(hap2_snp)
    sample.dip_cand = f(list(set(hap1_snp + hap2_snp)))

    # Remove extra *1 alleles.
    def f(l,target_locus):
        if target_locus=="cyp2c19":
            wt_star="*38"
        else:
            wt_star="*1"
        if len(l) == 1:
            return
        for i in reversed(range(len(l))):
            if l[i].name == wt_star:
                del l[i]

    for sample in persons:
        f(sample.hap[0].cand,gene)
        f(sample.hap[1].cand,gene)
        f(sample.dip_cand,gene)

    # Order the haplotypes.
    for sample in persons:
        if not sample.gt:
            continue
        if sort_star_names([sample.hap[0].cand[0].name, sample.hap[1].cand[0].name])[0] == sample.hap[1].cand[0].name:
            sample.hap[0], sample.hap[1] = sample.hap[1], sample.hap[0]

    # Predict the phenotype.
    for person in persons:
        print(person.hap[0].cand[0].name,person.hap[1].cand[0].name,person.pt)


    # for k, v in star_dict.items():
    #     for y in v.core:
    #         print (v.name, y.pos, y.wt, y.var)


        