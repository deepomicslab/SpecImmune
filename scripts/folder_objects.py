import os
import re
import numpy as np
import sys

class My_folder():

    def __init__(self, args):
        root_dir = args["o"]
        ID = args["n"]

        self.root_dir = root_dir
        self.outdir = root_dir + "/" + ID
        self.sample_prefix = self.outdir + "/" + ID
        self.genes_dir = self.outdir + "/Genes/"
        self.step2_genes_dir = self.outdir + "/Genes_step2/"
        self.sequence_dir = self.outdir + "/Sequences/"
        self.reads_dir = self.outdir + "/Reads/"
        self.for_viz_dir = self.outdir + "/For_viz/"
        self.visualization_dir = self.outdir + "/Visualization/"  # for visualization
        self.gene_sample_prefix = self.genes_dir + "/" + ID

    def make_dir(self):
        if not os.path.exists(self.root_dir):
            os.system("mkdir %s"%(self.root_dir))
        if not os.path.exists(self.outdir):
            os.system("mkdir %s"%(self.outdir))
        if not os.path.exists(self.genes_dir):
            os.system("mkdir %s"%(self.genes_dir))
        if not os.path.exists(self.sequence_dir):
            os.system("mkdir %s"%(self.sequence_dir))
        if not os.path.exists(self.reads_dir):
            os.system("mkdir %s"%(self.reads_dir))
        if not os.path.exists(self.step2_genes_dir):
            os.system("mkdir %s"%(self.step2_genes_dir))
        if not os.path.exists(self.for_viz_dir):
            os.system("mkdir %s"%(self.for_viz_dir))
        if not os.path.exists(self.visualization_dir):
            os.system("mkdir %s"%(self.visualization_dir))
            

