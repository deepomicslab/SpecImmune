import os
import re
import numpy as np
import sys

class My_folder():

    def __init__(self, args):
        self.root_dir = args["o"]
        self.outdir = args["o"] + "/" + args["n"]
        self.sample_prefix = self.outdir + "/" + args["n"]
        
    def make_dir(self):
        if not os.path.exists(self.root_dir):
            os.system("mkdir %s"%(self.root_dir))
        if not os.path.exists(self.outdir):
            os.system("mkdir %s"%(self.outdir))
