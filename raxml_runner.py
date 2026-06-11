import os
import sys
import random
import argparse
from pathlib import Path
import shutil
from ete3 import Tree
from Bio import Phylo
import itertools
import util



class raxmlParams():
  def __init__(self):
    self.raxmlbin = "/usr/bin/raxml-ng"  #full path to raxml-ng binary
    self.alignment_files = []  #list of full paths to alignment files
    self.output_dir = ""

    self.model = "GTR+G"
    self.fast_bootstrap = False
    self.bootstrap_replicates = 100
    self.raxml_other_args = ""
    self.show_command = True
    self.log_command_in_file = None
    self.start_from_random_trees = None
    
    
    
def run_raxml(params):
    
    util.make_dir(params.output_dir)
    
    treesdir = os.path.join(params.output_dir, "trees")
       
    for afile in params.alignment_files:
        command = f"raxml-ng --all --msa {afile} --model {params.model} --threads auto --prefix {params.output_dir}/{Path(afile).name} --redo "
        
        if params.start_from_random_trees:
            command += " --tree rand{" + str(params.start_from_random_trees) + "}" 
            
        if params.fast_bootstrap:
            command += f" --bs-trees {params.bootstrap_replicates}"

      
        if params.show_command:
            print(f"Running {command}")


    
        if not params.log_command_in_file is None:
            os.system(' echo "' + command + '" >> ' + params.log_command_in_file + '"')

        os.system(command)
     




