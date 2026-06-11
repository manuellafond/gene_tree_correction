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



class ParGenesParams():
  def __init__(self):
    self.pargenesbin = "/home/manuel/git/ParGenes/pargenes/pargenes.py"  #unless your name is also manuel, you probably need to change this
    self.alignment_files = []  #list of full paths to alignment files
    self.output_dir = ""
    self.raxmlbin = "/usr/bin/raxml-ng"  #full path to raxml-ng binary
    self.model = "GTR+G"
    self.fast_bootstrap = False
    self.bootstrap_replicates = 100
    self.raxml_other_args = ""
    self.show_command = True
    self.log_command_in_file = None
    self.cores = 5    #how many cores will pargenes use?
    self.start_from_random_trees = None
    
    
    
def run_pargenes(pargenes_params):
    
    util.make_dir(pargenes_params.output_dir)
    
    aldir = os.path.join(pargenes_params.output_dir, "alignments")
    treesdir = os.path.join(pargenes_params.output_dir, "trees")
    
    util.make_dir(aldir)
    
    #pargenes needs all alignment files in one place
    for afile in pargenes_params.alignment_files:
       if os.path.exists(afile):
           afile_path = Path(afile)
           shutil.copy(afile, os.path.join(aldir,  afile_path.name))
       
    command = f"python {pargenes_params.pargenesbin} -a {aldir} -o {treesdir} -c {pargenes_params.cores} -s 0 -p 1 -d nt --raxml-binary {pargenes_params.raxmlbin}"
    
    
    #--model GTR+G --blopt nr_safe --redo
    raxml_args = f"--model {pargenes_params.model} --blopt nr_safe --redo"
    
    if pargenes_params.fast_bootstrap:
        command += f" -b {pargenes_params.bootstrap_replicates} "
        #raxml_args += " -f a "


    if not pargenes_params.start_from_random_trees is None:
        raxml_args += " --tree rand{" + str(pargenes_params.start_from_random_trees) + "} "
    
    command += f" -R '{raxml_args}'"

    
    if pargenes_params.show_command:
        print(f"Running {command}")
    
    if not pargenes_params.log_command_in_file is None:
        os.system(' echo "' + command + '" >> commands.txt')

    os.system(command)
     

