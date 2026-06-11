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

parser = argparse.ArgumentParser(description="Script that calls speciesrax on the given family file and computes RF with the true tree.")


parser.add_argument(
    "-d", "--dir",
    help="Directory"
)


args = parser.parse_args()

outdir = "tmp"
util.make_dir(outdir)


#1: RUN PARGENES

gene_tree_files = []
alignment_files = []
for i in range(1, 101):
    alignment_file = os.path.join(args.dir, "1", f"dataset_{i:03d}_TRUE.phy")
    alignment_files.append(alignment_file)

#INFER GENE TREES WITH PARGENES
for i in range(1, 101):
    treefile = os.path.join(outdir, "pargenes/trees/mlsearch_run/results", f"dataset_{i:03d}_TRUE_phy", f"dataset_{i:03d}_TRUE_phy.raxml.bestTree")
    gene_tree_files.append( treefile )


raxml_args = " --model GTR+G4 --redo --extra bs-start-rand "

util.run_pargenes(alignment_files, os.path.join(outdir, "pargenes"), "/usr/bin/raxml-ng", raxml_args)




#2: RUN SPECIESRAX

family_file = os.path.join(outdir, "speciesrax_family.txt")
true_family_file = os.path.join(outdir, "speciesrax_family_true_trees.txt")
true_tree_files = []
for i in range(1, 101):
    true_tree_files.append( os.path.join(args.dir, "1", f"g_trees{i:03d}.trees") )

#make both true and inferred to test
util.make_generax_family_file_from_trees(family_file, gene_tree_files)    
util.make_generax_family_file_from_trees(true_family_file, true_tree_files)


generax_bin = "~/git/GeneRax/build/bin/generax"


generax_outdir = os.path.join(outdir, "speciesrax")

command = f"{generax_bin} --si-strategy HYBRID  -s MiniNJ -f {family_file} -p {generax_outdir} --rec-model UndatedDTL "
command += " --strategy SKIP "

#command += " --prune-species-tree "


#tests to see if it improves
command += " --per-family-rates"
#command += " --si-spr-radius 5  --si-small-root-radius 10  --si-big-root-radius 20 "  

print(f"Running: {command}")


os.system(command)


simphy_species_tree_file = os.path.join(args.dir, "1", "s_tree.trees")
simphy_sptree_file_cleaned = simphy_species_tree_file + ".cleaned.newick"

generax_sptreefile = os.path.join(generax_outdir, "species_trees", "inferred_species_tree.newick")
generax_sptreefile_cleaned = generax_sptreefile + ".cleaned.newick"


util.clean_simphy_species_tree(simphy_species_tree_file, simphy_sptree_file_cleaned)

util.clean_simphy_species_tree(generax_sptreefile, generax_sptreefile_cleaned)

rf = util.get_RF(simphy_sptree_file_cleaned, generax_sptreefile_cleaned)

print(f"RF={rf}")



