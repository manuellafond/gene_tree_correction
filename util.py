import os
import sys
import generate_families_with_simphy as GFWS
import random
import argparse
from pathlib import Path
import shutil
from ete3 import Tree
from Bio import Phylo
import itertools
import re
from collections import defaultdict


def run_raxml(seqfile, raxmlbin, returned_collapsed_treefile = False):
    
    
    outfile = seqfile + ".raxml.bestTree"
    if returned_collapsed_treefile:
        outfile = seqfile + ".raxml.bestTreeCollapsed"
    
    
    command = f"{raxmlbin} --msa {seq_file} --model GTR+G --redo"

    print(f"Running {command}")

    if not args.skipphylo:
        os.system(command)
    
    return outfile
    
    
    
def run_iqtree(seqfile, bootstrap_replicates = 1000):
    command = f"iqtree -s {seq_file} -B {bootstrap_replicates} -nt AUTO -redo"

    print(f"Running {command}")

    os.system(command)




#calls astral-pro.  The script concatenates all tree_files into work_file
#The gene names are converted to species names, eg 25_10_5 becomes 25
def run_apro_from_symphy(tree_files, work_file, out_file):
    
    strout = ""
    
    for tfile in tree_files:
        if os.path.exists(tfile):
            with open(tfile, 'r') as file:
                content = file.read()
                tree = Tree(content)
            
                for leaf in tree.iter_leaves():
                    leaf.name = leaf.name.split("_")[0]
                
                strout += tree.write(format=1) + "\n"
        else:
            print("WARNING: " + tfile + " does not exist")

    with open(work_file, "w") as f:
        f.write(strout)
            

    command = f'java -D"java.library.path=/home/manuel/git/A-pro/ASTRAL-MP/lib" -jar /home/manuel/git/A-pro/ASTRAL-MP/astral.1.1.6.jar -i {work_file} -o {out_file}'
    print("Running " + command)
    
    os.system(command)







def make_generax_family_file_from_trees(out_filename, tree_files, alignment_files = []):
    

    out = open(out_filename, 'w')

    out.write("[FAMILIES]\n")

    for i in range(0, len(tree_files)):
        
        if os.path.exists(tree_files[i]):
            out.write(f"- family_{i + 1}\n")
            out.write(f"starting_gene_tree = {tree_files[i]}\n")
            
            if len(alignment_files) > i and os.path.exists(alignment_files[i]): 
                out.write(f"alignment = {alignment_files[i]}\n")
                out.write("subst_model = GTR+G\n")
        
        
    out.close()





def run_speciesrax(family_file, generax_outdir, generax_bin, use_strategy_spr = False):

    '''command = f"{generax_bin} --si-strategy HYBRID  -s MiniNJ -f {family_file} -p {generax_outdir} --rec-model UndatedDTL "
    
    
    if not use_strategy_spr:
        command += " --strategy SKIP "
    else:
        command += " --strategy SPR "
    
    #recommended in docs
    #command += " --prune-species-tree "
    

    
     
    #tests to see if it improves
    command += " --per-family-rates"
    #command += " --skip-family-filtering "
    #command += " --si-spr-radius 3 "  
    '''
    
    command = f"{generax_bin} -f {family_file} -p {generax_outdir} -s MiniNJ --si-strategy HYBRID --rec-model UndatedDTL --strategy SKIP --per-family-rates --skip-family-filtering"
    
    print(f"Running: {command}")
    
    os.system(' echo "' + command + '" >> commands.txt')
    
    os.system(command)







def write_to_file(filename, contents):
    with open(filename, 'w') as f:
        f.write(contents)




def make_dir(path, clear_if_exists = True):
    dir_path = Path(path)

    if clear_if_exists and dir_path.exists():
        shutil.rmtree(dir_path)

    dir_path.mkdir(parents=True, exist_ok=True)






#bins is the bin dict to fill
def count_support_bins(newick_file, bins = None, bin_size=10, all_supports_list = None):
    # Load the tree
    tree = Tree(newick_file, format = 1)

    if bins is None:
        bins = defaultdict(int)

    for node in tree.traverse():
        # Skip leaves (they usually don't have support values)
        if node.is_leaf():
            continue
        
        
        if node.name == "":
            continue
            
        support = float(node.name)    #node.support


        # Skip nodes with no support (optional: adjust if needed)
        if support is None:
            continue
            
        if not all_supports_list is None:
            all_supports_list.append(support)
            

        if support == 0:
            bins["0"] += 1
            bin_label = "0"
        else:
            # Determine bin (e.g., 83 -> 80–89)
            bin_start = int(support // bin_size) * bin_size
            if bin_start == 0:
                bin_start = 1
            bin_end = bin_start + bin_size - 1
            bin_label = f"{bin_start}-{bin_end}"

        

        bins[bin_label] += 1



def get_avg_support(treefile):
    
    supps = []
    count_support_bins(treefile, all_supports_list = supps)
    
    #if len(supps) == 0:
    #    return 0
    
    return float(sum(supps)) / float(len(supps))





#takes a tree in newick, makes nodes with 0 bootstrap to 1
def fix_zero_bootstrap(input_file, output_file):

    t = Tree(input_file, format=1)

    for node in t.traverse():
        if not node.is_leaf():

            if node.name != "" and float(node.name) == 0:
                node.name = "1"
    
    
    t.write(outfile=output_file, format=1)

    #newick_str = t.write(format=1)
    # replace scientific notation with fixed-point
    # disclaimer: chatgpt's solution
    #def sci_to_fixed(match):
    #    return "{:.6f}".format(float(match.group()))

    # regex matches numbers like 1e-06
    #newick_fixed = re.sub(r"\d+\.?\d*e[+-]?\d+", sci_to_fixed, newick_str, flags=re.IGNORECASE)

    # save to file
    #with open(output_file, "w") as f:
    #    f.write(newick_fixed)







def get_RF(treefile1, treefile2, normalize = True, unrooted = True):
    t1 = Tree(treefile1)
    t2 = Tree(treefile2)
    rf, max_rf, *_ = t1.robinson_foulds(t2, unrooted_trees=unrooted)

    if (normalize):
        return rf / max_rf
    return rf
    
    
    
def ete3_rf(treefile1, treefile2, unrooted = True, get_all_return = False, format = 1):

    tree1 = Tree(treefile1, format = 1)
    tree2 = Tree(treefile2, format = 1)
    #if (len(tree2.children) == 3):
    #  tree2.set_outgroup(tree2.children[0])
    #if (len(tree1.children) == 3):
    #  tree1.set_outgroup(tree1.children[0])

    leaves1 = {leaf.name for leaf in tree1.iter_leaves()}
    leaves2 = {leaf.name for leaf in tree2.iter_leaves()}

    if leaves1 != leaves2:
        print("WARNING: trees have different leaves")


    res = tree1.robinson_foulds(tree2, unrooted_trees=unrooted)   #, skip_large_polytomies = True, correct_by_polytomy_size = True)
    
    if get_all_return:
        return res
    else:
        return float(res[0]) / float(res[1])
    






def get_average_branch_length(tree_file, include_root=False):
    """
    Return the average branch length of all branches in a Newick tree file.
    Disclaimer: this is chatgpt

    Parameters
    ----------
    tree_file : str
        Path to the Newick tree file.
    include_root : bool, default False
        Whether to include the root branch length (if present).
        In most rooted trees, the root has dist = 0 and is usually excluded.

    Returns
    -------
    float
        Mean branch length. Returns 0.0 if no branches are found.
    """
    tree = Tree(tree_file, format=1)

    branch_lengths = []
    for node in tree.traverse():
        # Skip the root unless explicitly requested
        if node.is_root() and not include_root:
            continue
        branch_lengths.append(node.dist)

    if not branch_lengths:
        return 0.0

    return sum(branch_lengths) / len(branch_lengths)



def get_tree_height(tree_file):
    """
    Return the tree height: the maximum root-to-leaf distance
    (sum of branch lengths along the path).

    Parameters
    ----------
    tree_file : str
        Path to the Newick tree file.

    Returns
    -------
    float
        Longest root-to-leaf path length.
    """
    tree = Tree(tree_file, format = 1)

    # get_distance() automatically sums branch lengths
    return max(tree.get_distance(leaf) for leaf in tree.iter_leaves())







def count_leaves(tree_file):
    tree = Tree(tree_file, format=1)
    return len(tree)



def first_word_from_file(filename):
    with open(filename, "r", encoding="utf-8") as f:
        words = f.readline().split()
    return words[0] if words else None



def run_iqtree_on_all(alignment_files, output_dir, nb_bootstrap = 1000, skip_existing = True):

    #no outgroup??


    make_dir(output_dir)
    
    aldir = os.path.join(output_dir, "alignments")
    
    
    make_dir(aldir)
    
    for afile in alignment_files:
       afile_path = Path(afile)
       afilecopy_filename = os.path.join(aldir,  afile_path.name)
       shutil.copy(afile, afilecopy_filename)
       
       
       command = f"iqtree -s {afilecopy_filename} -nt 4 -m JC -B {nb_bootstrap} --redo"
    
    
       treefilename = afilecopy_filename + ".treefile"
       
       if skip_existing and os.path.exists(treefilename):
           continue
           
       os.system(command)
       

    
    



def run_pargenes(alignment_files, output_dir, raxmlbin, raxml_args = '--model GTR+G --blopt nr_safe --redo', use_model_test = False):
    
    make_dir(output_dir)
    
    aldir = os.path.join(output_dir, "alignments")
    treesdir = os.path.join(output_dir, "trees")
    
    
    
    make_dir(aldir)
    
    for afile in alignment_files:
       afile_path = Path(afile)
       shutil.copy(afile, os.path.join(aldir,  afile_path.name))
       
    if use_model_test:
        command = f"python /home/manuel/git/ParGenes/pargenes/pargenes.py -a {aldir} -o {treesdir} -c 10 -d nt --use-modeltest -R '{raxml_args}' --raxml-binary {raxmlbin}"
    else: 
       #command = f"python /home/manuel/git/ParGenes/pargenes/pargenes.py -a {aldir} -o {treesdir} -c 6 -d nt -R '{raxml_args}' --raxml-binary {raxmlbin}"
        command = f"python /home/manuel/git/ParGenes/pargenes/pargenes.py -a {aldir} -o {treesdir} -b 0 -c 5 -s 0 -p 1 -d nt -R '{raxml_args}' --raxml-binary {raxmlbin}"
    print(f"Running {command}")
    
    os.system(' echo "' + command + '" >> commands.txt')

    os.system(command)



def run_pargenes_on_all(simphy_dir, output_dir, raxmlbin, pattern = "TRUE.phy"):
    source_dir = Path(simphy_dir)

    alignment_files = []


    for file_path in source_dir.iterdir():
        if file_path.is_file() and "TRUE.phy" in file_path.name:
            alignment_files.append(file_path)
            
    run_pargenes(alignment_files, output_dir, raxmlbin)






def clean_simphy_species_tree(input_file, output_file):
    # Read the tree
    tree = Phylo.read(input_file, "newick")

    # Remove single quotes from the entire file content
    with open(input_file, "r") as f:
        newick_str = f.read().replace("'", "")

    # Remove internal node names
    for clade in tree.find_clades():
        #clade.branch_length = None
        #clade.confidence = None
        if not clade.is_terminal():  # If it's an internal node
            clade.name = None  # Remove the name

    # Write back the cleaned tree
    Phylo.write(tree, output_file, "newick")

    # Ensure the single quotes are removed by re-processing the output file
    # ML doesn't know why this is there...
    with open(output_file, "r") as f:
        cleaned_str = f.read().replace("'", "")
        #cleaned_str = cleaned_str.replace(":0.00000", "")


    with open(output_file, "w") as f:
        f.write(cleaned_str)


def is_bad_simphy_dir(simphy_dir):
    simphy_species_tree_file = os.path.join(simphy_dir, "1", "s_tree.trees")
    
    if not os.path.exists(simphy_species_tree_file):
        return True
       
       
 
    if os.path.getsize(simphy_species_tree_file) < 1:
        return True


    return False







def parse_sim_params_from_dir(dirname):
    """
    Extract X, A, B, C, P from strings of the form:

    ssim_dtl_s25_f100_sitesX_GTR_bl1.0_dA_lB_tC_gc0.0_p0.0_popP_ms0.0_mf0.0_seed3001

    where X, A, B, C, and P are strings that do not contain underscores.

    Returns
    -------
    tuple[str, str, str, str, str]
        (X, A, B, C, P)
    """
    pattern = (
        r"sites([^_]+)"   # X
        r"_GTR_bl1\.0"
        r"_d([^_]+)"      # A
        r"_l([^_]+)"      # B
        r"_t([^_]+)"      # C
        r"_gc0\.0_p0\.0"
        r"_pop([^_]+)"    # P
        r"_"
    )

    m = re.search(pattern, dirname)
    if not m:
        raise ValueError(f"Could not parse: {dirname}")

    X, A, B, C, P = m.groups()
    
    
    #get seed cuz I'm bad with regex
    seed = dirname.split("_")[-1].replace("seed", "")
    
    return X, A, B, C, P, seed











