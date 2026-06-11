import os
import sys
from pathlib import Path
import re
import util
from collections import Counter
from ete3 import Tree
import argparse
import math


#Example usage
# python make_genetree_csv.py -d=my_exp_p1e7_dl2_t1_s300/ -o test.csv

parser = argparse.ArgumentParser(description="Goes through a directory created with run_exp.py, and outputs a csv file with all the info you might need.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)



parser.add_argument(
    "-d", "--dir",
    help="The directory created by run_exp.py (that you provided with -o)."
)



parser.add_argument(
    "-o", "--outfile",
    help="Name of output csv file."
)





parser.add_argument(
    "--append",
    action="store_true",
    help="Set to append rows to the specified csv file.  Otherwise will overwrite.  Note, if append is set, header row will not be appended."
)



args = parser.parse_args()







def get_rf_iqtree(basedir, treenum):
    
    filename = os.path.join(basedir, "rf", "iqtree_" + str(treenum) + ".rf")
    
    if not os.path.isfile(filename):
        return -1
    
    with open(filename, "r") as f:
        rf = float(f.read().strip())
        
    return rf



def get_rf_eccetera(basedir, treenum):
    
    #TODO: copypasta from previous function
    filename = os.path.join(basedir, "rf", "eccetera_" + str(treenum) + ".rf")
    
    if not os.path.isfile(filename):
        return -1
    
    with open(filename, "r") as f:
        rf = float(f.read().strip())
        
    return rf




def get_gene_tree_copy_numbers(filename):
    #DISCLAIMER: this is chatgpt
    """
    Read a Newick tree file with ete3 and analyze leaf names of the form X_Y_Z.
    DISCLAIMER: this is mostly chatgpt code.

    Parameters
    ----------
    filename : str
        Path to a Newick-formatted tree file.

    Returns
    -------
    max_count : int
        Maximum number of occurrences of any single X value.
    mean_count : float
        Average number of occurrences across all distinct X values.
    counts : dict
        Dictionary mapping each X value to its number of occurrences.

    Example
    -------
    If leaf names are:
        A_1_x, A_2_y, B_1_z, B_2_w, B_3_q

    Then:
        counts = {'A': 2, 'B': 3}
        max_count = 3
        mean_count = 2.5
    """
    # Load the tree
    tree = Tree(filename, format=1)

    # Extract X from each leaf name (X_Y_Z)
    x_values = []
    for leaf in tree.iter_leaves():
        parts = leaf.name.split("_")
        if len(parts) < 3:
            raise ValueError(
                f"Leaf name '{leaf.name}' does not match expected format X_Y_Z"
            )
        x_values.append(parts[0])

    # Count occurrences of each X
    counts = Counter(x_values)

    # Compute statistics
    max_count = max(counts.values()) if counts else 0
    mean_count = sum(counts.values()) / len(counts) if counts else 0.0

    return max_count, mean_count












#returns a dict where key = tree index, val = indelible_rate
def extract_indelible_rates(control_file: Path):
    #DISCLAIMER: this is mostly chatgpt code
    results = {}
    current_model = None

    with control_file.open() as f:
        for line in f:
            line = line.strip()

            # Detect model section
            if line.startswith("[MODEL]"):
                # e.g. "[MODEL] modelA_3"
                current_model = line.split()[1].split("_")[1]

            # Detect rates line inside a model
            elif line.startswith("[rates]") and current_model:
                parts = line.split()
                # eg parts = ['[rates]', '0', '2.139077', '0']
                second_number = float(parts[2])
                results[current_model] = second_number

    return results






directory = args.dir
if not os.path.isdir(directory):
    print(f"Error: {directory} is not a valid directory")
    sys.exit(1)


write_mode = "a" if args.append else "w"

outfile = open(args.outfile, write_mode)



header = "seed,gtnum,sites,duprate,lossrate,transferrate,pop,indelible_rate,avg_bs_iqtree,avg_bs_iqtree_bin,nb_leaves,avg_brlen,height,max_copy,mean_copy,rf_eccetera,rf_iqtree"
#print(header)

if not args.append:
    outfile.write(header + "\n")


for childdir in Path(directory).iterdir():
    #at this point childdir should be in a simphy dir, of the form ssim_dtl_...

    control_file = childdir / "1" / "control.txt"

    if not control_file.is_file():
        print(str(control_file) + " does not exist, skipping")
        continue        
    
    rates_per_sim = extract_indelible_rates(control_file)

    for i in range(1, 101):
   
        simphy_gtreefile = os.path.join(str(childdir), "1", f"g_trees{i:03d}.trees")
   
        #################################
        # SIMULATION PARAMS
        #################################
        sites, duprate, lossrate, transferrate, pop, seed = util.parse_sim_params_from_dir(os.path.basename(str(childdir)))
   
        #################################
        #INDELIBLE RATE
        #################################
        indelible_rate = rates_per_sim[str(i)]

        #################################
        # AVERAGE BOOTSTRAP
        #################################
        '''
        # This block obtains bootstrap values for pargenes, or falls back to raxml if not found
        avg_bs = -1
   
        bstree = os.path.join(str(childdir), "pargenes/trees/supports_run/results", f"dataset_{i:03d}_TRUE_phy.support.raxml.support")
   
        if not os.path.exists(bstree):
            bstree = os.path.join(str(childdir), "raxml", f"dataset_{i:03d}_TRUE.phy.raxml.support")
   
   
        if os.path.exists(bstree):
            supps = []
            util.count_support_bins(bstree, all_supports_list = supps)


            avg_bs = float(sum(supps)) / float(len(supps))
        '''
        
        #bootstrap from iqtree
        avg_bs = -1
        bs_iqtree_filename = os.path.join(str(childdir), "iqtree", "alignments", f"dataset_{i:03d}_TRUE.phy.treefile")
        
        if os.path.exists(bs_iqtree_filename):
            supps = []
            util.count_support_bins(bs_iqtree_filename, all_supports_list = supps)

            avg_bs = float(sum(supps)) / float(len(supps))
        else:
            print(f"Bootstrap tree {bs_iqtree_filename} does not exist")

        avg_bs_bin = -1
        if avg_bs > -1:
            avg_bs_bin = math.floor(avg_bs / 10) * 10

        #################################
        # NB LEAVES, AVG BR LEN, HEIGHT
        #################################
        nb_leaves = util.count_leaves(simphy_gtreefile)
   
        avg_brlen = util.get_average_branch_length(simphy_gtreefile)
   
        height = util.get_tree_height(simphy_gtreefile)
   
   
   
        #################################
        # COPY NUMBERS
        #################################
        max_copy, mean_copy = get_gene_tree_copy_numbers(simphy_gtreefile)
   
   
        #################################
        # RF VALS
        #################################
        rf_ecce = get_rf_eccetera(str(childdir), i)
        rf_iqtree = get_rf_iqtree(str(childdir), i)
   
   
   
   
        #################################
        # DONE, WRITE THE LINE
        #################################
        line = f"{seed},{i},{sites},{duprate},{lossrate},{transferrate},{pop},{indelible_rate},{avg_bs},{avg_bs_bin},{nb_leaves},{avg_brlen},{height},{max_copy},{mean_copy},{rf_ecce},{rf_iqtree}"
   
        #print(line)
        outfile.write(line + "\n")
   
    print("Done with " + str(childdir))
       
outfile.close()



















