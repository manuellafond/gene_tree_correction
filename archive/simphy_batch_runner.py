import os
import sys
from ete3 import Tree
from Bio import Phylo
import generate_families_with_simphy as GFWS



def clean_simphy_species_tree(input_file, output_file):
    # Read the tree
    tree = Phylo.read(input_file, "newick")

    # Remove single quotes from the entire file content
    with open(input_file, "r") as f:
        newick_str = f.read().replace("'", "")

    # Remove internal node names
    for clade in tree.find_clades():
        if not clade.is_terminal():  # If it's an internal node
            clade.name = None  # Remove the name

    # Write back the cleaned tree
    Phylo.write(tree, output_file, "newick")

    # Ensure the single quotes are removed by re-processing the output file
    # ML doesn't know why this is there...
    with open(output_file, "r") as f:
        cleaned_str = f.read().replace("'", "")

    with open(output_file, "w") as f:
        f.write(cleaned_str)





out_rootdir = "sims"
params = GFWS.SimphyParameters()

params.seed = 42   #don't forget to generate random seeds
params.families_number = 100


GFWS.generate_from_parameters(params, out_rootdir)

output_dir = GFWS.get_output_dir(params, out_rootdir)


simphy_species_tree_file = os.path.join(output_dir, "1", "s_tree.trees")

for i in range(1, params.families_number + 1):
    #run iqtree
    seq_file = os.path.join(output_dir, "1", f"dataset_{i:03d}_TRUE.phy")	#:03d means length 3 padded with 0s
    command = f"iqtree -s {seq_file} -B 1000 -nt AUTO -redo"

    print(f"Running {command}")

    os.system(command)

    simphy_tree_file = os.path.join(output_dir, "1", f"g_trees{i:03d}.trees")
    iqtree_tree_file= f"{seq_file}.treefile"  
    
    t1 = Tree(simphy_tree_file)
    t2 = Tree(iqtree_tree_file)

    # Compute RF distance (unrooted)
    rf, max_rf, *_ = t1.robinson_foulds(
        t2,
        unrooted_trees=True
    )

    print("RF distance:", rf)
    print("Max RF:", max_rf)
    print("Normalized RF:", rf / max_rf)
    
    
    
    
    #run ecceTERA    
    python clean_species_tree.py "${dir}/sim_${i}/s_tree.newick" "${dir}/sim_${i}/s_tree_cleaned.newick"
    clean_simphy_species_tree(

    threshold = 75
    dup_cost = 2
    command = f"ecceTERA_linux64 species.file={simphy_species_tree_file} \
        gene.file={iqtree_tree_file} \
        dated=0 compute.T=false collapse.threshold={threshold} collapse.mode=1 \
        dupli.cost={dup_cost} loss.cost=1 print.newick=true resolve.trees=1"
        
    print(f"Running {command}"

    
    sys.exit()




