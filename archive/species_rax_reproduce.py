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

import util

parser = argparse.ArgumentParser(description="Script that calls generate_families_with_simphy with the given arguments.")


parser.add_argument(
    "-o", "--outrootdir",
    help="Output directory for simphy files", 
    default="simphy_out"
)

parser.add_argument(
    "-r", "--rep",
    type=int,
    help="Number of replicates", 
    default=50
)



parser.add_argument(
    "--vary",
    type=str,
    help="Which parameter to vary.  Must be DL or DTL", 
    default="DL"
)



parser.add_argument(
    "--method",
    type=str,
    help="Which method to use.  Must be speciesrax or apro", 
    default="speciesrax"
)





parser.add_argument(
    "-g", "--genes",
    help="Number of gene families", 
    default=100
)


#parser.add_argument(
#    "--seed",
#    help="Seed given to simphy", 
#    default=random.randint(0, sys.maxsize)
#)


parser.add_argument(
    "--skipsimphy",
    action="store_true",
    help="Set to skip simphy simulations, and indelible if exists"
)



parser.add_argument(
    "--phylomethod",
    type=str,
    default="pargenes",
    help="Which method will infer gene trees.  Either raxml or pargenes"
)


parser.add_argument(
    "--skipphylo",
    action="store_true",
    help="Set to skip phylogeny reconstruction (won't call raxml nor pargenes)"
)




parser.add_argument(
    "--speciesraxbin",
    help="Full path to the SpeciesRax binary", 
    default="generax"
)


parser.add_argument(
    "--raxmlbin",
    help="Full path to the RAxML binary", 
    default="/usr/bin/raxml-ng"
)



parser.add_argument(
    "--workdir",
    help="Directory for temp files", 
    default="work"
)



use_strategy_spr = False  #todo make arg


args = parser.parse_args()



#TODO : do not hardcode this
generax_bin = "~/git/GeneRax/build/bin/generax"


#long string that will contain the csv file output, todo output after each iteration instead
stats_str = ""

dtl_rates = [0.5, 1, 2, 3, 5]
dl_rates = [1, 0, 0.5, 5, 2]

reps = range(1, args.rep + 1) 



if args.vary == "DL":
    parameter = dl_rates
elif args.vary == "DTL":
    parameter = dtl_rates
else:
    print("Invalid vary parameter")
    sys.exit()


#for (dtl_rate, rep) in itertools.product(dtl_rates, reps):

for (pval, rep) in itertools.product(parameter, reps):
    
    
    params = GFWS.SimphyParameters()
    #params.seed = args.seed   
    
    #TODO: what if we encounter a previous seed?  Unlikely but possible
    #params.seed = random.randint(0, 2**32)
    params.seed = 3000 + rep
    params.families_number = args.genes
    

    if args.vary == "DL":
        params.dup_rate = pval
        params.loss_rate = pval   #must be equal to dup_rate else simphy fails
        params.transfer_rate = 1.0  #1.0 ??
    elif args.vary == "DTL":
        params.dup_rate = pval
        params.loss_rate = pval   #must be equal to dup_rate else simphy fails
        params.transfer_rate = pval
    
    
    output_dir = GFWS.get_output_dir(params, args.outrootdir)
    simphy_species_tree_file = os.path.join(output_dir, "1", "s_tree.trees")

    #simulates both sequences and indelible
    if not args.skipsimphy:
        try:
            GFWS.generate_from_parameters(params, args.outrootdir)
        except Exception as e:
            print(e)
            continue
    else:
        #if we skip, check that dir exists, otherwise continue  
        if util.is_bad_simphy_dir(output_dir):
            continue
            
    
        #or maybe simphy did its job previously, but indelible failed, and so we try it again here        
        if not os.path.exists( os.path.join(output_dir, "1", "dataset_001_TRUE.phy"  ) ):
            try:
                GFWS.generate_from_parameters(params, args.outrootdir, skip_simphy = True)
            except Exception as e:
                print(e)
                continue

    



    

    gene_tree_files = []

    alignment_files = []
    for i in range(1, params.families_number + 1):
        alignment_file = os.path.join(output_dir, "1", f"dataset_{i:03d}_TRUE.phy")
        alignment_files.append(alignment_file)


    if args.phylomethod == "pargenes":
        #INFER GENE TREES WITH PARGENES
        for i in range(1, params.families_number + 1):
            treefile = os.path.join(output_dir, "pargenes/trees/mlsearch_run/results", f"dataset_{i:03d}_TRUE_phy", f"dataset_{i:03d}_TRUE_phy.raxml.bestTree")
            gene_tree_files.append( treefile )

        
        if not args.skipphylo or not os.path.exists(gene_tree_files[0]):
            #util.run_pargenes_on_all(os.path.join(output_dir, "1"), os.path.join(output_dir, "pargenes"), args.raxmlbin)
            
            #raxml_args = " --model GTR+G4 --redo --extra bs-start-rand "
            util.run_pargenes(alignment_files, os.path.join(output_dir, "pargenes"), args.raxmlbin)
    
        
    elif args.phylomethod == "raxml":
        #INFER GENE TREES WITH RAxML
        for i in range(1, params.families_number + 1):
            seq_file = os.path.join(output_dir, "1", f"dataset_{i:03d}_TRUE.phy")    #:03d means length 3 padded with 0s
            raxml_treefile = util.run_raxml(seq_file, args.raxmlbin, returned_collapsed_treefile = False)   
            gene_tree_files.append(raxml_treefile)
    else:
        print("Skipping gene tree reconstruction")

    


    if args.method == "speciesrax":
    
        family_filename = "speciesrax_family.txt"
        if use_strategy_spr:
            family_filename = "speciesrax_family_withalign.txt"
    
        family_file = os.path.join(output_dir, family_filename)
    
        util.make_generax_family_file_from_trees(family_file, gene_tree_files, alignment_files)

        generax_outdir = os.path.join(output_dir, "speciesrax")
    
        #calls GeneRax in the SpeciesRax setting, in particular initial SpeciesTree is the MiniNJ one
        #command = f"{generax_bin} --strategy SKIP --prune-species-tree --si-strategy HYBRID  -s MiniNJ -f {family_file} -p {generax_outdir}"
        #command += " --si-spr-radius 3 "
        #print(f"Running {command}")
        
        util.run_speciesrax(family_file, generax_outdir, generax_bin, use_strategy_spr = use_strategy_spr)
        

        generax_sptreefile = os.path.join(generax_outdir, "species_trees", "inferred_species_tree.newick")
        simphy_sptree_file_cleaned = simphy_species_tree_file + ".cleaned.newick"
        
        print(f"Cleaning {simphy_species_tree_file}")
        util.clean_simphy_species_tree(simphy_species_tree_file, simphy_sptree_file_cleaned)
    
        generax_sptreefile_cleaned = generax_sptreefile + ".cleaned.newick"
        util.clean_simphy_species_tree(generax_sptreefile, generax_sptreefile_cleaned)

        rf = util.get_RF(simphy_sptree_file_cleaned, generax_sptreefile_cleaned)
    
    
    elif args.method == "apro":
        
        apro_sptree_file = os.path.join(output_dir, "apro_sptree.newick")
        util.run_apro_from_symphy(gene_tree_files, os.path.join(output_dir, "apro_trees.txt"), apro_sptree_file)
        
        simphy_sptree_file_cleaned = simphy_species_tree_file + ".cleaned.newick"
        util.clean_simphy_species_tree(simphy_species_tree_file, simphy_sptree_file_cleaned)
        
        apro_sptreefile_cleaned = apro_sptree_file + ".cleaned.newick"
        util.clean_simphy_species_tree(apro_sptree_file, apro_sptreefile_cleaned)
        rf = util.get_RF(simphy_sptree_file_cleaned, apro_sptreefile_cleaned)

    
    
    with open(os.path.join(output_dir, f"{args.method}_rf.txt"), "w") as f:
        f.write(str(rf))

    print(f"RF={rf}")
    
    stats_str += f"{args.vary}={pval},rep={rep},{rf}\n"


with open(f"stats_param_{args.vary}_{args.method}.csv", 'w') as f:
    f.write(stats_str)



