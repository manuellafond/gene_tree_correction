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
import pargenes_runner
import util
import raxml_runner
from collections import defaultdict

parser = argparse.ArgumentParser(description="Script that generates a specified number of simphy+indelible runs.   For each of them, gene trees are reconstructed from the alignments and corrected with eccetera.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)


#example using iqtree for 25 simphyruns, 300 sites, x2 dl rate, x1 transfer rate, pop = 1e7
#python run_exp.py --phylomethod=iqtree --skipexisting -r 25 --sites=300 --dlrate=2 --trate=1 --pop=10000000 -o my_exp_p1e7_dl2_t1_s300


parser.add_argument(
    "-o", "--outrootdir",
    help="Output directory for simphy+alignment+gene trees+correction files.  If directory exists, does not destroy existing files.", 
    default="run_exp_out"
)

parser.add_argument(
    "-r", "--rep",
    type=int,
    help="Number of simphy replicates", 
    default=50
)



parser.add_argument(
    "--sites",
    type=int,
    help="Number of average sites", 
    default=300
)



parser.add_argument(
    "-p", "--phylomethod",
    help="Method to infer trees, one of [raxml, pargenes, iqtree]", 
    default="raxml"
)




parser.add_argument(
    "--pop",
    type=int,
    help="Population for simphy", 
    default=10
)




parser.add_argument(
    "-g", "--genes",
    help="Number of gene families in simphy.  Must contain three digits, i.e., between 100 and 999.", 
    default=100
)





parser.add_argument(
    "--skipeccetera",
    action="store_true",
    help="Set to skip eccetera"
)



parser.add_argument(
    "--eccetera_bs",
    type=float,
    default=70,
    help="Threshold for eccetera to collapse branches"
)


parser.add_argument(
    "--skipexisting",
    action="store_true",
    help="Set to skip any command if the relevant file exists (not set by default)"
)


parser.add_argument(
    "--startseed",
    help="Starting seed.  First simphy run will use that, then next one that seed +1, and so on.", 
    type=int,
    default=3001
)




parser.add_argument(
    "--raxmlbin",
    help="Full path to the RAxML-ng binary", 
    default="/usr/bin/raxml-ng"
)




parser.add_argument(
    "--simphydir",
    help="Full path to the directory containing the Simphy bin dir", 
    default="/home/manuel/SimPhy/"
)

parser.add_argument(
    "--simphybin",
    help="Full path to the directory containing the Simphy bin dir", 
    default="/home/manuel/SimPhy/bin/simphy"
)


parser.add_argument(
    "--pargenesbin",
    help="Full path to pargenes", 
    default="/home/manuel/git/ParGenes/pargenes/pargenes.py"
)


parser.add_argument(
    "--iqtreebin",
    help="Full path to iqtree", 
    default="iqtree"
)



parser.add_argument(
    "--ecceterabin",
    help="Full path to eccetera", 
    default="/home/manuel/git/ecceTERA/bin/ecceTERA"
)




parser.add_argument(
    "--generaxbin",
    help="Full path to generax", 
    default="~/git/GeneRax/build/bin/generax"
)


parser.add_argument(
    "--workdir",
    help="Directory for temp files", 
    default="work"
)



parser.add_argument(
    "--dlrate",
    help="Dup loss rate multiplier",
    type=float, 
    default=1
)

parser.add_argument(
    "--trate",
    help="Transfer rate multiplier",
    type=float, 
    default=1
)





######################################
# ASTRAL-PRO arguments
######################################
parser.add_argument(
    "--apro_mode",
    action="store_true",
    help="Set to run astral pro.  Species tree will be in simphy dir, named s_tree.trees.apro."  
         "eccetera will use that tree, *in addition to* the simphy one."
)


parser.add_argument(
    "--apro_libpath",
    help="Directory that contains the Astral-Pro binaries", 
    default="/home/manuel/git/A-pro/ASTRAL-MP/lib"
)


parser.add_argument(
    "--apro_binname",
    help="Name of the binary file of Astral-Pro (the java jar file)", 
    default="/home/manuel/git/A-pro/ASTRAL-MP/astral.1.1.6.jar"
)





######################################
# Random NNIs arguments
######################################
parser.add_argument(
    "--nni_mode",
    action="store_true",
    help="Set to run eccetera but after applying --nni_k NNIs on the species tree. "
         "Species tree will be in simphy dir, named s_tree.trees.nni_{k}.  "
         "eccetera will use that tree, *in addition to* the simphy one."
)


#list of the number of random nnis to try
#e.g., --nni_k_list '1 2 3'
parser.add_argument(
  "--nni_k_list",  
  nargs="*",  
  type=int,
  default=[3]
)




######################################
# Clustered NNIs arguments
######################################
parser.add_argument(
    "--nnicluster_mode",
    action="store_true",
    help="Set to run eccetera but after applying 3 NNIs on the species tree. "
         "These NNIs are 'clustered', meaning that they are applied consecutively around the same branch."
         "(i.e., we choose a branch and the sibling branch (which changes) is relocated thrice."
         "The idea is ti simulate a series of problems occurring in a small portion of the tree."
)


#list of the number of random nnis to try
#e.g., --nniclusters_list '1 2 3'
parser.add_argument(
  "--nniclusters_list",  
  nargs="*",  
  type=int,
  default=[3]
)







######################################
# End of arguments
######################################





args = parser.parse_args()


generax_bin = args.generaxbin

reps = range(1, args.rep + 1) 


############################################
# Helper functions
#returns suffix of eccetera files depending on parameters
def get_eccetera_suffix(bs_threshold, sptree_method):
    suffix = f"_{bs_threshold}"

    if sptree_method != "simphy":
        suffix += f"_{sptree_method}"
    
    return suffix



#for nni, expected format is nnirand_{k}
def get_sptree_filename_by_method(sptree_method, output_dir, extra_args = []):
    if sptree_method == "simphy":
        return os.path.join(output_dir, "1", "s_tree.trees")
    elif sptree_method == "apro":
        return os.path.join(output_dir, "1", "s_tree.trees.apro")
    elif sptree_method.startswith("nnirand_"):
        k = sptree_method.replace("nnirand_", "")
        return os.path.join(output_dir, "1", f"s_tree.trees.nni_{k}")
    elif sptree_method.startswith("nnicluster_"):
        nbclusters = sptree_method.replace("nnicluster_", "")
        return os.path.join(output_dir, "1", f"s_tree.trees.nnicluster_{nbclusters}")

############################################



for rep in range(1, args.rep+1):
    
    
    GFWS.simphy_path = args.simphydir
    GFWS.simphy_bin = args.simphybin
    
    params = GFWS.SimphyParameters()   #GFWS = BMorel's simphy script

    ##################################################################################    
    # step 1: run simphy
    ##################################################################################
    params.seed = (args.startseed - 1) + rep
    params.families_number = args.genes
    params.sites = args.sites    
    params.population = args.pop

    params.dup_rate = args.dlrate
    params.loss_rate = args.dlrate   #must be equal to dup_rate else simphy fails
    params.transfer_rate = args.trate

    
    output_dir = GFWS.get_output_dir(params, args.outrootdir)
    
    
    simphy_species_tree_file = get_sptree_filename_by_method("simphy", output_dir)
    simphy_gene_tree_files = [""] * args.genes
    for i in range(1, args.genes + 1):
        simphy_gene_tree_files[i-1] = os.path.join(output_dir, "1", f"g_trees{i:03d}.trees")
    
    #our check for "simphy was already run" = the species tree exists
    skip = False
    if args.skipexisting and os.path.exists(simphy_species_tree_file):
        skip = True
    
    if not skip:
        GFWS.generate_from_parameters(params, args.outrootdir)



    ##################################################################################
    # step 1.1: run Astral-Pro to infer a species tree, if needed
    ##################################################################################
    if args.apro_mode:
        apro_workfile = os.path.join(output_dir, "1", "apro.trees")
        apro_outfile = get_sptree_filename_by_method("apro", output_dir)
        util.run_apro_from_symphy(simphy_gene_tree_files, apro_workfile, apro_outfile, apro_lib_path = args.apro_libpath, apro_bin_path = args.apro_binname)
        
        
        
    ##################################################################################
    # step 1.2: perform NNIs on the species tree if command line asked for it
    ##################################################################################
    if args.nni_mode:
        
        for k in args.nni_k_list:
            nni_outfile = get_sptree_filename_by_method(f"nnirand_{k}", output_dir)
            print(nni_outfile)

            util.random_nni(simphy_species_tree_file, k, nni_outfile)


    ##################################################################################
    # step 1.3: perform clustered NNIs on the species tree if command line asked for it
    ##################################################################################
    if args.nnicluster_mode:
        
        for nbclusters in args.nniclusters_list:
            nnicluster_outfile = get_sptree_filename_by_method(f"nnicluster_{nbclusters}", output_dir)
            print(nnicluster_outfile)

            
            util.cluster_nni(simphy_species_tree_file, nbclusters, nnicluster_outfile)


    
    ##################################################################################    
    # step 2: run gene tree inference method (pargenes/raxml/iqtree)
    ##################################################################################
    if args.phylomethod == "pargenes":
        
        alignment_files = []
        for i in range(1, args.genes + 1):
            alignment_file = os.path.join(output_dir, "1", f"dataset_{i:03d}_TRUE.phy")
            alignment_files.append(alignment_file)

        

        pgparams = pargenes_runner.ParGenesParams()    
    
        pgparams.alignment_files = alignment_files
        pgparams.output_dir = os.path.join(output_dir, "pargenes")
        pgparams.fast_bootstrap = True
        pgparams.bootstrap_replicates = 100
        #pgparams.start_from_random_trees = 1
    
        gene_tree_files = [""] * args.genes    
        for i in range(1, args.genes + 1):
            treefile = os.path.join(output_dir, "pargenes/trees/mlsearch_run/results", f"dataset_{i:03d}_TRUE_phy", f"dataset_{i:03d}_TRUE_phy.raxml.bestTree")
            gene_tree_files[i-1] = treefile
    
        skip = False
        if os.path.exists(gene_tree_files[0]) and args.skipexisting:
            skip = True
    
        if not skip:
            pargenes_runner.run_pargenes(pgparams)
            
            
        

        gene_tree_files_with_bs = [""] * args.genes      #bs means bootstrap here (not something else)
        gene_tree_files_with_bs_nozero = [""] * args.genes  #some software (eccetera I think) refuse tree with bootstrap 0, so we replace 0 with 1
        for i in range(1, args.genes + 1):
        
            treefile_with_bs = os.path.join(output_dir, "pargenes/trees/supports_run/results", f"dataset_{i:03d}_TRUE_phy.support.raxml.support")
            gene_tree_files_with_bs[i-1] = treefile_with_bs
            
            nozero_file = treefile_with_bs + ".nozero"
            if os.path.exists(treefile_with_bs): 
                util.fix_zero_bootstrap(treefile_with_bs, nozero_file)
            
            gene_tree_files_with_bs_nozero[i-1] = nozero_file
            
            
        '''
        #just a sanity check - the support trees should be identical to best trees
        for i in range(1, args.genes + 1):
            treefile = raxml_gene_tree_files[i-1]
            treefile_with_bs = raxml_gene_tree_files_with_bs[i-1]
            if os.path.exists(treefile):
                print(f"Comparing {treefile} and {treefile_with_bs}")
                rf = util.ete3_rf(treefile, treefile_with_bs)
                print(f"RF={rf}")
                if rf > 0:
                    print(f"{treefile} and {treefile_with_bs} are different")
                    sys.exit()
        '''        
                
    ##################################################################################    
    # step 2.ALT: run raxml
    ##################################################################################
    if args.phylomethod == "raxml":

        alignment_files = []
        for i in range(1, args.genes + 1):
            alignment_file = os.path.join(output_dir, "1", f"dataset_{i:03d}_TRUE.phy")
            alignment_files.append(alignment_file)


        raxmlparams = raxml_runner.raxmlParams()    
    
        raxmlparams.alignment_files = alignment_files
        raxmlparams.output_dir = os.path.join(output_dir, "raxml")
        raxmlparams.fast_bootstrap = True
        raxmlparams.bootstrap_replicates = 100
        raxmlparams.start_from_random_trees = 10

        
        #TODO: lots of copy-pasting from pargenes step here
        gene_tree_files = [""] * args.genes
        gene_tree_files_with_bs = [""] * args.genes      #bs means bootstrap here, not something else
        gene_tree_files_with_bs_nozero = [""] * args.genes        
        
        for i in range(1, args.genes + 1):
            treefile = os.path.join(output_dir, "raxml", f"dataset_{i:03d}_TRUE.phy.raxml.bestTree")
            gene_tree_files[i-1] = treefile
        
            treefile_with_bs = os.path.join(output_dir, "raxml", f"dataset_{i:03d}_TRUE.phy.raxml.support")
            gene_tree_files_with_bs[i-1] = treefile_with_bs
        
        
    
        skip = False
        if os.path.exists(gene_tree_files[0]) and args.skipexisting:
            skip = True
    
        if not skip:
            raxml_runner.run_raxml(raxmlparams)
            
        

                    
        for i in range(1, args.genes + 1):
            
            
            nozero_file = gene_tree_files_with_bs[i-1] + ".nozero"
            if os.path.exists(treefile_with_bs): 
                util.fix_zero_bootstrap(treefile_with_bs, nozero_file)
            
            gene_tree_files_with_bs_nozero[i-1] = nozero_file
            
    
    
    
    ##################################################################################    
    # step 2.ALT: run iqtree
    ##################################################################################
    if args.phylomethod == "iqtree":
    
        def get_iqtree_treefile(i):
           return os.path.join(output_dir, "iqtree", "alignments", f"dataset_{i:03d}_TRUE.phy.treefile")
    

        alignment_files = []
        for i in range(1, args.genes + 1):
            alignment_file = os.path.join(output_dir, "1", f"dataset_{i:03d}_TRUE.phy")
            
            treefile = get_iqtree_treefile(i)
            
            #don't consider alignment file if treefile exists already
            if os.path.exists(treefile) and args.skipexisting:
                continue
             
            
            if os.path.exists(alignment_file):
                nb_seqs = util.first_word_from_file(alignment_file)
            else:
                nb_seqs = None
            
            #do not bother with trees having 3 or less leaves, iqtree won't like them
            if not nb_seqs is None and int(nb_seqs) >= 4:
                alignment_files.append(alignment_file)

        iqtree_dir = os.path.join(output_dir, "iqtree")
    
    
        if len(alignment_files) > 0:
            util.run_iqtree_on_all(alignment_files, iqtree_dir, skip_existing = args.skipexisting, iqtree_bin = args.iqtreebin)
            
        
        gene_tree_files = [""] * args.genes
        gene_tree_files_with_bs = [""] * args.genes      #bs means bootstrap here, not something else
        gene_tree_files_with_bs_nozero = [""] * args.genes
                    
        for i in range(1, args.genes + 1):
            treefile = get_iqtree_treefile(i)
            gene_tree_files[i-1] = treefile
        
            treefile_with_bs = treefile  #in iqtree, inferred tree has support on it directly
            gene_tree_files_with_bs[i-1] = treefile_with_bs
            
            nozero_file = treefile_with_bs + ".nozero"
            if os.path.exists(treefile_with_bs): # and not os.path.exists(nozero_file):
                util.fix_zero_bootstrap(treefile_with_bs, nozero_file)
            else:
                pass
                #print(treefile_with_bs + " does not exist")
                #sys.exit()
            
            gene_tree_files_with_bs_nozero[i-1] = nozero_file
            
    
    ##################################################################################
    # step 3: run eccetera on the bootstrap trees
    # NOTE: if args.apro_mode is set, will *also* use stree.trees.apro, not simphy's tree
    ##################################################################################

    
    eccetera_corrected_trees = {}

    error_log = "eccetera_errors.log"
    
    possible_bootstraps = [50, 70]
    
    possible_sptree_methods = ["simphy"]
    if args.apro_mode:
        possible_sptree_methods.append("apro")
        
    if args.nni_mode:
        for k in args.nni_k_list:
            possible_sptree_methods.append(f"nnirand_{k}")
            
    if args.nnicluster_mode:        
        for nbclusters in args.nniclusters_list:            
            possible_sptree_methods.append(f"nnicluster_{nbclusters}")
            

    

    for bs_threshold in possible_bootstraps:
        for sptree_method in possible_sptree_methods:

            eccetera_suffix = get_eccetera_suffix(bs_threshold, sptree_method)

            eccetera_corrected_trees[(bs_threshold, sptree_method)] = [""] * args.genes

            for i in range(1, args.genes + 1):
                gtreefile = gene_tree_files_with_bs_nozero[i-1]

                eccetera_gfilename = gtreefile + ".eccetera" + eccetera_suffix
                eccetera_corrected_trees[(bs_threshold, sptree_method)][i-1] = eccetera_gfilename

                sptree_for_eccetera = get_sptree_filename_by_method(sptree_method, output_dir) 

                if os.path.exists(gtreefile):
                    command = f"{args.ecceterabin} species.file={sptree_for_eccetera} gene.file={gtreefile} dated=0 compute.T=false "
                    command += f"collapse.threshold={bs_threshold} collapse.mode=1 resolve.trees=1 verbose=true amalgamate=false "    
                    command += f"print.newick=true print.newick.gene.tree.file={eccetera_gfilename} degree.limit=12 "

                    print(f"Executing\n{command}")

                    skip = args.skipeccetera
                    if args.skipexisting and os.path.exists(eccetera_gfilename):
                        skip = True

                    if not skip:
                        os.system(command)

                        if not os.path.exists(eccetera_gfilename):
                            with open(error_log, "a") as err:
                                err.write(f"ERROR: ecceTERA did not create {eccetera_gfilename} "        
                                          f"(gene tree: {gtreefile}, threshold={bs_threshold})\n")

                else:
                    print(f"{gtreefile} does not exist, skipping ecceTERA")
                    with open(error_log, "a") as err:
                        err.write(f"ERROR: {gtreefile} does not exist "    
                                  f"(threshold={bs_threshold})\n")
                
            
    
    
    ##################################################################################
    # step 4: compute *unrooted* RF values
    # NOTE: for eccetera RF, uses eccetera_suffix from above
    ##################################################################################
    rfdir = os.path.join(output_dir, "rf")
    util.make_dir(rfdir)

    for i in range(1, args.genes + 1):
        gfile = gene_tree_files_with_bs[i-1]
        simphy_gfile = simphy_gene_tree_files[i-1]

        if os.path.exists(gfile):
            print(f"comparing\n{gfile}\n{simphy_gfile}")
            urf_phylomethod = util.ete3_rf(gfile, simphy_gfile)

            util.write_to_file(
                os.path.join(rfdir, f"{args.phylomethod}_{i}.rf"),
                str(urf_phylomethod)
            )

            for bs_threshold in possible_bootstraps:
                for sptree_method in possible_sptree_methods:

                    eccetera_suffix = get_eccetera_suffix(bs_threshold, sptree_method)
            
                    eccetera_gfile = eccetera_corrected_trees[(bs_threshold, sptree_method)][i-1]

                    if os.path.exists(eccetera_gfile):
                        print(f"comparing\n{eccetera_gfile}\n{simphy_gfile}")
                        urf_eccetera = util.ete3_rf(eccetera_gfile, simphy_gfile)

                        util.write_to_file(
                            os.path.join(rfdir, f"eccetera_{i}{eccetera_suffix}.rf"),
                            str(urf_eccetera)
                        )

                        #outstr = f"urf_phylomethod={urf_phylomethod}\nurf_eccetera={urf_eccetera}"
                        #print(outstr) 
            




