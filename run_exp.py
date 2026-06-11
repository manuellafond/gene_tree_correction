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
    "--simphybin",
    help="Full path to the directory containing the Simphy bin dir", 
    default="/home/manuel/SimPhy/"
)


parser.add_argument(
    "--pargenesbin",
    help="Full path to pargenes", 
    default="/home/manuel/git/ParGenes/pargenes/pargenes.py"
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



args = parser.parse_args()



#TODO : do not hardcode this
generax_bin = args.generaxbin


reps = range(1, args.rep + 1) 




for rep in range(1, args.rep+1):
    
    
    GFWS.simphy_path = args.simphybin
    
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
    
    
    simphy_species_tree_file = os.path.join(output_dir, "1", "s_tree.trees")
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
            util.run_iqtree_on_all(alignment_files, iqtree_dir, skip_existing = args.skipexisting)
            
        
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
    ##################################################################################
    
    eccetera_corrected_trees = [""] * args.genes
    
    bs_threshold = args.eccetera_bs
    for i in range(1, args.genes + 1):
        gtreefile = gene_tree_files_with_bs_nozero[i-1]
        
        eccetera_gfilename = gtreefile + ".eccetera"
        eccetera_corrected_trees[i-1] = eccetera_gfilename

        #example
        #/home/manuel/git/ecceTERA/bin/ecceTERA species.file=run_exp_out/ssim_dtl_s25_f100_sites100_GTR_bl1.0_d1_l1_t1.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed3001/1/s_tree.trees gene.file=run_exp_out/ssim_dtl_s25_f100_sites100_GTR_bl1.0_d1_l1_t1.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed3001/pargenes/trees/supports_run/results/dataset_002_TRUE_phy.support.raxml.support dated=0 compute.T=false collapse.threshold=70 collapse.mode=1 resolve.trees=1 verbose=true print.newick=true print.newick.gene.tree.file=run_exp_out/ssim_dtl_s25_f100_sites100_GTR_bl1.0_d1_l1_t1.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed3001/pargenes/trees/supports_run/results/dataset_002_TRUE_phy.support.raxml.support.eccetera
        
        if os.path.exists(gtreefile):
            command = f"{args.ecceterabin} species.file={simphy_species_tree_file} gene.file={gtreefile} dated=0 compute.T=false "
            command += f"collapse.threshold={bs_threshold} collapse.mode=1 resolve.trees=1 verbose=true amalgamate=false "
            command += f"print.newick=true print.newick.gene.tree.file={eccetera_gfilename} degree.limit=12 "

            print(f"Executing\n{command}")
            
            skip = args.skipeccetera
            if args.skipexisting and os.path.exists(eccetera_gfilename):
                skip = True
            
            if not skip:     
                os.system(command)

        else:
            print(f"{gtreefile} does not exist, skipping eccetera")            



    ##################################################################################    
    # step 4: compute *unrooted* RF values
    ##################################################################################
    rfdir = os.path.join(output_dir, "rf")
    util.make_dir( rfdir )
    for i in range(1, args.genes + 1):
        gfile = gene_tree_files_with_bs[i-1]
        eccetera_gfile = eccetera_corrected_trees[i-1]
        simphy_gfile = simphy_gene_tree_files[i-1]
        
        if os.path.exists(gfile) and os.path.exists(eccetera_gfile):
            print(f"comparing\n{gfile}\n{simphy_gfile}")
            urf_phylomethod = util.ete3_rf(gfile, simphy_gfile)
            
            print(f"comparing\n{eccetera_gfile}\n{simphy_gfile}")
            urf_eccetera = util.ete3_rf(eccetera_gfile, simphy_gfile)
            
            
            util.write_to_file( os.path.join(rfdir, f"{args.phylomethod}_{i}.rf"), str(urf_phylomethod) )
            util.write_to_file( os.path.join(rfdir, f"eccetera_{i}.rf"), str(urf_eccetera) )
            
            #outstr = f"urf_phylomethod={urf_phylomethod}\nurf_eccetera={urf_eccetera}"
            #print(outstr)
            
            

