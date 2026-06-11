import argparse
from pathlib import Path
import sys
from ete3 import Tree
import util

parser = argparse.ArgumentParser(description="Read rf.txt from each subdirectory in dir.")
parser.add_argument("dir", help="Path to parent directory")

parser.add_argument("method", help="Method, either 'speciesrax' (default), 'apro', or 'mininj', or 'speciesraxbm' or 'mininjbm', which uses the structure of Benoit's experiments", default='speciesrax')


args = parser.parse_args()
parent_dir = Path(args.dir)


unrooted_rf = True





def ete3_rf(treefile1, treefile2, unrooted = True, get_all_return = False):
    #print(treefile1)
    #print(treefile2)
    tree1 = Tree(treefile1, format = 1)
    tree2 = Tree(treefile2, format = 1)
    if (len(tree2.children) == 3):
      tree2.set_outgroup(tree2.children[0])
    if (len(tree1.children) == 3):
      tree1.set_outgroup(tree1.children[0])

    leaves1 = {leaf.name for leaf in tree1.iter_leaves()}
    leaves2 = {leaf.name for leaf in tree2.iter_leaves()}

    if leaves1 != leaves2:
        print("WARNING: trees have different leaves")


    res = tree1.robinson_foulds(tree2, unrooted_trees=unrooted, skip_large_polytomies = True, correct_by_polytomy_size = True)
    
    if get_all_return:
        return res
    else:
        return float(res[0]) / float(res[1])



def get_RF(treefile1, treefile2, normalize = True, unrooted = True):
    t1 = Tree(treefile1, format=1)
    t2 = Tree(treefile2, format=1)
    rf, max_rf, *_ = t1.robinson_foulds(t2, unrooted_trees=unrooted)

    if (normalize):
        return rf / max_rf
    return rf


if not parent_dir.is_dir():
    print("Error: Provided path is not a directory.")
    sys.exit(1)


rf_per_dl = dict()



for d in parent_dir.iterdir():
    if d.is_dir():
        
        #extract duprate from folder name (chatgpt's solution tbh)

        left, right = d.name.split("_l", 1)   # split only at first _l
        
        duprate = right.split("_", 1)[0]        # duprate ends at next underscore

        ''' that's crap, ignore
        if args.method == "speciesraxtest":
            simphy_file = d / "1" / "s_tree.trees.cleaned.newick"
            inferred_file = d / "speciesrax" / "species_trees" / "inferred_species_tree.newick.cleaned.newick"
            inferred_file_llr = d / "speciesrax" / "species_trees" / "species_tree_llr.newick"
            inferred_file_llr_cleaned = d / "speciesrax" / "species_trees" / "species_tree_llr.newick.cleaned.newick"
            
            inferred_file_supp = d / "speciesrax" / "species_trees" / "species_tree_root_support.newick"
            inferred_file_supp_cleaned = d / "speciesrax" / "species_trees" / "species_tree_root_support.newick.cleaned.newick"
            
            if inferred_file_llr.exists():
                util.clean_simphy_species_tree(str(inferred_file_llr), str(inferred_file_llr_cleaned))
                
            if inferred_file_supp.exists():
                util.clean_simphy_species_tree(str(inferred_file_supp), str(inferred_file_supp_cleaned))
            
            if inferred_file.exists():
                print(inferred_file)
                rf1 = ete3_rf(str(inferred_file), str(simphy_file), unrooted = unrooted_rf, get_all_return = True)
                
                rf2 = ete3_rf(str(inferred_file_supp_cleaned), str(simphy_file), unrooted = unrooted_rf, get_all_return = True)
                
                print(f"rfinf={rf1}")
                print(f"rfllr={rf2}")
        
            continue
        '''
        
        simphy_file = d / "1" / "s_tree.trees"
        #simphy_file = d / "1" / "s_tree.trees.cleaned.newick"

        rf_file = None

        if args.method == "speciesrax":
            
            inferred_file = d / "speciesrax" / "species_trees" / "inferred_species_tree.newick"
            #inferred_file = d / "speciesrax" / "species_trees" / "inferred_species_tree.newick.cleaned.newick"
        elif args.method == "apro":
            inferred_file = d / "apro_sptree.newick"
        elif args.method == "mininj":
            #inferred_file_base = d / "speciesrax" / "species_trees" / "starting_species_tree.newick"
            inferred_file = d / "speciesrax" / "species_trees" / "starting_species_tree.newick"
            
            #if not inferred_file.exists() and inferred_file_base.exists():
            #    util.clean_simphy_species_tree(str(inferred_file_base), str(inferred_file))

        elif args.method == "speciesraxbm" or args.method == "mininjbm":
            simphy_file = d / "1" / "s_tree.trees"
            inferred_file_name = "/home/manuel/git/phd_experiments/results/SpeciesRax/" + d.name
            inferred_file_name += "/MiniNJ_start_raxml-ng/run_--rec-model_UndatedDTL_--si-strategy_HYBRID_--per-family-rates_--skip-family-filtering_--run_speciesrax-dtl-raxml-ng-perfam-HYBRID.GTR_2/generax/species_trees/"
            
            if args.method == "speciesraxbm":
                inferred_file_name += "inferred_species_tree.newick"
            elif args.method == "mininjbm":
                inferred_file_name += "starting_species_tree.newick"

            
            if not Path(inferred_file_name).exists():
                inferred_file_name = inferred_file_name.replace("raxml-ng-perfam-HYBRID.GTR_2", "raxml-ng-perfam-HYBRID.GTR_1")
            
            if not Path(inferred_file_name).exists():
                inferred_file_name = inferred_file_name.replace("raxml-ng-perfam-HYBRID.GTR_1", "raxml-ng-perfam-HYBRID.GTR_0")
                
            if not Path(inferred_file_name).exists():
                print(inferred_file_name + " does not exist")
                sys.exit()
            
            inferred_file = Path(inferred_file_name)
            
            if not inferred_file.exists():
                print(f"WARNING: {inferred_file} does not exist")
                
            rf_file = Path(d / "rfinferred.txt")
            if unrooted_rf:
                rf_file = Path(d / "urfinferred.txt")
                
            if args.method == "mininjbm":
                rf_file = Path(str(rf_file).replace("rfinferred.txt", "rfstarting.txt"))
            
            if not rf_file.exists():
                print(f"WARNING: {rf_file} does not exist")

        if inferred_file.exists():
            
            
            rf = ete3_rf(str(inferred_file), str(simphy_file), unrooted = unrooted_rf)
            xrf = ete3_rf(str(inferred_file), str(simphy_file), unrooted = (not unrooted_rf))
            

            
            #test code, trying to find a directory where mininj is better than speciesrax
            #found: 
            #outsprax_DL/ssim_dtl_s25_f100_sites100_GTR_bl1.0_d1_l1_t1.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed3056
            #outsprax_DL/ssim_dtl_s25_f100_sites100_GTR_bl1.0_d1_l1_t1.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed3099
            '''
            if args.method == "mininj":  # and not str(d).endswith("3056"):
                sprax_spfile = d / "speciesrax" / "species_trees" / "inferred_species_tree.newick.cleaned.newick"
                rfsp = ete3_rf(str(sprax_spfile), str(simphy_file), unrooted = unrooted_rf)
                
                if rfsp > rf:
                    print(d)
                    print(f"njrf={rf} spraxrf={rfsp}")
                    
                    ret = ete3_rf(str(sprax_spfile), str(simphy_file), unrooted = unrooted_rf)
                    print(ret)
                    
                    sys.exit()
            '''

            
            if rf_file != None:
                with open(str(rf_file), 'r') as file:
                    file_rf = file.read()
                    #rf = file_rf
                    if float(rf) != float(file_rf):
                        print(d)
                        print(f"RF difference: mine = {rf}, theirs = {file_rf}")
                        sys.exit()
            
            
            if duprate not in rf_per_dl:
                rf_per_dl[duprate] = []
            rf_per_dl[duprate].append(float(rf))
                
                
print(f"Unrooted={unrooted_rf}")

for duprate, rfs in rf_per_dl.items():
    print(f"paramval={duprate}  nbsims={len(rf_per_dl[duprate])}")

for duprate, rfs in rf_per_dl.items():
    avg = sum(rfs) / len(rfs)
    print(f"paramval={duprate},{avg}")
                
                
                
                
                
                
                
                
