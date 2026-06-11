import argparse
from pathlib import Path
import sys
from ete3 import Tree
import regex as re
import os
from collections import defaultdict
import util




filename_suffix = "_300sites"


datasets = { #"dtl2" : "ssim_dtlsim_s25_f100_sites100_GTR_bl1.0_d2.0_l2.0_t2.0_p0.0_pop10_mu1.0_theta0.0_seed", 
             #"dtl1" : "run_exp_out/ssim_dtl_s25_f100_sites100_GTR_bl1.0_d1_l1_t1.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed", 
             "dtl1" : "out_run_exp_300/ssim_dtl_s25_f100_sites300_GTR_bl1.0_d1_l1_t1.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed"
             #"dtl3" : "ssim_dtlsim_s25_f100_sites100_GTR_bl1.0_d3.0_l3.0_t3.0_p0.0_pop10_mu1.0_theta0.0_seed", 
             #"dtl0.5" : "run_exp_out/ssim_dtl_s25_f100_sites100_GTR_bl1.0_d0.5_l0.5_t1.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed3001", 
           }



support_histo = defaultdict(int)
support_vals = list()
support_avg_per_tree = list()



for ds in datasets:

    for seed in range(3001, 3026):
        dirname = f"{datasets[ds]}{seed}"

        
    
        for i in range(1, 100+1):
            
            bstree = os.path.join(dirname, "pargenes/trees/supports_run/results", f"dataset_{i:03d}_TRUE_phy.support.raxml.support")
            
            
            print(bstree)
            
            if os.path.exists(bstree):
            
                tree_sups = list()
                util.count_support_bins(bstree, support_histo, bin_size=10, all_supports_list = tree_sups)

                for s in tree_sups:
                    support_vals.append(s)
                if len(tree_sups) > 0:
                    support_avg_per_tree.append( float(sum(tree_sups)) / float(len(tree_sups))  )

                    
                    
                if len(tree_sups) > 0:
                    avg = ( float(sum(tree_sups)) / float(len(tree_sups))  )
                    print("Avg support = " + str(avg))
                    
                    #if avg > 90 and len(tree_sups) > 15:
                    #if avg < 10 and nbzero > 10:
                    #    print(tree_sups)
                    #    sys.exit()
                    
                
            

totalcount = 0
for bin in sorted(support_histo.keys()):
    print(f"{bin}  {support_histo[bin]}")
    totalcount += support_histo[bin]
    
    
for bin in sorted(support_histo.keys()):
    print(f"{bin}  {float(support_histo[bin])/float(totalcount)}")


if not support_vals is None:
    strout = "\n".join( str(sup) for sup in support_vals)
    util.write_to_file("all_support" + filename_suffix + ".csv", strout)
    
if not support_avg_per_tree is None:
    strout = "\n".join( str(sup) for sup in support_avg_per_tree)
    util.write_to_file("all_avg_support" + filename_suffix + ".csv", strout)



                
                
                
