import os
import sys
from itertools import product


runs = 25

dlmult_list = [0, 0.5 ,1 ,2 ,5]

#transfer_mult = [0, 0.5,1 ,2 ,5]
transfermult_list = [1]  #for now just testing with baseline transfer rates

pop_list = [10, 1e7, 1e8, 1e9]


#sites = [100, 300, 500]
sites_list = [300]


#eccetera_thresholds = [50,70]
eccetera_thresholds = [50]



#for testing
runs = 5
dlmult_list = [1]
transfermult_list = [1]  #for now just testing with baseline transfer rates
pop_list = [10, 1e7]
sites_list = [300]
eccetera_thresholds = [50, 70]






#CHANGE PATH TO EXECUTABLES HERE

simphy_dir = "/home/manuel/git/SimPhy_1.0.2/"
simphy_bin = "/home/manuel/git/SimPhy_1.0.2/bin/simphy_lnx32"    #BMorel's script needs both simphy path and the bin

#iqtree version must be >= 3
iqtree_bin = "/home/manuel/git/iqtree-3.0.1-Linux/bin/iqtree3"


eccetera_bin = "/home/manuel/git/ecceTERA/bin/ecceTERA"

#INDELible must also be installed.  The executable must be in the system PATH, ie we need to be able to launch the command
#> indelible
#To change that, modify generate_families_with_simphy.py


#directory of Astral-Pro, and name of binary jar file
use_apro = True
apro_libpath = "/home/manuel/git/A-pro/ASTRAL-MP/lib"
apro_binname = "/home/manuel/git/A-pro/ASTRAL-MP/astral.1.1.6.jar"



use_nni = True
nni_k_list = "1 2 10"   #list of nb nnis to try



use_nni_clusters = True
nni_clusters_list = "1 2 3"






for dl, tr, pop, sites, ecce_threshold in product(dlmult_list, transfermult_list, pop_list, sites_list, eccetera_thresholds):

    #eg run: python run_exp.py --phylomethod=iqtree --skipexisting -r 25 --sites=300 --dlrate=2 --trate=1 --pop=10000000 --eccetera_bs=70 -o my_exp_p1e7_dl2_t1_s300
    
    pop = int(pop)  #otherwise it uses float and causes errors
    
    #outdir = f"allexp_p{pop}_dl{dl}_t{tr}_s{sites}"
    outdir = f"aprotest_p{pop}_dl{dl}_t{tr}_s{sites}"
    
    command = f"python run_exp.py --phylomethod=iqtree --skipexisting -r {runs} --sites={sites} --dlrate={dl} --trate={tr} "
    command += f" --eccetera_bs={ecce_threshold} --pop={pop} -o {outdir}"
    
    #set all binary paths
    command += f" --iqtreebin={iqtree_bin} --simphydir={simphy_dir} --simphybin={simphy_bin} --ecceterabin={eccetera_bin} "
    
    if use_apro:
        command += f" --apro_mode --apro_libpath={apro_libpath} --apro_binname={apro_binname}"

    if use_nni:        
        command += f" --nni_mode --nni_k_list {nni_k_list}"


    if use_nni_clusters:
        command += f" --nnicluster_mode --nniclusters_list {nni_clusters_list}"

    
    
    print("Running: " + command)
    
    os.system(command)




