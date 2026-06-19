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


#CHANGE PATH TO EXECUTABLES HERE

simphy_dir = "/home/manuel/git/SimPhy_1.0.2/"
simphy_bin = "/home/manuel/git/SimPhy_1.0.2/bin/simphy_lnx32"    #BMorel's script needs both simphy path and the bin

#iqtree version must be >= 3
iqtree_bin = "/home/manuel/git/iqtree-3.1.3-Linux/bin/iqtree3"


eccetera_bin = "/home/manuel/git/ecceTERA/bin/ecceTERA"

#INDELible must also be installed.  The executable must be in the system PATH, ie we need to be able to launch the command
#> indelible
#To change that, modify generate_families_with_simphy.py



for dl, tr, pop, sites in product(dlmult_list, transfermult_list, pop_list, sites_list):

    #eg run: python run_exp.py --phylomethod=iqtree --skipexisting -r 25 --sites=300 --dlrate=2 --trate=1 --pop=10000000 -o my_exp_p1e7_dl2_t1_s300
    
    pop = int(pop)  #otherwise it uses float and causes errors
    
    outdir = f"allexp_p{pop}_dl{dl}_t{tr}_s{sites}"
    
    command = f"python run_exp.py --phylomethod=iqtree --skipexisting -r {runs} --sites={sites} --dlrate={dl} --trate={tr} --pop={pop} -o {outdir}"
    
    #set all binary paths
    command += f" --iqtreebin={iqtree_bin} --simphydir={simphy_dir} --simphybin={simphy_bin} --ecceterabin={eccetera_bin} "
    
    print("Running: " + command)
    
    os.system(command)

