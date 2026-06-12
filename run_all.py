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




for dl, tr, pop, sites in product(dlmult_list, transfermult_list, pop_list, sites_list):

    #eg run: python run_exp.py --phylomethod=iqtree --skipexisting -r 25 --sites=300 --dlrate=2 --trate=1 --pop=10000000 -o my_exp_p1e7_dl2_t1_s300
    
    outdir = f"allexp_p{pop}_dl{dl}_t{tr}_s{sites}"
    
    command = f"python run_exp.py --phylomethod=iqtree --skipexisting -r {runs} --sites={sites} --dlrate={dl} --trate={tr} --pop={pop} -o {outdir}"
    
    print("Running: " + command)
    
    os.system(command)

