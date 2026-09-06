# gene_tree_correction

The main script to call is run_exp.py.  It launches a given number of simphy runs along with indelible alignments, reconstructs all gene trees with iqtree, corrects them with eccetera, and then compares the iqtree and eccetera trees with the true gene trees.  For each run, a direcotry is created.  Each such directory contains an "rf" directory containing the rf distances, one per file. 

Dependencies: iqtree (version >= 3), Simphy, INDELible, ecceTERA must be installed.  run_exp.py needs ete3 and Biopython.
You can specify the path to the executables in the command line, EXCEPT INDELible which must be in the PATH.  

Please use

> python run_exp.py --help

for details on the arguments.
Example usage:

> python run_exp.py --phylomethod=iqtree --skipexisting -r 25 --sites=300 --dlrate=2 --trate=1 --pop=10000000 -o my_exp_p1e7_dl2_t1_s300

Note, there are various arguments to specify the path to the required executables (iqtree, eccetera, simphy).

There is also the script

> python run_all.py

which will just run every combination of dl rate, tr rate, and population.



# Running with Astral-Pro

If you add --apro_mode to the list of command line arguments, run_exp.py will infer the Astral-Pro species tree in each simphy directory (that new species tree is added to the
simphy directory).
It will then run eccetera using this species tree --- in addition to running eccetera with the simphy species tree.

You need to specify where to find Astral-Pro using
> --apro_mode --apro_libpath=[your path to the lib directory of Astral_pro]
> --apro_binname=[your path to the Astral-Pro jar file]

See the maximal example below.



# Running with Random NNIs

If you add --nni_mode, run_exp.py will apply random NNIs on the simphy species tree.  You can specify a list of random NNIs to apply, and a species tree will be generated 
for each number in that list.  For example

> --nni_mode --nni_k_list 1 2 10 

will generate one species tree obtained from one NNI, then another species tree with 2 NNI, then another with 10 NNI.

eccetera wil lthen be run on all these species trees, in addition to the simphy one (and the apro one if specified).



If you add --nnicluster_mode, run_exp.py will apply random *cluster NNIs* on the simphy species tree. 
A cluster NNI applies three consecutive NNIs on adjacent edges, i.e., it chooses a directed path of length 3 and does an NNI on each edge of the path.  
This is to simulate a "local" portion of the tree that is erroneous.

You can specify a list of cluster NNIs to apply, and a species tree will be generated for each number in that list.  For example

> --nnicluster_mode --nniclusters_list 1 2 3

will generate three species trees.  The first has one modified cluster, the second has two, the third has three (so, that third tree undergoes 9 NNIs in total).

eccetera wil lthen be run on all these species trees, in addition to the simphy one (and the apro one if specified).


# Maximal example

> python run_exp.py --phylomethod=iqtree --skipexisting -r 5 --sites=300 --dlrate=1 --trate=1  --eccetera_bs=50 --pop=10 -o aprotest_p10_dl1_t1_s300 --iqtreebin=/home/manuel/git/iqtree-3.0.1-Linux/bin/iqtree3 --simphydir=/home/manuel/git/SimPhy_1.0.2/ --simphybin=/home/manuel/git/SimPhy_1.0.2/bin/simphy_lnx32 --ecceterabin=/home/manuel/git/ecceTERA/bin/ecceTERA  --apro_mode --apro_libpath=/home/manuel/git/A-pro/ASTRAL-MP/lib --apro_binname=/home/manuel/git/A-pro/ASTRAL-MP/astral.1.1.6.jar --nni_mode --nni_k_list 1 2 10 --nnicluster_mode --nniclusters_list 1 2 3




# Generating csv data

The data can be aggregated in a csv file using the script make_genetree_csv.py.  This script will go through all the simphy runs in a specified directory to extract stats on each single gene tree (provided iqtree reconstructed them, e.g., if a gene tree has two leaves iqtree will not reconstruct it). All values that could not be computed, for instance the rf of a non-existing tree, appear as -1 in the csv file.
Run

> python make_genetree_csv.py --help

for more info.
Example usage: 

> python make_genetree_csv.py -d my_exp_p1e7_dl2_t1_s300/ -o test.csv

If you have multpiple directories with data, add the flag "--append" to the above command to add lines to the csv (instead of overwriting it).


# Other notes

test.csv contains a csv of partial experiments with dup-loss rate x2, transfer rate x1, no ILS, 300 sites.
test_pivot.ods is the same, but there is a tab with a pivot table that analyzes the RF with respect to average bootstrap.

The python scripts that contain "speciesrax" in the name were there to test simphy generation and inference with speciesrax and others.  Also, the directory "archive" contains several older files used to validate the methodology.
