# gene_tree_correction

The main script to call is run_exp.py.  It launches a given number of simphy runs along with indelible alignments, reconstructs all gene trees with iqtree, corrects them with eccetera, and then compares the iqtree and eccetera trees with the true gene trees.  For each run, a direcotry is created.  Each such directory contains an "rf" directory containing the rf distances, one per file. 
Please use

> python run_exp.py --help

for details on the arguments.
Example usage:

> python run_exp.py --phylomethod=iqtree --skipexisting -r 25 --sites=300 --dlrate=2 --trate=1 --pop=10000000 -o my_exp_p1e7_dl2_t1_s300


#Generating csv data

The data can be aggregated in a csv file using the script make_genetree_csv.py.  This script will go through all the simphy runs in a specified directory to extract stats on each single gene tree (provided iqtree reconstructed them, e.g., if a gene tree has two leaves iqtree will not reconstruct it). All values that could not be computed, for instance the rf of a non-existing tree, appear as -1 in the csv file.
Run

> python make_genetree_csv.py --help

for more info.
Example usage: 

> python make_genetree_csv.py -d=my_exp_p1e7_dl2_t1_s300/ -o test.csv
