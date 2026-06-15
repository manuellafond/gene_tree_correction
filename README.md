# gene_tree_correction

The main script to call is run_exp.py.  It launches a given number of simphy runs along with indelible alignments, reconstructs all gene trees with iqtree, corrects them with eccetera, and then compares the iqtree and eccetera trees with the true gene trees.  For each run, a direcotry is created.  Each such directory contains an "rf" directory containing the rf distances, one per file. 

Dependencies: iqtree and Simphy must be installed.  run_exp.py needs ete3 and Biopython.

Please use

> python run_exp.py --help

for details on the arguments.
Example usage:

> python run_exp.py --phylomethod=iqtree --skipexisting -r 25 --sites=300 --dlrate=2 --trate=1 --pop=10000000 -o my_exp_p1e7_dl2_t1_s300

Note, there are various arguments to specify the path to the required executables (iqtree, eccetera, simphy).

There is also the script

> python run_all.py

which will just run every combination of dl rate, tr rate, and population.

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
