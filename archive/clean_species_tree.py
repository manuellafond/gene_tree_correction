from Bio import Phylo
import sys

def clean_newick(input_file, output_file):
    # Read the tree
    tree = Phylo.read(input_file, "newick")

    # Remove single quotes from the entire file content
    with open(input_file, "r") as f:
        newick_str = f.read().replace("'", "")

    # Remove internal node names
    for clade in tree.find_clades():
        if not clade.is_terminal():  # If it's an internal node
            clade.name = None  # Remove the name

    # Write back the cleaned tree
    Phylo.write(tree, output_file, "newick")

    # Ensure the single quotes are removed by re-processing the output file
    with open(output_file, "r") as f:
        cleaned_str = f.read().replace("'", "")

    with open(output_file, "w") as f:
        f.write(cleaned_str)

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    clean_newick(input_file, output_file)
