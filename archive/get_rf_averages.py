import os
import sys
from pathlib import Path
import re

def read_single_number(filepath):
    """Read a file containing a single number."""
    with open(filepath, "r") as f:
        line = f.readline().strip()
        return float(line)




def process_directory(directory):
    eccetera_values = []
    raxml_values = []
    count = 0

    # Iterate over all files in directory
    for filename in os.listdir(directory):
        if filename.startswith("eccetera_") and filename.endswith(".rf"):
            # Get index i
            i_part = filename[len("eccetera_"):-len(".rf")]
            ecc_file = os.path.join(directory, filename)
            raxml_file = os.path.join(directory, f"raxml_{i_part}.rf")

            # Check if the corresponding raxml file exists
            if os.path.isfile(raxml_file):
                try:
                    ecc_val = read_single_number(ecc_file)
                    rax_val = read_single_number(raxml_file)
                    eccetera_values.append(ecc_val)
                    raxml_values.append(rax_val)
                    count += 1
                except Exception as e:
                    print(f"Warning: could not read files for index {i_part}: {e}")

    if count == 0:
        print("No valid file pairs found.")
        return None

    avg_ecc = sum(eccetera_values) / count
    avg_rax = sum(raxml_values) / count

    #print(f"Number of file pairs read: {count}")
    #print(f"Average for eccetera: {avg_ecc}")
    #print(f"Average for raxml: {avg_rax}")
    
    return { "count" : count, "ecc_sum" : sum(eccetera_values), "rax_sum" : sum(raxml_values) }

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python get_rf_averages.py /path/to/directory")
        sys.exit(1)

    directory = sys.argv[1]
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory")
        sys.exit(1)
        
    phylomethod = "raxml"
    
    if "--iqtree" in sys.argv:
        phylomethod = "iqtree"
        
        
    rf_per_params = {}
        
    for childdir in Path(directory).iterdir():
    
        rf_dirname = "rf"
        if phylomethod == "iqtree":
            rf_dirname = "rf_iqtree"
    
        if childdir.is_dir() and Path(os.path.join(str(childdir), rf_dirname)).is_dir():
            rfparams = process_directory(os.path.join(str(childdir), rf_dirname))
            
            if rfparams is None:
                print(f"dir={str(childdir)} is None")
                #sys.exit()
                continue
                
            #extract duprate from folder name (chatgpt's solution tbh)
            match = re.search(r'_(d\d\.\d_l\d\.\d_t\d\.\d)_', str(childdir))
            if not match:
                match = re.search(r'_(d\d\_l\d\_t\d\.\d)_', str(childdir))
                if not match:
                    print("Regex doesn't match for dir " + str(childdir))
                    continue
            dirparams = match.group(1)
            
            if not dirparams in rf_per_params:
                rf_per_params[dirparams] = {"count" : 0, "ecc_sum" : 0, "rax_sum" : 0}
            rf_per_params[dirparams]["count"] += rfparams["count"]
            rf_per_params[dirparams]["ecc_sum"] += rfparams["ecc_sum"]
            rf_per_params[dirparams]["rax_sum"] += rfparams["rax_sum"]

    for param in sorted(rf_per_params):
        count = rf_per_params[param]["count"]
        ecc_sum = rf_per_params[param]["ecc_sum"]
        rax_sum = rf_per_params[param]["rax_sum"]
        print(param)
        
        
        print(f"ecc_avg = { float(ecc_sum)/float(count) }   count = {count}")
        print(f"{phylomethod}_avg = { float(rax_sum)/float(count) }   count = {count}")


