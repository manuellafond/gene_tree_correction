import os
import sys
from pathlib import Path
import re
import util



def extract_rates(control_file: Path):
    results = {}
    current_model = None

    with control_file.open() as f:
        for line in f:
            line = line.strip()

            # Detect model section
            if line.startswith("[MODEL]"):
                # e.g. "[MODEL] modelA_3"
                current_model = line.split()[1].split("_")[1]

            # Detect rates line inside a model
            elif line.startswith("[rates]") and current_model:
                parts = line.split()
                # eg parts = ['[rates]', '0', '2.139077', '0']
                second_number = float(parts[2])
                results[current_model] = second_number

    return results






if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python get_rf_averages.py /path/to/directory")
        sys.exit(1)

    directory = sys.argv[1]
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory")
        sys.exit(1)
        
    outfilename = "rate_vs_support.csv"
    outfile = open(outfilename, "w")
    
    outfile.write("rate,avgbs\n")
        
    rf_per_params = {}
        
    for childdir in Path(directory).iterdir():
        #at this point we should be in a simphy dir
        
        control_file = childdir / "1" / "control.txt"
        
        if not control_file.is_file():
            print(str(control_file) + " does not exist")
            continue        
            
        rates_per_sim = extract_rates(control_file)
        
        for i in range(1, 101):
           bstree = os.path.join(str(childdir), "pargenes/trees/supports_run/results", f"dataset_{i:03d}_TRUE_phy.support.raxml.support")
            
            
           print(bstree)
            
           if os.path.exists(bstree):
            
               supps = []
               util.count_support_bins(bstree, all_supports_list = supps)
    
    
               avg = float(sum(supps)) / float(len(supps))
            
               if len(supps) > 5:
                   outfile.write(f"{rates_per_sim[str(i)]},{avg}\n")
               
    outfile.close()
 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
