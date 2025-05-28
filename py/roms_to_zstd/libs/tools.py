import os
import netCDF4 as nc
import numpy as np
import re
    
def getlist(headpath, tail):
       
    files = []
    for file in os.listdir(headpath):
        
        if tail in file:
            files.append(headpath+file)
            
    return files




class ExtractVar:
    
    def __init__(self, namelistfile):
        """Initialize the class by reading the namelist file into memory."""
        with open(namelistfile, 'r') as file:
            self.namelist_lines = file.readlines()  
        
    def getvar(self, var):
        pattern = rf"^\s*{var}\s*=\s*(.*)"  
        for line in self.namelist_lines:
            match = re.search(pattern, line)
            if match:
                value_part = match.group(1).split("!")[0].strip()  
                num_match = re.search(r"[-+]?\d*\.?\d+[dDeE]?[+-]?\d*", value_part)  
                if num_match:
                    return float(num_match.group().replace('d0', '').replace('e0', ''))
                else:
                    print(f"Error: No numerical value found for variable {var}")
                    return None
        print(f"Error: Variable {var} not found in namelist file")
        return None