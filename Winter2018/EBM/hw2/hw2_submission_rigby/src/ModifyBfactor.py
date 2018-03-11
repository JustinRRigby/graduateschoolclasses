#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  
#
# Modified by NWV 18 Jan 2013
# Modified by EWJW 20 Jan 2014
# Modified by CT 04 Jan 2016 -> Python 3.5 compatible

import os, sys

# choose working directory and filenames before you execute
#base_directory = os.path.expanduser("/Users/tylernstarr/Dropbox/evol-biol-mol/homeworks")
base_directory = "/home/justin/Desktop/graduateschoolclasses/Winter_2018/Evolution_of_Biological_Molecules/"
pdb_file = "hw2/data/1GWRa.pdb"
new_bfact = "hw2/results/b-factorList.txt"
new_pdb_file = "hw2/results/1GWRa_New_B_Factor.pdb"

# Define a new function which will read in your new_bfact file

def loadDataFile(data_file,key_field=0,data_field=1):
   
    # Read in file
    f = open(data_file,'r')
    data = f.readlines()
    f.close()

    # Strip out blank lines and comments
    data = [d for d in data if d[0] != "#" and d.strip() != ""]

    # Make an empty dictionary to store your information
    data_dict = {}

    # for each line...
    for record in data:
        try:
            # Add the key (i.e. residue number) to the dictionary and replace its b-value
            field = record.split()
            key = field[key_field]
            data_dict.update([(key,float(field[data_field]))])

        # Write an error if your b-factor file is messed up            
        except IndexError:
            sys.stderr.write("Mangled data, skipping line:\n%s" % record)
            continue
        except ValueError:
            sys.stderr.write("Mangled data, skipping line:\n%s" % record)
            continue
    # Return the new b-factor dictionary       
    return data_dict


# Define a function to apply the b-factors to the old pdb file

def pdbBfactor(pdb_file,data_dict):
    
    # Goes through pdb line by line.  If residue is in dictionary data_dict,
    # the b-factor is replaced in output_file by the dictionary value.  If the
    # residue is not in the dictionary, the b-factor is given value 0.0.
    # Returns void.
    
    # Make an empty list
    out = []

    # For each line in the pdb file (which you read in in the main() function, see below
    for line in pdb_file:
        # If the line begins with ATOM
        if line[0:6] == "ATOM  ":
            # Characters 23:36 are the residue numbers
            resnum = line[23:26].strip()
            # If the residue number is in your dictionary from loadDataFile...
            if resnum in data_dict.keys():
                # Add the new B-factor to it
                out.append("%s%6.2F%s" % (line[:60],data_dict[resnum],
                                          line[66:]))
            # Else write an NA for that B-factor value
            else:
                out.append("%s%6s%s" % (line[:60],"NA",line[66:]))
        
        elif line[0:6] == "HETATM":
            out.append("%s%6s%s" % (line[:60],"NA",line[66:]))
        
        else:
            out.append(line)
    
    return out

####################
# Now the functions are defined, we need to execute them

# Change the working directory
os.chdir(base_directory)

# Open an output file for writing
outfile = open(os.path.join(base_directory,new_pdb_file), 'w')

# Read in pdb file
print("reading pdb_file", pdb_file)
f = open(os.path.join(base_directory,pdb_file),'r')
pdb = f.readlines()
f.close()

# Execute loadDataFile
print("reading data file", pdb_file)
data_dict = loadDataFile(os.path.join(base_directory,new_bfact))

# Execute pdbBfactor
print("replacing b-factors with substitute output from data_file"),
out = pdbBfactor(pdb,data_dict)

# For each line in the list "out", print it to your outfile
print("writing output to ", new_pdb_file)
for line in out:
    outfile.write(line)
# Close the outfile
outfile.close()

print("script complete")
