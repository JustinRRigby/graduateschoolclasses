#!python

###############################################################
# NWV 6 Jan 2013                                            
# This script will read in a multiple sequence alignment (MSA) in a 
# user-specified format and return a histogram of the numbers 
# of sites in the alignment with a specific number of         
# alternative amino acids / nucleotides.                      
###############################################################

#updated 9 Jan 2014 MS
#updated 31 Dec 2013 DAD
#updated 24 Nov 2013 TNS
#updated 03 Jan 2017 CGT <- update to 3.5


import os, re, sys
from Bio import AlignIO
from matplotlib import pyplot

# "chdir" = "change directory" Tell Python to change where it
# will look for files. For example, you'll probably be storing
# your evol-biol-mol files in a directory like
# C:/Users/Nicholas/Documents/EvolBiolMol/HW (for windows users)
# or /Users/Mo/Dropbox/evol-biol-mol/homeworks (mac/linux users)
# you can write the path above as ~/Mo/Dropbox/evol-biol-mol/homeworks

base_directory = "/home/justin/Desktop/graduateschoolclasses/Winter_2018/Evolution_of_Biological_Molecules/"
os.chdir(base_directory)

# Define a variable for each of your input and output files.
# If you have your directory structured correctly, you should not
# have to change these lines.

alignment_filename = "hw1/data/TEM1_aln_with_PSE_seq.fasta"
results_filename = "hw1/results/AlignmentAnalysis.txt"
histogram_filename = "hw1/results/SiteDiversityHist.png"

# Open the results file so we can write ('w') to it. This statment
# literally says that the variable 'results_file'  is the document
# "c:/Dropbox/evol-biol-mol/homeworks/hw1/results/AlignmentAnalysis.txt"
# (see what os.path.join does?) and enables that program to write to it.

##results_file = file(os.path.join(base_directory,results_filename),'w')
results_file = open(results_filename,'w')

# Use Biopython to read in the entire alignment. Read the Biopython
# HOWTOs for the AlignIO module at http://biopython.org/wiki/AlignIO
# This line essentially says "Use the AlignIO module to read the
# 'alignment_file'(which is our sequence alignment in fasta # format) and 
# store that information in the variable 'alignment',which our script
# can now go on and use.

alignment = AlignIO.read(alignment_filename, "fasta")

# Print a summary of the alignment to the screen.

print(alignment)

# Store the length of the alignment in the variable 'aln_length' and print
# to screen. len(alignment[0,:]) asks "How long is the first line
# in the alignment?". Similarly, we could ask how long the second line is
# using len(alignment[1,:]). You can think of the alignment, in this case,
# as an excel sheet where alignment[1,:] would be the entirety of row 2. In
# this way we can access any element of the alignment. I recommend reading
# the documentation at http://docs.python.org/2/tutorial/index.html. [x,y]
# are coordinates, [row, column]

aligned_seq_length = len(alignment[0,:])
print("Length of aligned sequences (including gaps): {}".format(aligned_seq_length))

# Store the number of sequences in your alignment and print it to your
# screen.

num_seqs = len(alignment)
print("Number of sequences in the alignment: {}".format(num_seqs))

# Print this information to our results file for safe keeping. 

results_file.write("Alignment length: {L}\n".format(L=aligned_seq_length))
results_file.write("Number of sequences in the alignment: {n}\n".format(n=num_seqs));

# Now to the meat of the script. The goal: Determine, for each column (site)
# in the alignment, how many _different_ amino acids are present (a measure
# of how conserved the site is).

# Result: Write, and plot, a histogram of the number of sites that contain a
# particular number of different AAs.

# Initialize an array with all zeros. See the Python tutorial on data structure
# to better understand what an array (vs. a dictionary for example) is.

histogram = [0]*21

# Position 1 of the histogram counts the sites with 1 amino acid, etc.

# In site_counts we'll collect a list of all the amino-acid counts at each site.

site_counts = []

# To identify particular sites and compare those to published studies, we need residue numbers
# Residue numbers in our alignment don't necessarily track those used in the literature,
# and we must not count gaps.

published_sequence_index = 0
published_site_number = 25 # starting residue number in the alignment

# Check to make sure the published sequence index is the right one
print(alignment[published_sequence_index].id)
assert alignment[published_sequence_index].id.startswith('TEM-1-imipenem')

# For each site...
for site_number in range(aligned_seq_length):
    # pull out the site data, and give it a convenient name
    site_aas = alignment[:,site_number]

    # increment the site number as used in publications,
    # but only if there's no gap.
    if site_aas[published_sequence_index] != '-':
        published_site_number += 1

    # debugging
    #if published_site_number in [69,70,182]:
    	# this should print M, S and M
    	#print site_aas[published_sequence_index]

    # if less than 20% of the sequences have gaps...
    fraction_gaps = site_aas.count('-')/float(num_seqs)
    if fraction_gaps < 0.2:
        # find the number of unique amino acids at the site. We use
        # a trick: put the amino acids in a set (which, like a
        # mathematical set, only has one of each item), then
        # ask how big (long) the set is.
        number_of_unique_aas = len(set(site_aas))
        if number_of_unique_aas == 1:
          print(published_site_number)
        # ...and add 1 to the histogram at that number.
        histogram[number_of_unique_aas] += 1
        site_counts += [number_of_unique_aas]

        # identify sites with only one unique AA
	if number_of_unique_aas == 1:
		print published_site_number, site_aas[published_sequence_index]

            
results_file.write("number.of.aas\tnumber.of.sites\n")
for ind in range(len(histogram)):
    results_file.write("%s\t%s\n" % (ind, histogram[ind]))
# We're done with this file, so close it.
results_file.close()
print("Wrote results to " + results_filename)

pyplot.hist(site_counts, bins=range(0,20,1), facecolor='blue')
pyplot.xlabel('Number of unique amino acids')
pyplot.ylabel('Count')
pyplot.title('Frequency of Sites with a Given\nNumber of Amino Acids')
pyplot.savefig(histogram_filename)
pyplot.show() # Show this histogram live
print("Plotted histogram in " + histogram_filename)
