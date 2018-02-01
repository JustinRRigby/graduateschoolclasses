#!python

###############################################################
# NWV 6 Jan 2013                                            
# This script will read in a multiple sequence alignment (MSA) in a 
# user-specified format and return a histogram of the numbers 
# of sites in the alignment with a specific number of         
# alternative amino acids / nucleotides.                      
###############################################################

#updated 04 Jan 2016 CT
#updated 20 Jan 2014 EWJW
#updated 31 Dec 2013 DAD
#updated 24 Nov 2013 TNS



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

##2 Change these filenames!
##2 also, results_file changed to hist_results_filename, 
##2 hist_plot_filename changed to hist_plot_filename
alignment_filename = "hw2/data/SR214-plus-1GWRa.fasta"
hist_results_filename = "hw2/results/AlignmentAnalysisHist.txt"
hist_plot_filename = "hw2/results/SiteDiversityHist.png"
##2 here is an ideal place to define a second results file

# Open the results file so we can write ('w') to it. This statment
# literally says that the variable 'hist_results_file'  is the document
# "c:/Dropbox/evol-biol-mol/homeworks/hw1/results/AlignmentAnalysis.txt"
# (see what os.path.join does?) and enables that program to write to it.

hist_results_file = open(os.path.join(base_directory,hist_results_filename),'w')
##2 Open your new results file for writing. And close it when you're done.

# Use Biopython to read in the entire alignment. Read the Biopython
# HOWTOs for the AlignIO module at http://biopython.org/wiki/AlignIO
# This line essentially says "Use the AlignIO module to read the
# 'alignment_file'(which is our MSA in fasta # format) and store that
# information in the variable 'alignment',which our script can now go on
# and use.

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
##2 Added comments to these lines, so scripts that expect a table do not get confused
hist_results_file.write("# Alignment length: {L}\n".format(L=aligned_seq_length))
hist_results_file.write("# Number of sequences in the alignment: {n}\n".format(n=num_seqs))

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
##2 replaced site_counts = [] with allocating a fixed-size list,
## because incrementing large variables can be slow and problematic
site_counts = [0] * aligned_seq_length

# To identify particular sites and compare those to published studies, we need residue numbers
# Residue numbers in our alignment don't necessarily track those used in the literature,
# and we must not count gaps.

##2 Change this to agree with residue numbering in 1GWRa.pdb
published_site_number = 313 # starting residue number in the alignment

# Check to make sure the published sequence index is the right one
##2 this is a crucial diagnostic! edit it!
published_sequence_index = 0
print(alignment[published_sequence_index].id)
assert alignment[published_sequence_index].id.startswith('1GWRa')

##2 If you are writing to the results file inside your for loop,
##2 this would be an ideal place to write a header line, if desired 

# For each site...
for site_number in range(aligned_seq_length):
    # pull out the site data, and give it a convenient name
    site_aas = alignment[:,site_number]
    
    ##2 this block moved outside the if fraction_gaps clause
    ##2 now, the number of unique aas is calculated at every site
    # find the number of unique amino acids at the site. We use
    # a trick: put the amino acids in a set (which, like a
    # mathematical set, only has one of each item), then
    # ask how big (long) the set is.
    number_of_unique_aas = len(set(site_aas))


    # increment the site number as used in publications,
    # but only if there's no gap.
    if site_aas[published_sequence_index] != '-':
        ##2 Here is an ideal place to write a line to your sitewise results file
	print published_site_number, " ", number_of_unique_aas
        published_site_number += 1 ##2 increment AFTER output
	

    # debugging 
    ##2 very useful for comparing to 1GWRa.pdb--change expected amino acids accordingly. 
    #if published_site_number in [69,70,182]:
    	# this should print M, S and M
    	##2 print statement was edited to make more informative.
	#print "published site", published_site_number, ", residue", site_aas[published_sequence_index]

    # if less than 20% of the sequences have gaps...
    fraction_gaps = site_aas.count('-')/float(num_seqs)
    if fraction_gaps < 0.2:
        # ...and add 1 to the histogram at that number.
        histogram[number_of_unique_aas] += 1

        # identify sites with only one unique AA
        #if number_of_unique_aas == 1:
        #    print published_site_number, site_aas[published_sequence_index]
    
    ##2 this block also moved outside fraction_gaps clause
    ##2 also replaced appending by indexing
    ##2 for HW1 this would have made the plotted and outputted histograms disagree
    #site_counts += [number_of_unique_aas]
    site_counts[site_number] = number_of_unique_aas


##2 Write number of amino acids per site results to hist_results_file
# write header with names of fields
hist_results_file.write("number.of.aas\tnumber.of.sites\n")
# populate those fields for each possible number of aas
for ind in range(len(histogram)):
    hist_results_file.write("%s\t%s\n" % (ind, histogram[ind]))
# We're done with this file, so close it.
hist_results_file.close()
print("Wrote site diversity histogram to " + hist_results_filename)

##2 plot histogram of number of aas per site
pyplot.hist(site_counts, bins=range(0,20,1), facecolor='blue')
##2 directly uses our constructed histogram data
pyplot.bar(left=range(len(histogram)), height=histogram)
pyplot.xlabel('Number of unique amino acids')
pyplot.ylabel('Count')
pyplot.title('Frequency of Sites with a Given\nNumber of Amino Acids')
pyplot.savefig(hist_plot_filename)
pyplot.show() # Show this histogram live
print("Plotted histogram in " + hist_plot_filename)

##2 plot bar chart of number of aas per site, by site; diagnostic output
pyplot.figure(2)
pyplot.bar(range(aligned_seq_length),site_counts,color='black')
pyplot.xlabel('site')
pyplot.ylabel('Number of unique amino acids')
pyplot.title('Frequency of Sites with a Given\nNumber of Amino Acids')
pyplot.show() # Show this figure live
