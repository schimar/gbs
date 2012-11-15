
from __future__ import division

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd

###

def import_raw_loci(filename):
    '''Retrieve sequencing data from the text file and store it in a numpy array'''
    return np.genfromtxt(filename, dtype=str, delimiter='\t')


def sort_loci(data):
    '''Strips the column headers as well as the first column from the input array and sorts the loci and individuals according to highest abundance of bases'''
    indivID = data[0, 1:]
    locus = data[1:, 0]
    data_strip_header = sp.delete(data, 0, 0)
    data_strip_loci_row = sp.delete(data_strip_header, 0, 1)
    data_not_N = data_strip_loci_row != 'N'
    row_sums = -1*(np.sum(data_not_N, axis=1))
    col_sums = -1*(np.sum(data_not_N, axis=0))
    row_order = row_sums.argsort()
    col_order = col_sums.argsort()
    indivID_sorted = indivID[col_order]
    locus_sorted = locus[row_order]
    data_sorted = data_strip_loci_row[row_order][:, col_order]
    return data_sorted

# split values in the alleles column (and probably handle that in sort_loci()

for allele in alleles:
    first = allele[0]
    second = allele[2]
    first_allele.append(first)
    second_allele.append(second)

data = np.insert(data, 0, first_allele, axis = 1)
data = np.insert(data, 1, second_allele, axis = 1)
data = sp.delete(data, 2, 1)

#######





########################################################################
# PLOTTING

# Now plot Locus distributions:
loci_per_ind = -1*(col_sums)
plt.hist(loci_per_ind, 60)
plt.xlabel('Number of loci Represented', fontsize=20)
plt.ylabel('Number of Individuals', fontsize= 20)
plt.savefig(pp, format='pdf')
#plt.show()#This will need to change to plt.figure for multiple graphs
# Go to http://matplotlib.org/gallery.html - try heat map. Also explore saving figures
#

# Next part works nicely on its own (when above plot removed) but fails when both together
# Also ask Ethan about savefile
# Now plot Individual distributions:
inds_per_locus = -1*(row_sums)
plt.hist(inds_per_locus, 140)# you can simply measure the length ofindivID or max of
# Inds_per_locus to make an appropriate number of bins!
plt.xlabel('Number of loci Represented', fontsize=20)
plt.ylabel('Number of Inds', fontsize= 20)
plt.savefig(pp, format='pdf')
#plt.show()

pp.close()
plot_data = data_sorted.copy()
#plot data is the sorted array(descening) converted to binary for plotting

for (x,y), value in np.ndenumerate(plot_data):
     if value == "N":
          plot_data[x,y] = 0
     else:
          plot_data[x,y] = 1
plot_data = np.array(plot_data, dtype=int)

# plot_data is now a simple binary array. But it seems hard to plot. Here is an alternative:
# Redo the above loop over plot data.sorted. ndenumerate returns the x,y coordinates (above
# I used it to loop through the array). But you could use it to create plot_data as
# x,y coordinates that can then feed into a scatter plot: see some great ideas at:
# www.prettygraph.com/blog/how-to-plot-a-scatter-plot-using-matplotlib/

# attempt at a start
#plot_data = np.array()
#for (x,y), value in np.ndenumerate(data_sorted):
#     if value == "N":
#          plot_data.append(index)
#
#print plot_data

if __name__ == "main":

    raw_data = import_raw_loci('hmp_play.txt')

    data = sort_loci(raw_data)
