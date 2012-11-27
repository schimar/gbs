
from __future__ import division

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd

###
def filter_single_col(base_list, count_list, threshold = 4):
    '''Returns a list of nucleotides filtered (threshold, default = 4) by number of occurence of a sequencing run at specific loci in a list of bases'''
    output = []
    for i, base in enumerate(base_list):
        count_1, count_2 = count_list[i].split('|')
        if base == 'N':
            value = 'N'
        else:
            if int(count_1) < threshold and int(count_2) < threshold:
                value = 'N'
            else:
                value = base
        output.append(value)
    return output

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



# you may have to mirror those changes in the 'order' files, too
# or: maybe return those (row_id and col_header) as well ?
# pdata = import_raw_loci('hmp_play.txt')
# pdata = sp.delete(pdata, [2,3,4,5,6,7,8,9,10],1)




if __name__ == "main":

    # filter

    hmp = pd.read_table('hmp_play.txt', index_col = 0, header = 0)
    hmc = pd.read_table('hmc_play02.txt', index_col = 0, header = 0)
    #######

    first_allele = []
    second_allele = []
    for allele in hmp.ix[:, 'alleles']:
        first_allele.append(allele.split('/')[0])
        second_allele.append(allele.split('/')[1])

    hmp = hmp.ix[:, 'FLFL04':'WWA30']
    hmp.insert(0, 'allele_1', first_allele)
    hmp.insert(1, 'allele_2', second_allele)
    ###
    hmp_trimmed = hmp.ix[:30,2:]
    results = []
    for i, col in enumerate(hmp_trimmed.columns):
        base_values = hmp_trimmed[col]
        count_values = hmc[hmc.columns[i]]
        results.append(gbs.filter_single_col(base_values, count_values))

    results.insert(0, list(hmp.allele_1))
    results.insert(1, list(hmp.allele_2))

    df = pd.DataFrame(zip(*results), index = hmp.index[:30], columns=hmp.columns, dtype = np.str)


    # sort loci
    raw_data = import_raw_loci('hmp_play.txt')

    data = sort_loci(raw_data)
