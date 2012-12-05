
from __future__ import division
import csv

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd

###
def filter_single_col(base_list, count_list, threshold = 4):
    '''Returns a list of nucleotides filtered (threshold, default = 4) by number of occurence of a sequencing run at specific loci in a list of bases'''
    output = []
    for i, base in enumerate(base_list):
        count_1, count_2 = map(int, count_list[i].split('|'))
        if base == 'N':
            value = 'N'
        else:
            if count_1 < threshold and count_2 < threshold:
                value = 'N'
            else:
                value = base
        output.append(value)
    return output

def get_loci_from_iupac_codes(base_list, count_list, allele_list, threshold = 4):
    '''Returns a list of nucleotides where ambiguity codes have been changed to their respective value based on a list of measured allele levels'''
    ambig = []
    value = ''
    for i, base in enumerate(base_list):
        count_1, count_2 = map(int, count_list[i].split('|'))
        allele_1, allele_2 = allele_list[i].split('/')
        if base == 'N':
            value = 'N'
        else:
            if count_1 >= threshold and count_2 < threshold:
                value =  str(allele_1 + '/' + allele_1)
            elif count_1 < threshold and count_2 >= threshold:
                value = str(allele_2 + '/' + allele_2)
            elif count_1 >= threshold and count_2 >= threshold:
                value = str(allele_1 + '/' + allele_2)
        ambig.append(value)
    return ambig


def import_raw_loci(filename):
    '''Retrieve sequencing data from the text file and store it in a numpy array'''
    return np.genfromtxt(filename, dtype=str, delimiter='\t')


def sort_loci_def(data, column_headers, row_headers):
    '''Strips the column headers as well as the first column from the input np.array and sorts the loci and individuals according to highest abundance of bases'''
    indivID = data[0, 1:]
    locus = data[1:, 0]
    data_strip_header = sp.delete(data, 0, 0)
    data_strip_loci_row = sp.delete(data_strip_header, 0, 1)
    data_not_N = data_strip_loci_row != 'N'
    row_sums = -1*(np.sum(data_not_N, axis=1))
    col_sums = -1*(np.sum(data_not_N, axis=0))
    row_order = row_sums.argsort()
    col_order = col_sums.argsort()
    indivID_sorted = indivID[row_order] # was changed (i.e. indivID[col_order])
    locus_sorted = locus[col_order]
    data_sorted = data_strip_loci_row[row_order][:, col_order]
    return data_sorted

def sort_loci_np(data, column_headers, row_headers):
    '''Strips the column headers as well as the first column from the input np.array and sorts the loci and individuals according to highest abundance of bases'''
    data_not_N = data != 'N'
    row_sums = -1*(np.sum(data_not_N, axis=1))
    col_sums = -1*(np.sum(data_not_N, axis=0))
    row_order = row_sums.argsort()
    col_order = col_sums.argsort()
    indivID_sorted = column_headers[col_order] # was changed (i.e. indivID[col_order])
    locus_sorted = row_headers[row_order]
    data_sorted = data[row_order][:, col_order]
    return data_sorted

def sort_loci_pdDF(data):
    '''Strips the column headers as well as the first column from the input pd.DataFrame and sorts the loci and individuals according to highest abundance of bases'''
    data_strip_loci = data.ix[:, 'FLFL04': 'WWA30']
    data_not_N = data_strip_loci != 'N'
    row_sums = -1*(np.sum(data_not_N, axis=1))
    col_sums = -1*(np.sum(data_not_N, axis=0))
    row_order = row_sums.argsort()
    col_order = col_sums.argsort()
    indivID_sorted = data.index[row_order]
    locus_sorted = data.columns[col_order]
    data_sorted = np.array(data_strip_loci)[row_order][:, col_order]
    data_sorted = pd.DataFrame(data_sorted, index = indivID_sorted, columns= locus_sorted, dtype = np.str)
    return data_sorted



if __name__ == "main":

    # filter

    hmp = pd.read_table('HapMap.hmp.txt', index_col = 0, header = 0)
    hmc = pd.read_table('HapMap.hmc.txt', index_col = 0, header = 0)

    #######
    hmp_trimmed = hmp.ix[:, 'FLFL04':'WWA30']

    filter_results = []
    for i, col in enumerate(hmp_trimmed.columns):
        base_values = hmp_trimmed[col]
        count_values = hmc[hmc.columns[i]]
        filter_results.append(filter_single_col(base_values, count_values))

    df = pd.DataFrame(zip(*filter_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)

    unambiguous_results = []
    for i, col in enumerate(hmp_trimmed.columns):
        base_list = hmp_trimmed[col]
        count_list = hmp_trimmed[hmc.columns[i]]
        unambiguous_results.append(get_loci_from_iupac_codes(base_list, count_list, hmp.alleles))

    df_unambiguous = pd.DataFrame(zip(*unambiguous_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)


    data_sorted4 = sort_loci_pdDF(df_unambiguous)

    data_sorted4.insert(0, 'alleles', hmp.ix[:, 'alleles'])

    df2 = hmp.ix[:, 'chrom': 'QCcode']
    df = data_sorted4.join(df2)

    # and now the ambig_filter... better before the insert/join


    # some plotting to see the distributions

    df.to_csv("data_sorted4.csv")



