'''Filter and sorting for GBS (Genotyping by Sequencing) HapMap files'''

from __future__ import division
import csv
import sys
import itertools
import re

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd


###
def drop_N_columns(data, drop_level = 0.9):
    '''Returns a pd.DataFrame, where all columns, which consist of more than 90 per cent (default) 'N's are being dropped'''
    data_dropped = data.copy()
    for i, col in enumerate(data_dropped.columns):
        base_series = data_dropped[col]
        N_to_length = np.sum(base_series == 'N')/len(base_series)
        if N_to_length > drop_level:
            data_dropped = data_dropped.drop(col, axis=1)
    return data_dropped

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

def get_alleles_zero(base_list, count_list, allele_list):
    '''Returns a list of nucleotides where ambiguity codes have been changed to their respective value based on a list of measured allele levels'''
    result = []
    value = ''
    for i, base in enumerate(base_list):
	count_1, count_2 = map(int, count_list[i].split('|'))
	allele_1, allele_2 = allele_list[i].split('/')
	if base == 'N':
	    value = 'NA'
	else:
	    if count_1 > 0 and count_2 == 0:
		value = str(allele_1 + '/' + allele_1)
	    elif count_1 == 0 and count_2 > 0:
		value = str(allele_2 + '/' + allele_2)
	    elif count_1 > 0 and count_2 > 0:
		value = str(allele_1 + '/' + allele_2)
	    else:
		value = 'NA'
	result.append(value)
    return result  

def get_alleles_4base(base_list, count_list, allele_list, threshold = 4):
    '''Returns a list of nucleotides where ambiguity codes have been changed to their respective value based on a list of measured allele levels and the read-depth (default = 4)'''
    ambig = []
    value = ''
    for i, base in enumerate(base_list):
        count_1, count_2 = map(int, count_list[i].split('|'))
        allele_1, allele_2 = allele_list[i].split('/')
        if base == 'N':
            value = 'NA'
        else:
            if count_1 >= threshold and count_2 < threshold:
                value =  str(allele_1 + '/' + allele_1)
            elif count_1 < threshold and count_2 >= threshold:
                value = str(allele_2 + '/' + allele_2)
            elif count_1 >= threshold and count_2 >= threshold:
                value = str(allele_1 + '/' + allele_2)
	    else:
		value = 'NA'
        ambig.append(value)
    return ambig


def get_alleles_adv(base_list, count_list, allele_list, threshold = 4):
    '''Returns a list of nucleotides where ambiguity codes have been changed to their respective value based on a list of measured allele levels. A read-depth threshold (default = 4) is applied. If one allele is over the threshold, but not the second one, it will be checked for having at least double the amount of the threshold, to qualify as heterozygote (if not, it'll be '?'  '''
    ambig = []
    value = ''
    for i, base in enumerate(base_list):
        count_1, count_2 = map(int, count_list[i].split('|'))
        allele_1, allele_2 = allele_list[i].split('/')
        if base == 'N':
            value = 'NA'
        else:
            if count_1 >= threshold and count_2 < threshold:
                if count_1 < 2*threshold:
		    value = str(allele_1 + '/' + '?')
		else: 
		    value =  str(allele_1 + '/' + allele_1)
	    elif count_1 < threshold and count_2 >= threshold:
		if count_2 < 2*threshold:
		    value = str(allele_2 + '/' + '?')
		else:
		    value = str(allele_2 + '/' + allele_2)
            elif count_1 >= threshold and count_2 >= threshold:
                value = str(allele_1 + '/' + allele_2)
	    else:
		value = 'NA' 
	ambig.append(value)
    return ambig

# here: MAF filter

def get_alleles_MAF(base_list, count_list, allele_list, MAF = 0.45):
    '''Returns a list of nucleotides where ambiguity codes have been changed to their respective value based on a list of measured allele levels. The MAF will be calculated and, if smaller than the given MAF (default = 0.45), it will be defined as homozygote, otherwise as heterozygous.'''
    result = []
    value = ''
    for i, base in enumerate(base_list):
        count_1, count_2 = map(int, count_list[i].split('|'))
        allele_1, allele_2 = allele_list[i].split('/')
	if base == 'N':
	    value = 'NA'
	else:
	    if count_1 > 0 and count_2 > 0:
		if count_1 > count_2:
		    if count_2/ (count_1 + count_2) >= MAF:
			value = str(allele_1 + '/' + allele_2)
		    else: 
			value = str(allele_1 + '/' + allele_1)
		elif count_2 > count_1:
		    if count_1/ (count_1 + count_2) >= MAF:
			value = str(allele_1 + '/' + allele_2)
		    else:
			value = str(allele_2 + '/' + allele_2)
		else: 
		    value = str(allele_1 + '/' + allele_2)
	    elif count_1 > 0 and count_2 == 0:
		value = str(allele_1 + '/' + allele_1) 
	    elif count_2 > 0 and count_1 == 0:
		value = str(allele_2 + '/' + allele_2)
	result.append(value)
    return result 

####

def sort_loci_pdDF(data):
    '''Strips the column headers as well as the first column from the input pd.DataFrame and sorts the loci and individuals according to highest abundance of bases'''
    #data_strip_loci = data.ix[:, 'FLFL04': 'WWA30']
    data_not_N = data != 'NA'
    row_sums = -1*(np.sum(data_not_N, axis=1))
    col_sums = -1*(np.sum(data_not_N, axis=0))
    row_order = row_sums.argsort()
    col_order = col_sums.argsort()
    indivID_sorted = data.index[row_order]
    #locus_sorted = data.columns[col_order]
    data_sorted = np.array(data)[row_order]#[:, col_order]
    data_sorted = pd.DataFrame(data_sorted, index = indivID_sorted, columns= data.columns)
    return data_sorted



if __name__ == "__main__":
    # that's old stuff, currently, don't use this from the command line, but instead use the gbs.py version !
    if len(sys.argv > 1):
        hmp = pd.read_table(sys.argv[2], index_col = 0, header = 0)
        hmc = pd.read_table(sys.argv[3], index_col = 0, header = 0)
    else:
	hmp = pd.read_table('HapMap.hmp.txt', index_col = 0, header = 0)
	hmc = pd.read_table('HapMap.hmc.txt', index_col = 0, header = 0)
	
	# only use the actual samples
	hmp_trimmed = hmp.ix[:, 'FLFL04':'WWA30']
	# drop all the columns with more than 90 (default) per cent 'N's (not yet needed)
	# hmp_trimmed = drop_N_columns(hmp_trimmed)
	
	# filter according to read depth count in hmc (default threshold = 4)
	#filter_results = []
	#for i, col in enumerate(hmp_trimmed.columns):
	#	base_values = hmp_trimmed[col]
	#	count_values = hmc[hmc.columns[i]]
	#	filter_results.append(filter_single_col(base_values, count_values))
	
	# df = pd.DataFrame(zip(*filter_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)
	
	# transform ambiguous iupac codes to unambiguous nucleotides
	unambiguous_results = []
	for i, col in enumerate(hmp_trimmed.columns):
		base_list = hmp_trimmed[col]
		count_list = hmc[hmc.columns[i]]
		unambiguous_results.append(get_loci_from_iupac_codes(base_list, count_list, hmp.alleles))
	
	data_unambiguous = pd.DataFrame(zip(*unambiguous_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)
	
	# transform dataset into genepop format (first step)
	numeric_alleles = []
	for col in data_unambiguous.columns:
	    allele_list = data_unambiguous[col]
	    numeric_alleles.append(get_genepop_codes(allele_list))
	
	data_numeric = pd.DataFrame(zip(*numeric_alleles), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)

	# sort the data according to abundance of bases (as opposed to 'N's)
	data_sorted4 = sort_loci_pdDF(data_unambiguous)
	
	# the alleles column gets to be inserted
	data_sorted4.insert(0, 'alleles', hmp.ix[:, 'alleles'])
	# and the additional information as well
	df2 = hmp.ix[:, 'chrom': 'QCcode']
	df = data_sorted4.join(df2)
	
	# some plotting to see the distributions
	# write the output to a csv file
	df.to_csv("data_sorted4.csv")

