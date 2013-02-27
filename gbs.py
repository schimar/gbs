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

def get_alleles_4base(base_list, count_list, allele_list, threshold = 4):
    '''Returns a list of nucleotides where ambiguity codes have been changed to their respective value based on a list of measured allele levels and the read-depth (default = 4)'''
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


def get_alleles_adv(base_list, count_list, allele_list, threshold = 4):
    '''Returns a list of nucleotides where ambiguity codes have been changed to their respective value based on a list of measured allele levels. A read-depth threshold (default = 4) is applied. If one allele is over the threshold, but not the second one, it will be checked for having at least double the amount of the threshold, to qualify as heterozygote (if not, it'll be '?'  '''
    ambig = []
    value = ''
    for i, base in enumerate(base_list):
        count_1, count_2 = map(int, count_list[i].split('|'))
        allele_1, allele_2 = allele_list[i].split('/')
        if base == 'N':
            value = 'N'
        else:
            if count_1 >= threshold and count_2 < threshold:
                if count_1 <= 2*threshold:
		    value = str(allele_1 + '/' + '?')
		elif count_1 > 2*threshold:
		    value =  str(allele_1 + '/' + allele_1)
	    elif count_1 < threshold and count_2 >= threshold:
		if count_2 <= 2*threshold:
		    value = str(allele_2 + '/' + '?')
		elif count_2 > 2*threshold:
		    value = str(allele_2 + '/' + allele_2)
            elif count_1 >= threshold and count_2 >= threshold:
                value = str(allele_1 + '/' + allele_2)
	    else:
		value = 'N' 
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
	    value = 'N'
	else:
	    if count_1 > 0 and count_2 > 0:
		if count_1 > count_2:
		    if ((count_1 - count_2)/count_1) >= MAF:
			value = str(allele_1 + '/' + allele_2)
		    elif ((count_1 - count_2)/count_1) < MAF:
			value = str(allele_1 + '/' + allele_1)
		elif count_2 > count_1:
		    if ((count_2 - count_1)/count_2) >= MAF:
			value = str(allele_1 + '/' + allele_2)
		    elif ((count_2 - count_1)/count_2) < MAF:
			value = str(allele_2 + '/' + allele_2)
		elif count_1 > 0 and count_2 == 0:
		    value = str(allele_1 + '/' + allele_1) 
		elif count_2 > 0 and count_1 == 0:
		    value = str(allele_2 + '/' + allele_2)
	result.append(value)
    return result 
####

def get_homo_prob(base_list, count_list, allele_list):
    '''Related to the filter functions ('get_alleles_4base', 'get_alleles_adv' and 'get_alleles_MAF'), this func calculates the probability of being homozygous for the given loci from the respective output data frame.'''
    result = []
    value = ''
    for i, base in enumerate(base_list):
        count_1, count_2 = map(int, count_list[i].split('|'))
        allele_1, allele_2 = allele_list[i].split('/')
        if base == 'N':
            value = 'N'
        else:
	    if count_1 > count_2:
		value = (count_1-count_2)/count_1
	    elif count_2 > count_1:
		value = (count_2-count_1)/count_2
	    else:
		value = (count_1-count_2)/count_1
	result.append(value)
    return result

def get_pooled_loci_no_N(data):
    '''Append all the cells that are not "N" to one big list'''
    prob_homo_values = []
    for column in data.columns:
        for value in data.ix[:, column]:
            if value == 'N':
                continue
            else:
                prob_homo_values.append(value)
    return prob_homo_values


def get_zygosity_types(pooled_data):
    '''loops over a pooled list of all the values from the matrix (see "backtrack the count_values") and creates a new list of the same length and assigns zygosity types (0 = 'N', 11 = homo, 12 = hetero and 13 = "N/?"'''
    allele_types = []
    for alleles in all_pooled_values:
	# nucleo = dict([['A', '1'], ['C', '1'], ['G', '1'], ['T', '1'], ['?', '2']])
	if alleles == 'N':
	    value = 0
	elif alleles == '':
	    value = 0
	else:
	    allele_1, allele_2 = alleles.split('/')
	    if allele_1 == allele_2:
		value = 11
	    elif allele_1 == '?' or allele_2 == '?':
		value = 13
	    else:
		value = 12
	    #nucleo.get(allele_1) + nucleo.get(allele_2)
	allele_types.append(value)
    return allele_types


def import_raw_loci(filename):
    '''Retrieve sequencing data from the text file and store it in a numpy array'''
    return np.genfromtxt(filename, dtype=str, delimiter='\t')


def sort_loci_pdDF(data):
    '''Strips the column headers as well as the first column from the input pd.DataFrame and sorts the loci and individuals according to highest abundance of bases'''
    #data_strip_loci = data.ix[:, 'FLFL04': 'WWA30']
    data_not_N = data != 'N'
    row_sums = -1*(np.sum(data_not_N, axis=1))
    col_sums = -1*(np.sum(data_not_N, axis=0))
    row_order = row_sums.argsort()
    col_order = col_sums.argsort()
    indivID_sorted = data.index[row_order]
    #locus_sorted = data.columns[col_order]
    data_sorted = np.array(data)[row_order]#[:, col_order]
    data_sorted = pd.DataFrame(data_sorted, index = indivID_sorted, columns= data.columns)
    return data_sorted


def get_single_counts(count_list):
    '''Take a list or column (in pandas.DataFrame) of form 'no.|no.' and return the respective numbers in two seperate columns, appended in a list'''
    co_list = []
    for count in count_list:
        if count == '0|0':
            continue
        else:
            count_1, count_2 = count.split('|')
        co_list.append([count_1, count_2])
    return co_list


def get_OLD_allele_types(data):
    type_list = []
    for(x,y) in data:
        if x == 0 or y == 0:
            value = 0
        elif x >= 4 and y >= 4:
            if x >= 14 or y >= 14:
                value = 2
            elif x < 14 or y < 14:
                value = 1
        elif x >= 4 and y < 4 or x < 4 and y >= 4:
            if x - y == abs(1) or y - x == abs(1):
                value = 3
            elif x - y == abs(2) or y - x == abs(2):
                value = 4
            elif x - y == abs(3) or y - x == abs(3):
                value = 5
            elif x - y > abs(3) or y - x > abs(3):
                value = 6
        type_list.append(value)
    return type_list

def get_genepop_codes(allele_list):
    '''Transforms the alleles (in the form of e.g. 'A/A') into numeric type (where 01 = A, 02 = C, 03 = G, 04 = T)'''
    output = []
    nucleo = dict([['A', '01'], ['C', '02'], ['G', '03'], ['T', '04'], ['?', '00']])
    for alleles in allele_list:
	if alleles == 'N':
	    value = '0000'
	else:
	    allele_1, allele_2 = alleles.split('/')
	    value = nucleo.get(allele_1) + nucleo.get(allele_2)
	output.append(value)
    return output
    

#    if sys.argv > 1:
#        hmp = pd.read_table(sys.argv[2], index_col = 0, header = 0)
#        hmc = pd.read_table(sys.argv[3], index_col = 0, header = 0)
#    else:
hmp = pd.read_table('HapMap.hmp.txt', index_col = 0, header = 0)
hmc = pd.read_table('HapMap.hmc.txt', index_col = 0, header = 0)

# only use the actual samples
hmp_trimmed = hmp.ix[:, 'FLFL04':'WWA30']
#################################################################
#################################################################
# drop all the columns with more than 90 (default) per cent 'N's (not yet needed)
# hmp_trimmed = drop_N_columns(hmp_trimmed)


#################################################################
##### FILTER 

# there are 3 versions (MAF not yet finished)
# you're basically looping over the column and call the filter functions for each column. Then, you create a pd.dataFrame.
#################################################################


#################################################################
# simple filter with threshold 4 (for different threshold, change it in the get_alleles_4base(..., threshold= <value>)
#################################################################

alleles_4base_results = []
for i, col in enumerate(hmp_trimmed.columns):
    base_list = hmp_trimmed[col]
    count_list = hmc[hmc.columns[i]]
    alleles_4base_results.append(get_alleles_4base(base_list, count_list, hmp.alleles))

data_4base = pd.DataFrame(zip(*alleles_4base_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)
###

#################################################################
# advanced filter with threshold (default = 4) and '?' 
# (where threshold of 2nd allele is < 4, if 1st allele < 2* threshold)
#################################################################

adv_fil_results = []
for i, col in enumerate(hmp_trimmed.columns):
    base_list = hmp_trimmed[col]
    count_list = hmc[hmc.columns[i]]
    adv_fil_results.append(get_alleles_adv(base_list, count_list, hmp.alleles))

data_adv = pd.DataFrame(zip(*adv_fil_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)

#################################################################
# MAF filter 
#################################################################

MAF_results = []
for i, col in enumerate(hmp_trimmed.columns):
    base_list = hmp_trimmed[col]
    count_list = hmc[hmc.columns[i]]
    MAF_results.append(get_alleles_MAF(base_list, count_list, hmp.alleles))

data_MAF = pd.DataFrame(zip(*MAF_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)


#################################################################
#################################################################
# get the PROBABILITY of being HOMOZYGOUS

# (based on the cells of the input dataFrame (here it's the initial (post-UNEAK) dataset), it looks up the values in the hmc count-file and calculates the values (with 0 = 'fully' heterozygous and 1 = 'fully' homozygous))
#################################################################

homo_prob_results = []
for i, col in enumerate(hmp_trimmed.columns):
    base_list = hmp_trimmed[col]
    count_list = hmc[hmc.columns[i]]
    homo_prob_results.append(get_homo_prob(base_list, count_list, hmp.alleles))

data_prob = pd.DataFrame(zip(*homo_prob_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)

prob_homo_values = pd.Series(get_pooled_loci_no_N(data_prob))

# NOTE: from here on, I wrote the pd.Series to a csv, to work in R (to plot & the reshape package has a nice counting function)

#################################################################
# backtrack the count values based on the input filter (not yet finished, see BGS_import_pandas, line 308f.)
#################################################################

# create one big list:

all_pooled_values = []
for column in data.columns:
    for value in data.ix[:, column]:
        all_pooled_values.append(value)

zyg_types = pd.Series(get_zygosity_types(all_pooled_values))






#################################################################
# sort the data according to abundance of bases (as opposed to 'N's)
#################################################################
data_sorted4 = sort_loci_pdDF(data_adv)



#################################################################
# transform dataset into genepop format (first step)
#################################################################

numeric_alleles = []
for col in data_unambiguous.columns:
    allele_list = data_unambiguous[col]
    numeric_alleles.append(get_genepop_codes(allele_list))

data_numeric = pd.DataFrame(zip(*numeric_alleles), index = hmp_trimmed.ix[2:,:].index, columns=hmp_trimmed.ix[2:,:].columns, dtype = np.str)

# NOTE: from here on, the bases are numerically encoded
# you still need to tweak the output file, refer to genepop format for further details

#################################################################
# write the output file:
#################################################################

data.to_csv(<output_file_name>.csv)


#################################################################
# if you want to append the initial columns to the dataFrame (with whichever filter)
#################################################################
# the alleles column gets to be inserted
data_sorted4.insert(0, 'alleles', hmp.ix[:, 'alleles'])
# and the additional information as well
df2 = hmp.ix[:, 'chrom': 'QCcode']
df = data_sorted4.join(df2)
# some plotting to see the distributions
# write the output to a csv file
df.to_csv("data_sorted4.csv")

