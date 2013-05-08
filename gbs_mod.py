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

import Hardy_Weinberg_Equilibrium_exact_test_user_Kantale as hwe

###
def drop_N_individuals(data, drop_level = 0.9):
    '''Returns a pd.DataFrame, where all columns, which consist of more than 90 per cent (default) 'N's are being dropped'''
    data_dropped = data.copy()
    for i, col in enumerate(data_dropped.columns):
	base_series = data_dropped[col]
        N_ratio = np.sum(base_series == 'N')/len(base_series)
	if N_ratio > drop_level:
	    del data_dropped[col]
    return data_dropped

def drop_N_loci(hmp, hmc, drop_level = 0.9):
    '''Returns a pd.DataFrame, where all rows (i.e. loci), which consist of more than 90 per cent (default) 'N's are being dropped. hmp and hmc have to be specified, in order to keep the same dimensions. Additionally, a drop_list is being returned with the rows to be discarded.'''
    drop_list = []
    for i, locus in enumerate(hmp.index):
	base_series = hmp.xs(locus)
	N_ratio = np.sum(base_series == 'N')/len(base_series)
	if N_ratio > drop_level:
	    drop_list.append(locus)
    df_hmp = hmp.drop(drop_list)
    df_hmc = hmc.drop(drop_list)
    return df_hmp, df_hmc, drop_list

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


def get_alleles_zero(base_list, count_list, allele_list, allele_sep='/', NA= 'N'):
    '''Returns a list of nucleotides where ambiguity codes have been changed to their respective value based on a list of measured allele levels. Seperator between alleles (default= '/') can be specified with allele_sep= '<seperator>'. Missing values can be specified with NA= '<NA>' '''
    result = []
    value = ''
    for i, base in enumerate(base_list):
	count_1, count_2 = map(int, count_list[i].split('|'))
	allele_1, allele_2 = allele_list[i].split('/')
	if base == NA:
	    value = NA
	else:
	    if count_1 > 0 and count_2 == 0:
		value = str(allele_1 + allele_sep + allele_1)
	    elif count_1 == 0 and count_2 > 0:
		value = str(allele_2 + allele_sep + allele_2)
	    elif count_1 > 0 and count_2 > 0:
		value = str(allele_1 + allele_sep + allele_2)
	    else:
		value = NA
	result.append(value)
    return result  

def get_alleles_4base(base_list, count_list, allele_list, threshold = 4, allele_sep= '/', NA= 'N'):
    '''Returns a list of nucleotides where ambiguity codes have been changed to their respective value based on a list of measured allele levels and the read-depth (default = 4). Seperator between alleles (default= '/') can be specified with allele_sep= '<seperator>'. Missing values can be specified with NA= '<NA>' '''
    ambig = []
    value = ''
    for i, base in enumerate(base_list):
        count_1, count_2 = map(int, count_list[i].split('|'))
        allele_1, allele_2 = allele_list[i].split('/')
        if base == NA:
            value = NA
        else:
            if count_1 >= threshold and count_2 < threshold:
                value =  str(allele_1 + allele_sep + allele_1)
            elif count_1 < threshold and count_2 >= threshold:
                value = str(allele_2 + allele_sep + allele_2)
            elif count_1 >= threshold and count_2 >= threshold:
                value = str(allele_1 + allele_sep + allele_2)
	    else:
		value = NA
        ambig.append(value)
    return ambig


def get_alleles_adv(base_list, count_list, allele_list, threshold = 4, allele_sep= '/', NA= 'N'):
    '''Returns a list of nucleotides where ambiguity codes have been changed to their respective value based on a list of measured allele levels. A read-depth threshold (default = 4) is applied. If one allele is over the threshold, but not the second one, it will be checked for having at least double the amount of the threshold, to qualify as heterozygote (if not, it'll be '?'. Seperator between alleles (default= '/') can be specified with allele_sep= '<seperator>'. Missing values can be specified with NA= '<NA>' '''
    ambig = []
    value = ''
    for i, base in enumerate(base_list):
        count_1, count_2 = map(int, count_list[i].split('|'))
        allele_1, allele_2 = allele_list[i].split('/')
        if base == NA:
            value = NA
        else:
            if count_1 >= threshold and count_2 < threshold:
                if count_1 < 2*threshold:
		    value = str(allele_1 + allele_sep + '?')
		else: 
		    value =  str(allele_1 + allele_sep + allele_1)
	    elif count_1 < threshold and count_2 >= threshold:
		if count_2 < 2*threshold:
		    value = str(allele_2 + allele_sep + '?')
		else:
		    value = str(allele_2 + allele_sep + allele_2)
            elif count_1 >= threshold and count_2 >= threshold:
                value = str(allele_1 + allele_sep + allele_2)
	    else:
		value = NA
	ambig.append(value)
    return ambig

# here: MAF filter

def get_alleles_MAF(base_list, count_list, allele_list, MAF = 0.45, allele_sep= '/', NA= 'N'):
    '''Returns a list of nucleotides where ambiguity codes have been changed to their respective value based on a list of measured allele levels. The MAF will be calculated and, if smaller than the given MAF (default = 0.45), it will be defined as homozygote, otherwise as heterozygous. Seperator between alleles (default= '/') can be specified with allele_sep= '<seperator>'. Missing values can be specified with NA= '<NA>' '''
    result = []
    value = ''
    for i, base in enumerate(base_list):
        count_1, count_2 = map(int, count_list[i].split('|'))
        allele_1, allele_2 = allele_list[i].split('/')
	if base == NA:
	    value = NA
	else:
	    if count_1 > 0 and count_2 > 0:
		if count_1 > count_2:
		    if count_2/ (count_1 + count_2) >= MAF:
			value = str(allele_1 + allele_sep + allele_2)
		    else: 
			value = str(allele_1 + allele_sep + allele_1)
		elif count_2 > count_1:
		    if count_1/ (count_1 + count_2) >= MAF:
			value = str(allele_1 + allele_sep + allele_2)
		    else:
			value = str(allele_2 + allele_sep + allele_2)
		else: 
		    value = str(allele_1 + allele_sep + allele_2)
	    elif count_1 > 0 and count_2 == 0:
		value = str(allele_1 + allele_sep + allele_1) 
	    elif count_2 > 0 and count_1 == 0:
		value = str(allele_2 + allele_sep + allele_2)
	    else:
		value = NA
	result.append(value)
    return result 


####

def sort_loci_pdDF(data, NA= 'N'):
    '''Strips the column headers as well as the first column from the input pd.DataFrame and sorts the loci and individuals according to highest abundance of bases'''
    #data_strip_loci = data.ix[:, 'FLFL04': 'WWA30']
    data_not_N = data != NA
    row_sums = -1*(np.sum(data_not_N, axis=1))
    col_sums = -1*(np.sum(data_not_N, axis=0))
    row_order = row_sums.argsort()
    col_order = col_sums.argsort()
    indivID_sorted = data.index[row_order]
    #locus_sorted = data.columns[col_order]
    data_sorted = np.array(data)[row_order]#[:, col_order]
    data_sorted = pd.DataFrame(data_sorted, index = indivID_sorted, columns= data.columns)
    return data_sorted

####
def get_hwe_exact(locus_pop_subset, pot_alleles, allele_sep = '/', NA = 'N'):
    '''Returns a pd.Series of exact Hardy-Weinberg-tests for a given set of SNPs for one population (calling the function  of Wigginton et al. 2005). A pd.Series of given alleles, found at each locus has to specified (pot_alleles). If a locus has only Ns, 'N' (default allele_sep) will be returned. If there's only one allele to be found at a given locus, 'NA' will be returned (not to be confused with the allele_sep; if allele_sep = 'NA' is specified, then you won't be able to distinguish the monomorphic loci and the ones where all individuals are 'N').'''
    obs_het= 0
    obs_hom1= 0
    obs_hom2= 0
    is_N = locus_pop_subset == NA
    if is_N.all():
	return NA
    else:
	for i, val in enumerate(locus_pop_subset):
	    pot_allele_1, pot_allele_2 = pot_alleles.split(allele_sep)
	    if val == NA:
		continue
	    else:
		pop_allele_1, pop_allele_2 = val.split(allele_sep)
		if pop_allele_1 == '?' or pop_allele_2 == '?':
		    continue
		elif pop_allele_1 != pop_allele_2:
		    obs_het += 1
		elif pop_allele_1 == pot_allele_1:
		    obs_hom1 += 1
		else:
		    obs_hom2 += 1
	if obs_hom2 == 0 and obs_het == 0:
	    return 'NA'
	elif obs_hom1 == 0 and obs_het == 0:
	    return 'NA'
	else:
	    return hwe.Hardy_Weinberg_Equilibrium_exact_test_user_Kantale(obs_het, obs_hom1, obs_hom2)
####
def get_expected_heterozygosity(genotype_array):
    total = 0; p = 0
    p_allele = 'x' # p-allele is the first real (non-N) allele encountered and is used to caclulate p
    for genotype in genotype_array:
        for allele in genotype:
            if p_allele == 'x': # if the first allele not yet encountered
                if (allele != 'N') and (allele != '/') and (allele != '?'):
                    total +=1; p_allele = allele; p+=1
            else: #p-allele now set to an allele in the list
                if allele == p_allele:
                    total +=1; p+=1
                elif (allele != 'N') and (allele != '/') and (allele != '?'): # don't count 'N' or '/' in total
                    total +=1
    if total == 0:
        exp_het = 0
    else:
        exp_het = 2*(float(p)/total)*(1-(float(p)/total))
    return exp_het, total
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

##
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
    for alleles in pooled_data:
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

def get_hmp_replica_summ(repl_1, repl_2):
    '''Gives a summary list for the comparison of 2 replicate samples'''
    NN, Nn, equal_not_N, not_equal, ambig = (0,0,0,0,0)
    replica = []
    N_1 = repl_1 == 'N'
    N_2 = repl_2 == 'N'
    compare_repl = repl_1 == repl_2
    nucleo = dict([['K', ['G','T']], ['M', ['A', 'C']], ['R', ['A', 'G']], ['S', ['G', 'C']], ['W', ['A', 'T']], ['Y', ['T', 'C']]])
    for i, val_1 in enumerate(N_1):
	val_2 = N_2[i]
	if val_1 and val_2:
	    NN += 1
	elif val_1 or val_2:
	    Nn += 1
	else:
	    equal = compare_repl[i]
	    if equal:
		equal_not_N += 1
	    else:
		base_1 = repl_1[i]
		base_2 = repl_2[i]
		if nucleo.get(base_2):
		    if base_1 == nucleo.get(base_2)[0] or nucleo.get(base_2)[1]:
			ambig += 1
		elif nucleo.get(base_1):
		    if base_2 == nucleo.get(base_1)[0] or nucleo.get(base_1)[1]:
			ambig += 1
		else:
		    not_equal += 1
    replica.append([equal_not_N, not_equal, ambig, equal_not_N+not_equal+ambig, NN, Nn])
    return replica

def get_filter_replica_summ(repl_1, repl_2):
    '''Gives a summary list for the comparison of 2 replicate samples of the filtered output'''
    NN, Nn, equal_not_N, not_equal, ambig = (0,0,0,0,0)
    replica = []
    N_1 = repl_1 == 'N'
    N_2 = repl_2 == 'N'
    compare_repl = repl_1 == repl_2
    for i, val_1 in enumerate(N_1):
	val_2 = N_2[i]
	if val_1 and val_2:
	    NN += 1
	elif val_1 or val_2:
	    Nn += 1
	else:
	    rep_1_allele_1, rep_1_allele_2 = repl_1[i].split('/')
	    rep_2_allele_1, rep_2_allele_2 = repl_2[i].split('/')
	    equal = compare_repl[i]
	    if equal:
		equal_not_N += 1
	    elif rep_1_allele_1 == rep_2_allele_1 or rep_1_allele_2 == rep_2_allele_2:
		ambig += 1
	    else:
		not_equal += 1
    replica.append([equal_not_N, not_equal, ambig, equal_not_N+not_equal+ambig, NN, Nn])
    return replica



def import_raw_loci(filename):
    '''Retrieve sequencing data from the text file and store it in a numpy array'''
    return np.genfromtxt(filename, dtype=str, delimiter='\t')




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

def get_structure_format(allele_list, allele_sep= ' ', NA= 'N'):
    '''Transforms the alleles (in the form of e.g. 'A/A') into numeric type (where 01 = A, 02 = C, 03 = G, 04 = T)'''
    output = []
    nucleo = dict([['A', '1'], ['C', '2'], ['G', '3'], ['T', '4'], ['?', '-9']])
    for alleles in allele_list:
	if alleles == NA:
	    value = NA
	else:
	    allele_1, allele_2 = alleles.split(allele_sep)
	    value = nucleo.get(allele_1) + allele_sep + nucleo.get(allele_2)
	output.append(value)
    return output

