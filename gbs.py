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

from gbs_mod import *
import Hardy_Weinberg_Equilibrium_exact_test_user_Kantale as hwe

###


#    if len(sys.argv) > 1:
#        hmp = pd.read_table(sys.argv[2], index_col = 0, header = 0)
#        hmc = pd.read_table(sys.argv[3], index_col = 0, header = 0)
#    else:
hmp = pd.read_table('HapMap.hmp.txt', index_col = 0, header = 0)
hmc = pd.read_table('HapMap.hmc.txt', index_col = 0, header = 0)

# only use the actual samples
hmp_trimmed = hmp.ix[:, 'FLFL04':'WWA30']

#################################################################
#################################################################
# drop all the loci (rows) with more than 90 (default) per cent 'N's 
# then columns and now filter:

# 1) drop_N_loci (on hmp_trimmed)
# 2) drop_N_individuals (on 1)
# 3) filter (4base, adv, or MAF) 
# 4) drop_N_loci 

hmp_trim_drop1, hmc_drop1, drop_list = drop_N_loci(hmp_trimmed, hmc)

alleles_drop = hmp.alleles.drop(drop_list)

hmp_trim_drop2 = drop_N_individuals(hmp_trim_drop1)

#################################################################
##### FILTER 

# there are 4 versions 
# you're basically always looping over the column and calling the filter functions for each column. Then, you create a pd.dataFrame.
#################################################################

#################################################################
# SIMPLE FILTER for the initial dataset (without a threshold !!)
#################################################################
data_hmp = hmp_trim_drop2.copy()
data_hmc = hmc_drop1.copy()

alleles_zero_results = []
for i, col in enumerate(data_hmp.columns):
    base_list = data_hmp[col]
    count_list = data_hmc[data_hmc.columns[i]]
    alleles_zero_results.append(get_alleles_zero(base_list, count_list, alleles_drop))

data_zero = pd.DataFrame(zip(*alleles_zero_results), index = data_hmp.index, columns=data_hmp.columns, dtype = np.str)

data_zero_drop, hmc_zero_drop, drop_list_zero = drop_N_loci(data_zero, hmc_drop1)
#################################################################
# SIMPLE FILTER with threshold 4 (for different threshold, change it in the get_alleles_4base(..., threshold= <value> (same goes for allele_sep and )
#################################################################
data_hmp = hmp_trim_drop2.copy()
data_hmc = hmc_drop1.copy()

alleles_4base_results = []
for i, col in enumerate(data_hmp.columns):
    base_list = data_hmp[col]
    count_list = data_hmc[data_hmc.columns[i]]
    alleles_4base_results.append(get_alleles_4base(base_list, count_list, alleles_drop))

data_4base = pd.DataFrame(zip(*alleles_4base_results), index = data_hmp.index, columns=data_hmp.columns, dtype = np.str)

data_4base_drop, hmc_4base_drop2, drop_list_4base = drop_N_loci(data_4base, hmc_drop1)
###

#################################################################
# ADVANCED FILTER with threshold (default = 4) and '?' 
# (where threshold of 2nd allele is < 4, if 1st allele < 2* threshold)
#################################################################
data_hmp = hmp_trim_drop2.copy()
data_hmc = hmc_drop1.copy()

adv_fil_results = []
for i, col in enumerate(data_hmp.columns):
    base_list = data_hmp[col]
    count_list = data_hmc[data_hmc.columns[i]]
    adv_fil_results.append(get_alleles_adv(base_list, count_list, alleles_drop))

data_adv = pd.DataFrame(zip(*adv_fil_results), index = data_hmp.index, columns= data_hmp.columns, dtype = np.str)

data_adv_drop, hmc_drop_adv, drop_list_adv = drop_N_loci(data_adv, hmc_drop1)
#################################################################
# MAF filter 
#################################################################
data_hmp = hmp_trim_drop2.copy()
data_hmc = hmc_drop1.copy()

MAF_results = []
for i, col in enumerate(data_hmp.columns):
    base_list = data_hmp[col]
    count_list = data_hmc[data_hmc.columns[i]]
    MAF_results.append(get_alleles_MAF(base_list, count_list, alleles_drop))

data_MAF = pd.DataFrame(zip(*MAF_results), index = data_hmp.index, columns= data_hmp.columns, dtype = np.str)

data_MAF_drop, hmc_drop_MAF, drop_list_MAF = drop_N_loci(data_MAF, hmc_drop1)

# Note: for combining 4base and MAF, change 'hmp_trimmed' to 'data_4base' (you need to do the 4base-filter first) 
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
# backtrack the count values based on the input filter 
#################################################################

# create one big list:

all_pooled_values = []
for column in data_zero.columns:
    for value in data_zero.ix[:, column]:
        all_pooled_values.append(value)

zyg_types = pd.Series(get_zygosity_types(all_pooled_values))

zyg_types.to_csv("allele_types2_4base_zero.csv")
# NOTE: from here on, I wrote the pd.Series to a csv, to work in R (to plot & the reshape package has a nice counting function)

#################################################################
# sort the data according to abundance of bases (as opposed to 'N's)

#################################################################
data_sorted = sort_loci_pdDF(data_4base)


# to get the population names:
pops = []
for val in data_sorted.columns:
    pop = re.findall('([A-Z]+)', val)
    pops.append(pop)


tdata = data_sorted.transpose()
pops = pd.Series(pops, index= tdata.index)
tdata.insert(0, 'population', pops)

tdata.to_csv("4base_MAF_sorted.csv")

#################################################################
# Calculate HWE exact test for each population
#################################################################
# change <data_4base...> and <drop_list...> to respective filter output

data = data_zero_drop.copy() 
##
alleles = hmp.alleles.drop(drop_list_zero)

grouped = data.groupby(lambda x: re.match("[A-Z]+", x).group(), axis=1)
#
results= []
for name, group in grouped:
    hwe_list = []
    for i, locus in enumerate(group.index):
	hwe_list.append(get_hwe_exact(group.xs(locus), alleles[i]))
    results.append(hwe_list)

hwe_df = pd.DataFrame(zip(*results), columns= ['FLFL', 'HSPQ', 'KFO', 'MI', 'SFQ', 'WWA'], index = data.index, dtype= np.float64)

hwe_df.to_csv("hwe_MAF.csv")

# loop into one list (w/ pops in second column)

#################################################################
# get expected heterozygosity
#################################################################
data = data_adv_drop.copy() 
##

grouped = data.groupby(lambda x: re.match("[A-Z]+", x).group(), axis=1)
#
results= []
for name, group in grouped:
    exp_het_list = []
    for i, locus in enumerate(group.index):
	exp_het_list.append(get_expected_heterozygosity(group.xs(locus)))
    results.append(exp_het_list)


hwe_df = pd.DataFrame(zip(*results), columns= ['FLFL', 'HSPQ', 'KFO', 'MI', 'SFQ', 'WWA'], index = data.index, dtype= np.float64)



#################################################################
# transform dataset into genepop format (first step)
#################################################################
data = data_sorted.copy()

genepop_alleles = []
for col in data.columns:
    allele_list = data[col]
    genepop_alleles.append(get_genepop_codes(allele_list))

data_numeric = pd.DataFrame(genepop_alleles, index = data.columns, columns=data.index, dtype = np.str)

# NOTE: from here on, the bases are numerically encoded
# you still need to tweak the output file, refer to genepop format for further details


#################################################################
# write output file:
#################################################################

data_numeric.to_csv(<output_file_name>.csv)


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

