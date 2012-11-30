from __future__ import division

import numpy as np
import pandas as pd
import scipy as sp

###############

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

##############

hmp = pd.read_table('hmp_play.txt', index_col = 0, header = 0)
hmc = pd.read_table('hmc_play02.txt', index_col = 0, header = 0)

e = hmp.alleles[:30].copy()

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
    results.append(filter_single_col(base_values, count_values)) # gbs...
###

results.insert(0, list(hmp.allele_1))
results.insert(1, list(hmp.allele_2))


df = pd.DataFrame(zip(*results), index = hmp.index[:30], columns=hmp.columns, dtype = np.str)
#####################################

iupac = pd.DataFrame([['N', 'G', 'A', 'T', 'C'], ['G', 'G',None,None,None], ['A',None ,'A',None,None], ['T',None ,None,'T',None], ['C',None ,None,None,'C'], ['B', 'G', None, 'T', 'C'], ['D', 'G', 'A', 'T',None], ['H',None ,'A', 'T', 'C'], ['K', 'G',None,'T',None], ['M',None ,'A', None,'C'], ['R', 'G', 'A',None,None], ['S', 'G',None,None,'C'],['V', 'G', 'A',None,'C'], ['W',None ,'A', 'T',None], ['Y',None ,None,'T', 'C']], columns = ('amb_code', 'G', 'A', 'T', 'C'), dtype = np.str)

# target:           5                               10                            15
['A', 'T', 'C', 'A', 'A', 'N', 'N', 'N', 'A', 'G', 'A', 'G', 'T', 'C', 'T',
 'G', 'T', 'G', 'T', 'A',      'C', 'C', 'A', 'T', 'G', 'N', 'N', 'A', 'G', 'N']
# amb --- c     for testing.... (TP26 and TP28: amb_codes are wrong...)
c = pd.Series(
['M', 'W', 'H', 'D', 'R', 'N', 'N', 'N', 'W', 'R', 'H', 'V', 'K', 'M', 'Y',
 'A', 'G', 'G', 'T', 'A',      'C', 'N', 'W', 'Y', 'Y', 'T', 'S', 'W', 'R', 'N'], index = hmp.index[:30])

e.ix[12] = 'G/T'
e.ix[15] = 'C/G'
e.ix[16] = 'G/T'

d = pd.Series(    #            5                                 10
['5|1', '0|7', '1|6', '4|0', '7|0', '0|3', '2|2', '1|3', '4|2', '0|10',
 '4|1', '0|5', '0|7', '0|9', '3|11', '1|4', '0|7', '1|4', '2|6', '7|0',
 '5|0', '0|8', '7|1', '0|4', '1|12', '0|1', '1|3', '10|2', '0|4', '1|3'], index = hmp.index[:30])

#####################################################
ambig = []
for i, base in enumerate(c):
    count_1, count_2 = d[i].split('|')
    allele_1, allele_2 = e[i].split('/')
    for j, ambi in enumerate(iupac.ix[5:, 0]): # starting from the 5th row; we don't need 'N's and single nucleotides.
        if base == 'N':
            value = 'N'
        elif base == ambi:
            if count_1 >= 4 and count_2 < 4: # then check with the allele_1 & 2 and, if
                value = allele_1
                break
            elif count_1 < 4 and count_2 >= 4:
                value = allele_2
                break
            elif count_1 >= 4 and count_2 >= 4:
                value = str(allele_1 + '/' + allele_2)
          # elif base == ambi and count_1 >= 4 and count_2 >= 4:
           #     value = 'what now?'
           #     break
        else:
            continue
    ambig.append(value)

