"""" Very important: Without the "index_col=0", to get the row names to act as the index you need to remove the word "locus" (or whatever the column 1 name is) and move it to the second row all on its own. So then the real data starts on row 3. By including the index at import you can take a regular spreadsheet"""
import numpy as np
import pandas as pd
import scipy as sp

########################################################################

########################################################################
########################################################################
from gbs_cli import *

hmp = pd.read_table('HapMap.hmp.txt', index_col = 0, header = 0)
hmc = pd.read_table('HapMap.hmc.txt', index_col = 0, header = 0)

# only use the actual samples
hmp_trimmed = hmp.ix[:, 'FLFL04':'WWA30']
# drop all the columns with more than 90 (default) per cent 'N's
hmp_trimmed = drop_N_columns(hmp_trimmed)

# filter according to read depth count in hmc (default threshold = 4)
filter_results = []
for i, col in enumerate(hmp_trimmed.columns):
    base_values = hmp_trimmed[col]
    count_values = hmc[hmc.columns[i]]
    filter_results.append(filter_single_col(base_values, count_values))

df = pd.DataFrame(zip(*filter_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)

# transform ambiguous iupac codes to unambiguous nucleotides and (2) alleles 
unambiguous_results = []
for i, col in enumerate(hmp_trimmed.columns):
    base_list = hmp_trimmed[col]
    count_list = hmc[hmc.columns[i]]
    unambiguous_results.append(get_loci_from_iupac_codes(base_list, count_list, hmp.alleles))

data_unambiguous = pd.DataFrame(zip(*unambiguous_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)

# alternatively:
# data = np.array(zip(*results))

# sort the data according to abundance of bases (as opposed to 'N's)
data_sorted4 = sort_loci_pdDF(data_unambiguous)

# get_genepop_codes:

numeric_alleles = []
for col in data_unambiguous.columns:
    allele_list = data_unambiguous[col]
    numeric_alleles.append(get_genepop_codes(allele_list))



##### do we need this?
# the alleles column gets to be inserted
data_sorted4.insert(0, 'alleles', hmp.ix[:, 'alleles'])
# and the additional information as well
df2 = hmp.ix[:, 'chrom': 'QCcode']
df = data_sorted4.join(df2)
# some plotting to see the distributions
# write the output to a csv file
df.to_csv("data_sorted4.csv")

########################################################################
################################# later...
first_allele = []
second_allele = []
for allele in hmp.ix[:, 'alleles']:
    first_allele.append(allele.split('/')[0])
    second_allele.append(allele.split('/')[1])

hmp = hmp.ix[:, 'FLFL04':'WWA30']
hmp.insert(0, 'allele_1', first_allele)
hmp.insert(1, 'allele_2', second_allele)

#######

#!head -n 10 hmc_play


#####################################################
###################### test w/ single columns
a = hmc.ix[:,1].copy()
b = hmp.ix[:,3].copy()[:30]
a[0] = '5|0'
a[1] = '0|7'
b[0] = 'T'
b[1] = 'C'

############# simplified version (w/ just one series/list) working, now extend the whole thing to the whole pie

single = []
for i, base in enumerate(b):
    count_1, count_2 = a[i].split('|')
    if base == 'N':
        value = 'N'
    else:
        if int(count_1) < 4 and int(count_2) < 4:
            value = 'N'
        elif int(count_1) < 4 and int(count_2) >= 4:
            value = base # note: check allele_1 and allele_2 and ambiguity codes!
        elif int(count_1) >= 4 and int(count_2) < 4:
            value = base
        else:
            value = base # maybe here: check for ambiguity and potential alleles
    single.append(value)

########################################################################
df = pd.DataFrame([['A', 'A', 'C', 'N'], ['N', 'T', 'T', 'N'], ['N', 'T', 'A', 'N'], ['N', 'G', 'T', 'N'], ['N', 'A', 'C', 'N'], ['N', 'A', 'T', 'N'], ['N', 'C', 'G', 'N']])


def drop_N_columns(data, drop_level = 0.9):
    data_dropped = data.copy()
    for i, col in enumerate(data_dropped.columns):
        base_series = data_dropped[col]
        N_to_length = np.sum(base_series == 'N')/len(base_series)
        if N_to_length > drop_level:
            data_dropped = data_dropped.drop(col, axis=1)
    return data_dropped



def drop_N_columns(base_list, drop_level = 0.9):
    for base in base_list:
        sum_none = 0
        if base == 'N':
            sum_none += 1
            if sum_none/len(base_list) < drop_level:
                return True
            else:
                return False


data_dropped = total_data2.drop(['open', 'close'], axis=1)


########################################################################
# ambiguity codes
# - if ambiguity code doesn't make sense, flag that position (how?)
# - so compare allele_1 & allele_2 with actual value



iupac = pd.DataFrame([['N', 'G', 'A', 'T', 'C'], ['G', 'G',None,None,None], ['A',None ,'A',None,None], ['T',None ,None,'T',None], ['C',None ,None,None,'C'], ['B', 'G', None, 'T', 'C'], ['D', 'G', 'A', 'T',None], ['H',None ,'A', 'T', 'C'], ['K', 'G',None,'T',None], ['M',None ,'A', None,'C'], ['R', 'G', 'A',None,None], ['S', 'G',None,None,'C'],['V', 'G', 'A',None,'C'], ['W',None ,'A', 'T',None], ['Y',None ,None,'T', 'C']], columns = ('amb_code', 'G', 'A', 'T', 'C'), dtype = np.str)


########################################################################
# ambiguity codes
# - if ambiguity code doesn't make sense, flag that position (how?)
# - so compare allele_1 & allele_2 with actual value
labels = ['PT1', 'PT2', 'PT3', 'PT4', 'PT5', 'PT6', 'PT7', 'PT8', 'PT9', 'PT10', 'PT11', 'PT12', 'PT13', 'PT14', 'PT15']
c = pd.Series(['A', 'T', 'C', 'G', 'W', 'B', 'D', 'H', 'K', 'M', 'N', 'R', 'S', 'V', 'Y'], index = labels)
# d = a[:30]

# al_1 = hmp.allele_1[:15]
# al_2 = hmp.allele_2[:15]
e = hmp.alleles[:30].copy() # old hmp with 'x/y' values in .alleles
###
# c = df.FLFL23[:30].copy()
d = hmc.FLFL23[:30].copy()

# target:             5                       10                        15
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


########################################################################

co = hmc.ix[:,'FLFL04': 'WWA30']

# loop over all the columns and get the values in one big list
co_values = []
for column in co.columns:
    for value in co.ix[:, column]:
        co_values.append(value)

# now call the func with the formerly created list
co_val = np.array(get_single_counts(co_values), dtype = np.int)
####

sub_val = co_val[:30].copy()
sub_val[22] = [4,3]
sub_val[23] = [3,4]
sub_val[24] = [5,3]
sub_val[0] = [3,6]
sub_val[1] = [6,3]
sub_val[2] = [3, 7]
sub_val[4] = [7,3]
sub_val[5] = [14, 22]
####


def get_allele_types(data):
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
####
allele_types = get_allele_types(co_val), dtype = np.int
co_val = pd.DataFrame(co_val)
co_val.insert(2, 'allele_type', allele_types)


###
plt.hist(type_list, bins = 7)
plt.axis([0, 6, 0, 130000])
plt.grid(b=None, which='major', axis='both')

####
plt.scatter(co_val[:,0], co_val[:,1])

# subset with just the ones
sub_co_val = co_val[co_val.allele_type > 3].copy()

########################################################################
###################################################################
###################################################################

# convert pd.DataFrame to Genepop

data = pd.read_csv('data_sorted4.csv', header=0, index_col=0)

gbs = data.ix[:, 'MI17':'KFO4replicate']

###

# first loop over those and create the genotypes
#  01:A 02:C 03:G 04:T    (as strings)

data = data_unambiguous.ix[:30,:].copy()
a = data.MI17.copy()
b = data.ix[:,140].copy()

data = pd.read_csv("subset_unambiguous_4.csv", header = 0, index_col = 0)

##################################################################
### get_genepop_codes:

numeric_alleles = []
for col in data.columns:
    allele_list = data[col]
    numeric_alleles.append(get_genepop_codes(allele_list))
# I could possibly go on, without the zip(*_) 

#hmpp = hmp_trimmed.ix[:30,:]

data_numeric = pd.DataFrame(numeric_alleles, index = hmp.columns, columns= hmp.index)

###
pops = []
for pop in data.columns:
    pops.append(re.findall('([A-Z]+)', pop))

np.unique(np.array(pops))
###

###
prev_pop = ''
index = []
for i, row in enumerate(data_numeric.index):
    population = re.findall('([A-Z]+)', row)
    if population == prev_pop:
	index.append(row)
    else:
	index.append('POP')
	index.append(row)
    prev_pop = population

output = pd.DataFrame(index = index, columns = data_numeric.columns, dtype=str)
##
df1.join(df2, how='outer')

concat(output, data_numeric, axis=0, join='outer', join_axes=None, ignore_index=False)


# looping over rows seems to be inefficient, so maybe zip(* back and loop over columns ?
data_numeric.xs('FLFL04')



# if no new population:
data_numeric.ix[i, :]

# if new population:
insert line w/ POP





################################################
####### now with columns...
data = pd.DataFrame(zip(*numeric_alleles), index = data_sorted4.index, columns=data_sorted4.columns, dtype = np.str)
###
data = data.ix[:3000, :].copy()

###
prev_pop = ''
header = []
for i, pop in enumerate(data.columns):
    population = re.findall('([A-Z]+)', pop)
    if population == prev_pop:
	header.append(pop)
    else:
	header.append('POP')
	header.append(pop)
    prev_pop = population

###
empty_list = []
for val in range(len(data)):
    empty_list.append([])

new = []
for i, sample in enumerate(header):
    if sample == 'POP':
	new.append(empty_list)
    else:
	new.append(list(data[sample]))


loci = str(' ')
for i in range(len(data.index)+1)[1:]:
    loci += ('Loc' + str(i) + ', ')

zipped = zip(*new)
header.insert(0, '')
zipped.insert(0, header)
new_unzipped = zip(*zipped)

new.insert(0, loci)

###
# how to insert that shit without problems at writing...

w = csv.writer(open('output.gen','w'))

w.writerow(new_unzipped)

outfile = open('output.gen', 'w')

writer = csv.writer(outfile)


for row in new_unzipped:
    writer.writerow(row)
outfile.close()

# the rows & cols aren't right...

arr = pd.DataFrame(new, index = header, columns = data.index)



###
