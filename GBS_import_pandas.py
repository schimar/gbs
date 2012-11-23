"""" Very important: Without the "index_col=0", to get the row names to act as the index you need to remove the word "locus" (or whatever the column 1 name is) and move it to the second row all on its own. So then the real data starts on row 3. By including the index at import you can take a regular spreadsheet"""
import numpy as np
import pandas as pd
import scipy as sp
#data = pandas.read_csv("play.csv",index_col=0)
# file path if needed: /Users/paulwolf/Documents/Manuscripts_and_Projects/Aspen_GBS/UNEAK_8_Aug _12/aspen_hapmap/
#output_file = open("locus_dist_plot", 'w')
#df=pandas.DataFrame(data)
#not_N = df.values != 'N'
#col_sums = scipy.sum(not_N, axis=0)
#row_sums = scipy.sum(not_N, axis=1)



#for row in df.values:
#    inds_per_locus = 0
#    for ind in locus:
#        if ind != 'N':
#            inds_per_locus += 1
#    new_row = str(row) + "\n"

 #   output_file.write(new_row)
#output_file.close()

#You can see a particular column by indexing its name:
#print df.p8
#print""
#You can see rows by indexing the range. Note that you need to use numpy to view the array. Just the index works but it does not show under a print statement.
#print "Four rows of DataFrame:", df['TP2':'TP5']
#row_set=np.array(df['TP2':'TP5'])
#print""
#print "Now converted to np.array:"
#print row_set
#start=data[:6]
#print start
#df = DataFrame(data)
#print data['TP2', 'TP3']
########################################################################

# hmp = pd.read_csv("hmp_play.txt", index_col = 0, header = 1)
# hmc = pd.read_csv("hmc_play", index_col = 0, header = 1)

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

#######

#!head -n 10 hmc_play
#######
# not_N = hmp.values != 'N'
# print not_N[i]


#not_N = df.values != 'N'

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
######################
######################
def filter_single_row(base_list, count_list, threshold = 4):
    '''Returns a list of nucleotides filtered (threshold, default = 4) by number of occurence of a sequencing run at specific loci in a list of bases'''
    output = []
    for i, base in enumerate(base_list):
        count_1, count_2 = count_list[i].split('|')
        if base == 'N':
            value = 'N'
        else:
            if int(count_1) < threshold and int(count_2) < threshold:
                value = 'N'
            elif int(count_1) < threshold and int(count_2) >= threshold:
                value = base # note: check allele_1 and allele_2 and ambiguity codes!
            elif int(count_1) >= threshold and int(count_2) < threshold:
                value = base
            else:
                value = base # maybe here: check for ambiguity and potential alleles
        output.append(value)
    return output
#######

hmc.columns[:152] # all the meaningful columns

def get_count_values(count_input):
    '''Returns a list of values from a pandas.DataFrame with the respective number of occurences of a sequencing run at specific loci'''
    count_col_list = []
    for count_col in count_input[count_input.columns[:152]]:
        for count in count_input[count_col]:
            count_col_list.append(count)
    return count_col_list
#######

def get_base_values(base_input):
    '''Returns a list of values from the pandas.DataFrame consisting of nucleotide positions at specific loci'''
    base_col_list = []
    for base_col in base_input[base_input.columns[2:]]:
        for base in base_input[base_col]:
            base_col_list.append(base)
    return base_col_list
#######

base_values = get_base_values(hmp.ix[:30,:])
count_values = get_count_values(hmc)

base_filter = filter_single_row(base_values, count_values)

###

for col in base_filter:



# ambiguity codes

# unused columns have to be cut out in order to just process the real stuff
# probably append again after this is finished (??)
# ----> maybe use df.filter to restrict to columns based on a certain (regex ??) pattern



##############################################
## def my_awesome_gbs_function(hmp., hmc, threshold = 4): # threshold 4 per default

# 1) read in both input files (hmc and hmp)

# 2) a dataFrame similar to hmp (input !!!)
df = pd.DataFrame({}, index = hmp.index[:30], columns=hmp.columns, dtype = np.str)

# 3) probably extract columns allele_1 and allele_2 into seperate list/series ?

# 4) ...the loop...

# 5) insert N's and bases into new df (in the loop or outside?

## --->> probably at the start, specify the potential bases allele_1 & allele_2

