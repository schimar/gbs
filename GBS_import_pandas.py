"""" Very important: Without the "index_col=0", to get the row names to act as the index you need to remove the word "locus" (or whatever the column 1 name is) and move it to the second row all on its own. So then the real data starts on row 3. By including the index at import you can take a regular spreadsheet"""
import numpy as np
import pandas as pd
import scipy as sp
data = pandas.read_csv("play.csv",index_col=0)
# file path if needed: /Users/paulwolf/Documents/Manuscripts_and_Projects/Aspen_GBS/UNEAK_8_Aug _12/aspen_hapmap/
output_file = open("locus_dist_plot", 'w')
df=pandas.DataFrame(data)
not_N = df.values != 'N'
col_sums = scipy.sum(not_N, axis=0)
row_sums = scipy.sum(not_N, axis=1)



#for row in df.values:
#    inds_per_locus = 0
#    for ind in locus:
#        if ind != 'N':
#            inds_per_locus += 1
    new_row = str(row) + "\n"
    output_file.write(new_row)
output_file.close()

#You can see a particular column by indexing its name:
#print df.p8
#print""
#You can see rows by indexing the range. Note that you need to use numpy to view the array. Just the index works but it does not show under a print statement.
#print "Four rows of DataFrame:", df['TP2':'TP5']
row_set=np.array(df['TP2':'TP5'])
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
hmc = pd.read_table('hmc_play', index_col = 0, header = 0)
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

!head -n 10 hmc_play

### hmp ###
for i in range(2,53):
    for val in hmp.ix[:, i]:
        if val == 'N':
            continue
        if val != 'N':
            print val


not_N = hmp.values != 'N'
print not_N[i]




not_N = df.values != 'N'


for val in hmp.ix[:, 10]:
    if val == 'N':
        print 'juhu'   # how to jump to the next value in the list? or do nothing??
        continue
    if val != 'N':
        print val


for name in hmp.columns:
    print name

### hmc ###

### ----->>>> first column is probably shifted (??)

for val in hmc.ix[:, 1]:
    print val.split('|')

###################### test w/ single columns
a = hmc.ix[:,1].copy()
b = hmp.ix[:,3].copy()[:30]

single = []
for base in b:
    counts = (val.split('|') for val in a)
    for j, [count_1, count_2] in enumerate(counts):
        if base == 'N':
            value = 'N'  # or write 'N' to the new df (hmp)!!
        if base != 'N':
            if int(count_1) < 4: # drop the base in new df
                value = 'N'
            if int(count_2) < 4:
                value = 'N'
            if int(count_1) >= 4: # write the according base into the new df
                value = base
            if int(count_2) >= 4:
                value = base
    single.append(value)

# --->>> fix the logic:
#                       if value gets to be N when count_1 is below 0, then it is being overwritten with the secong allele...


# --->>> if both count_1 and count_2 are >= 4    (!!!)
# ambiguity codes

# unused columns have to be cut out in order to just process the real stuff
# probably append again after this is finished (??)



##############################################
## def my_awesome_gbs_function(hmp., hmc, threshold = 4): # threshold 4 per default

# 1) read in both input files (hmc and hmp)

# 2) a dataFrame similar to hmp (input !!!)
df = pd.DataFrame({}, index = hmp.index[:30], columns=hmp.columns, dtype = np.str)

# 3) probably extract columns allele_1 and allele_2 into seperate list/series ?

# 4) ...the loop...

# 5) insert N's and bases into new df (in the loop or outside?

## --->> probably at the start, specify the potential

