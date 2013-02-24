###

from gbs_cli import *

hmp = pd.read_table('HapMap.hmp.txt', index_col = 0, header = 0)
hmc = pd.read_table('HapMap.hmc.txt', index_col = 0, header = 0)

hmp_trimmed = hmp.ix[:, 'FLFL04':'WWA30']
# hmp = hmp.ix[:30,:]
# hmc = hmc.ix[:30,:]
###
p = hmp.ix[:30,:]
c = hmc.ix[:30,:]

###

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

###

adv_fil_results = []
for i, col in enumerate(hmp_trimmed.columns):
    base_list = hmp_trimmed[col]
    count_list = hmc[hmc.columns[i]]
    adv_fil_results.append(get_alleles_adv(base_list, count_list, hmp.alleles))

data_adv = pd.DataFrame(zip(*adv_fil_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)
##################################
##################################
def get_homo_prob(base_list, count_list, allele_list):
    '''Related to the filter functions ('get_alleles_4base', 'get_alleles_adv' and 'get_alleles_MAF'), this func calculates the probability of being homozygous for the given loci from the respective output data frame or the initial data set.'''
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
	    elif count_1 == count_2:
		value = 0
	result.append(value)
    return result

###

homo_prob_results = []
for i, col in enumerate(hmp_trimmed.columns):
    base_list = hmp_trimmed[col]
    count_list = hmc[hmc.columns[i]]
    homo_prob_results.append(get_homo_prob(base_list, count_list, hmp.alleles))

data_prob = pd.DataFrame(zip(*homo_prob_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)

# pool everything in one list
	
def get_pooled_loci(data):
    '''Append all the cells that are not "N" to one big list'''
    prob_homo_values = []
    for column in data.columns:
        for value in data.ix[:, column]:
            if value == 'N':
                continue
            else:
                prob_homo_values.append(value)
    return prob_homo_values

### 
prob_homo_values = pd.Series(get_pooled_loci(data_prob))


##
plt.hist(prob_homo_values)

####################
# MAF filter
####################

# MAF >= 0.45 = heterozygote
# --->  think about the division by zero !!!!!
def get_alleles_MAF(base_list, count_list, allele_list, threshold = 4, MAF = 0.45):
    '''Returns a list of nucleotides where ambiguity codes have been changed to their respective value based on a list of measured allele levels. A read-depth threshold (default = 4) is applied. If one allele of a locus is over the threshold and, however, not the second one, then the MAF will be calculated and, if smaller than the given MAF (default = 0.45), it will be defined as homozygote.'''
    result = []
    value = ''
    for i, base in enumerate(base_list):
        count_1, count_2 = map(int, count_list[i].split('|'))
        allele_1, allele_2 = allele_list[i].split('/')
        if base == 'N':
            value = 'N'
        else:
            if count_1 >= threshold and count_2 < threshold:
		if count_2/count_1 < MAF:
		    value =  str(allele_1 + '/' + allele_1)
		elif count_2/count_1 >= MAF:
		    value =  str(allele_1 + '/' + '?')
	    elif count_1 < threshold and count_2 >= threshold:
		if count_1/count_2 < MAF:
		    value = str(allele_2 + '/' + allele_2)
		elif count_1/count_2 >= MAF:
		    value = str(allele_2 + '/' + '?')
            elif count_1 >= threshold and count_2 >= threshold:
		if count_1 < count_2:
		    if count_1/count_2 < MAF:
			value = str(allele_2 + '/' + allele_2)
		    elif count_1/count_2 >= MAF:
			value = str(allele_1 + '/' + allele_2)
		elif count_2 < count_1:
		    if count_2/count_1 < MAF:
			value = str(allele_1 + '/' + allele_1)
		    elif count_2/count_1 >= MAF:
			value = str(allele_1 + '/' + allele_2)
	    else:
		value = 'N' 
        result.append(value)
    return result
###

MAF_results = []
for i, col in enumerate(hmp_trimmed.columns):
    base_list = hmp_trimmed[col]
    count_list = hmc[hmc.columns[i]]
    MAF_results.append(get_alleles_MAF(base_list, count_list, hmp.alleles))

data_MAF = pd.DataFrame(zip(*MAF_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)

###


# pool the content into one big list.

unambig_values = []
for column in hmp_trimmed.columns:
    for value in hmp_trimmed.ix[:, column]:
        unambig_values.append(value)
###

# not finished yet!
def get_single_alleles(allele_list):
    '''group allele types (homo, hetero, unclear (e.g. A/?), N)'''
    all_list = []
    for count in count_list:
        if count == 'N':
            
        else:
            count_1, count_2 = count.split('|')
        co_list.append([count_1, count_2])
    return all_list
#######################################################
#######################################################

# REPLICATES TEST

#######################################################

hmc_trim = hmc.ix[:30, 'FLFL04':'WWA30']
hmp_trim = hmp.ix[:30, 'FLFL04':'WWA30']
####

pops = []
for pop in hmp_trim.columns:
    if re.findall('([A-Z]+[0-9]+[a-z]+)', pop):
	pops.append(re.findall('([A-Z]+[0-9]+)', pop))
	pops.append(re.findall('([A-Z]+[0-9]+[a-z]+)', pop))

rep_samples = np.unique(np.array(pops))
###

hmp.columns[rep_samples]

for sample in rep_samples:
    print sample



for i, alleles in hmp_trim.FLFL08:
    replica = hmp_trim.FLFL08replicate[i]
    if alleles == replica:
	print 'yes'
    
a = hmc_trim.FLFL08 == hmc_trim.FLFL08replicate
b = hmc_trim.FLFL08 != hmc_trim.FLFL08replicate

np.sum(a) # gives the number of True's

np.sum(b)/len(hmp_trim)
# now that's for the counts, but we're also interested in the alleles!
