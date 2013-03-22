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

######################################################################
######################################################################

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

############################################
### NOTE: there are 9 3rd replicates !!!
############################################

### course of action: 

# take the post-UNEAK replicates and filter those


 # hmp 
# 1) loop over columns and find the pairs of replicates
# 2) compare them and count the number of True and False
# 3) create a new dataFrame, with True and False counts


# hmc

# 1) loop over columns and find the pairs of replicates
# 2) take count_1 and count_2 of both and compare them.


#######
# first, assign pairs of replicates (with re's?)
# sth like:
def get_hmp_replica_ratio(repl_1, repl_2):
    true = 0
    false = 0
    replica_ratio = []
    compare_repl = repl_1 == repl_2
    for equal in compare_repl:
	if equal:
	    true += 1
	else:
	    false += 1
    replica_ratio.append([true/len(repl_1), false/len(repl_1)])
    return replica_ratio

def get_hmp_replica_summ(repl_1, repl_2):
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



###
replica_list = []
for pop in hmp_trimmed.columns:
    if re.findall('([A-Z]+[0-9]+[a-z]+[0-9]+)', pop):
	continue
    elif re.findall('([A-Z]+[0-9]+[a-z]+)', pop):
	replica_1 = re.findall('([A-Z]+[0-9]+[a-z]+)', pop)[0]
	replica_2 = re.findall('([A-Z]+[0-9]+)', pop)[0]
	replica_list.append(get_hmp_replica_ratio(hmp_trimmed[replica_1], hmp_trimmed[replica_2]))

# now, I omitted the 'replicate2' ...

col_names_replica = []
for pop in hmp_trimmed.columns:
    if re.findall('([A-Z]+[0-9]+[a-z]+[0-9]+)', pop):
	continue
    elif re.findall('([A-Z]+[0-9]+[a-z]+)', pop):
	col_names_replica.append(re.findall('([A-Z]+[0-9]+)', pop)[0])


# 
for val in replica_list:
    a.append(val[0])


repl_ratio = pd.DataFrame(zip(*a), index= ['true', 'false'], columns= col_names_replica)


repl_ratio.columns = pd.Index(col_names_replica) # unhashable type: list (why's that?)



####################################


hmp = pd.read_table('HapMap.hmp.txt', index_col = 0, header = 0)
hmc = pd.read_table('HapMap.hmc.txt', index_col = 0, header = 0)

# only use the actual samples
hmp_trimmed = hmp.ix[:, 'FLFL04':'WWA30']


alleles_MAF_results = []
for i, col in enumerate(hmp_trimmed.columns):
    base_list = hmp_trimmed[col]
    count_list = hmc[hmc.columns[i]]
    alleles_MAF_results.append(get_alleles_MAF(base_list, count_list, hmp.alleles))

data_MAF = pd.DataFrame(zip(*alleles_MAF_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)


all_pooled_values = []
for column in data_MAF.columns:
    for value in data_MAF.ix[:, column]:
        all_pooled_values.append(value)

zyg_types = pd.Series(get_zygosity_types(all_pooled_values))

zyg_types.to_csv("allele_types2_MAF.csv")
###########################################

# adv filter + MAF

adv_fil_results = []
for i, col in enumerate(hmp_trimmed.columns):
    base_list = hmp_trimmed[col]
    count_list = hmc[hmc.columns[i]]
    adv_fil_results.append(get_alleles_adv(base_list, count_list, hmp.alleles))

data_adv = pd.DataFrame(zip(*adv_fil_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)
#
alleles_adv_MAF_results = []
for i, col in enumerate(data_adv.columns):
    base_list = data_adv[col]
    count_list = hmc[hmc.columns[i]]
    alleles_adv_MAF_results.append(get_alleles_MAF(base_list, count_list, hmp.alleles))

data_adv_MAF = pd.DataFrame(zip(*alleles_adv_MAF_results), index = hmp_trimmed.index, columns=hmp_trimmed.columns, dtype = np.str)

adv_MAF_pooled_values = []
for column in data_adv_MAF.columns:
    for value in data_adv_MAF.ix[:, column]:
        adv_MAF_pooled_values.append(value)

zyg_types = pd.Series(get_zygosity_types(adv_MAF_pooled_values))

zyg_types.to_csv("allele_types2_simple_MAF.csv")
