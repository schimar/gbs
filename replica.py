from gbs_mod import *


#################################################################
# REPLICATE FUNC for initial data  -   hmp
#################################################################

replica_list = []
for pop in hmp_trimmed.columns:
    if re.findall('([A-Z]+[0-9]+[a-z]+[0-9]+)', pop):
	continue
    elif re.findall('([A-Z]+[0-9]+[a-z]+)', pop):
	replica_1 = re.findall('([A-Z]+[0-9]+[a-z]+)', pop)[0]
	replica_2 = re.findall('([A-Z]+[0-9]+)', pop)[0]
	replica_list.append(get_hmp_replica_summ(hmp_trimmed[replica_1], hmp_trimmed[replica_2]))


####
# get the sample names 

col_names_replica = []
for pop in hmp_trimmed.columns:
    if re.findall('([A-Z]+[0-9]+[a-z]+[0-9]+)', pop):
	continue
    elif re.findall('([A-Z]+[0-9]+[a-z]+)', pop):
	col_names_replica.append(re.findall('([A-Z]+[0-9]+)', pop)[0])

index = ['equal_not_N', 'not_equal', 'ambig', 'total_not_N', 'NN', 'Nn']
####
a = []
for val in replica_list:
    a.append(val[0])

##
repl_summ = pd.DataFrame(zip(*a), index= index, columns= col_names_replica)
##
repl_summ.to_csv("adv_replicate_summary.csv")

#################################################################
# REPLICATE FUNC for filtered data
#################################################################
# here with the simple 4base filter; just change data...

data = data_4base.copy()

replica_list = []
for pop in data.columns:
    if re.findall('([A-Z]+[0-9]+[a-z]+[0-9]+)', pop):
	continue
    elif re.findall('([A-Z]+[0-9]+[a-z]+)', pop):
	replica_1 = re.findall('([A-Z]+[0-9]+[a-z]+)', pop)[0]
	replica_2 = re.findall('([A-Z]+[0-9]+)', pop)[0]
	replica_list.append(get_filter_replica_summ(data[replica_1], data[replica_2]))

####
# get the sample names 

col_names_replica = []
for pop in data.columns:
    if re.findall('([A-Z]+[0-9]+[a-z]+[0-9]+)', pop):
	continue
    elif re.findall('([A-Z]+[0-9]+[a-z]+)', pop):
	col_names_replica.append(re.findall('([A-Z]+[0-9]+)', pop)[0])

index = ['equal_not_N', 'single_mism', 'double_mism', 'ambig', 'total_not_N', 'NN', 'Nn']

####
proper = []
for val in replica_list:
    proper.append(val[0])


repl_summ_data = pd.DataFrame(zip(*proper), index= index, columns= col_names_replica)

repl_summ_data.to_csv("MAF_replicate_summary.csv")





