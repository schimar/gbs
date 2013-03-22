from gbs_cli import *


#################################################################
# REPLICATE FUNCs
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



repl_summ = pd.DataFrame(zip(*a), index= index, columns= col_names_replica)

repl_summ.to_csv("replicate_summary.csv")
