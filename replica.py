# Copyright (c) <2013>, <Martin P Schilling>
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

####
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


#################################################################
# MISMATCH / READ-DEPTH ANALYSIS
#################################################################
data_hmp = hmp_trimmed.copy()
data_hmc = hmc[data_hmp]

alleles_zero_results = []
for i, col in enumerate(data_hmp.columns):
    base_list = data_hmp[col]
    count_list = data_hmc[data_hmc.columns[i]]
    alleles_zero_results.append(get_alleles_zero(base_list, count_list, hmp.alleles)) # instead of alleles_drop, hmp.alleles to get all loci

data_zero = pd.DataFrame(zip(*alleles_zero_results), index = data_hmp.index, columns=data_hmp.columns, dtype = np.str)
####

data = data_zero.copy()

match_list_1_total, match_list_2_total, match_homozyg_total, mismatch_list_1_total, mismatch_list_2_total, mismatch_homozyg_total, mismatch_type_total = ([], [], [], [], [], [], [])
for pop in data.columns:
    if re.findall('([A-Z]+[0-9]+[a-z]+[0-9]+)', pop):
	continue
    elif re.findall('([A-Z]+[0-9]+[a-z]+)', pop):
	replica_1 = re.findall('([A-Z]+[0-9]+[a-z]+)', pop)[0]
	replica_2 = re.findall('([A-Z]+[0-9]+)', pop)[0]
	match_list_1, match_list_2, match_homozyg, mismatch_list_1, mismatch_list_2, mismatch_homozyg, mismatch_type = get_replica_read_depth_minmax(data, data_hmc, replica_1, replica_2)
	match_list_1_total.append(match_list_1)
	match_list_2_total.append(match_list_2)
	match_homozyg_total.append(match_homozyg)
	mismatch_list_1_total.append(mismatch_list_1)
	mismatch_list_2_total.append(mismatch_list_2)
	mismatch_homozyg_total.append(mismatch_homozyg)
	mismatch_type_total.append(mismatch_type)
	

# matches
match_1 = get_pooled_list(match_list_1_total)
match_2 = get_pooled_list(match_list_2_total)
match_homo = get_pooled_list(match_homozyg_total)

match = zip(*[match_1, match_2, match_homo])
df_match = pd.DataFrame(match, columns = ['match_1', 'match_2', 'homozygote'])

# match
select_min_list, select_max_list = ([], [])
for i, val_1 in enumerate(df_match.match_1): 
    val_2 = df_match.match_2[i]
    val_list = [val_1[0], val_1[1], val_2[0], val_2[1]]
    if df_match.homozygote[i] == 1:
	select_min = int(np.min(val_list)/2)
	select_max = int(np.max(val_list)/2)
    else:
	select_min = int(np.min(val_list))
	select_max = int(np.max(val_list))
    select_min_list.append(select_min)
    select_max_list.append(select_max)

df_match['min'] = pd.Series(select_min_list)
df_match['max'] = pd.Series(select_max_list)
# df_match.to_csv{"match.csv")

####
# mismatches
mismatch_1 = get_pooled_list(mismatch_list_1_total)
mismatch_2 = get_pooled_list(mismatch_list_2_total)
mismatch_homo = get_pooled_list(mismatch_homozyg_total)
mismatch_type = get_pooled_list(mismatch_type_total)

mismatch = zip(*[mismatch_1, mismatch_2, mismatch_homo, mismatch_type])
df_mismatch = pd.DataFrame(mismatch, columns = ['mismatch_1', 'mismatch_2', 'homozygote', 'type'])


select_min_list, select_max_list = ([], [])
for i, val_1 in enumerate(df_mismatch.mismatch_1): 
    val_2 = df_mismatch.mismatch_2[i]
    val_list = [val_1[0], val_1[1], val_2[0], val_2[1]]
    if df_mismatch.homozygote[i] == 1:
	val_list[0] = val_list[0]/2
	val_list[1] = val_list[1]/2
	select_min = np.float(np.min(val_list))
	select_max = np.float(np.max(val_list))
    elif df_mismatch.homozygote[i] == 2:
	val_list[2] = val_list[2]/2
	val_list[3] = val_list[3]/2
	select_min = np.float(np.min(val_list))
	select_max = np.float(np.max(val_list))
    else:
	select_min = np.min(val_list)
	select_max = np.max(val_list)
    select_min_list.append(select_min)
    select_max_list.append(select_max)

df_mismatch['min'] = pd.Series(select_min_list)
df_mismatch['max'] = pd.Series(select_max_list)
# df_mismatch.to_csv("mismatch.csv")
