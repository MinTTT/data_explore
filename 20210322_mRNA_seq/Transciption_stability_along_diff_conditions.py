from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
from scipy.stats import rankdata
from statsmodels.stats.multitest import fdrcorrection
import pandas as pd
import numpy as np
import os
import sciplot as splt
# %%  Import all data to the memory.
DIR = r"./20210322_mRNA_seq"
growth_df_ps = r"./20210322_mRNA_seq/sample_growth_rate.xlsx"
expression_files_ps = [file.name for file in os.scandir(DIR)
                       if file.is_file() and file.name.split('.')[-1] == 'csv' and file.name.split('_')[0] == '1321']
seq_id_list = ['_'.join(file.split('_')[0:2]) for file in expression_files_ps]
growth_pf = pd.read_excel(growth_df_ps, usecols=range(5))
# info concat
all_pf = {}
strains_dict = {'NH3.23': [], 'NH3.24': []}

for index, expression_pf in enumerate(expression_files_ps):
    seq_id = seq_id_list[index]
    growth_rate = growth_pf[growth_pf['RNA_seq_ID'] == seq_id]['Growth_rate'].values[0]
    medium = growth_pf[growth_pf['RNA_seq_ID'] == seq_id]['Medium'].values[0]
    df = pd.read_csv(os.path.join(DIR, expression_pf))
    sample_id = expression_pf.split('.')[0]
    print(sample_id)
    if int(sample_id.split('_')[1]) <= 5:
        strains_dict['NH3.23'].append(sample_id)
    else:
        strains_dict['NH3.24'].append(sample_id)
    table_length = len(df)
    df['Sample_ID'] = [sample_id] * table_length
    df['Medium'] = [medium] * table_length
    df['Growth_rate'] = [growth_rate] * table_length
    df = df[df.apply(lambda row: int(row.values[0]) != 1430, axis=1)]  # remove lacI in genome
    df['TPM'] = 1e6 * df['counts'] / df['length'] / \
                        np.sum(df['counts'] / df['length'])
    all_pf[sample_id] = df
# compute average of samples

all_group_name = list(set(['_'.join(sample.split('_')[:-1]) for sample in list(all_pf.keys())]))
all_group_name.sort(key=lambda x: int(x.split('_')[-1]))  # rank the sample name
all_sample_name = [samp for samp in list(all_pf.keys())]
# get average coverage_FPKM data:
all_group_data = {}
for group in all_group_name:
    sample_name = [sample for sample in all_sample_name if '_'.join(sample.split('_')[:-1]) == group]
    data1, data2 = all_pf[sample_name[0]], all_pf[sample_name[1]]
    data_avg = pd.DataFrame(data={'sample_1': data1['TPM'],
                                  'sample_2': data2['TPM'],
                                  })
    data_avg['data_mean'] = np.mean(data_avg[['sample_1', 'sample_2']], axis=1)
    data_avg = pd.concat([data1.iloc[:, :6], data_avg], axis=1)
    all_group_data[group] = data_avg


all_group_data = {}
for group in all_group_name:
    sample_name = [sample for sample in all_sample_name if '_'.join(sample.split('_')[:-1]) == group]
    data1, data2 = all_pf[sample_name[0]], all_pf[sample_name[1]]
    data_avg = pd.DataFrame(data={'sample_1': data1['TPM'],
                                  'sample_2': data2['TPM'],
                                  })
    data_avg['data_mean'] = np.mean(data_avg[['sample_1', 'sample_2']], axis=1)
    data_avg = pd.concat([data1.iloc[:, :6], data_avg], axis=1)
    all_group_data[group] = data_avg


expression_level_col_names = ['sample_1', 'sample_2']
data_1, data_2 = all_group_data['1321_4'], all_group_data['1321_9']

data_1_expression = data_1[expression_level_col_names]
data_2_expression = data_2[expression_level_col_names]

gene_list = list(set.union(set(data_1['gene'].tolist()), set(data_2['gene'].tolist())))

mwu_pvalue_list = []
log_fc_list = []
for i, gene in enumerate(gene_list):
    exp_level_1, exp_level_2 = np.squeeze(data_1_expression[data_1['gene'] == gene].values), \
                               np.squeeze(data_2_expression[data_2['gene'] == gene].values)
    log_fc_list.append(np.log(exp_level_1.mean() / exp_level_2.mean()))
    try:
        _, mwu_pvalue = mannwhitneyu(exp_level_1, exp_level_2)
    except ValueError:
        mwu_pvalue = np.nan
    mwu_pvalue_list.append(mwu_pvalue)

mwu_pvalue_list = np.array(mwu_pvalue_list)
nan_filter = ~np.isnan(mwu_pvalue_list)
mwu_pvalue_list = mwu_pvalue_list[nan_filter]

log_fc_list = np.array(log_fc_list)
log_fc_list = log_fc_list[nan_filter]

qvalues = fdrcorrection(mwu_pvalue_list)
filtered_gene_list = np.array(gene_list)[nan_filter]

# volcano plot for GED
fig7, ax7 = plt.subplots(1, 1)  # volcano plot

ax7.scatter(log_fc_list, -np.log10(mwu_pvalue_list))
ax7.scatter(log_fc_list[qvalues[1] < .3], -np.log10(mwu_pvalue_list)[qvalues[1] < .3], color='r')

fig7.show()

sig_gene_list = filtered_gene_list[np.logical_and(qvalues[1] < .21, log_fc_list < np.log(.1))]

# print('\n'.join(filtered_gene_list[ np.logical_and(qvalues[1] < .21, log_fc_list > np.log(2))]))

fig8, ax8 = plt.subplots(1, 1, figsize=(8, 8))

gene_table_filter = [True if gene in sig_gene_list else False for gene in data_1['gene']]

ax8.scatter(data_1['data_mean'], data_2['data_mean'], c='#B2BABB', alpha=.5)
ax8.scatter(data_1[gene_table_filter]['data_mean'], data_2[gene_table_filter]['data_mean'], c='r')
# ax8.scatter(data_1[data_1['gene']=='mcaS']['data_mean'], data_2[data_2['gene']=='mcaS']['data_mean'], color='r')
ax8.set_xscale('log')
ax8.set_yscale('log')
ax8.set_xlim(2e-2, 2e4)
ax8.set_ylim(2e-2, 2e4)
ax8.set_xlabel('NH3.23 (LacI/GFP)')
ax8.set_ylabel('NH3.24 (TetR/RFP)')
ax8.set_xticks([1, 100, 10000])
splt.aspect_ratio(1)
fig8.show()
print('\n'.join(sig_gene_list))
# Plot the transcriptome comparison
fig5, ax5 = plt.subplots(5, 5, figsize=(8 * 5, 8 * 5))

corss_x, corss_y = np.meshgrid(range(5), range(5, 10))
corss_x = corss_x.flatten()
corss_y = corss_y.flatten()
index_x, index_y = np.meshgrid(range(5), range(5))
for en_i, i in enumerate(zip(index_x.flatten(), index_y.flatten())):
    sample_x_index = corss_x[en_i]
    sample_y_index = corss_y[en_i]
    print(sample_x_index, sample_y_index)
    ax5[i[0], i[1]].scatter(all_group_data[all_group_name[sample_x_index]]['data_mean'],
                            all_group_data[all_group_name[sample_y_index]]['data_mean'], c='#B2BABB')
    circuits_filter = all_group_data[all_group_name[i[0]]]['locus_tag'].apply(
        lambda tag: tag.split('_')[0] == 'ts' or tag.split('_')[0] == 'ls')
    ax5[i[0], i[1]].scatter(all_group_data[all_group_name[sample_x_index]][circuits_filter]['data_mean'],
                            all_group_data[all_group_name[sample_y_index]][circuits_filter]['data_mean'], c='#82E0AA')
    ax5[i[0], i[1]].scatter(all_group_data[all_group_name[sample_x_index]][gene_table_filter]['data_mean'],
                            all_group_data[all_group_name[sample_y_index]][gene_table_filter]['data_mean'], c='r')
    ax5[i[0], i[1]].set_xlabel(all_group_name[sample_x_index])
    ax5[i[0], i[1]].set_ylabel(all_group_name[sample_y_index])

    ax5[i[0], i[1]].set_xscale('log')
    ax5[i[0], i[1]].set_yscale('log')
    ax5[i[0], i[1]].set_xlim(2e-2, 2e4)
    ax5[i[0], i[1]].set_ylim(2e-2, 2e4)
    ax5[i[0], i[1]].set_xticks([1, 1E2, 1E4])
    splt.aspect_ratio(1, ax5[i[0], i[1]])
fig5.show()
#%%
mean_exp = []
std_exp = []
cv_exp = []
expression_level = []
for gene in filtered_gene_list:
    exp_level = []

    for key in list(all_group_data.keys()):
        data = all_group_data[key]
        exp_level.append(data[data['gene'] == gene]['data_mean'].values[0])
    expression_level.append(exp_level)
    mean = np.mean(exp_level)
    std = np.std(exp_level)
    mean_exp.append(mean)
    std_exp.append(std)
    cv_exp.append(std / mean)

gene_expression_statics = pd.DataFrame(data={'gene': filtered_gene_list,
                                             'expression': expression_level,
                                             'mean': mean_exp,
                                             'std': std_exp,
                                             'cv': cv_exp})

fig9, ax9 = plt.subplots(1, 1, figsize=(15, 15))

ax9.scatter(gene_expression_statics['mean'], gene_expression_statics['cv'], alpha=.5, c='#B2BABB' )
ax9.set_xscale('log')
ax9.set_yscale('log')
# annotate_list = ['rho', 'cyoA', 'ispA', 'fabB']  # qPCR ref gene candidates
annotate_list = ['ynaE', 'ydgJ', 'yiiE', 'sthA', 'HHKABJCC_00549', 'aslA', 'yaiT']  # qPCR ref gene candidates


for gene in annotate_list:
    gene_filter = gene_expression_statics['gene'] == gene
    ax9.scatter(gene_expression_statics[gene_filter]['mean'], gene_expression_statics[gene_filter]['cv'], label=gene)


splt.aspect_ratio(1)
ax9.set_xlabel('Average transcription level \n(a.u., n=10)')
ax9.set_ylabel('CV')

ax9.legend()
fig9.show()