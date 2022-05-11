import os
import sys
# […]

# Libs
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np  # Or any other
import time
from sciplot import whitegrid, aspect_ratio

whitegrid()
# -*- coding: utf-8 -*-

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
# %% create the dataframes of growth rate vs coverage_FPKM
ls_df = []
for sample_id in strains_dict['NH3.23']:
    ls_df.append(all_pf[sample_id])
genes_list = ls_df[0]['gene'].values.tolist()
locus_tag_list = ls_df[0]['locus_tag'].values.tolist()
ls_df = pd.concat(ls_df)
ls_df.index = pd.Series(range(len(ls_df)))
gr_set = list(set(ls_df['Growth_rate'].values))
ls_vs_gr_avg = {}
ls_vs_gr_std = {}
gr_set.sort()
for gr in gr_set:
    partial_ls = ls_df[ls_df['Growth_rate'] == gr]
    sample_set = list(set(partial_ls['Sample_ID']))
    cvg_fpkm = []
    for sample in sample_set:
        cvg_fpkm.append(partial_ls[partial_ls['Sample_ID'] == sample].loc[:, 'TPM'].values.reshape(-1, 1))
    avg = np.mean(cvg_fpkm, axis=0)
    std = np.std(cvg_fpkm, axis=0)
    ls_vs_gr_avg[gr] = avg.flatten().tolist()
    ls_vs_gr_std[gr] = std.flatten().tolist()

ls_vs_gr_avg['gene'] = genes_list
ls_vs_gr_avg['locus_tag'] = locus_tag_list

ls_vs_gr_std['gene'] = genes_list
ls_vs_gr_std['locus_tag'] = locus_tag_list

ls_vs_gr_avg = pd.DataFrame(data=ls_vs_gr_avg)
ls_vs_gr_std = pd.DataFrame(data=ls_vs_gr_std)
# ls_vs_gr_avg.to_csv(os.path.join(DIR, 'NH3.24_vs_gr.csv'))


ts_df = []
for sample_id in strains_dict['NH3.24']:
    ts_df.append(all_pf[sample_id])
genes_list = ts_df[0]['gene'].values.tolist()
locus_tag_list = ts_df[0]['locus_tag'].values.tolist()
ts_df = pd.concat(ts_df)
ts_df.index = pd.Series(range(len(ts_df)))
gr_set = list(set(ts_df['Growth_rate'].values))
ts_vs_gr_avg = {}
ts_vs_gr_std = {}
gr_set.sort()
for gr in gr_set:
    partial_ts = ts_df[ts_df['Growth_rate'] == gr]
    sample_set = list(set(partial_ts['Sample_ID']))
    cvg_fpkm = []
    for sample in sample_set:
        cvg_fpkm.append(partial_ts[partial_ts['Sample_ID'] == sample].loc[:, 'TPM'].values.reshape(-1, 1))
    avg = np.mean(cvg_fpkm, axis=0)
    std = np.std(cvg_fpkm, axis=0)
    ts_vs_gr_avg[gr] = avg.flatten().tolist()
    ts_vs_gr_std[gr] = std.flatten().tolist()

ts_vs_gr_avg['gene'] = genes_list
ts_vs_gr_avg['locus_tag'] = locus_tag_list

ts_vs_gr_std['gene'] = genes_list
ts_vs_gr_std['locus_tag'] = locus_tag_list

ts_vs_gr_avg = pd.DataFrame(data=ts_vs_gr_avg)
ts_vs_gr_std = pd.DataFrame(data=ts_vs_gr_std)

# %% draw RNAP flux pattern
fld_ts_vs_gr = ts_vs_gr_avg[ts_vs_gr_avg.apply(lambda row: 0 not in row.values[:-1].astype(np.float),
                                               axis=1)]  # filter 0 contained rows

gr_vs_ts = pd.DataFrame(data=fld_ts_vs_gr.values[:, :-1].astype(np.float).T,
                        index=fld_ts_vs_gr.columns[:-1], columns=fld_ts_vs_gr.values[:, -1])
gr_vs_ts_corr = gr_vs_ts.corr()  # correlation matrix

import seaborn as sns
from scipy.cluster import hierarchy

distence_mat = hierarchy.distance.pdist(gr_vs_ts_corr.values)
link = hierarchy.linkage(distence_mat, method='complete')
index = hierarchy.fcluster(link, 0.5 * distence_mat.max(), 'distance')
columns = [gr_vs_ts_corr.index.tolist()[i] for i in np.argsort(index)]
sort_ts_corr = gr_vs_ts_corr.reindex(index=columns, columns=columns)
#
fig2, ax2 = plt.subplots(1, 1, figsize=(20, 20))
ax2.imshow(sort_ts_corr, origin='lower', cmap='coolwarm')
fig2.show()
#
sort_gr_ts = gr_vs_ts.reindex(columns=columns)
fig3, ax3 = plt.subplots(1, 1, figsize=(5, 100))
ax3.imshow(sort_gr_ts.T / np.max(sort_gr_ts.T, axis=1).values.reshape(-1, 1), cmap='coolwarm', aspect='auto')
ax3.set_yticks(np.arange(len(sort_gr_ts.T)))
ax3.set_yticklabels(sort_gr_ts.T.index.to_list(), fontsize=1)
ax3.set_xticks(np.arange(len(sort_gr_ts.T.columns)))
ax3.set_xticklabels(sort_gr_ts.T.columns.to_list(), fontsize=1)
ax3.tick_params(width=0.5, length=1)
for key, spine in ax3.spines.items():
    spine.set_linewidth(0.5)
fig3.show()
fig3.savefig(os.path.join(DIR, 'NH3.24_RNAP_flux_pattern.svg'))
fig3.savefig(os.path.join(DIR, 'NH3.24_RNAP_flux_pattern.pdf'))
sort_gr_ts.T.to_csv(os.path.join(DIR, "NH3.24_sorted_RNAP_flux_pattern.csv"))
# %%
gr_avg_dic = {"NH3.24": ts_vs_gr_avg, "NH3.23": ls_vs_gr_avg}
part_genes = ['tetR-sfgfp', 'torI', 'phnO', 'phnC', 'waaJ', 'rclR', 'fis', 'comR', 'aspS', 'bamA', 'bamB', 'hflK',
              'hisS', 'minD', 'cpdA', 'surA']
strain = 'NH3.24'
gr_avg = gr_avg_dic[strain]

fig1, ax1 = plt.subplots(4, 4, figsize=(50, 50))
ax_list = ax1.flatten().tolist()
for i, gene_name in enumerate(part_genes):
    ax = ax_list[i]
    gene_expression_level = gr_avg[gr_avg['gene'] == gene_name].iloc[:, :-2].values.flatten()
    norm_exp = gene_expression_level / gene_expression_level.max()
    ax.plot(gr_set, norm_exp, label=gene_name, ls=':')
    ax.scatter(gr_set, norm_exp)
    ax.legend()
    ax.set_ylim(0, 1)
# ax1.set_yscale('log')
# ax1.set_ylim(-1000, 16000)
fig1.show()

# %% find RNA flux of genes related to RNA
func_ribosome = lambda gene_name: gene_name[:3] == 'rpl' or gene_name[:3] == 'rps' \
                                  or gene_name[:3] == 'sra' or gene_name[:3] == 'rpm'

gr_vs_ribo = {}

for strain, data in gr_avg_dic.items():
    gr_vs_ribo[strain] = data[data['gene'].apply(func_ribosome)]

for strain, data in gr_vs_ribo.items():
    data.to_csv(os.path.join(DIR, f'{strain}_gr_vs_ribosome_proteins.csv'))

for strain, data in gr_avg_dic.items():
    fpkm_ribo = data.iloc[:, :5]  # type: pd.DataFrame
    gr_rate = fpkm_ribo.columns.to_list()

    fpkm_ribo_data = fpkm_ribo.values.T
    norm_fpkm_ribo_data = fpkm_ribo_data.sum(axis=1)

# %%

gr_avg_dic = {"NH3.24": ts_vs_gr_avg, "NH3.23": ls_vs_gr_avg}


def locus_filter(tag): return tag.split('_')[0] in ['ts', 'ls']


gr = []
ratio = []
strains = []
for name, data in gr_avg_dic.items():
    filter = data['locus_tag'].apply(locus_filter).tolist()
    circuits_FPKM = np.sum(data.iloc[filter, :5], axis=0)
    total_FPKM = np.sum(data.iloc[:, :5], axis=0)
    circuits_FPKM_ratio = circuits_FPKM / total_FPKM
    gr += circuits_FPKM_ratio.index.tolist()
    ratio += circuits_FPKM_ratio.values.tolist()
    strains += [name] * len(circuits_FPKM_ratio)

gr_circuits_ratio = dict(growth_rate=gr, circuits_FPKM_ratio=ratio, strains=strains)
df_gr_circuits_ratio = pd.DataFrame(data=gr_circuits_ratio)

fig4, ax4 = plt.subplots(1, 1, figsize=(16, 16))

ax4.scatter(df_gr_circuits_ratio[df_gr_circuits_ratio['strains'] == 'NH3.24']['growth_rate'],
            df_gr_circuits_ratio[df_gr_circuits_ratio['strains'] == 'NH3.24']['circuits_FPKM_ratio'])
ax4.scatter(df_gr_circuits_ratio[df_gr_circuits_ratio['strains'] == 'NH3.23']['growth_rate'],
            df_gr_circuits_ratio[df_gr_circuits_ratio['strains'] == 'NH3.23']['circuits_FPKM_ratio'])
ax4.set_xlabel('Growth rate ($h^{-1}$)')
ax4.set_ylabel('RNP flux')

fig4.show()

# %% Show normalization quality

fig6, ax6 = plt.subplots(1, 1, figsize=(30, 8))
box_data = []
for i, key in enumerate(list(all_pf.keys())):
    box_data.append(all_pf[key]['TPM'])
    ax6.text(i + 0.5, 3e4, '_'.join(key.split('_')[1:]))
ax6.boxplot(box_data, 1)
ax6.set_yscale('log')
ax6.set_ylim(0.1, 1e5)
ax6.set_xlabel('Sample #')
ax6.set_ylabel('Transcription level \n(a.u., coverage_FPKM)')

fig6.show()
# %% Compare the samples
import sciplot as splt

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

    ax5[i[0], i[1]].set_xscale('log')
    ax5[i[0], i[1]].set_yscale('log')
    ax5[i[0], i[1]].set_xlim(2e-2, 2e4)
    ax5[i[0], i[1]].set_ylim(2e-2, 2e4)

    splt.aspect_ratio(1, ax5[i[0], i[1]])
fig5.show()

# data1, data2 = all_pf['1321_4_1'], all_pf['1321_9_1']
# circuits_filter = data1['locus_tag'].apply(lambda tag: tag.split('_')[0] == 'ts' or tag.split('_')[0] == 'ls')
# ax5.scatter(data1['coverage_FPKM'], data2['coverage_FPKM'], c='#B2BABB')
# ax5.scatter(data1[circuits_filter]['coverage_FPKM'], data2[circuits_filter]['coverage_FPKM'], c='#82E0AA')
# splt.aspect_ratio(1)
# ax5.set_xscale('log')
# ax5.set_yscale('log')
#
# fig5.show()

# %% GDE test

from scipy.stats import mannwhitneyu
from scipy.stats import rankdata
from statsmodels.stats.multitest import fdrcorrection

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

# %%
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

fig9, ax9 = plt.subplots(1, 1, figsize=(10, 10))

ax9.scatter(gene_expression_statics['mean'], gene_expression_statics['cv'], alpha=.5, c='#B2BABB' )
ax9.set_xscale('log')
ax9.set_yscale('log')
rho_filter = gene_expression_statics['gene'] == 'rho'
cyoa_filter = gene_expression_statics['gene'] == 'cyoA'
ispA_filter = gene_expression_statics['gene'] == 'ispA'
fabB_filter = gene_expression_statics['gene'] == 'fabB'

ax9.scatter(gene_expression_statics[rho_filter]['mean'], gene_expression_statics[rho_filter]['cv'], label='rho')
ax9.scatter(gene_expression_statics[cyoa_filter]['mean'], gene_expression_statics[cyoa_filter]['cv'], label='cyoA')
ax9.scatter(gene_expression_statics[ispA_filter]['mean'], gene_expression_statics[ispA_filter]['cv'], label='ispA')
ax9.scatter(gene_expression_statics[fabB_filter]['mean'], gene_expression_statics[fabB_filter]['cv'], label='fabB')

splt.aspect_ratio(1)
ax9.set_xlabel('Average transcription level \n(a.u., n=10)')
ax9.set_ylabel('CV')

ax9.legend()
fig9.show()





# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
from random import shuffle
import sys


def calculate_p(df):
    df['symbol'] = df['target'].apply(
        lambda x: x.split('_')[0] if (x.split('_')[-3].endswith('-P1P2')) else x.split('_')[0] + '_' +
                                                                               x.split('_')[-3].split('-')[-1])
    genes = [i for i in list(df['symbol'].unique()) if 'non-targeting' not in i]
    num_of_genes = len(genes)
    ntc_epsilons = list(df[df['symbol'].str.contains('non-targeting')]['epsilon'])
    gene_eps_p = {}
    for gene in genes:
        df_gene = df[df['symbol'] == gene]
        epsilons = list(df_gene['epsilon'])
        avr_epsilon = np.mean(sorted(epsilons, key=abs)[-3:])
        x, pvalue = mannwhitneyu(epsilons, ntc_epsilons, alternative='two-sided')
        gene_eps_p[gene] = [avr_epsilon, pvalue]
    for j in range(num_of_genes):
        shuffle(ntc_epsilons)
        ntc_selected = ntc_epsilons[:5]
        avr_epsilon = np.mean(sorted(ntc_selected, key=abs)[-3:])
        x, pvalue = mannwhitneyu(ntc_selected, ntc_epsilons, alternative='two-sided')
        gene_eps_p['NTC_' + str(j)] = [avr_epsilon, pvalue]
    df = pd.DataFrame(gene_eps_p).T
    df.columns = ['epsilon', 'pvalue']
    df['product'] = df['epsilon'] * (-np.log10(df['pvalue']))
    df.sort_values('product', ascending=False)
    df.reset_index(inplace=True)
    return df


def product_threshold_fdr(df, fdr=0.05):
    maxi = abs(df['product']).max()
    for pro in np.arange(0, maxi, 0.1):
        df_thres = df[abs(df['product']) > pro]
        if (1.0 * len(df_thres[df_thres['index'].str.contains('NTC')]) / len(df_thres)) < fdr:
            break
    return pro, df_thres


def volcano_plot(df, pheno='epsilon', p='pvalue', product='product', savefig=True, figure_name='volcano_plot',
                 product_thres=20, genes_to_label=[]):
    genes = list(df['index'])
    phenotype = list(df[pheno])
    p_value = list(-np.log10(df[p]))
    products = list(df[product])
    plt.figure(figsize=[10, 8])
    texts = []
    hits_epsilon = []
    hits_p = []
    if len(genes_to_label) != 0:
        clusters = [i[:-1] for i in genes_to_label if i[-1] == '*']
        for i in clusters:
            for gene in genes:
                if i in gene:
                    genes_to_label.append(gene)
        for x, y, s in zip(phenotype, p_value, genes):
            if s in genes_to_label:
                texts.append(s)
                hits_epsilon.append(x)
                hits_p.append(y)

    plt.scatter(df[pheno], -np.log10(df[p]), c='red', s=8, label='other genes')
    plt.scatter(hits_epsilon, hits_p, c='r', s=20)
    df_ntc = df[df['index'].str.contains('NTC')]
    plt.scatter(df_ntc[pheno], -np.log10(df_ntc[p]), c='grey', s=8, label="NTC")
    for i in range(len(texts)):
        plt.annotate(texts[i], (hits_epsilon[i], hits_p[i]))
    plt.xlabel('Phenotype', fontsize=14)
    plt.ylabel('-log10 P', fontsize=14)

    plt.ylim(0, max(-np.log10(df[p])) + 0.5)
    plt.xlim(min(df[pheno]) - 0.5, max(df[pheno]) + 1)
    plt.title(figure_name, fontsize=18)
    plt.legend(loc=3, fontsize='large', fancybox=True)
    if savefig == True:
        plt.savefig(output_path + '/' + figure_name + '.pdf')
    plt.show()


if len(sys.argv) != 3:
    sys.stderr.write("USAGE: %s <l2es_path> <output_path>\n" % sys.argv[0])
    sys.exit(1)
l2es_path = sys.argv[1]
output_path = sys.argv[2]

print('counts threshold: ')
counts_threshold = float(raw_input('-->'))
print('FDR: ')
fdr = float(raw_input('-->'))
print('Figure name: ')
figure_name = raw_input('-->')
print('path of gene list: ')
gene_list_path = raw_input('-->')

if len(gene_list_path) == 0:
    genes_to_label = []
else:
    with open(gene_list_path.strip(), 'r') as f:
        genes_to_label = [i.strip() for i in f.readlines()]
for i in genes_to_label:
    print
    i

df = pd.read_table(l2es_path, names=['target', 'epsilon', 'count1', 'count2'])
df = df[~((df['count1'] < counts_threshold) & (df['count2'] < counts_threshold))]
df_eps_p = calculate_p(df)
product_thres, df_hits = product_threshold_fdr(df_eps_p, fdr=fdr)
volcano_plot(df_eps_p, figure_name=figure_name, savefig=True, product_thres=product_thres,
             genes_to_label=genes_to_label)
df_hits.sort_values('product', ascending=False).to_csv(
    output_path + '/' + figure_name + '_fdr0.05_product%s_hits.csv' % product_thres, index=False)
df_eps_p.sort_values('product', ascending=False).to_csv(output_path + '/' + figure_name + '_all_gene.csv', index=False)
