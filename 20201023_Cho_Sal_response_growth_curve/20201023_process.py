# -*- coding: utf-8 -*-

"""
{Description}
{License_info}
"""

# Built-in/Generic Imports
#%%
import utiilitys as ut
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sciplot as splt
splt.whitegrid()
#%%
dir = r'./20201023_Cho_Sal_response_growth_curve'
plate3 = ut.PlateReader('palte_3')
plate3.file_load(r'./20201023_Cho_Sal_response_growth_curve/20201023_Cho_Sal_response_od.csv', method='nl', unit='hour')
plate3.find_samples()
plate3.correct_data(method='noblank')
plate3.compute_growth_rare(trimed=True)
plate3.plot_parameters(paras='rate', unit='\mathrm{h}^{-1}', dir=dir)
plate3.plot_parameters(paras='doubling_time', unit='\mathrm{h}', dir=dir)
plate3.plot_all_curve(dir=dir)
plate3.growth_para_export()


#%%
cyto_data = pd.read_csv(r'./20201023_Cho_Sal_response_growth_curve/20201023_Cho_Sal_response_cytometry.csv')
cyto_gr_link = pd.read_csv(r'./20201023_Cho_Sal_response_growth_curve/20201023_Cho_Sal_response_sample_link.csv')
cyto_gr_link = pd.merge(cyto_gr_link, cyto_data, on='Tube Name:')
gr_cyto_data = pd.merge(plate3.growth_parameters, cyto_gr_link, left_on='#_well', right_on='#_well_od')
gr_cyto_data['Inducer'] = [name.split('-')[0] for name in list(gr_cyto_data['sample_x'])]
sample_summary = [gr_cyto_data[gr_cyto_data['sample_x'] == s].describe() for s in plate3.samples_name]
gr_cyt_summ_col_name = list(sample_summary[0].columns + '_mean') + list(sample_summary[0].columns + '_std')
gr_cyto_summ = pd.DataFrame(data=None, columns=gr_cyt_summ_col_name)
gr_cyto_summ['sample'] = plate3.samples_name
gr_cyto_summ['Inducer'] = [name.split('-')[0] for name in list(gr_cyto_summ['sample'])]
for index, name in enumerate(plate3.samples_name):
    gr_cyto_summ.iloc[index, 0:-2] = sample_summary[index].loc[['mean', 'std']].values.flatten()

gr_cyto_summ[gr_cyto_summ.columns[0:-2]] = gr_cyto_summ[gr_cyto_summ.columns[0:-2]].astype(float)

gr_cyto_summ.to_csv(dir+'/20201023_Cho_Sal_response_growth_curve_mean_std.csv')

#%%

data = gr_cyto_summ[gr_cyto_summ['Inducer'] == 'Sal']
x = data['concentration(mM)_mean']
y = data['P2 Mean Green-H_mean']
yerr = data['P2 Mean Green-H_std']
fig1, ax1 = plt.subplots()
ax1.errorbar(x, y, yerr=yerr, marker='o', ls='none')
ax1.set_xscale('log')
ax1.set_xlim(1e-5, 100)
ax1.set_ylim(-10, 650)
ax1.grid(False)
ax1.set_xlabel('Inducer Conc. (mM)')
ax1.set_ylabel('Output (RPU)')
ax1.set_title('Response Function of NahR')
fig1.show()

data = gr_cyto_summ[gr_cyto_summ['Inducer'] == 'Cho']
x = data['concentration(mM)_mean']
y = data['P2 Mean Green-H_mean']
yerr = data['P2 Mean Green-H_std']
fig2, ax2 = plt.subplots()
ax2.errorbar(x, y, yerr=yerr, marker='o', ls='none')
ax2.set_xscale('log')
ax2.set_xlim(2e-5, 200)
ax2.set_ylim(-20, 650)
ax2.grid(False)
ax2.set_xlabel('Inducer Conc. (mM)')
ax2.set_ylabel('Output (RPU)')
ax2.set_title('Response Function of BetI')
fig2.show()

#%%
fig3, ax3 = plt.subplots()
sns.scatterplot(data=gr_cyto_data, x='concentration(mM)', y='rate', hue='Inducer', ax=ax3)
ax3.set_xscale('log')
ax3.set_xlim(1e-6, 30)
fig3.show()
