#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sciplot as splt
from scipy import stats
#%%
pf = pd.read_csv(r'./20200108/20200109_cm_data_sum.csv')
pf2 = pd.read_csv(r'./20191228/expression_sum.csv')
#%%
media = ['0 $\mathrm{\mu M}$', '2 $\mathrm{\mu M}$',
         '4 $\mathrm{\mu M}$']
media_index = [0, 2, 4]
strain = ['pK4-Blank', 'J23101', 'PLtetO-1', 'Ptrc']
pf['media'] = [media[media_index.index(int(name[1]))] for name in pf['sample']]
pf['strain'] = [strain[int(name[-1]) - 1] for name in pf['sample']]
pf['FITC/FSC_ALL'] = pf['P1 Mean FITC-A'] / pf['P1 Mean FSC-A'] * 10**4   # scale * 10000
pf['FITC/FSC_EXP'] = pf['P3 Mean FITC-A'] / pf['P3 Mean FSC-A'] * 10**4

media2 = ['RDM_Glucose', 'RDM_Glycerol', 'MOPS_CAA_Glucose',
         'MOPS_CAA_Glycerol', 'MOPS_Glucose']
strain2 = ['pK4-Blank', 'ptrc-c8', 'ptrc-b11', 'j23101-b11']
pf2['media'] = [media2[int(name[0]) - 1] for name in pf2['sample']]
pf2['strain'] = [strain2[int(name[-1]) - 1] for name in pf2['sample']]
pf2['rate'] *= 60
pf2['doubling_time'] /= 60

#%% combine two data
pf2.rename(columns={'FITC': 'P1 Mean FITC-A', 'FSC': 'P1 Mean FSC-A', 'FITC/FSC': 'FITC/FSC_ALL'}, inplace=True)
#%%
pf_all = pd.concat([pf, pf2], sort=False)
pf_all.to_csv(r'./20200108/all_data.csv')
#%%
norm_expression = []
for index in range(len(pf)):
    strain_na = pf['strain'][index]
    expression = pf['expression'][index]
    max_exp = pf['expression'].loc[pf['strain'] == strain_na].max()
    if max_exp:
        norm_expression.append(expression/max_exp)
    else:
        norm_expression.append(0)
pf['norm_exp'] = norm_expression
#%%
norm_fitc = []
for index in range(len(pf)):
    strain_na = pf['strain'][index]
    expression = pf['FITC'][index]
    max_exp = pf['FITC'].loc[pf['strain'] == strain_na].max()
    if max_exp:
        norm_fitc.append(expression/max_exp)
    else:
        norm_fitc.append(0)
pf['norm_fitc'] = norm_fitc

#%%
splt.whitegrid()
fig1, ax = plt.subplots(1, 1, figsize=(10, 9))
sns.scatterplot(x='rate', y='expression', hue='strain', data=pf)
sns.lineplot(x='rate', y='expression', hue='strain', data=pf, legend=False)
# ax.set_yscale('log')
ax.set_xlim(0, pf['rate'].max() + 0.1)
ax.set_xlabel('Growth Rate ($\mathrm{h}^{-1}$)')
ax.set_ylabel('mVenus Level (a.u./plate reader)')
fig1.show()

#%%
splt.whitegrid()
fig1, ax = plt.subplots(1, 1, figsize=(10, 9))
splt.twoaxiserrorbar(x='rate', y='expression', hues=['strain', 'media'], data=pf_all)
ax.set_xlim(0, pf['rate'].max() + 0.1)
ax.set_xlabel('Growth Rate ($\mathrm{h}^{-1}$)')
ax.set_ylabel('mVenus Level (a.u./plate reader)')
fig1.show()
#%%
y_list = ['expression', 'FITC/FSC_ALL', 'FITC/FSC_EXP', 'P3 Mean FSC-A', 'P3 Mean FITC-A']
r, c = [3, 2]
fig2, ax = plt.subplots(r, c, figsize=(20, 28))
for i, j in enumerate(y_list):
    y_ind = i // c
    x_ind = (i - y_ind * c) % r
    plt.sca(ax[y_ind, x_ind])
    splt.twoaxiserrorbar(x='rate', y=j, hues=['strain', 'media'], data=pf)
    cax = plt.gca()
    cax.set_xlabel('Growth Rate ($\mathrm{h}^{-1}$)')
    cax.set_ylabel(j)
fig2.show()


#%% FITC/FSC
fig2, ax = plt.subplots(1, 1, figsize=(9, 8))
sns.scatterplot(x='rate', y='FITC/FSC_ALL', hue='strain', data=pf)
sns.lineplot(x='rate', y='FITC/FSC_ALL', hue='strain', data=pf, legend=False)
# ax.set_yscale('log')
ax.set_xlim(0, pf['rate'].max() + 0.1)
ax.set_xlabel('Growth Rate ($\mathrm{h}^{-1}$)')
ax.set_ylabel('mVenus Level (a.u./FITC/FSC)')
fig2.show()

#%% FITC/FSC _P3
fig2, ax = plt.subplots(1, 1, figsize=(9, 8))
sns.scatterplot(x='rate', y='FITC/FSC_EXP', hue='strain', data=pf)
sns.lineplot(x='rate', y='FITC/FSC_EXP', hue='strain', data=pf, legend=False)
# ax.set_yscale('log')
ax.set_xlim(0, pf['rate'].max() + 0.1)
ax.set_xlabel('Growth Rate ($\mathrm{h}^{-1}$)')
ax.set_ylabel('mVenus Level (a.u./FITC/FSC_EXP)')
fig2.show()


#%% FITC_lim
fig3, ax = plt.subplots(1, 1, figsize=(9, 8))
sns.scatterplot(x='rate', y='FITC-lim', hue='strain', data=pf)
sns.lineplot(x='rate', y='FITC-lim', hue='strain', data=pf, legend=False)
# ax.set_yscale('log')
ax.set_xlim(0, pf['rate'].max() + 0.1)
ax.set_xlabel('Growth Rate ($\mathrm{h}^{-1}$)')
ax.set_ylabel('mVenus Level (a.u./FITC_Lim)')
fig3.show()

#%% FITC/FSC VS Plate reader
fig4, ax = plt.subplots(1, 1, figsize=(9, 8))
sns.regplot(x='expression', y='FITC/FSC_ALL', data=pf)
sns.scatterplot(x='expression', y='FITC/FSC_ALL', hue='strain', data=pf)
sns.scatterplot(x='expression', y='FITC/FSC', hue='strain', data=pf2)
ax.set_xlabel('mVenus Level (a.u./plate reader)')
ax.set_ylabel('mVenus Level (a.u./FITC/FSC)')
fig4.show()

#%% ALL data plot in one figure
x, y = ['expression', 'P1 Mean FITC-A']
fig5, ax = plt.subplots(1, 1, figsize=(9, 8))
sns.regplot(x=x, y=y, data=pf_all)
sns.scatterplot(x=x, y=y, hue='strain', data=pf_all)
_, _, r_sqr, _, _ = stats.linregress(pf_all[x], pf_all[y])
# ax.text(x=pf_all[x].mean()+9000, y=pf_all[y].max(), text=(r'$\mathrm{R^{2}}$ =  %f' % r_sqr), s=2)
ax.set_xlabel('mVenus Level (a.u./plate reader)')
ax.set_ylabel('mVenus Level (a.u./FITC/FSC_ALL)')
fig5.show()

#%%
fig6, ax = plt.subplots(1, 1, figsize=(9, 8))
fig1, ax = plt.subplots(1, 1, figsize=(9, 8))
sns.scatterplot(x='rate', y='norm_exp', hue='strain', data=pf)
sns.lineplot(x='rate', y='norm_exp', hue='strain', data=pf, legend=False)
# ax.set_yscale('log')
ax.set_xlim(0, pf['rate'].max() + 0.1)
ax.set_xlabel('Growth Rate ($\mathrm{h}^{-1}$)')
ax.set_ylabel('Normalized_mVenus Level (a.u./plate reader)')
fig1.show()

#%%
fig7, ax = plt.subplots(1, 1, figsize=(9, 8))
sns.scatterplot(x='rate', y='norm_fitc', hue='strain', data=pf)
sns.lineplot(x='rate', y='norm_fitc', hue='strain', data=pf, legend=False)
# ax.set_yscale('log')
ax.set_xlim(0, pf['rate'].max() + 0.1)
ax.set_xlabel('Growth Rate ($\mathrm{h}^{-1}$)')
ax.set_ylabel('mVenus Level (a.u./norm_fitc)')
fig7.show()

#%%
fig7, ax = plt.subplots(1, 1, figsize=(9, 8))
sns.scatterplot(x='rate', y='P3 Mean FSC-A', hue='strain', data=pf)
sns.lineplot(x='rate', y='P3 Mean FSC-A', hue='strain', data=pf, legend=False)
# ax.set_yscale('log')
ax.set_xlim(0, pf['rate'].max() + 0.1)
ax.set_xlabel('Growth Rate ($\mathrm{h}^{-1}$)')
ax.set_ylabel('FSC_EXP')
fig7.show()

#%% SAVE all data
pf_all.to_csv(r'./20200108/data_all.csv')