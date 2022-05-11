#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sciplot as splt

#%%
pf = pd.read_csv(r'./20191228/expression_sum.csv')
#%%
media = ['RDM_Glucose', 'RDM_Glycerol', 'MOPS_CAA_Glucose',
         'MOPS_CAA_Glycerol', 'MOPS_Glucose']
strain = ['pK4-Blank', 'ptrc-c8', 'ptrc-b11', 'j23101-b11']
pf['media'] = [media[int(name[0]) - 1] for name in pf['sample']]
pf['strain'] = [strain[int(name[-1]) - 1] for name in pf['sample']]
pf['rate'] *=60
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
pf['Translation_Activity'] = pf['rate'] * pf['FITC/FSC']
pf.to_csv(r'./20191228/constitutive_expression_in_different_media_statistic_summary.csv')
#%%
splt.whitegrid()
fig1, ax = plt.subplots(1, 1, figsize=(9, 8))
sns.scatterplot(x='rate', y='expression', hue='strain', data=pf)
sns.lineplot(x='rate', y='expression', hue='strain', data=pf, legend=False)
# ax.set_yscale('log')
ax.set_xlim(0, pf['rate'].max() + 0.1)
ax.set_xlabel('Growth Rate ($\mathrm{h}^{-1}$)')
ax.set_ylabel('mVenus Level (a.u./plate reader)')
fig1.show()


#%% FITC/FSC
fig2, ax = plt.subplots(1, 1, figsize=(9, 8))
sns.scatterplot(x='rate', y='FITC/FSC', hue='strain', data=pf)
sns.lineplot(x='rate', y='FITC/FSC', hue='strain', data=pf, legend=False)
# ax.set_yscale('log')
ax.set_xlim(0, pf['rate'].max() + 0.1)
ax.set_xlabel('Growth Rate ($\mathrm{h}^{-1}$)')
ax.set_ylabel('mVenus Level (a.u./FITC/FSC)')
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
sns.regplot(x='expression', y='FITC/FSC', data=pf)
sns.scatterplot(x='expression', y='FITC/FSC', hue='strain', data=pf)
ax.set_xlabel('mVenus Level (a.u./plate reader)')
ax.set_ylabel('mVenus Level (a.u./FITC/FSC)')
fig4.show()

#%%
fig5, ax = plt.subplots(1, 1, figsize=(9, 8))
sns.regplot(x='expression', y='FITC', data=pf)
sns.scatterplot(x='expression', y='FITC', hue='strain', data=pf)
ax.set_xlabel('mVenus Level (a.u./plate reader)')
ax.set_ylabel('mVenus Level (a.u./FITC)')
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
sns.scatterplot(x='rate', y='FSC', hue='strain', data=pf)
sns.lineplot(x='rate', y='FSC', hue='strain', data=pf, legend=False)
# ax.set_yscale('log')
ax.set_xlim(0, pf['rate'].max() + 0.1)
ax.set_xlabel('Growth Rate ($\mathrm{h}^{-1}$)')
ax.set_ylabel('FSC')
fig7.show()


#%%  PA vs lambda
fig8, ax = plt.subplots(1, 1, figsize=(9, 8))
sns.scatterplot(x='rate', y='Translation_Activity', hue='strain', data=pf)
sns.lineplot(x='rate', y='Translation_Activity', hue='strain', data=pf, legend=False)
# ax.set_yscale('log')
ax.set_xlim(0, pf['rate'].max() + 0.1)
ax.set_xlabel('Growth Rate ($\mathrm{h}^{-1}$)')
ax.set_ylabel('Translation_Activity')
fig8.show()