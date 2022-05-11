#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import sciplot as sci

#%%
sci.whitegrid()
df = pd.read_csv(r'./20200629_DIFFERENT_ORIGN_ORIENTATION\DIFFERENT_OIGIN.csv')

#%%
# show growth rate vs copy number

fig1, ax1 = plt.subplots(1, 1)
sns.regplot(x='Growth_Rate', y='Copy_Number', data=df, ax=ax1)
sns.scatterplot(x='Growth_Rate', y='Copy_Number', data=df, ax=ax1, hue='Name')
# ax1.set_yscale('log', basey=2)
ax1.set_xlabel('Growth Rate ($\mathrm{h^{-1}}$)')
ax1.set_ylabel('Copy Number (Reltiv. P/G)')
ax1.legend(edgecolor=(0, 0, 0))
fig1.show()


#%%
fig2, ax2 = plt.subplots(1, 2, figsize=(18, 8))

a1 = sns.barplot(x='Name', y='Copy_Number', data=df, yerr=df['Copy_Number_STD'][:7], hue='State', ax=ax2[0])
ax2[0].tick_params(axis='x', rotation=45)
sns.barplot(x='Name', y='Growth_Rate', data=df, hue='State', ax=ax2[1])

fig2.show()


#%%
fig3, ax3 = plt.subplots(1, 2, figsize=(28, 10))
sns.scatterplot(x='Growth_Rate', y='Green-H', data=df, hue='Name', ax=ax3[0])
sns.scatterplot(x='Growth_Rate', y='Red-H_Green', data=df, hue='Name', ax=ax3[1])
ax3[0].set_yscale('log')
ax3[1].set_yscale('log')
fig3.show()

#%%
fig3, ax3 = plt.subplots(1, 2, figsize=(28, 10))
sns.scatterplot(x='Copy_Number', y='Green-H', data=df[df['State'] == 'Green'], hue='Name', ax=ax3[0])
sns.scatterplot(x='Copy_Number', y='Red-H', data=df[df['State'] == 'Red'], hue='Name', ax=ax3[1])
# ax3[0].set_yscale('log')
# ax3[1].set_yscale('log')
fig3.show()
#%%
df['Rel_Green'] = df['Green-H'] / df['Green-H'].max()
df['Rel_Red'] = df['Red-H'] / df['Red-H'].max()
fig4, ax4 = plt.subplots(1, 1, figsize=(8, 8))
sns.regplot(x='Copy_Number', y='Rel_Green', data=df[df['State'] == 'Green'], label='G state')
sns.regplot(x='Copy_Number', y='Rel_Red', data=df[df['State'] == 'Red'], label='R state')
# sns.lmplot(x='Copy_Number_Green', y='Red-H_Green', data=df[df['State'] == 'Red'], ax=ax4[1])
# ax3[0].set_yscale('log')
# ax3[1].set_yscale('log')
ax4.set_xlabel('Copy Number (Reltiv. P/G)')
ax4.set_ylabel('Fluorescence Intensity (Reltiv.)')
ax4.legend()
fig4.show()
