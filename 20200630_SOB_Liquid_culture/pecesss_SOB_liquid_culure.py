#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import sciplot as sci

#%%
try:
    if len(df):
        print(f'{df.shape} size data was deleted.')
        del df
except NameError:
    print('data is not imported.')

sci.whitegrid()
df = pd.read_csv(r'.\20200630_SOB_Liquid_culture\Exp_20200623_SOB_Liqid_statistic.csv',
                 na_values='####')

#%%
df =  df[df['Mask'] == 1]
#%%

df.insert(len(df.columns), column='Normalized_Green', value=(df['P2 Mean Green-H'] / df['P1 Mean Green-H']))
df.insert(len(df.columns), column='Normalized_Red', value=(df['P2 Mean Red-H'] / df['P1 Mean Red-H']))

#%%
df = df.sort_values(axis=0, by='Sample ID:')
#%%

fig, ax = plt.subplots(1 , 1)
ax.scatter(df['Sample ID:'], df['Normalized_Red'], c='r')
ax.plot(df['Sample ID:'], df['Normalized_Red'], '-r', label='mCherry')
ax.set_xlabel('Time (min)')
ax.set_ylabel('Relative Intensity (mCherry)', color='r')

ax.grid(False)
ax.set_ylim(0, 1.5)
ax2 = ax.twinx()
ax2.scatter(df['Sample ID:'], df['Normalized_Green'], c='g')
ax2.plot(df['Sample ID:'], df['Normalized_Green'], '-g', label='GFP')
ax2.set_ylabel('Relative Intensity (GFP)', color='g')
fig.tight_layout()
ax2.set_ylim(0, 3)


fig.show()
