#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import sciplot as sci
from scipy.optimize import leastsq
def convert_ratio(ratio):
    if ratio is np.nan:
        return ratio
    else:
        return np.float(ratio[:-1])/100.
#%%
try:
    if len(df):
        print(f'{df.shape} size data was deleted.')
        del df
except NameError:
    print('data is not imported.')
df = pd.read_csv(r'./20200622/Exp_20200620_1_statistic.csv', na_values='####')
sample_name = [name.split('-')[-1] for name in df['Tube Name:']]
df.insert(loc=1, column='Sample Name', value=sample_name)

#%%
generation_dic = dict(seed=0, G1=np.log2(147.5/2.5),
                  G2=np.log2(147.5/2.5*150),
                  G3=np.log2(147.5/2.5*150*150),
                  G4=np.log2(147.5/2.5*150*150*150),
                  G5=np.log2(147.5/2.5*150**4))

generation = [generation_dic[name.split('-')[0]] for name in df['Tube Name:']]
df.insert(2, 'Generation', generation)
#%%
green_ratio = np.array([convert_ratio(ratio) for ratio in df['H2-LR % Parent']])
green_red_ratio = np.array([convert_ratio(ratio) for ratio in df['H2-UR % Parent']])
red_ratio = np.array([convert_ratio(ratio) for ratio in df['H2-UL % Parent']])
total = green_ratio + green_red_ratio + red_ratio
green_ratio = green_ratio / total
green_red_ratio = green_red_ratio / total
red_ratio = red_ratio / total
df.insert(len(df.columns), 'green_ratio', green_ratio)
df.insert(len(df.columns), 'green_red_ratio', green_red_ratio)
df.insert(len(df.columns), 'red_ratio', red_ratio)
#%%
sci.whitegrid()
fig, ax = plt.subplots(3, 1, figsize=(8, 26))
sns.lineplot(x='Generation', y='green_red_ratio', data=df.loc[df['Sample Name'] != 'M5'],
            hue='Sample Name', ax=ax[0])
sns.scatterplot(x='Generation', y='green_red_ratio', data=df.loc[df['Sample Name'] != 'M5'],
            hue='Sample Name', ax=ax[0], s=400, legend=False)
sns.lineplot(x='Generation', y='red_ratio', data=df.loc[df['Sample Name'] != 'M5'],
            hue='Sample Name', ax=ax[1], legend=False)
sns.scatterplot(x='Generation', y='red_ratio', data=df.loc[df['Sample Name'] != 'M5'],
            hue='Sample Name', ax=ax[1], legend=False, s=400)
sns.lineplot(x='Generation', y='green_ratio', data=df.loc[df['Sample Name'] != 'M5'],
            hue='Sample Name', ax=ax[2], legend=False)
sns.scatterplot(x='Generation', y='green_ratio', data=df.loc[df['Sample Name'] != 'M5'],
            hue='Sample Name', ax=ax[2], legend=False, s=400)
fig.show()
#%%
fig1, ax1 = plt.subplots(1, 1, figsize=(8, 8))
sns.scatterplot(x='green_ratio', y='H1-LR Mean Cherry-H',
                data=df.loc[np.all([df['Generation'] == generation_dic['G5'], df['Sample Name'] != 'M5'], axis=0)],
                hue='Sample Name', ax=ax1, s=400)
fig1.show()

#%%
df.to_csv(r'./20200622/Exp_20200620_1_statistic_modified.csv')


