#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.style as mtps
from scipy.stats import linregress

mtps.use('ggplot')


#%%
def conver_time(time, factor):
    time = sum([i*j for i, j in zip(map(int, time.split(':')), factor)])
    return time


def growth_curve_fitting(time, od):
    lgod = np.log(od)
    lambda_, n0, r_sqr, p_val, std_err = linregress(time, lgod)
    dt = np.log(2) / lambda_
    list = [lambda_, dt, np.exp(n0), r_sqr, p_val, std_err]
    return list


#%%
df = pd.read_csv(r'20191023/venus_NB_maturation_plus.csv')
sample_prefix = ['C', 'B']
#%% growth curve
growth_mask = [False if pd.isnull(i) else True for i in df['OD']]
growth_df = df.loc[growth_mask].copy()
sample = [i[0] for i in growth_df['Tube Name:']]
growth_df['Sample'] = sample
time = [conver_time(i, (60, 1)) for i in growth_df['Time']]
time = [i - time[0] for i in time]
growth_df['Time'] = time
growth_paras = [growth_curve_fitting(growth_df['Time'].loc[growth_df['Sample'] == i],
                                     growth_df['OD'].loc[growth_df['Sample'] == i]) for i in sample_prefix]
fig_1, ax1 = plt.subplots(1, 1)
sns.scatterplot(x='Time', y='OD', data=growth_df, hue='Sample', ax=ax1, s=99)
x = np.linspace(growth_df['Time'].min(), growth_df['Time'].max(), 50)
y1 = np.exp(x * growth_paras[0][0] + np.log(growth_paras[0][2]))
y2 = np.exp(x * growth_paras[1][0] + np.log(growth_paras[1][2]))
ax1.plot(x, y1, '--')
ax1.plot(x, y2, '--')
ax1.set_yscale('log')
fig_1.show()
#%% fluorescence
treat_dic = {'R': 'Room Temperature', 'E': 'Ice', 'I': 'Incubation'}
treat = [treat_dic[i[1]] if i[1] in treat_dic else 'growth' for i in df['Tube Name:']]
sample = [i[0] for i in df['Tube Name:']]
flu_df = df.copy()
flu_df['Treat'] = treat
flu_df['Sample'] = sample
time = [conver_time(i, (60, 1)) for i in flu_df['Time']]
time = [i - time[0] for i in time]
flu_df['Time'] = time
fig_2, ax = plt.subplots(2, 1, figsize=(20, 10))
sns.scatterplot(x='Time', y='P3 Mean FITC-A', data=flu_df, hue='Treat', style='Sample', ax=ax[0], s=99)
cv = [float(i[:-1]) for i in flu_df['P1 CV FITC-A']]
flu_df['P1 CV FITC-A'] = cv
sns.scatterplot(x='Time', y='P3 Mean FSC-Width', data=flu_df, hue='Treat', style='Sample', ax=ax[1], s=99)
ax[0].set_xlim(0, 400)
ax[1].set_xlim(0, 400)
fig_2.show()