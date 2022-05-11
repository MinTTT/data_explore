#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import sciplot as sci
from scipy.optimize import leastsq

def time_process(time, unit):
    day, time = time.split()
    if unit == 's':
        factor = [24*60*60, 60*60, 60, 1]
    if unit == 'min':
        factor = [24*60, 60, 1, 1./60]
    if unit == 'hr':
        factor = [24, 1, 1./60, 1./(60*60)]
    day = [day.split('-')[-1]]
    time = day + time.split(':')
    return sum([i*j for i, j in zip(map(int, time), factor)])


def convert_ratio(ratio):
    if ratio is np.nan:
        return ratio
    else:
        return np.float(ratio[:-1])/100.


#%%
od_df = pd.read_excel(r'./20200327/20200327-ratio.xlsx', na_values='####',
                      )
od_df = od_df.iloc[:, :4]
# cyto_df = pd.read_csv(r'./20200327/2020032_ratio_cytometry.csv', na_values='####',
#                       skiprows=2)
# import data from JWZ
cyto_df = pd.read_csv(r'./20200327/Exp_20200327__MOPS_EZ_GLY_JWZ.csv', na_values='####', skiprows=2)

#%%
od_df.index = od_df['ID']
cyto_df.index = cyto_df['Tube Name:']
data_raw = pd.concat([od_df, cyto_df], axis=1)
data_raw = data_raw.drop(data_raw.columns[4], axis=1)
#%%
data_raw['Time'] = [time_process(t, unit='hr') for t in data_raw['Time']]
data_raw = data_raw.sort_values(by=['Time'])
data_raw['Time'] = data_raw['Time'] - data_raw['Time'].min()

#%%
data_raw['GFP_ratio'] = [convert_ratio(ratio) for ratio in data_raw['GFP % Parent']]
data_raw['mCherry_ratio'] = [convert_ratio(ratio) for ratio in data_raw['mCherry % Parent']]
data_raw['GR_ratio'] = [convert_ratio(ratio) for ratio in data_raw['GR % Parent']]
data_raw['GFP_ratio'] = data_raw['GFP_ratio'] / (data_raw['GFP_ratio'] + data_raw['mCherry_ratio'] + data_raw['GR_ratio'])
data_raw['mCherry_ratio'] = data_raw['mCherry_ratio'] / (data_raw['GFP_ratio'] + data_raw['mCherry_ratio'] + data_raw['GR_ratio'])
data_raw['GR_ratio'] = data_raw['GR_ratio'] / (data_raw['GFP_ratio'] + data_raw['mCherry_ratio'] + data_raw['GR_ratio'])

data_raw.to_csv(r'./20200327/20200327_ratio_sum_from_JWZ.csv')
#%%
# data = data_raw.loc[data_raw['Sample'] == '#4-1', :]
data = data_raw
sci.whitegrid()
fig1, ax = plt.subplots(1, 2, figsize=(18, 8))
sns.scatterplot(x='Time', y='mCherry_ratio', hue='Sample', data=data, ax=ax[0], legend=False)
sns.lineplot(x='Time', y='mCherry_ratio', hue='Sample', data=data, ax=ax[0], legend=False)
sns.scatterplot(x='Time', y='OD', hue='Sample', data=data, ax=ax[1])
ax[1].set_yscale('log')
ax[1].set_ylim(0.002, 0.25)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
fig1.show()
#%% perform data fitting  3 parameters
def g_func(pars, x, y):
    a, b, c = pars
    return np.log((a/y - 1.)/(b - 1.)) - c * x


data = data_raw.loc[data_raw['Sample'] == '#1-1', :]
t = data['Time'][4:]
y = data['GFP_ratio'][4:]

pars_opt, pars_cov = leastsq(g_func, (-0.0001, 0.8, -0.0001), args=(t, y))

phi_sst, phi_0, kml = pars_opt[0], pars_opt[0]/pars_opt[1], pars_opt[2]
delta_lambda = kml / phi_sst
k = kml + delta_lambda
print(delta_lambda, k)

t_sim = np.arange(t.min(), t.max())
y_sim = 1. / ((1./phi_0 - 1./phi_sst) * np.exp(kml * t_sim) + 1./phi_sst)

fig2, ax = plt.subplots(1, 1)
sns.scatterplot(x=t, y=y)
sns.lineplot(x=t_sim, y=y_sim)
fig2.show()

#%% two parameters
def g2_func(paras, x, y, phi_0):
    a, b = paras
    return np.log((a/y - 1.)/(a/phi_0 - 1.)) - b * x


data = data_raw.loc[data_raw['Sample'] == '#1-1', :]
t = data['Time'][3:]
y = data['mCherry_ratio'][3:]

pars_opt, pars_cov = leastsq(g2_func, (-0.01, -0.01), args=(t, y, y[0]))

phi_sst, kml = pars_opt
delta_lambda = kml / phi_sst
k = kml + delta_lambda
print(delta_lambda, k)

t_sim = np.arange(t.min(), t.max())
y_sim = 1. / ((1./y[0] - 1./phi_sst) * np.exp(kml * t_sim) + 1./phi_sst)

fig2, ax = plt.subplots(1, 1)
sns.scatterplot(x=t, y=y)
sns.lineplot(x=t_sim, y=y_sim)
fig2.show()