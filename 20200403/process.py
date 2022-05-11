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
od_df = pd.read_excel(r'./20200403/20200331_ratio_EZ_MOPS_GLY_2.xlsx', na_values='NA',
                      )
od_df = od_df.iloc[:, :5]
# full the time
od_time = []
for inx, time in enumerate(od_df['Time']):
    if time is not np.nan:
        print(time)
        od_time.append(time)
    else:
        od_time.append(od_df['Time'][inx-1])
od_df['Time'] = od_time
# cyto_df = pd.read_csv(r'./20200327/2020032_ratio_cytometry.csv', na_values='####',
#                       skiprows=2)
# import data from JWZ
cyto_df = pd.read_csv(r'./20200403/20200331_MOPS_EZ_gly.csv', na_values='####', skiprows=2)
# trime unnecessary row
drop_mask = [True if tu_name[0] == '#' else False for tu_name in cyto_df['Tube Name:']]
cyto_df = cyto_df.loc[drop_mask, :]
cyto_df.index = np.arange(len(cyto_df))
data = np.empty((len(cyto_df), 5)).tolist()
for inx3, sample in enumerate(cyto_df['Tube Name:']):
    sample_no = sample.split('-')[0]
    data[inx3][:] = od_df.loc[od_df['No.'] == sample_no, :].values[0]
data = pd.DataFrame(data=data, columns=od_df.columns)

data_raw = pd.concat([data, cyto_df], axis=1)
data_raw['sample'] = data_raw['sample'].replace('_',  '-', regex=True)

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

data_raw.to_csv(r'./20200403/20200331_ratio_EZ_MOPS_GLY_2.xlsx_sum_from_JWZ.csv')
#%%
# data = data_raw.loc[data_raw['Sample'] == '#4-1', :]
data = data_raw
sci.whitegrid()
fig1, ax = plt.subplots(1, 1, figsize=(12, 8))
sns.lineplot(x='Time', y='mCherry_ratio', hue="sample", style="aTc conc. (ng/mL)",
             data=data, ax=ax, legend=False, markers=True)
# sns.scatterplot(x='Time', y='OD', hue='sample', data=data, ax=ax[1])
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
fig1.show()

#%% GR ratio vs time
# data = data_raw.loc[data_raw['Sample'] == '#4-1', :]
data = data_raw
sci.whitegrid()
fig1, ax = plt.subplots(1, 1, figsize=(12, 8))
sns.lineplot(x='Time', y='GR_ratio', hue="sample", style="aTc conc. (ng/mL)",
             data=data, ax=ax, legend=False, markers=True)
# sns.scatterplot(x='Time', y='OD', hue='sample', data=data, ax=ax[1])
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
fig1.show()
#%% GFP ratio vs time
# data = data_raw.loc[data_raw['Sample'] == '#4-1', :]
data = data_raw
sci.whitegrid()
fig1, ax = plt.subplots(1, 1, figsize=(12, 8))
sns.lineplot(x='Time', y='GFP_ratio', hue="sample", style="aTc conc. (ng/mL)",
             data=data, ax=ax, legend=False, markers=True)
# sns.scatterplot(x='Time', y='OD', hue='sample', data=data, ax=ax[1])
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
fig1.show()
#%%
sample = '2-2'
mask = data['sample'] == sample
x = data['Time'].loc[mask]
y_data = [data['GFP_ratio'].loc[mask],
     data['GR_ratio'].loc[mask],
     data['mCherry_ratio'].loc[mask]]
sci.whitegrid()
fig1, ax = plt.subplots(1, 1, figsize=(12, 8))
plt.stackplot(x, y_data[0], y_data[1], y_data[2], labels=['GFP_ratio',
                                                          'GR_ratio',
                                                          'mCherry_ratio'],
              colors=['#70F42F', '#F4B42F', '#F4512F'])
# sns.scatterplot(x='Time', y='OD', hue='sample', data=data, ax=ax[1])
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
ax.set_xlim((x.min(), x.max()))
ax.set_ylim((0, 1))
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