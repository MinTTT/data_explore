#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sciplot as splt



#%%
rdf = pd.read_csv(r'./20191222/mRFP Q66C charac.csv')




#%%
def convert_time(time):
    factor = [60., 1.]
    day, time = time.split(r' ')
    min = np.sum([h*m for h, m in zip(map(float, time.split(r':')), factor)])
    day_min = float(day.split(r'/')[-1]) * 24 * 60.
    min = min + day_min
    return min

time_all = np.asarray([convert_time(i) for i in rdf['Record Time:']])
time_all_relatv = pd.Series(data=time_all - time_all.min(), name='Time')
rdf['Time'] = time_all_relatv



#%%
splt.whitegrid()
fig1, ax = plt.subplots(1, 1, figsize=(8, 8))
sns.scatterplot(x='Time', y='PE-A', hue='Sample', data=rdf.loc[rdf['Condition'] == 37.], axes=ax)
ax.set_yscale('log')
ax.set_xlabel('Time (min)')
fig1.show()
#%%
splt.whitegrid()
fig1, ax = plt.subplots(1, 1, figsize=(8, 8))
sns.scatterplot(x='Time', y='PE-A', hue='Sample', data=rdf.loc[rdf['Condition'] == 0.], axes=ax)
ax.set_yscale('log')
ax.set_xlabel('Time (min)')
fig1.show()
#%%
splt.whitegrid()
fig1, ax = plt.subplots(1, 1, figsize=(8, 8))
sns.scatterplot(x='Time', y='PE-A', hue='Sample', data=rdf.loc[rdf['Condition'] == 0.], axes=ax)
ax.set_yscale('log')
ax.set_xlim((980, 1100))
ax.set_ylim((-1000, 16000))
ax.set_xlabel('Time (min)')
fig1.show()
