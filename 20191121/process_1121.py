#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.style as mtps


mtps.use('ggplot')


#%%
def conver_time(time, factor):
    time = sum([i*j for i, j in zip(map(int, time.split(':')), factor)])
    return time

def conver_ratio(num):
    if isinstance(num, str):
        if num[-1] == '%':
            denum = num.split('%')
            denu = float(denum[0]) / 100.
            return denu
        else:
            return num
    else:
        return num

time_factor = [60., 1., 1./60]  # min

df = pd.read_csv(r'./20191121/20191120.csv')
medium_list = ['#1', '#2', '#3', '#4', '#5']
df['Time'] = [conver_time(i, time_factor) for i in df['Time']]
df_values = df.values
func_convert = np.frompyfunc(conver_ratio, 1, 1)
df_convert = func_convert(df_values)
df = pd.DataFrame(data=df_convert, columns=df.columns, index=df.index)


#%%
fig1, ax1 = plt.subplots(1, 1)
df['Time'] = df['Time'].astype(np.float)
df['H1-LR % Parent'] = df['H1-LR % Parent'].astype(np.float)
sns.lineplot(x='Time', y='H1-LR % Parent', hue='Medium', data=df, ax=ax1)
plt.show()