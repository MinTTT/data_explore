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
fig_1, ax1 = plt.subplots(1, 1)
sns.scatterplot(x='Time', y='OD', data=growth_df, hue='Sample', ax=ax1)
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
#%%
df = pd.read_csv(r'20191007/venus_NB_library.csv')
# data = df.values[0:11, 0:3]
# depart_columns = ['P1 Mean FITC-A', 'P1 Mean ECD-A']
depart_columns = ['P1 Mean FITC-A']
repeat_columns = df.loc[:, [i not in set(depart_columns) for i in df.columns]]
column_name = repeat_columns.columns
column_name = column_name.append(pd.Index(['RI', 'channel']))
repeat_columns = repeat_columns.values
len_column = len(repeat_columns)
# 3 dimension data
data = np.array([np.hstack((repeat_columns, df.loc[:, i].values.reshape(-1, 1),
                            np.resize(np.array([i]), (len_column, 1)))) for i in depart_columns])

# stack in data frame
data = data.reshape(1, -1, data.shape[-1])
data = data[0, :, :]
rectruct_df = pd.DataFrame(data=data, columns=column_name)

#%%
# samples_name = data[:, 0].reshape(-1, 1)
# # reconstruct data,
# samples_names = np.hstack((samples_name, samples_name)).reshape(-1, 1).astype('str')
# channels_value = data[:, 1:3].reshape(-1, 1)
# channels_name = ['Green Fluorescence', 'Red Fluorescence']
# channels_names = np.resize(np.array(channels_name), (22, 1))
# rectruct_data = np.hstack((samples_names, channels_value, channels_names))
# rectruct_df = pd.DataFrame(data=rectruct_data, columns=['sample_name', 'RI', 'channel'])
# rectruct_df['RI'] = rectruct_df['RI'].astype('float64')
# rectruct_df['sample_name'] = rectruct_df['sample_name'].astype('category')
# rectruct_df['channel'] = rectruct_df['channel'].astype('category')
#%%
colors = ['#57D1C9', '#ED5485', '#FFFBCB']
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
g = sns.catplot(x='Tag',
                y='RI', hue='channel', kind='bar', data=rectruct_df, ax=ax, palette=colors[0:2])
xticks = ax.get_xticklabels()
ax.set_xticklabels(xticks, rotation=45)
# ax.set_yscale('log')
fig.show()
#%%
df = rectruct_df.loc[rectruct_df['strain_number'] == 1]
df = df.sort_values('RI', ascending=False)
fig2, ax2 = plt.subplots(1, 1, figsize=(28, 12))
sns.catplot(x='Tube Name:', y='RI', data=df, ax=ax2, kind='bar')
ax2.set_xlabel('Tube Name')
xticks = ax2.get_xticklabels()
ax2.set_yscale('log')
ax2.set_xticklabels(xticks, rotation=90)
fig2.savefig(r'./20191007/venus_NB_library.pdf')
#%%
export_data = df[['Tube Name:', 'RI']]
export_data.to_csv(r'./20191007/venus_NB_library_sorted.csv')