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
od_df = pd.read_csv(r'./20200709_SOB_M1_4_5_7_pECJ3/20200708_SOB_OD.csv')
flu_df = pd.read_csv(r'./20200709_SOB_M1_4_5_7_pECJ3/20200708_SOB_cyto.csv')

#%%
tube_name = flu_df['Tube Name:']
mask = [True if ((name.split('-')[0] != '1') and (name.split('-')[0] != '15')) else False for name in tube_name ]
flu_df = flu_df.loc[mask, :]
flu_df.index = pd.Index(range(len(flu_df)))
#%%
od_dic = {0: 'M1', 1: 'M4', 2: 'M5', 3: 'M7', 4: 'pECJ3'}
od_data = []
for row in range(len(od_df)):
    row_content = od_df.iloc[row, :]
    time = row_content[0]
    ods = row_content[1:]
    od_index = np.nonzero(~np.isnan(ods).values)[0][0]
    od = ods[od_index]
    sample_name = od_dic[od_index]
    od_data.append([time, od, sample_name])

od_data = pd.DataFrame(data=od_data, columns=['Time', 'OD', 'Sample'])

#%%
data_all = flu_df.join(od_data)

#%%
fig1, ax1 = plt.subplots(1, 1)
sns.lineplot(x='Time', y='OD', hue='Sample', data=data_all, ax=ax1)
sns.scatterplot(x='Time', y='OD', hue='Sample', data=data_all, ax=ax1, legend=False)
ax1.set_yscale('log')
fig1.show()
#%%
fig2, ax2 = plt.subplots(2, 2, figsize=(10, 10))
sns.lineplot(x='Time', y='P2 Mean GREEN-H', hue='Sample', data=data_all, ax=ax2[0, 0])
sns.scatterplot(x='Time', y='P2 Mean GREEN-H', hue='Sample', data=data_all, ax=ax2[0, 0], legend=False)
sns.lineplot(x='Time', y='P2 Mean RED-H', hue='Sample', data=data_all, ax=ax2[0, 1])
sns.scatterplot(x='Time', y='P2 Mean RED-H', hue='Sample', data=data_all, ax=ax2[0, 1], legend=False)
sns.lineplot(x='Time', y='P3 Mean GREEN-H', hue='Sample', data=data_all, ax=ax2[1, 0])
sns.scatterplot(x='Time', y='P3 Mean GREEN-H', hue='Sample', data=data_all, ax=ax2[1, 0], legend=False)
sns.lineplot(x='Time', y='P3 Mean RED-H', hue='Sample', data=data_all, ax=ax2[1, 1])
sns.scatterplot(x='Time', y='P3 Mean RED-H', hue='Sample', data=data_all, ax=ax2[1, 1], legend=False)
fig2.show()
#%%
sample_list = list(od_dic.values())
channel = 'P3 Mean RED-H'
channel_list = ['P2 Mean RED-H', 'P2 Mean GREEN-H',
                'P3 Mean RED-H', 'P3 Mean GREEN-H']
list_of_data = []

for sample in sample_list:
    sample_data = data_all[data_all['Sample'] == sample]
    time = sample_data['Time']
    od = sample_data['OD']
    growth_rate = np.diff(np.log(od)) / np.diff(time) * 60.
    # flu = sample_data[channel]
    # P_OD = flu * od
    # DIFF_P_OD = np.diff(P_OD) / np.diff(time)
    OD_meid = od[0:-1] + np.diff(od)/2.
    # pa = DIFF_P_OD / OD_meid
    list_of_pa = []

    for channel in channel_list:
        flu = sample_data[channel]
        P_OD = flu * od
        DIFF_P_OD = np.diff(P_OD) / np.diff(time)
        pa = DIFF_P_OD / OD_meid
        list_of_pa.append(pa.to_list())

    sample_name_list = [sample] * len(pa)
    time = time[:-1] + np.diff(time)
    data = pd.DataFrame(data=np.array(([time.to_list(), list(growth_rate), sample_name_list]
                                      + list_of_pa)).T,
                        columns=(['Time', 'Lambda', 'Sample']+['PA (%s)' % channel for channel in channel_list]))
    list_of_data.append(data)


pa_df = pd.concat(list_of_data, ignore_index=True)
pa_df['Time'] = pa_df['Time'].astype(np.float)
channels_list = ['PA (%s)' % channel for channel in channel_list]
for channel in channels_list:
    pa_df[channel] = pa_df[channel].astype(np.float)
pa_df['Lambda'] = pa_df['Lambda'].astype(np.float)

pa_df.to_csv(r'./20200709_SOB_M1_4_5_7_pECJ3/20200708_SOB_pECJ3_M1_4_7_Promoter_activity_summary.csv')
data_all.to_csv(r'./20200709_SOB_M1_4_5_7_pECJ3/20200708_SOB_pECJ3_M1_4_7_cytometry_OD_summary.csv')
#%%
fig3, ax3 = plt.subplots(2, 2, figsize=(10, 10))
xx, yy = np.meshgrid(range(2), range(2))
for index, (i, j) in enumerate(zip(xx.flatten(), yy.flatten())):
    sns.lineplot(x='Time', y=channels_list[index], hue='Sample', data=pa_df, ax=ax3[i, j])
fig3.show()


