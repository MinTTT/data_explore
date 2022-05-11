# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sciplot as splt
import os
import flow_api as flowcyto
import seaborn as sns
import matplotlib.animation as animation
import matplotlib

# matplotlib.use('Agg')

splt.whitegrid()
get_sec = lambda x: x.item() / 1e9


def cal_gr(time, od, **kwargs):
    '''

    :param time:
    :param od:
    :param kwargs: df, flu
    :return:
    '''
    time_index = list(np.argsort(time))
    sorted_od = np.array(od)[time_index]
    sorted_time = np.array(time)[time_index]
    del_t = np.diff(sorted_time)
    del_seconds = np.array([get_sec(x) for x in del_t])
    if 'df' in kwargs:
        df = kwargs['df']
        df = [1. if np.isnan(f) else f for f in np.array(np.array(df)[time_index])]
        df = [np.prod(df[0:i+1]) for i in range(len(df))]
        accu_od = sorted_od * df
        gr = np.diff(np.log(accu_od)) / (del_seconds / 3600)
        gener_time = np.log2(accu_od / accu_od[0])
    else:
        gr = np.diff(np.log(sorted_od)) / (del_seconds / 3600)
        gener_time = np.log2(sorted_od / sorted_od[0])
    if 'flu' in kwargs:
        flu = kwargs['flu']
        sorted_flu = np.array(flu)[time_index]
        # print(sorted_flu)
        pa = np.divide(np.divide(np.diff(sorted_flu * accu_od), (del_seconds / 3600)), accu_od[1:])
        pa = np.insert(pa, 0, np.nan, 0)
        return np.insert(gr, 0, np.nan, 0), gener_time, pa
    else:
        pass
    return np.insert(gr, 0, np.nan, 0), gener_time


# %% set data dir
data_raw_path = r'./20201004_mdeium_shift/20201004_medium_down_shift.csv'
fcs_dir = r'E:\Fu_lab_Data\cytometer data\20201004_down_shift'
flowcyto_stat = pd.read_csv(r'./20201004_mdeium_shift/20201004_medium_down_shift_flowcyto.csv', skiprows=2,
                            na_values='####')
# %% data process: do statistic
data_raw = pd.read_csv(data_raw_path)
fcs_list = [f.name for f in os.scandir(fcs_dir) if f.name.split('.')[-1] == 'fcs']
fcs_list_num = [f.name.split('.')[0] for f in os.scandir(fcs_dir) if f.name.split('.')[-1] == 'fcs']
flowcyto_stat.index = flowcyto_stat['Tube Name:']
data_raw.index = data_raw['Number']
data_raw['Time'] = pd.to_datetime(data_raw['Time'])
data_raw['Time_h'] = [time.value / 1e9 / 3600 for time in (data_raw['Time'] - data_raw['Time'].min())]
sample_list = list(set(data_raw['Sample']))
data_summary = pd.concat([data_raw, flowcyto_stat], axis=1, sort=False)
data_summary['FCS'] = [True if str(num) in fcs_list_num else False for num in data_summary['Number']]
data_summary['Growth Rate'] = [np.nan] * len(data_summary)
# calculate growth rate
for index, samp in enumerate(sample_list):
    data_select = data_summary[data_summary['Sample'] == samp]
    gr, genr_num = cal_gr(data_select['Time'], data_select['OD600'], df=data_select['Dilution_Factor'])
    _, _, green_pa = cal_gr(data_select['Time'], data_select['OD600'], df=data_select['Dilution_Factor'],
                            flu=data_select['P1 AND NOT(H2-LL) Mean Green-H'])
    _, _, red_pa = cal_gr(data_select['Time'], data_select['OD600'], df=data_select['Dilution_Factor'],
                            flu=data_select['P1 AND NOT(H2-LL) Mean Red-H'])
    data_summary.loc[data_summary['Sample'] == samp, 'Growth Rate'] = gr
    data_summary.loc[data_summary['Sample'] == samp, 'Generation'] = genr_num
    data_summary.loc[data_summary['Sample'] == samp, 'Red Translation Efficiency'] = red_pa
    data_summary.loc[data_summary['Sample'] == samp, 'Green Translation Efficiency'] = green_pa




data_summary.to_csv(data_raw_path + '.summary.csv')
# gr = cal_gr(data_raw[data_raw['Sample'] == 'M5']['Time'], data_raw[data_raw['Sample'] == 'M5']['OD600'])


# %% plot gr
n_col = 2
n_row = len(sample_list) // n_col + 1 * len(sample_list) % n_col
splt.whitegrid()
fig2, ax2 = plt.subplots(n_row, n_col, figsize=(8 * n_col, 8 * n_row))

for index, samp in enumerate(sample_list):
    x_index = index % n_col
    y_index = index // n_row
    data_select = data_summary[data_summary['Sample'] == samp]
    gr = data_select['Growth Rate']
    gr_index = gr >= 0
    diluation_mask = [True if ~np.isnan(factor) else False for factor in data_select['Dilution_Factor']]
    ax2[x_index].scatter(data_select['Time_h'][gr_index], gr[gr_index], marker='D', facecolors='None', edgecolors='r')
    ax2[x_index].scatter(data_select['Time_h'][diluation_mask], [0] * len(data_select['Time_h'][diluation_mask]),
                         marker='^', facecolors='k', edgecolors='k')
    # ax2[y_index, x_index].legend()
    ax2[x_index].grid(False)
    ax2[x_index].set_xlim(0, data_summary['Time_h'].max() * 1.1)
    ax2[x_index].set_ylim(-0.12, 2.2)
    ax2[x_index].set_title(samp)
    ax2[x_index].set_xlabel('Time (h)')
    ax2[x_index].set_ylabel('Growth Rate ($h^{-1}$)')
fig2.show()

# %% plot flu
n_col = 2
n_row = len(sample_list) // n_col + 1 * len(sample_list) % n_col
splt.whitegrid()
fig3, ax3 = plt.subplots(n_row, n_col, figsize=(8 * n_col * 1.1, 8 * n_row))

for index, ax in enumerate(ax3):
    samp = sample_list[index]
    data_select = data_summary[data_summary['Sample'] == samp]
    gr = data_select['Growth Rate']
    gr_index = gr >= 0.001
    diluation_mask = [True if ~np.isnan(factor) else False for factor in data_select['Dilution_Factor']]
    ax.scatter(data_select['Time_h'], data_select['P1 AND NOT(H2-LL) Mean Red-H'], marker='D',
               facecolors='None', edgecolors='r')
    ax.scatter(data_select['Time_h'][diluation_mask], [0] * len(data_select['Time_h'][diluation_mask]),
               marker='^', facecolors='k', edgecolors='k')
    ax3_y = ax.twinx()
    ax3_y.scatter(data_select['Time_h'], data_select['P1 AND NOT(H2-LL) Mean Green-H'], marker='D',
                  facecolors='None', edgecolors='g')

    ax.grid(False)
    ax.set_ylabel('Red Channel (a. u.)', color='red')
    ax3_y.set_ylabel('Green Channel (a. u.)', color='green')
    ax.set_xlim(0.001, data_summary['Time_h'].max() * 1.1)
    ax.set_ylim(-25, 450)
    ax.set_xlabel('Time (h)')
    ax.set_title(samp)
    ax3_y.set_ylim(-40, 10000)
    ax3_y.grid(False)
fig3.show()


# %% plot translation activity
n_col = 2
n_row = len(sample_list) // n_col + 1 * len(sample_list) % n_col
splt.whitegrid()
fig3, ax3 = plt.subplots(n_row, n_col, figsize=(8 * n_col * 1.1, 8 * n_row))

for index, ax in enumerate(ax3):
    samp = sample_list[index]
    data_select = data_summary[data_summary['Sample'] == samp]
    gr = data_select['Growth Rate']
    gr_index = gr >= 0.001
    diluation_mask = [True if ~np.isnan(factor) else False for factor in data_select['Dilution_Factor']]
    ax.scatter(data_select['Time_h'], data_select['Red Translation Efficiency'], marker='D',
               facecolors='None', edgecolors='r')
    ax.scatter(data_select['Time_h'][diluation_mask], [0] * len(data_select['Time_h'][diluation_mask]),
               marker='^', facecolors='k', edgecolors='k')
    ax3_y = ax.twinx()
    ax3_y.scatter(data_select['Time_h'], data_select['Green Translation Efficiency'], marker='D',
                  facecolors='None', edgecolors='g')

    ax.grid(False)
    ax.set_ylabel('Red Channel Translational Efficiency (a. u.)', color='red')
    ax3_y.set_ylabel('Green Channel Translational Efficiency (a. u.)', color='green')
    ax.set_xlim(0.001, data_summary['Time_h'].max() * 1.1)
    ax.set_ylim(-50, 50)
    ax.set_xlabel('Time (h)')
    ax.set_title(samp)
    ax.yaxis.set_ticks_position('left')
    ax3_y.set_ylim(-500, 2000)
    ax3_y.grid(False)
    # ax.axes.tick_params(direction='in', length=20)
fig3.show()

# %% plot translation activity normalized and growth rate
n_col = 2
n_row = len(sample_list) // n_col + 1 * len(sample_list) % n_col
splt.whitegrid()
fig3, ax3 = plt.subplots(n_row, n_col, figsize=(8 * n_col * 1.1, 8 * n_row))

for index, ax in enumerate(ax3):
    samp = sample_list[index]
    data_select = data_summary[data_summary['Sample'] == samp]
    gr = data_select['Growth Rate']
    gr_index = gr >= 0.001
    dilution_mask = [True if ~np.isnan(factor) else False for factor in data_select['Dilution_Factor']]
    ax.scatter(data_select['Time_h'], data_select['Red Translation Efficiency']/data_select['Red Translation Efficiency'].max(),
               marker='D',
               facecolors='None', edgecolors='r')
    ax.scatter(data_select['Time_h'][dilution_mask], [0] * len(data_select['Time_h'][dilution_mask]),
               marker='^', facecolors='k', edgecolors='k')
    ax.scatter(data_select['Time_h'], data_select['Green Translation Efficiency']/data_select['Green Translation Efficiency'].max(),
               marker='D', facecolors='None', edgecolors='g')
    ax3_y = ax.twinx()
    ax3_y.scatter(data_select['Time_h'][gr_index], gr[gr_index], marker='D', facecolors='None', edgecolors='k')

    ax.grid(False)
    ax.set_ylabel('Normalized Translational Efficiency (a. u.)')
    ax3_y.set_ylabel('Growth Rate (a. u.)')
    ax.set_xlim(0.001, data_summary['Time_h'].max() * 1.1)
    ax.set_ylim(-.5, 1.1)
    ax.set_xlabel('Time (h)')
    ax.set_title(samp)
    ax3_y.set_ylim(-0.1, 6)
    ax3_y.grid(False)
    # ax.axes.tick_params(direction='in', length=20)
fig3.show()

# %% animation

def update_frame(index, ax, data):
    file_name = str(data['Number'][index]) + '.fcs'
    time = data['Time_h'][index]
    try:
        flow_data = flowcyto.FCSParser(fcs_dir + '\\' + file_name)
        red_h = flow_data.dataframe['ECD-H'] / flow_data.dataframe['FSC-H'] * 1000 + .0001
        green_h = flow_data.dataframe['FITC-H'] / flow_data.dataframe['FSC-H'] * 1000 + .0001
        ax.cla()
        sns.histplot(x=np.log(red_h), y=np.log(green_h), cmap='coolwarm', ax=ax)
        ax.grid(False)
        ax.set_xlim(0, np.log(9000))
        ax.set_ylim(-1, np.log(90000))
        ax.set_title(samp + ': %.2f h' % time)
    except FileNotFoundError:
        pass


# ax4.grid(False)
# ax4.set_xlim(0, np.log(5000))
# ax4.set_ylim(-0.5, np.log(9000))

for samp in sample_list:
    data_sclc = data_summary[np.logical_and(data_summary['Sample'] == samp, data_summary['FCS'] == True)]
    data_sclc = data_sclc.sort_values(by='Time_h')
    data_sclc.index = range(len(data_sclc))
    fig4, ax4 = plt.subplots(1, 1, figsize=(8, 8))
    ani = animation.FuncAnimation(fig4, update_frame, frames=len(data_sclc),
                                  fargs=(ax4, data_sclc), cache_frame_data=False)
    ani.save(f'./20201004_mdeium_shift/animation_{samp}.avi', writer='ffmpeg', fps=4)
# fig4.show()


# %% plot OD
fig1, ax1 = plt.subplots(1, 1)
ax1.grid(False)
for samp in sample_list:
    data_sclc = data_raw[data_raw['Sample'] == samp]
    ax1.scatter(data_sclc['Time_h'], data_sclc['OD600'])
ax1.set_yscale('log')
fig1.show()

# %% plot fcs file

flow_data = flowcyto.FCSParser(fcs_dir + '\\' + fcs_list[2])
red_h = flow_data.dataframe['ECD-H'] / flow_data.dataframe['FSC-H'] * 1000
green_h = flow_data.dataframe['FITC-H'] / flow_data.dataframe['FSC-H'] * 1000

fig3, ax3 = plt.subplots(1, 1, figsize=(8, 8))
# ax3.hist2d(np.log(red_h), np.log(green_h), bins=1000)
sns.histplot(x=np.log10(red_h), y=np.log10(green_h), ax=ax3, cmap='coolwarm')
ax3.grid(False)
ax3.set_xlim(0, np.log10(5000))
ax3.set_ylim(0, np.log10(9000))
# ax3.set_xscale('log')

fig3.show()
