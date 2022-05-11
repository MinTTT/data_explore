# -*- coding: utf-8 -*-

"""

@author: Pan M. CHU
@Email: pan_chu@outlook.com
"""

# Built-in/Generic Imports
import os
import sys
# […]

# Libs
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np  # Or any other
import time
from sciplot import whitegrid, aspect_ratio

whitegrid()


# […]

# Own modules

def cal_acc_od(od: list, dilute_facter: list) -> np.ndarray:
    acc_od = []
    odfactor = 1.
    for i in range(len(od)):
        if np.isnan(dilute_facter[i]):
            pass
        else:
            odfactor *= dilute_facter[i]
        acc_od.append(od[i] * odfactor)
    return np.array(acc_od)


def second2hour(seconds: float) -> float:
    hors, seconds = divmod(seconds, 3600)
    mins, seconds = divmod(seconds, 60)
    return hors + mins / 60 + seconds / 3600


def cal_alpha(ps, ods, ts):
    od1, od2 = ods
    p1, p2 = ps
    t1, t2 = ts
    od_prime = np.exp(0.5 * (np.log(od1) + np.log(od2)))
    t_prime = 0.5 * (t1 + t2)
    alpha = (p2 * od2 - p1 * od1) / (t2 - t1) / od_prime
    gr = np.log(od2 / od1) / (t2 - t1)
    return [alpha, t_prime, gr]


fmt2sec = lambda timefmt: time.mktime(time.strptime(timefmt, "%Y-%m-%d %H:%M:%S"))

sample_list = ["NH3.23", "NH3.24"]
sample_data_dict = dict()
for name in sample_list:
    od_df = pd.read_excel(r"./20210708_down_up_shift/20210630_downshift_od_data.xlsx", sheet_name=name,
                          usecols=range(6))
    od_df.columns.values[-1] = 'Tube Name'
    flu_df = pd.read_excel(r"./20210708_down_up_shift/20210630_downshift_cyto_data.xlsx", sheet_name=name,
                           usecols=range(4))
    data_all = pd.merge(od_df, flu_df, on='Tube Name', how='outer')
    acc_od = cal_acc_od(data_all['OD600'], data_all['Dilution_Factor'])
    seconds = np.array([fmt2sec(tm) for tm in data_all['Time']])
    hours = np.array([second2hour(secd) for secd in seconds])
    data_all['Accumulate_od'] = acc_od
    data_all['Seconds'] = seconds
    data_all['Hours'] = hours
    sample_data_dict[name] = data_all

alpa_time_dict = {}
for name in sample_list:
    data_frame = sample_data_dict[name]
    del data_frame['Dilution_Factor']
    data_mask = pd.isna(data_frame)
    data_mask = [False if data_mask.values[i, :].sum() >= 1 else True for i in range(len(data_mask))]
    data_frame = sample_data_dict[name].loc[data_mask]
    od_list = data_frame["OD600"]
    flu_list = data_frame["Green-H"]
    time_list = data_frame["Hours"]
    alph_time = [cal_alpha(flu_list[i:i + 2], od_list[i:i + 2], time_list[i:i + 2]) for i in range(len(flu_list) - 1)]
    alpa_time_dict[name] = np.array(alph_time)

colors = {"NH3.23": "#088408", "NH3.24": "#FF6000"}
pars = dict(s=700, marker='v', facecolors='none', lw=5)
fig1, ax = plt.subplots(1, 3, figsize=(25, 8))
for i, name in enumerate(sample_list):
    data = alpa_time_dict[name]
    data_mask = data[:, 2] > 0
    data = data[data_mask, :]
    ax[0].scatter(data[:, -1], data[:, 0], edgecolors=colors[name], **pars)
    ax[0].plot(data[:, -1], data[:, 0], '--', lw=1)

    ax[1].scatter(data[:, -1], data[:, 0] / data[:, -1], edgecolors=colors[name], **pars)
    ax[2].scatter(data[:, -1], data[:, -1] / data[:, 0], edgecolors=colors[name], **pars)
ax[0].set_ylabel('Relative $\\alpha_{R,G}$')
ax[1].set_ylabel('Relative $\\widetilde{\\alpha}_{R,G}$')
ax[2].set_ylabel('Relative $\\widetilde{K}_{DR,DG}$')
ax[0].set_xlim(0, 2)
ax[1].set_xlim(0, 2)
ax[2].set_xlim(0, 2)
ax[0].set_xlabel('Growth rate ($h^{-1}$)')
ax[1].set_xlabel('Growth rate ($h^{-1}$)')
ax[2].set_xlabel('Growth rate ($h^{-1}$)')
fig1.show()

fig1.savefig(f"{os.path.join(os.path.split(r'./20210708_down_up_shift/20210630_downshift_od_data.xlsx')[0], 'alpha_gr.svg')}", transparent=True)

colors = {"NH3.23": "#088408", "NH3.24": "#FF6000"}
pars = dict(s=700, marker='v', facecolors='none', lw=5)
fig2, ax = plt.subplots(1, 1, figsize=(16, 8))
for i, name in enumerate(sample_list):
    data = alpa_time_dict[name]
    data_mask = data[:, 2] > 0
    data = data[data_mask, :]
    ax.scatter(data[:, 1]-data[:, 1].min(), data[:, 0]/data[:, 0].max(), edgecolors=colors[name], **pars)

ax.set_ylim(-0.1, 1.2)
aspect_ratio(0.5)
ax.set_ylabel('Relative $\\alpha_{R,G}$')
ax.set_xlabel('Time ($h$)')

fig2.show()

fig2.savefig(f"{os.path.join(os.path.split(r'./20210708_down_up_shift/20210630_downshift_od_data.xlsx')[0], 'alpha_time.svg')}", transparent=True)

