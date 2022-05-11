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
# […]
from matplotlib import colors, cm
# Own modules
import sciplot
import sciplot as splt

splt.whitegrid()

# data_ps = r"./20210630_downshift/data/L3_green_red_flu_generation.xlsx"
data_ps = r'./20210630_downshift/data/downshift_20210630_initial_R.xlsx'
sheet_name = ['L3_R', 'M5_R', 'L4_R']
data_dict = {}
for name in sheet_name:
    data_frame = pd.read_excel(data_ps, usecols=range(5), sheet_name=name)

    data_mask = pd.isna(data_frame)
    data_mask = [False if data_mask.values[i, :].sum() >= 1 else True for i in range(len(data_mask))]
    data_frame = data_frame.loc[data_mask]
    data_dict[name] = data_frame

new_coolwarm = cm.get_cmap('Blues')
new_coolwarm = new_coolwarm(np.linspace(0, 1, 256))[50:]  # [-1::-1]  # [60:-60]
new_coolwarm = colors.ListedColormap(new_coolwarm)
norm = colors.Normalize(vmin=0, vmax=22)
#%%
fig1, axs = plt.subplots(1, 4, figsize=(30, 8))
for i, ax1 in enumerate(axs[:3]):
    data_frame = data_dict[sheet_name[i]]
    label_set = list(set(data_frame['Label']))
    for label in label_set:
        part_data = data_frame[data_frame['Label'] == label]
        if label == 1:

            scatter_pars = dict(s=500, c=part_data['Generation'])
            sct = ax1.scatter(part_data['Red-H']/1000, part_data['Green-H']/800,
                        cmap=new_coolwarm, norm=norm, **scatter_pars)
        else:
            scatter_pars = dict(s=500, facecolors='none', edgecolors=new_coolwarm(norm(part_data['Generation'])),
                                lw=6)
            ax1.scatter(part_data['Red-H']/1000, part_data['Green-H']/800, **scatter_pars)
        # # down shift
        # down_point = part_data[part_data['Shift'] == 'Down']
        # ax1.scatter(down_point['Red-H']/1000, down_point['Green-H']/800,
        #             facecolors='none', lw=5, edgecolors='#E91E63', s=500)
        # # up shift
        # down_point = part_data[part_data['Shift'] == 'Up']
        # ax1.scatter(down_point['Red-H']/1000, down_point['Green-H']/800,
        #             facecolors='none', lw=5, edgecolors='#FF9800', s=500)
    ax1.set_title(sheet_name[i], pad=20)
    ax1.set_xlim(1e-3, 2)
    ax1.set_ylim(1e-3, 2)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_aspect('equal')
    ax1.set_xticks([1e-2, 1])
    ax1.set_yticks([1e-2, 1])
    ax1.set_yticks([], minor=True)
    ax1.set_xticks([], minor=True)
    ax1.tick_params(axis='x', colors='#FF6000')
    ax1.tick_params(axis='y', colors='#008000')

# cax = plt.axes([0.75, 0.1, 0.025, 0.6])
cp = plt.colorbar(sct, ax=axs[-1], orientation='horizontal')
cp.ax.xaxis.set_ticks_position('top')

axs[-1].axis('off')
cp.set_ticks([0, 10, 20])

cp.set_label('Generation')
fig1.show()
fig1.savefig(os.path.join(os.path.split(data_ps)[0], f"{os.path.split(data_ps)[-1]}_image.svg"), transparent=True)

#%%

growth_pd = pd.read_excel(r"./20210630_downshift/data/time_growth_rate_generation.xlsx",
                          usecols=range(4))
data_mask = pd.isna(growth_pd)
data_mask = [False if data_mask.values[i, :].sum() >= 1 else True for i in range(len(data_mask))]
data_frame = growth_pd.loc[data_mask]
fig2, ax2 = plt.subplots(1, 1, figsize=(10, 14))
ax2.plot(data_frame['Time'], data_frame['Growth_rate'], '-', c='k', lw=7, zorder=0)

ax2.scatter(data_frame['Time'].iloc[::2], data_frame['Growth_rate'].iloc[::2],
            edgecolors=new_coolwarm(norm(data_frame['Generation'].iloc[::2])), s=500, facecolors='none',
            marker='s', lw=8)
# # down shift
# down_point = data_frame[data_frame['Shift'] == 'Down']
# ax2.scatter(down_point['Time'], down_point['Growth_rate']+.1,
#             facecolors='none', lw=5, edgecolors='#E91E63', s=500, marker='$\downarrow$')
# # up shift
# down_point = data_frame[data_frame['Shift'] == 'Up']
# # ax2.scatter(down_point['Time'], down_point['Growth_rate']+.1,
# #             edgecolor='none', lw=5, facecolor='#FF9800', s=500, marker='$\downarrow$')
# ax2.scatter(down_point['Time'], down_point['Growth_rate']+.1,
#             c='#FF9800', s=1000, marker='$\downarrow$')

ax2.set_ylim(-0.1, 2.1)
sciplot.aspect_ratio(0.7)
ax2.set_xticks([0, 12, 24, 36])
# ax2.set_yticks([0, 0.5, 1.0, 1.5])

cax = plt.axes([0.1, 0.8, 0.8, 0.04])
cp = plt.colorbar(sct, ax=ax2, cax=cax, orientation="horizontal")
cp.ax.xaxis.set_ticks_position('bottom')
axs[-1].axis('off')
cp.set_ticks([0, 10, 20])
cp.set_label('Generation', labelpad=5)

# ax2.set_aspect(0.7*data_frame['Time'].max()/data_frame['Growth_rate'].max())
fig2.show()
fig2.savefig(f"{os.path.join(os.path.split(data_ps)[0], 'growth_generation.svg')}", transparent=True)


