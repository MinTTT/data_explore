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
import pandas as pd
import numpy as np  # Or any other
# […]

# Own modules


pf_dxs = pd.read_excel(r'./20220407_plasmid_copy_number/20220407_plasmid_copy_number_for_steady_growth.xlsx', sheet_name ='dxs')
pf_cole1 = pd.read_excel(r'./20220407_plasmid_copy_number/20220407_plasmid_copy_number_for_steady_growth.xlsx', sheet_name ='ColE1')


sample_list = list(set(pf_dxs[pf_dxs['Sample type'] == 'Unknown']['Sample name']))
std_list = list(set(pf_dxs[pf_dxs['Sample type'] == 'Standard']['Sample name']))

pf_cole1['Sample group'] = ['_'.join(name.split('_')[:-1]) for name in pf_cole1['Sample name']]
pf_dxs['Sample group'] = ['_'.join(name.split('_')[:-1]) for name in pf_dxs['Sample name']]

sample_wo_repli = list(set(pf_cole1['Sample group'][pf_cole1['Sample group'] != '']))



# get sum and std
avg_dxs = []
std_dxs = []
avg_cole1 = []
std_cole1 = []
ct_dxs = []
ct_cole1 = []
for sample in sample_list:
    data_temp = pf_dxs[pf_dxs['Sample name'] == sample]
    avg_dxs.append(data_temp['Conc.'].mean())
    std_dxs.append(data_temp['Conc.'].std())
    ct_dxs.append(data_temp['Ct'].mean())

    data_temp = pf_cole1[pf_cole1['Sample name'] == sample]
    avg_cole1.append(data_temp['Conc.'].mean())
    std_cole1.append(data_temp['Conc.'].std())
    ct_cole1.append(data_temp['Ct'].mean())

std_ct_dxs = []
std_ct_cole1 = []

for std in std_list:
    data_temp = pf_dxs[pf_dxs['Sample name'] == std]
    std_ct_dxs.append(data_temp['Ct'].mean())

    data_temp = pf_cole1[pf_cole1['Sample name'] == std]
    std_ct_cole1.append(data_temp['Ct'].mean())

std_summary = pd.DataFrame(data={'Sample name': std_list+std_list,
                                 'Ct': std_ct_dxs+std_ct_cole1,
                                 'Gene': ['dxs']*len(std_list) + ['ColE1']*len(std_list)})

ddct_summary = pd.DataFrame(data={'Sample name': std_list,
                                  'dCt': np.array(std_ct_cole1) - np.array(std_ct_dxs)})

df_summary = pd.DataFrame(data={'Sample name': sample_list+sample_list,
                                'Conc.': avg_dxs + avg_cole1,
                                'Std Conc.': std_dxs + std_cole1,
                                'Ct': ct_dxs + ct_cole1,
                                'Gene': ['dxs']*len(sample_list) + ['ColE1']*len(sample_list)})
mean_ddCt = ddct_summary['dCt'].mean()
copy_number = []
dCt = []
for sample in sample_list:
    sample_temp = df_summary[df_summary['Sample name'] == sample]
    copy_number.append((sample_temp[sample_temp['Gene'] == 'ColE1']['Conc.'].values / sample_temp[sample_temp['Gene'] == 'dxs']['Conc.'].values)[0])
    dCt.append((sample_temp[sample_temp['Gene'] == 'ColE1']['Ct'].values - sample_temp[sample_temp['Gene'] == 'dxs']['Ct'].values)[0])
df_copy_number = pd.DataFrame(data={'Sample name': sample_list,
                                    'Copy number': copy_number,
                                    'dCt': dCt,
                                    'ddCt_fold_change': 2**-(np.array(dCt) - mean_ddCt)})

df_copy_number['Sample group'] = ['_'.join(name.split('_')[:-1]) for name in df_copy_number['Sample name']]

copy_avg = []
copy_std = []
copy_ddct_mean = []
copy_ddct_std = []

for sample in sample_wo_repli:
    data_temp = df_copy_number[df_copy_number['Sample group'] == sample]
    copy_std.append(data_temp['Copy number'].std())
    copy_avg.append(data_temp['Copy number'].mean())
    copy_ddct_mean.append(data_temp['ddCt_fold_change'].mean())
    copy_ddct_std.append(data_temp['ddCt_fold_change'].std())


# statistical data summary
df_copy_number_summary = pd.DataFrame(data={'Sample group': sample_wo_repli,
                                            'Copy number': copy_avg,
                                            'Copy number std': copy_std,
                                            'Copy number ddCt mean': copy_ddct_mean,
                                            'Copy number ddCt std': copy_ddct_std})