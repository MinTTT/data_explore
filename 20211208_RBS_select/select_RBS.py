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
#%%

rbs_file_dir = r'./'
rbs_files = os.listdir(rbs_file_dir)
rbs_files = [file for file in rbs_files if file.split('.')[-1] == 'csv']

# dir = os.path.dirname(rbs_file_ps)
for name in rbs_files:
    df = pd.read_csv(os.path.join(rbs_file_dir, name))
    RBS_range = slice(-21, -1)
    rbs_s_list = df['tir  [Translation Initiation Rate (au)]']
    max_s, min_s = rbs_s_list.max(), rbs_s_list.min()
    ln_s_range = np.ptp(np.log(rbs_s_list))
    # select 6 RBSs.

    n = 6
    step = ln_s_range / (n - 1)
    log_s_select = [np.log(min_s) + i * step for i in range(n)]
    s_select = np.exp(log_s_select)
    s_select_index = [np.argmin(np.abs(s - rbs_s_list)) for s in s_select]
    seclected_rbs_df = df.iloc[s_select_index, :]
    seclected_rbs_df['RBS_seq'] = [rbs[RBS_range]
                                   for rbs in seclected_rbs_df['RBS_sequence  [Sequence of the RBS variant]'].tolist()]

    seclected_rbs_df.to_csv(os.path.join(rbs_file_dir, name.strip('.csv') + 'selected.csv'))



