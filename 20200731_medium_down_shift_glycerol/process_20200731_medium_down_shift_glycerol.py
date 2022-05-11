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

#%%
growth_df = pd.read_csv(r'./20200731_medium_down_shift_glycerol/down_shift_OD_and_Time_statistic_summary.csv')
qpcr_df = pd.read_csv(r'./20200731_medium_down_shift_glycerol/down_shift_qPCR_mRNA_and_copy_number_statistic_summary.csv')
flu_df = pd.read_csv(r'./20200731_medium_down_shift_glycerol/down_shift_fluorescence_statistic_summary.csv')
#%%
all_data = pd.merge(growth_df, qpcr_df, on='Sample_Name', how='outer')
all_data = pd.merge(all_data, flu_df, on='Sample_Name', how='outer')
time_series = np.array([time_process(time, 'hr') for time in all_data['Time']])

all_data['Time'] = time_series - time_series[0]

#%%
all_data.to_csv(r'./20200731_medium_down_shift_glycerol/down_shift_all_data_in_one.csv')


