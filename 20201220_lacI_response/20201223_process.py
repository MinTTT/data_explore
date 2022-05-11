# -*- coding: utf-8 -*-

"""
{Description}
{License_info}
"""

# Built-in/Generic Imports
# %%
import utiilitys as ut
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sciplot as splt

splt.whitegrid()
# %%
dir = r'./20201220_lacI_response'
plate3 = ut.PlateReader('palte_3')
plate3.file_load(r'./20201220_lacI_response/20201220_LACI_RESPONSE_od.csv', method='nl', unit='hour')
plate3.find_samples()
plate3.correct_data(method='noblank')
plate3.compute_growth_rare(trimed=True)
plate3.plot_parameters(paras='rate', unit='\mathrm{h}^{-1}', dir=dir)
plate3.plot_parameters(paras='doubling_time', unit='\mathrm{h}', dir=dir)
plate3.plot_all_curve(dir=dir)
plate3.growth_para_export()

# %%
cyto_data = pd.read_csv(r'./20201220_lacI_response/20201220_lacI_response_cytometry.csv', skiprows=2)
cyto_gr_link = pd.read_csv(r'./20201220_lacI_response/20201220_LACI_RESPONSE_sample_link.csv')
cyto_gr_link = pd.merge(cyto_gr_link, cyto_data, on='Tube Name:')
gr_cyto_data = pd.merge(plate3.growth_parameters, cyto_gr_link, left_on='#_well', right_on='#_well_od')
gr_cyto_data['Inducer'] = [name.split('-')[0] for name in list(gr_cyto_data['sample_x'])]
sample_summary = [gr_cyto_data[gr_cyto_data['sample_x'] == s].describe() for s in plate3.samples_name]
gr_cyt_summ_col_name = list(sample_summary[0].columns + '_mean') + list(sample_summary[0].columns + '_std')
gr_cyto_summ = pd.DataFrame(data=None, columns=gr_cyt_summ_col_name)
gr_cyto_summ['sample'] = plate3.samples_name
gr_cyto_summ['Inducer'] = [name.split('-')[0] for name in list(gr_cyto_summ['sample'])]
for index, name in enumerate(plate3.samples_name):
    gr_cyto_summ.iloc[index, 0:-2] = sample_summary[index].loc[['mean', 'std']].values.flatten()

gr_cyto_summ[gr_cyto_summ.columns[0:-2]] = gr_cyto_summ[gr_cyto_summ.columns[0:-2]].astype(float)

gr_cyto_summ.to_csv(dir + '/20201220_LACI_RESPONSE_mean_std.csv')

# %% This part used to optimize the out put
from scipy.optimize import leastsq


def obj_func(pars, ip, op):
    y_min, y_max, n, k = pars
    y_prime = y_min + (y_max - y_min) * (k ** n / (ip ** n + k ** n))
    return op - y_prime

def opt_func(pars, x):
    y_min, y_max, n, k = pars
    return y_min + (y_max - y_min) * (k ** n / (x ** n + k ** n))


in_and_o = pd.read_csv(r'./20201220_lacI_response/20201220_lacI_Input_output.csv')
ip = in_and_o['Input']
op = in_and_o['Output']
pars_init = np.array([30, 800, 2, 1.5])
pars_opmz = leastsq(obj_func, pars_init, args=(ip, op,))


x = np.arange(0, 17, 0.1)
y = opt_func(pars_opmz[0], x)

fitting_curve = pd.DataFrame(data=dict(x=x, y=y))
fitting_curve.to_csv(r'./20201220_lacI_response/20201220_curve_fitting.csv')