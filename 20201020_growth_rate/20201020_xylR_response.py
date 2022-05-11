# -*- coding: utf-8 -*-

"""
{Description}
{License_info}
"""
#%%
import utiilitys as ut
import pandas as pd
#%%
dir = r'./20201020_growth_rate'
plate3 = ut.PlateReader('palte_3')
plate3.file_load(r'./20201020_growth_rate/growth rate_20201020_xylR.csv', method='nl', unit='hour')
plate3.find_samples()
plate3.correct_data(method='noblank')
plate3.compute_growth_rare(trimed=True)
plate3.plot_parameters(paras='rate', unit='\mathrm{h}^{-1}', dir=dir)
plate3.plot_parameters(paras='doubling_time', unit='\mathrm{h}', dir=dir)
plate3.plot_all_curve(dir=dir)
plate3.growth_para_export()

