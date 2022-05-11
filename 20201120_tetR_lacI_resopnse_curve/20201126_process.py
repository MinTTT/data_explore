# -*- coding: utf-8 -*-

"""
{Description}
{License_info}
"""

# Built-in/Generic Imports
import os
import sys
# […]

# Libs
import pandas as pd  # Or any other
# […]

# Own modules
import numpy as np
import pandas as pd

# %%
df = pd.read_csv(r'./20201120_tetR_lacI_resopnse_curve/20201122_tetRr_lacI_resopnse_curve_cyto_statistic.csv',
                 skiprows=2)

conc_dic = dict(A=50,
                B=10,
                C=2,
                D=0.4,
                E=0.08,
                F=0.016,
                G=0.0032,
                H=0
                )
sample_dic = {'2': 'C4', '5': 'B11', '8': 'TetR', '11': 'LacI'}

location = df['Tube Name:']

concentrations = [conc_dic[i[0]] for i in location]
samples = [sample_dic[i[1:]] for i in location]
df['Sample_Name'] = samples
df['Inducer_Conc'] = concentrations
df.to_csv(r'./20201120_tetR_lacI_resopnse_curve/20201122_tetRr_lacI_resopnse_curve_cyto_statistic_summary.csv')
