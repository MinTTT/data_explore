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
import numpy
import pandas as pd
import numpy as np  # Or any other
# […]

# Own modules
file_ps = r'.\20220413_mRNA_level\20220413_mRNA_qPCR_2.xlsx'

raw_data = pd.read_excel(file_ps)
# Index(['Well', 'Sample name', 'Sample type', 'Gene', 'Ct', 'Mean Ct',
#        'Conc. Std.', 'Mean Conc.', 'Std.Dev. Ct', 'Std.Dev. Mean Conc.'],
#       dtype='object')

sample_names = list(set([name for name in raw_data['Sample name'] if name.split('_')[0] in ['23', '24']]))

target_gene_list = ['lacI-1', 'tetR']
reference_gene = '16s'



laci_names = [name for name in sample_names if name.split('_')[0] in ['23']]
tetr_names = [name for name in sample_names if name.split('_')[0] in ['24']]
gene_names = list(set(raw_data['Gene']))

gene_order = dict()

for gene in gene_names:
    mean_conc = []
    std_conc = []
    name_columns = []
    for name in sample_names:
        row_columns = raw_data[np.logical_and(raw_data['Sample name'] == name, raw_data['Gene'] == gene)]
        if row_columns.empty:
            pass
        else:
            name_columns.append(name)
            mean_conc.append(row_columns['Mean Conc.'].values[0])
            std_conc.append(row_columns['Std.Dev. Mean Conc.'].values[0])
    gene_order[gene] = pd.DataFrame(data={'Sample name': name_columns,
                                          f'Mean Conc._{gene}': mean_conc,
                                          f'Conc. Std._{gene}': std_conc})

# DO Normalization (ref. 16s)
normalized_data = {}
for gene in target_gene_list:
    target_data = gene_order[gene]
    ref_data = gene_order[reference_gene]
    data = pd.merge(target_data, ref_data, on='Sample name')
    data['Target/Ref'] = data[f'Mean Conc._{gene}']/data[f'Mean Conc._{reference_gene}']
    normalized_data[gene] = data


