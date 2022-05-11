#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import sciplot as sci

#%%

data = pd.read_csv(r'./20200313/20200312.csv')


sci.whitegrid()
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
sns.lineplot(x='Time:', y='P1 Mean FSC-A', ax=ax, hue='Tube Name:', data=data)
fig.show()
