#%%
import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt
from scipy import interpolate
from tqdm import tqdm
import seaborn as sns
import sciplot as splt
import pandas as pd

#%%
def ratio_pop(args, t):
    delta_l, k, phi_0 = args
    phi_sst = 1. - k / delta_l
    if delta_l == 0:
        y = phi_0 * np.exp(-k * t)
        return y

    if k == delta_l:
        y = 1. / (delta_l * t + 1. / phi_0)
    else:
        y = 1. / ((1. / phi_0 - 1. / phi_sst) * np.exp((k - delta_l) * t) + 1. / phi_sst)
    return y