# %%
import sys
sys.path.append(r'./')
import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt
from tqdm import tqdm
import seaborn as sns
import sciplot as splt
import pandas as pd
from scipy.integrate import odeint


class ToggleSwitch:
    def __init__(self, paras=None):
        # initial the model parameters
        self.p_t = 40.
        self.p_l = 58
        self.Alpha2 = 0.5
        self.k_t = 4 * self.Alpha2
        self.k_l = 6 * self.Alpha2
        self.n_t = 2.
        self.n_l = 4.
        self.lambda_0 = 1.3
        self.theta_0 = 0.29  # from Klummp 2013, hr^-1
        self.beta = - self.lambda_0 / (self.lambda_0 + self.theta_0)
        self.Alpha = 0.00001
        self.alpha_t = 0.3 * self.Alpha
        self.alpha_l = 0.7 * self.Alpha
        self.paras_list = [self.p_t,
                           self.p_l,
                           self.Alpha2,
                           self.k_t,
                           self.k_l,
                           self.n_t,
                           self.n_l,
                           self.lambda_0,
                           self.beta,
                           self.Alpha,
                           self.alpha_t,
                           self.alpha_l]
        self.lac_cons = None
        self.tet_cons = None
        self.dev_laci_list = None
        self.gr_list = None
        self.dev_laci = None
        self.sst_index = None
        self.sst_laci_conc = None
        self.sst_tetr_conc = None
        self.sst_gr = None
        self.bistable = None
        self.laci_potential = None
        if paras is not None:
            self.paras_list = paras

    def tune_Alpha(self, alpha_):
        self.Alpha = alpha_
        self.alpha_t = 0.4 * self.Alpha
        self.alpha_l = 0.6 * self.Alpha

    def tune_lambda_0(self, lambda_):
        self.lambda_0 = lambda_
        self.beta = - self.lambda_0 / (self.lambda_0 + self.theta_0)

    def eq_func(self, tetr, laci):
        a = self.p_l / (1. + (laci / self.k_l) ** self.n_l)
        gr = self.lambda_0 / self.beta * (1. - tetr / self.p_l * (1. + (laci / self.k_l) ** self.n_l))
        b = gr * (1. - self.beta * gr / self.lambda_0)
        c = self.p_t / (1. + (tetr / self.k_t) ** self.n_t)
        return gr - self.lambda_0 * (1. - self.alpha_t * b * a - self.alpha_l * b * c)

    def dev_laci_func(self, gr, tetr, laci):
        w = gr * (1. - self.beta * gr / self.lambda_0)
        return w * self.p_t / (1. + (tetr / self.k_t) ** self.n_t) - gr * laci

    def compute_dev_laci(self, laci_):
        root_ = newton(self.eq_func, 1, args=(laci_,))  # root_ == tetR conc.
        lambda_t_ = self.lambda_0 / self.beta * (1. - root_ / self.p_l * (1. + (laci_ / self.k_l) ** self.n_l))
        dev_laci_ = self.dev_laci_func(lambda_t_, root_, laci_)
        return [root_, lambda_t_, dev_laci_]

    def potential(self, u, laci_):
        root_ = newton(self.eq_func, 1, args=(laci_,))  # root_ == tetR conc.
        lambda_t_ = self.lambda_0 / self.beta * (1. - root_ / self.p_l * (1. + (laci_ / self.k_l) ** self.n_l))
        dev_laci_ = self.dev_laci_func(lambda_t_, root_, laci_)
        return - dev_laci_

    def comput_potenital_laci(self):
        laci_potential = odeint(self.potential, 0, self.lac_cons)
        self.laci_potential = laci_potential.flatten()

    def solve_landscape(self, lac_cocs):
        self.lac_cons = lac_cocs
        self.dev_laci_list = np.array([self.compute_dev_laci(i) for i in self.lac_cons])
        self.tet_cons = self.dev_laci_list[:, 0]
        self.gr_list = self.dev_laci_list[:, 1]
        self.dev_laci = self.dev_laci_list[:, 2]

    def sovle_sst(self):
        sign_d_laci = np.sign(self.dev_laci)
        root = np.diff(sign_d_laci)
        self.sst_index = np.nonzero(root)[0]
        self.sst_laci_conc = self.lac_cons[self.sst_index]
        self.sst_tetr_conc = self.tet_cons[self.sst_index]
        self.sst_gr = self.gr_list[self.sst_index]
        if len(self.sst_index) > 1:
            self.bistable = True
        else:
            self.bistable = False

#%%
if __name__ == '__main__':
#%%
    toggle = ToggleSwitch()
    lac_cons = np.linspace(0, 10, 800)
    toggle.solve_landscape(lac_cons)
    # dev_laci_inter = interpolate.CubicSpline(x=toggle.lac_cons, y=toggle.dev_laci ** 2, )
    # dd_laci_inter = dev_laci_inter.derivative(1)
    # dd_laci_inter = dd_laci_inter(toggle.lac_cons)
    toggle.sovle_sst()

    promoter_landscape = np.linspace(0.05, 80)
    p_l, p_t = np.meshgrid(promoter_landscape, promoter_landscape)
    bistable_matrix = np.ones(p_l.shape) * 0
    for i in tqdm(range(p_l.shape[0])):
        for j in range(p_l.shape[1]):
            toggle.p_l = p_l[i, j]
            toggle.p_t = p_t[i, j]
            toggle.solve_landscape(lac_cons)
            toggle.sovle_sst()
            bistable_matrix[i, j] = toggle.bistable

#%% bistable or monostable ?
    fig1, ax = plt.subplots()
    im_matrix = bistable_matrix[-1::-1, :]
    pf = pd.DataFrame(data=im_matrix, index=np.around(promoter_landscape[-1::-1]),
                      columns=np.around(promoter_landscape))
    sns.heatmap(pf, ax=ax, cmap='YlGnBu',
                xticklabels=5,
                yticklabels=5)
    fig1.show()

#%%
    toggle = ToggleSwitch()
    lac_cons = np.linspace(0, 10, 800)
    toggle.lac_cons = lac_cons
    toggle.comput_potenital_laci()
    fig2, ax = plt.subplots()
    splt.whitegrid()
    sns.lineplot(x=toggle.lac_cons, y=toggle.laci_potential)
    fig2.show()
#%%
    p_t = 40.
    p_l = 58.
    Alpha2 = 0.5
    k_t = 6 * Alpha2
    k_l = 8 * Alpha2
    n_t = 2.
    n_l = 4.
    lambda_0 = 1.3
    beta = 0.8
    Alpha = 0.01
    alpha_t = 0.4 * Alpha
    alpha_l = 0.6 * Alpha

    # set lacI conc. = 1.

    def eq_func(tetr, laci):
        a = p_l / (1. + (laci / k_l) ** n_l)
        gr = lambda_0 / beta * (1. - tetr / p_l * (1. + (laci / k_l) ** n_l))
        b = gr * (1. - beta * gr / lambda_0)
        c = p_t / (1. + (tetr / k_t) ** n_t)
        return gr - lambda_0 * (1. - alpha_t * b * a - alpha_l * b * c)


    def dev_laci_func(gr, tetr, laci):
        w = gr * (1. - beta * gr / lambda_0)
        return w * p_t / (1. + (tetr / k_t) ** n_t) - gr * laci


    def compute_dev_laci(laci_):
        root_ = newton(eq_func, 1, args=(laci_,))  # root_ == tetR conc.
        lambda_t_ = lambda_0 / beta * (1. - root_ / p_l * (1. + (laci_ / k_l) ** n_l))
        dev_laci_ = dev_laci_func(lambda_t_, root_, laci_)
        return [root_, lambda_t_, dev_laci_]


    lacis = np.linspace(0, 15, 1000)
    dev_laci_list = np.array([compute_dev_laci(i) for i in lacis])
    dev_lacis = dev_laci_list[:, 2]

    fig, ax = plt.subplots(1, 3)
    ax[0].plot(lacis, dev_lacis)
    ax[0].plot(lacis, np.ones(shape=lacis.shape) * 0)
    # ax[0].set_ylim((-.1, .1))
    ax[1].plot(lacis, dev_laci_list[:, 0])
    ax[2].plot(lacis, dev_laci_list[:, 1])
    fig.show()
