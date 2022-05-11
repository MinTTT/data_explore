#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
import statsmodels.api as sm
import seaborn as sns
from scipy import fftpack
import sciplot as sp
from tqdm import tqdm
#%%
def framean(x=np.array([1, 2, 3]), length=2):
    lengthofrm = x.size - length + 1
    Index = np.arange(lengthofrm)
    frameanlist = [np.sum(x[i:(i+length-1)])/length for i in Index]
    return frameanlist


def od_fitting(time, od, range=(0.025, 0.5), ws=8, trimed=True, Fouier=True):
    '''
    predict growth parameters
    :param time: ndarray
    :param od: ndarray
    :param range: tuple
    :param ws: growth parameter fitting data range length
    :return: predicted od, predicted time, growth parameters [fitting r_square, constant, slope, ]
    '''
    # trime data, od range [0.15, 0.35]
    if trimed:
        mask = [True if (range[0] <= i <= range[1]) else False for i in od]
        tri_od = od[mask]
        tri_time = time[mask]
    else:
        tri_od = od
        tri_time = time
    # smooth via Fouier
    if len(tri_od) >= 5 and Fouier:
        F_od = fftpack.fft(tri_od)
        f_od = fftpack.fftfreq(len(tri_od), 1. / 10000)
        mask_f = np.where(f_od > 0)
        F_od_filtered = F_od * (abs(f_od) < 5000)
        od_filtered = fftpack.ifft(F_od_filtered)
        # interpolate, cubic interploate
        inter_od = interpolate.pchip(tri_time, od_filtered)
        num_interp = 100
        xs = np.linspace(tri_time[0], tri_time[-1], num=num_interp, endpoint=True)
        od_interp = inter_od(xs)
        # gradient
        ln_od_interp = np.log(od_interp)
        ln_od_interp_prime = np.gradient(ln_od_interp, xs)
        # dd
        d_ln_od_prime = np.gradient(ln_od_interp_prime, xs)
        # delete outlier
        mask = np.nonzero(abs(d_ln_od_prime) < 3.6)
        xs = xs[mask]  # for diauxic growth
        ln_od_interp_prime = ln_od_interp_prime[mask]
        # determine the compute frame
        max_od = ln_od_interp_prime.max()
        delta_od = tri_od - np.exp(max_od)
        linear_growth_end = np.abs(delta_od).argmin()
        window_len = ws
        sp = linear_growth_end - window_len + 1
        if sp >= 0:
            pass
        else:
            sp = 0
        window_index = np.arange(start=sp,
                                 stop=linear_growth_end,
                                 step=1, dtype=int)
        # linear regression
        x_pred = tri_time[window_index]
        y_pred = tri_od[window_index]
    else:
        # print(len(tri_od))
        x_pred = tri_time
        y_pred = tri_od
    Xs = sm.add_constant(x_pred)
    Model_Fit = sm.OLS(np.log(y_pred), Xs)
    Fit_Results = Model_Fit.fit()
    dia_growth = None
    print(Fit_Results.params)
    doubleing_time = np.log(2)/Fit_Results.params[1]
    fitting_params = [Fit_Results.rsquared, Fit_Results.params[0], Fit_Results.params[1], doubleing_time, dia_growth]
    od_pred = np.exp(Fit_Results.fittedvalues)
    return od_pred, x_pred, fitting_params


def dataprocess(time, od600, frame_size=10):
    """
    :param time: a list, time
    :param od600: dataframe, OD600
    :return: OD Mean, OD std, OD fitting,
    Fitting Parameters[r_squared, cont, slope, Doubling_Time],
    interpolate time, interpolate OD
    """
    if len(od600.shape) != 1:
        OD_Mean = od600.mean(axis=1)
        OD_Std = od600.std(axis=1)
    else:
        OD_Mean = od600
        OD_Std = None
    IndexofRC = np.where((OD_Mean > 0.02) & (OD_Mean < 0.4))
    x = time[IndexofRC]
    y = OD_Mean.iloc[IndexofRC[0]]
    #yer = OD_Std[IndexofRC[0]]
    xs = np.linspace(x[0], x[-1], 100)
    if type(y) == 'pandas.core.frame.DataFrame':
        curve = interpolate.pchip(x, y.values)
    else:
        curve = interpolate.pchip(x, y)
    ys = curve(xs)
    logys = np.log10(ys)
    dlogys = np.gradient(logys, xs)
    frameandlogys = framean(dlogys, frame_size)
    frameanxs = xs[0:(xs.size - frame_size + 1)]
    bottom_index = np.where(frameandlogys > frameandlogys[0] * 0.96)
    bottom_time = xs[(bottom_index[0][-1] + 100 / x.size * 0.3 * frame_size).astype(int)]
    exp_time = x[np.where(x < bottom_time)]
    exp_growth = y.iloc[np.where(x < bottom_time)[0]]

    if np.max(frameandlogys > frameandlogys[0] * 1.03):
        print("Diauxic growth!")
        dia_growth = True
    else:
        dia_growth = False

    X = exp_time
    X = sm.add_constant(X)
    Y = np.log10(exp_growth)
    Model_Fit = sm.OLS(Y, X)
    Fit_Results = Model_Fit.fit()
    r_squared = Fit_Results.rsquared
    cont = Fit_Results.params[0]

    slope = Fit_Results.params[1]
    y_predict = Fit_Results.fittedvalues
    OD_predict = 10.0**y_predict
    Doubling_Time = np.log(2)/slope
    # print(Fit_Results.summary())
    # print("Doubling time = %.*f min." % (2, Doubling_Time))
    Fitting_Para = [r_squared, cont, slope, Doubling_Time, dia_growth]
    return OD_Mean, OD_Std, OD_predict, exp_time, Fitting_Para, xs, ys, frameanxs, frameandlogys


def time_covert(time, unit='min'):
    '''
    :param time:
    :param unit: 's', 'min', or 'hour'
    :return:
    '''
    if unit == 'min':
        scale = [60., 1., 1/60.]
    elif unit == 'hour':
        scale = [1., 1/60., 1/(60 * 60.)]
    elif unit == 's':
        scale = [60 * 60., 60., 1.]
    sum_of_time = sum([i*j for i, j in zip(map(int, time.split(':')), scale)])
    return sum_of_time

#%%
class PlateReader:
    def __init__(self, name):
        self.name = name
        self.origin_data = None
        self.convered_data = None
        self.blanks_name = []
        self.samples_name = None
        self.corrected_samples = None
        self.samples_data = None
        self.blanks_data = None
        self.time_serial = None
        self.file_path = None
        self.plate_data = None
        self.growth_parameters = pd.DataFrame(columns=['r_square', 'const', 'rate', 'doubling_time', 'diauxic_growth',
                                                       'sample', '#_well'])
        self.predict_growth = pd.DataFrame(columns=['OD', 'time', 'sample', '#_well'])

    def file_load(self, path_of_file, method='l', **kwargs):
        self.file_path = path_of_file
        self.plate_data = pd.read_csv(self.file_path, header=[0, 1])
        convered_data = pd.DataFrame(columns=['OD', 'sample', '#_well'])
        if method == 'nl':
            self.time_serial = self.plate_data['Unnamed: 0_level_0', 'time'].values
            self.time_serial = np.array([time_covert(i, unit=kwargs['unit']) for i in self.time_serial])
        for samples in self.plate_data:
            lenofsample = len(self.plate_data[samples])
            if samples[1] != 'time':
                convered_data = convered_data.append(pd.DataFrame({'OD': list(self.plate_data[samples]),
                                                        'sample': [samples[0]]*lenofsample,
                                                        '#_well': [samples[1]]*lenofsample}),
                                                        ignore_index=True)
            self.origin_data = self.plate_data
            self.convered_data = convered_data

    def find_blank(self, prefix='Blank'):
        sample_set = set(self.convered_data['sample'])
        self.blanks_name = [i for i in list(sample_set) if i[0:len(prefix)] == prefix]

    def find_samples(self):
        sample_name = [i for i in list(set(self.convered_data['sample']))
                       if (i not in self.blanks_name) and (i[:len('Unnamed')] != 'Unnamed')]
        self.samples_name = sample_name

    def correct_data(self, method='map'):
        # samples_name = self.samples_name
        # convered_data = self.convered_data

        sample_mask = [True if i1 in self.samples_name else False for i1 in self.convered_data['sample']]
        self.samples_data = self.convered_data.loc[sample_mask]
        self.wells_name = list(set(self.samples_data['#_well']))
        self.corrected_data = self.samples_data.copy()
        list_of_well = list(set(self.samples_data['#_well']))
        if method == 'map':
            blank_mask = [True if i2 in self.blanks_name else False for i2 in self.convered_data['sample']]
            blanks_data = self.convered_data.loc[blank_mask]
            self.blanks_data = blanks_data

            # find blank in same row
            list_of_row = [i3[0] for i3 in list_of_well]
            mask_of_blank = [[True if j[0] == i4 else False for j in list(blanks_data['#_well'])]
                             for i4 in list_of_row]
            for samples in list_of_well:
                index_of_row = list_of_well.index(samples)
                mapped_blank = blanks_data.loc[mask_of_blank[index_of_row]]
                setofmappedblank = list(set(mapped_blank['#_well']))
                blank_data = np.array([np.array(mapped_blank['OD'].loc[mapped_blank['#_well'] == blank]).reshape(1, -1)
                                       for blank in setofmappedblank])
                mean_blank = blank_data.mean(axis=0).flatten()
                self.corrected_data['OD'].loc[self.corrected_data['#_well'] == samples] = self.corrected_data['OD'].loc[self.corrected_data['#_well'] == samples] - mean_blank
        if method == 'noblank':
            blank_data = [self.corrected_data['OD'].loc[self.corrected_data['#_well'] == i5].min() for i5 in list_of_well]
            for i6, j in enumerate(list_of_well):
                blank = blank_data[i6]
                self.corrected_data['OD'].loc[self.corrected_data['#_well'] == j] -= np.float(blank)

    def compute_growth_rare(self, time_interval=None, **kwargs):
        list_of_well = list(set(self.corrected_data['#_well']))
        times = pd.Series(name='time')
        if time_interval:
            for wells in list_of_well:
                print('compute %s' % wells)
                od = self.corrected_data['OD'].loc[self.corrected_data['#_well'] == wells]
                len_of_od = len(od)
                sample_name = list(set(self.corrected_data['sample'].loc[self.corrected_data['#_well'] == wells]))[0]
                time = np.arange(0, (len_of_od - 1 + 0.01) * time_interval, time_interval)
                OD_predict, exp_time, Fitting_Para = od_fitting(time, od.values, trimed=kwargs['trimed'])
                time_index = od.index
                time = pd.Series(data=time, index=time_index, name='time')
                times = times.append(time)
                self.predict_growth = self.predict_growth.append(pd.DataFrame({
                    'OD': OD_predict,
                    'time': exp_time,
                    'sample': [sample_name] * len(OD_predict),
                    '#_well': [wells] * len(OD_predict)
                }), ignore_index=True)
                self.growth_parameters = self.growth_parameters.append(pd.Series(Fitting_Para + [sample_name, wells],
                                                                                 index=['r_square', 'const', 'rate',
                                                                                        'doubling_time',
                                                                                        'diauxic_growth',
                                                                                        'sample', '#_well']), ignore_index=True)
                self.corrected_data['time'] = times
        else:
            for wells in list_of_well:
                print('compute %s' % wells)
                od = self.corrected_data['OD'].loc[self.corrected_data['#_well'] == wells]
                sample_name = list(set(self.corrected_data['sample'].loc[self.corrected_data['#_well'] == wells]))[0]
                time = self.time_serial
                OD_predict, exp_time, Fitting_Para = od_fitting(time, od.values, trimed=kwargs['trimed'])
                time_index = od.index
                time = pd.Series(data=time, index=time_index, name='time')
                times = times.append(time)
                self.predict_growth = self.predict_growth.append(pd.DataFrame({
                    'OD': OD_predict,
                    'time': exp_time,
                    'sample': [sample_name] * len(OD_predict),
                    '#_well': [wells] * len(OD_predict)
                }), ignore_index=True)
                self.growth_parameters = self.growth_parameters.append(pd.Series(Fitting_Para + [sample_name, wells],
                                                                                 index=['r_square', 'const', 'rate',
                                                                                        'doubling_time', 'diauxic_growth',
                                                                                        'sample', '#_well']), ignore_index=True)
                self.corrected_data['time'] = times
        self.growth_parameters = self.growth_parameters.sort_values(by='doubling_time')  # ascending by doubling time

    def plot_all_curve(self, c_n=4, dir=None):
        sp.ggplot()  # update plot parameters.
        num_col = c_n
        num_row = int(np.ceil(len(self.wells_name) / num_col))
        fig2, ax2 = plt.subplots(num_row, num_col, figsize=(8.2 * num_col, 8.2 * num_row))
        for i in range(len(self.wells_name)):
            y = int(np.floor(i / num_col))
            x = int(i % num_col)
            well_name = self.wells_name[i]
            time = self.predict_growth['time'].loc[self.predict_growth['#_well'] == well_name]
            od = self.predict_growth['OD'].loc[self.predict_growth['#_well'] == well_name]
            sample_name = self.predict_growth['sample'].loc[self.predict_growth['#_well'] == well_name].values[
                0]
            dt = \
            self.growth_parameters['doubling_time'].loc[self.growth_parameters['#_well'] == well_name].values[0]
            rs = self.growth_parameters['r_square'].loc[self.growth_parameters['#_well'] == well_name].values[0]
            if \
            self.growth_parameters['diauxic_growth'].loc[self.growth_parameters['#_well'] == well_name].values[
                0]:
                label_tag = sample_name + '-' + well_name + '\nDoubling time: %.3f' % (dt) + '\nDiauxic Growth!'
            else:
                label_tag = sample_name + '-' + well_name + '\nDoubling time: %.3f' % (dt) + '\nr_square: %.3f ' % (
                    rs)
            ax2[y, x].scatter(self.corrected_data['time'].loc[self.corrected_data['#_well'] == well_name].values,
                              self.corrected_data['OD'].loc[self.corrected_data['#_well'] == well_name].values,
                              c='#38C0E1', marker='^')
            ax2[y, x].plot(time, od, label=label_tag, c='#E15938', alpha=0.6)
            ax2[y, x].set_yscale('log')
            ax2[y, x].set_ylim(min(od.values), max((od.values)) * 1.1)
            ax2[y, x].set_xlim(min(time.values), max(time.values))
            ax2[y, x].legend()
        if dir:
            fig2.savefig(dir + '/all_curves.pdf')
        else:
            pass
        fig2.show()

    def plot_parameters(self, paras='rate', unit='\mathrm{min}^{-1}', dir=None):
        '''
        :param paras: rate or duobling_time
        :param unit: s, min, h, or its LaTex style
        :return:
        '''
        sp.ggplot()
        fig3, ax3 = plt.subplots(1, 2, figsize=(1.5*len(self.wells_name)*2 + 2, 7.5))
        sns.boxplot(x='sample', y=paras, data=self.growth_parameters, ax=ax3[0])
        sns.barplot(x='sample', y=paras, data=self.growth_parameters, ax=ax3[1],
                    capsize=.2, edgecolor='.2', linewidth=3)
        for i in ax3:
            if paras == 'doubling_time':
                i.set_ylabel('Doubling Time ($%s$)' % (unit))
            else:
                i.set_ylabel('Growth Rate ($%s$)' % (unit))
            i.set_xlabel('Sample Name')
        if dir:
            fig3.savefig(dir + f'/parameters_plot_{paras}.pdf')
        else:
            pass
        fig3.show()

    def growth_para_export(self,):
        self.growth_parameters.to_csv(self.file_path + '_' + 'growth_parameters.csv')

#%%

if __name__ == '__main__':

    plate1 = PlateReader('plate_1')
    plate1.file_load(r'./growth_data.csv')
    plate1.find_blank()
    plate1.find_samples()
    plate1.correct_data()
    plate1.compute_growth_rare(time_interval=6/60., trimed=True)

    # process file with time
    plate2 = PlateReader('plate_2')
    plate2.file_load(r'./20191028.csv', method='nl', unit='min')
    plate2.find_blank()
    plate2.find_samples()
    plate2.correct_data()
    plate2.compute_growth_rare(trimed=False)
    plate2.plot_parameters(paras='rate', unit='\mathrm{min}')
    plate2.plot_all_curve()
    #
    sp.ggplot()
    fig1, ax = plt.subplots(1, 2, figsize=(20, 7.5))
    sns.catplot(x='sample', y='rate', data=plate2.growth_parameters, kind='box', ax=ax[0])
    sns.catplot(x='sample', y='rate', data=plate2.growth_parameters, kind='bar', ax=ax[1], capsize=0.2,
                edgecolor='.2', linewidth=3)
    for i in range(len(ax)):
        ax[i].set_ylabel('Rate ($\mathrm{min}^{-1}$)')
        ax[i].set_xlabel('Sample Name')
    fig1.show()
    #
    fig3, ax3 = plt.subplots(1, 2, figsize=(20, 7.5))
    sns.catplot(x='sample', y='doubling_time', data=plate2.growth_parameters, kind='box', ax=ax3[0])
    sns.catplot(x='sample', y='doubling_time', data=plate2.growth_parameters, kind='bar', ax=ax3[1],
                capsize=.2, edgecolor='.2', linewidth=3)
    for i in ax3:
        i.set_ylabel('Doubling Time ($\mathrm{min}$)')
        i.set_xlabel('Sample Name')
    fig3.show()
    #
    num_col = 4
    num_row = int(np.ceil(len(plate2.wells_name)/num_col))
    fig2, ax2 = plt.subplots(num_row, num_col, figsize=(5.2 * num_col, 5.1 * num_row))
    for i in range(len(plate2.wells_name)):
        y = int(np.floor(i / num_col))
        x = int(i % num_col)
        well_name = plate2.wells_name[i]
        time = plate2.predict_growth['time'].loc[plate2.predict_growth['#_well'] == well_name]
        od = plate2.predict_growth['OD'].loc[plate2.predict_growth['#_well'] == well_name]
        sample_name = plate2.predict_growth['sample'].loc[plate2.predict_growth['#_well'] == well_name].values[0]
        dt = plate2.growth_parameters['doubling_time'].loc[plate2.growth_parameters['#_well'] == well_name].values[0]
        rs = plate2.growth_parameters['r_square'].loc[plate2.growth_parameters['#_well'] == well_name].values[0]
        if plate2.growth_parameters['diauxic_growth'].loc[plate2.growth_parameters['#_well'] == well_name].values[0]:
            label_tag = sample_name + '-' + well_name + '\nDoubling time: %.3f' % (dt) +'\nDiauxic Growth!'
        else:
            label_tag = sample_name + '-' + well_name + '\nDoubling time: %.3f' % (dt) + '\nr_square: %.3f ' % (rs)
        ax2[y, x].scatter(plate2.samples_data['time'].loc[plate2.samples_data['#_well'] == well_name].values,
                          plate2.samples_data['OD'].loc[plate2.samples_data['#_well'] == well_name].values,
                          c='#38C0E1', marker='v')
        ax2[y, x].plot(time, od, label=label_tag, c='#E15938', alpha=0.6)
        ax2[y, x].set_yscale('log')
        ax2[y, x].set_xlim(min(time.values), max(time.values))
        ax2[y, x].legend()

    fig2.show()