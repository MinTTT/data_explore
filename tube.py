#%%
import numpy as np
import pandas as pd
import sciplot as sp
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from utiilitys import od_fitting


def conver_data(name, col_name, df):
    col_dex = [col_name.index(name), col_name.index(name) + 1]
    data = df.iloc[1:, col_dex]
    data.columns = ['Time', 'OD']
    data['Sample'] = [name.split('-')[0]] * len(data)
    data['Tube'] = [name] * len(data)
    return data


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
class Flask_Measure:
    def __init__(self, data_path, method=None):
        self.Data_path = data_path
        self.Data_dir = '/'.join(self.Data_path.split('/')[:-1])
        if method == 'raw':
            self.raw_data = pd.read_excel(self.Data_path)
            self.data = self.raw_data.iloc[:, :3]
            self.Tube = list(set(self.data['Tube']))
            self.Tube.sort()
            self.sample_list = list(set([tube.split('-')[0] for tube in self.Tube]))
            self.sample_list.sort()
            self.replicates_dic = {sample: [tube2 for tube2 in self.Tube if tube2.split('-')[0] == sample]
                                   for sample in self.sample_list}
            self.data['Sample'] = [tube.split('-')[0] for tube in self.data['Tube']]
        else:
            self.raw_data = pd.read_excel(self.Data_path)
            self.raw_data_col_name = self.raw_data.columns.to_list()
            self.Tube = [name for name in self.raw_data_col_name if name[0:4] != 'Unna']
            self.Sample = list(set([tube.split('-')[0] for tube in self.Tube]))
            self.data = pd.concat([conver_data(name2, self.raw_data_col_name, self.raw_data) for name2 in self.Tube])

        self.time_unit = None
        self.predict_growth = None
        self.growth_parameters = None
        self.data_trime_range = None
        self.sample_parameters = None


    def datprocess(self, unit, trim_range=(0.02, 0.15)):
        self.data_trime_range = trim_range
        self.time_unit = unit
        self.data['Time'] = [time_process(time, self.time_unit) for time in self.data['Time']]
        self.data['OD'] = self.data['OD'].astype(np.float64)

        # detect transfer index
        self.data = self.data.sort_values(by='Time')
        # detect transfer index (we transfer culture to fresh media before its od reach to 0.2)
        # data['Transfer'] = np.ones(len(data))
        for tube3 in self.Tube:
            tube_data = self.data.loc[self.data['Tube'] == tube3, :]
            diff_1 = np.diff(tube_data['OD'])
            # filter noise delta_od <0.09
            diff_1 = np.array([0.001 if -0.09 < diff < 0 else diff for diff in diff_1])
            diff_1_index = np.nonzero(diff_1 - abs(diff_1))[0]
            transfer_index = np.ones(len(tube_data)) * len(diff_1_index)
            for count, index in enumerate(diff_1_index[-1::-1]):
                transfer_index[:index+1] = len(diff_1_index) - count - 1.
            self.data.loc[self.data['Tube'] == tube3, 'Transfer'] = transfer_index
        self.predict_growth = pd.DataFrame(data=[], columns=self.data.columns)

        gr_paras = []
        for tube3 in self.Tube:
            time = self.data['Time'].loc[self.data['Tube'] == tube3].values
            time = time - time.min()
            self.data['Time'].loc[self.data['Tube'] == tube3] = time.copy()
            # last transfer time and od
            last_fransfer_num = self.data['Transfer'].loc[self.data['Tube'] == tube3].max()
            data = self.data.loc[(self.data['Tube'] == tube3) & (self.data['Transfer'] == last_fransfer_num)]
            od = data['OD'].loc[(data['OD'] <= self.data_trime_range[1]) &
                                (data['OD'] >= self.data_trime_range[0])]
            time = data['Time'].loc[(data['OD'] <= self.data_trime_range[1]) &
                                    (data['OD'] >= self.data_trime_range[0])]
            OD_predict, exp_time, Fitting_Para = od_fitting(time, od, range=(0, 0.24), trimed=False, Fouier=False)
            mask = (self.data['Tube'] == tube3) & (data['OD'] <= self.data_trime_range[1]) & (data['OD'] >= self.data_trime_range[0])
            df_temp = pd.DataFrame(data=[])
            df_temp['OD'] = OD_predict
            df_temp['Time'] = exp_time
            df_temp['Tube'] = [tube3] * len(OD_predict)
            df_temp['Sample'] = [self.data['Sample'].loc[self.data['Tube'] == tube3].values[0]] * len(OD_predict)
            self.predict_growth = self.predict_growth.append([df_temp])
            gr_paras.append([tube3, tube3.split('-')[0]] + Fitting_Para)
        self.growth_parameters = pd.DataFrame(np.array(gr_paras))
        self.growth_parameters.columns = ['Tube', 'Sample', 'r_square', 'const', 'Growth Rate (1/%s)' % self.time_unit,
                                          'Doubling Time (%s)' % self.time_unit, 'dia_growth']
        self.sample_parameters = pd.DataFrame(columns=['Sample',
                                                                'Mean Growth Rate (1/%s)' % self.time_unit,
                                                                'STD of Mean Growth Rate (1/%s)' % self.time_unit,
                                                                'Mean Doubling Time (%s)' % self.time_unit,
                                                                'STD of Doubling Time (%s)' % self.time_unit]
                                              )
        for inx, sample in enumerate(self.sample_list):
            mask = self.growth_parameters['Sample'] == sample
            mean_gr = np.mean(self.growth_parameters['Growth Rate (1/%s)' % self.time_unit].loc[mask])
            std_gr = np.std(self.growth_parameters['Growth Rate (1/%s)' % self.time_unit].loc[mask])
            mean_dt = np.mean(self.growth_parameters['Doubling Time (%s)' % self.time_unit].loc[mask])
            std_dt = np.mean(self.growth_parameters['Doubling Time (%s)' % self.time_unit].loc[mask])
            row = [sample, mean_gr, std_gr, mean_dt, std_dt]
            self.sample_parameters.loc[inx] = row

    def plot_all_curve(self, c_n=4, dir=None):
        sp.ggplot()  # update plot parameters.
        num_col = c_n
        num_row = int(np.ceil(len(self.Tube) / num_col))
        fig2, ax2 = plt.subplots(num_row, num_col, figsize=(9 * num_col, 9 * num_row))
        ytick = np.logspace(np.log10(self.data['OD'].min() * 0.92), np.log10(0.2), num=6)
        for i in range(len(self.Tube)):
            paras_col = self.growth_parameters.columns.to_list()
            # print(paras_col)
            y = int(np.floor(i / num_col))
            x = int(i % num_col)
            well_name = self.Tube[i]
            last_fransfer_num = self.data['Transfer'].loc[self.data['Tube'] == well_name].max()
            time = self.predict_growth['Time'].loc[(self.predict_growth['Tube'] == well_name)]
            od = self.predict_growth['OD'].loc[(self.predict_growth['Tube'] == well_name)]
            sample_name = self.predict_growth['Sample'].loc[self.predict_growth['Tube'] == well_name].values[
                0]
            dt = \
            self.growth_parameters[paras_col[5]].loc[self.growth_parameters['Tube'] == well_name].values[0]
            rs = self.growth_parameters[paras_col[2]].loc[self.growth_parameters['Tube'] == well_name].values[0]
            if \
            self.growth_parameters[paras_col[-1]].loc[self.growth_parameters['Tube'] == well_name].values[
                0]:
                label_tag = sample_name + ':' + well_name + '\nDoubling time: %.3f %s' % (dt, self.time_unit) + '\nDiauxic Growth!'
            else:
                label_tag = sample_name + ':' + well_name + '\nDoubling time: %.3f %s' % (dt, self.time_unit) + '\nr_square: %.3f ' % (
                    rs)
            ax2[y, x].scatter(self.data['Time'].loc[(self.data['Tube'] == well_name) &
                                                   (self.data['Transfer'] == last_fransfer_num)].values,
                              self.data['OD'].loc[(self.data['Tube'] == well_name) &
                                                   (self.data['Transfer'] == last_fransfer_num)].values,
                              c='#38C0E1', marker='^')
            ax2[y, x].plot(time, od, label=label_tag, c='#E15938', alpha=0.6)
            ax2[y, x].set_yscale('log')
            ax2[y, x].set_yticks(ytick)
            ax2[y, x].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax2[y, x].set_ylim(ytick.min(), ytick.max() * 1.1)
            ax2[y, x].set_xlim(min(time.values) - (max(time.values) - min(time.values)) * 0.08, max(time.values) * 1.1)
            ax2[y, x].legend()
        if dir:
            fig2.savefig(dir+'/all_curves')
        else:
            fig2.show()

    def plot_parameters(self, paras='rate', dir=None):
        '''
        :param paras: rate or duobling_time
        :return:
        '''
        paras_col = self.growth_parameters.columns.to_list()
        if paras == 'rate':
            y_data = paras_col[4]
        elif paras == 'doubling_time':
            y_data = paras_col[5]
        else:
            print('paras = rate or doubling_time')
        sp.ggplot()
        fig3, ax3 = plt.subplots(1, 2, figsize=(20, 7.5))
        sns.boxplot(x='Sample', y=y_data, data=self.growth_parameters, ax=ax3[0])
        sns.barplot(x='Sample', y=y_data, data=self.growth_parameters, ax=ax3[1],
                    capsize=.2, edgecolor='.2', linewidth=3)
        for i in ax3:
            if paras == 'doubling_time':
                i.set_ylabel('Doubling Time ($\mathrm{%s}$)' % self.time_unit)
            else:
                i.set_ylabel('Growth Rate ($\mathrm{%s}^{-1}$)' % self.time_unit)
            i.set_xlabel('Sample Name')
        if dir:
            fig3.savefig(dir+'/gr_parameters_%s' % paras)
        else:
            fig3.show()

    def growth_para_export(self,):
        self.growth_parameters.to_csv(self.Data_path + 'growth_parameters.csv')
        self.sample_parameters.to_csv(self.Data_path + 'sample_parameters.csv')

    def growth_data_export(self):
        self.data.to_csv(self.Data_path + 'growth_data.csv')

    def do_all_process(self, unit):
        # compute, plot, export all in one!
        self.datprocess(unit)
        self.plot_parameters(paras='rate', dir=self.Data_dir)
        self.plot_parameters(paras='doubling_time', dir=self.Data_dir)
        self.plot_all_curve(dir=self.Data_dir)
        self.growth_para_export()
        self.growth_data_export()
        print('Successful! See files in directory: %s' % self.Data_dir)


if __name__ == '__main__':
    exp1 = Flask_Measure(r'./data/20200314/flask.xlsx')
    exp1.do_all_process(unit='min')

