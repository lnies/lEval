import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize
from datetime import datetime
import sys, os
import glob
import json#simplejson as json
import csv
import re
import ROOT
import lmfit
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)


class Freq_ratio():
    # class to calculate the frequency ratio for a mass measurement (MM), optimized for ISOLTRAP PI-ICR data
    # cite LMFIT via https://dx.doi.org/10.5281/zenodo.11813
    def __init__(self):
        self.degree = 3
        self.param_start = [1]*(self.degree+2)
        self.param_result = []
        self.y1df = {}
        self.y1names = []
        self.y1 = []
        self.y1err = []
        self.y1folder = ''      # cyc freq isotope 1
        self.y2df = {}
        self.y2names = []
        self.y2 = []
        self.y2err = []
        self.y2folder = ''      # cyc freq isotope 2
        self.x1df = {}
        self.x1helpdf = {}
        self.x1helpfolder = ''
        self.x1names = []
        self.x1 = []
        self.x1folder = ''      # time stamps isotope 1
        self.x2df = {}
        self.x2helpdf = {}
        self.x2helpfolder = ''
        self.x2names = []
        self.x2 = []
        self.x2folder = ''      # time stamps isotope 2
        self.isotopes = []
        self.dirname = ''
        self.run_folder = ''
        self.ratio = 0
        self.ratio_unc = 0
        self.red_chi_sq = 0
        self.covar = []
        self.correlation = 0
        self.ratio_dict = {}
        self.weighted_avg = 0
        self.total_unc = 0
        self.birge_ratio = 0
        self.merged_dict = {}
        self.lin_extrapolation = False
        self.mode = ''

        self.test = False


    def batch(self, upper_run_folder, isotopes, degree=3, z_classes=['Z1-1']):
        for z in z_classes:
            for run_folder in glob.glob(upper_run_folder+'{}/'.format(z)+'*/'):
                if 'run' in run_folder:
                    self.main(run_folder=run_folder, isotopes=isotopes, degree=degree, single_or_batch='batch', z_classes=z)
        # input_list = [['ratio', 'ratio_unc']]
        # for run in self.ratio_dict:
        #     input_list.append([self.ratio_dict[run]['ratio'], self.ratio_dict[run]['ratio_unc']])
        # self.weighted_avg, self.total_unc, self.birge_ratio = self.weighted_avg_calc(input_list)
        # print '\n\nWeighted average:      ', self.weighted_avg,'+/-', self.total_unc
        # self.ratio_dict['weighted_avg'] = self.weighted_avg
        # self.ratio_dict['weighted_avg_unc'] = self.total_unc
        # self.ratio_dict['birge_ratio'] = self.birge_ratio
        filename = '{}ratio_dict_{}_{}_deg{}_batch.json'.format(upper_run_folder,
                                                          isotopes[0], isotopes[1],
                                                          degree)
        with open(filename, 'w') as f:
            json.dump(self.ratio_dict, f, sort_keys=True, indent=4)


    def main(self, isotopes, run_folder, degree=3, single_or_batch='single', z_classes='Z1-1'):
        self.degree = degree
        self.get_start_param()
        if single_or_batch == 'single_easy':
            self.get_data_tof_icr(run_folder, isotopes)
            self.mode = 'single_easy'
        else:
            self.set_folders(run_folder, isotopes)
            self.get_data()

        if len(self.y2) > 1:    # make sure there are at least two points to fit
            self.fit()
            self.red_chi_sq_calc()
            print '\nRed.Chi.Sq.: ', self.red_chi_sq
            if self.covar is None:
                self.covar = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
            self.covar = self.covar * self.red_chi_sq
            self.param_result_err = np.array([np.sqrt(self.covar[i, i]) for i in range(len(self.param_result))])
            self.correlation = self.covar[len(self.param_result)-2, len(self.param_result)-1]       # calculates the corellation
        else:
            self.lin_extrapolation = True
            self.fit()
            self.param_result_err = np.sqrt(np.diag(self.covar))
        self.ratio_calc()
        print 'Frequency ratio for {} / {} :'.format(self.isotopes[0], self.isotopes[1]), self.ratio
        self.plotting()
        if z_classes not in self.ratio_dict.keys():
            self.ratio_dict[z_classes] = {}
        self.ratio_dict[z_classes]['{}'.format(self.dirname)] = {'red_chi_sq': self.red_chi_sq,
                                                                 'isotopes': self.isotopes,
                                                                 'ratio': self.ratio,
                                                                 'ratio_unc': self.ratio_unc,
                                                                 'degree': self.degree}
        if single_or_batch == 'single':
            self.save_update(z_classes)
        self.lin_extrapolation = False


    def fill_tree(self):
        f = ROOT.TFile.Open('{}{}-{}.root'.format(self.run_folder, self.isotopes[0], self.isotopes[1]), 'recreate')
        tree_name = '{}{}'.format(self.isotopes[0], self.isotopes[1])
        tree = ROOT.TTree('{}{}'.format(self.isotopes[0], self.isotopes[1]), '{}{}'.format(self.isotopes[0], self.isotopes[1]))
        tmp_data = np.empty((3,1), dtype=float)
        x1_list = self.x1
        y1_list = self.y1
        y1err_list = self.y1err
        # x2_list = self.x2
        # y2_list = self.y2
        # y2err_list = self.y2err
        tree.Branch('data', tmp_data, 'x1/D:y1/D:y1err/D')
        # tree.Branch('data', tmp_data, 'x2/D:y2/D:y2err/D')
        for i in range(len(x1_list)):
            tmp_data[0][0] = x1_list[i]
            tmp_data[1][0] = y1_list[i]
            tmp_data[2][0] = y1err_list[i]
            # tmp_data[0][0] = x2_list[i]
            # tmp_data[1][0] = y2_list[i]
            # tmp_data[2][0] = y2err_list[i]
            tree.Fill()
        tree.Print()
        tree.Write()
        f.Close()


    def merge_single(self, upper_z_folder, isotopes):
        # function to merge all results calculated in single mode
        z_classes = [i for i in os.listdir(upper_z_folder) if 'Z' in i and len(i)<5]
        z_classes.remove('Z0-0')
        for z in z_classes:
            self.merged_dict[z] = {}
            tmp_runs = [i for i in os.listdir(upper_z_folder+'{}/'.format(z)) if 'run' in i]
            for run in tmp_runs:
                filename = '{}{}/{}/ratio_dict_{}_{}.json'.format(upper_z_folder,
                                                                   z,run,isotopes[0],
                                                                   isotopes[1])
                print filename
                try:
                    with open(filename) as json_data:
                        tmp_dict = json.load(json_data)
                    self.merged_dict[z][run] = tmp_dict[z][run]
                except:
                    print('File note found')
        with open('{}/merged_ratio_dict_{}_{}.json'.format(upper_z_folder, isotopes[0], isotopes[1]), 'w') as f:
            json.dump(self.merged_dict, f, sort_keys=True, indent=4)


    def set_folders(self, run_folder, isotopes):
        self.isotopes = isotopes
        self.dirname = os.path.basename(os.path.normpath(run_folder))
        self.run_folder = run_folder
        self.x1folder = run_folder+isotopes[0]+'/p1p2/_time_info_p1p2.csv'
        self.x2folder = run_folder+isotopes[1]+'/p1p2/_time_info_p1p2.csv'
        self.x1helpfolder = run_folder+isotopes[0]+'/p1p2/'
        self.x2helpfolder = run_folder+isotopes[1]+'/p1p2/'
        self.y1folder = run_folder+isotopes[0]+'/piicr_cyc_freq_'+isotopes[0]+'.csv'
        self.y2folder = run_folder+isotopes[1]+'/piicr_cyc_freq_'+isotopes[1]+'.csv'#'_ground.csv'

    def get_start_param(self):
        if self.test:
            self.param_start = lmfit.Parameters()
            for i in range(self.degree+2):
                self.param_start.add('param{}'.format(i), value=1, max=1, min=-1)
            self.param_start['param1'].max = 5000000
            # self.param_start.append(1)
        else:
            self.param_start = [1]*(self.degree+2)

    def get_data_tof_icr(self, run_folder, isotopes):
        # bla
        self.isotopes = isotopes
        self.run_folder = run_folder

        self.y1df = pd.read_csv(run_folder+isotopes[0]+'.csv',
                                names=['frequency', 'unc', 'median_time'],
                                parse_dates=[2])
        self.y1 = self.y1df.iloc[:, 0].tolist()
        self.y1err = self.y1df.iloc[:, 1].tolist()

        self.x1df = self.y1df

        self.y2df = pd.read_csv(run_folder+isotopes[1]+'.csv',
                                names=['frequency', 'unc', 'median_time'],
                                parse_dates=[2])
        self.y2 = self.y2df.iloc[:, 0].tolist()
        self.y2err = self.y2df.iloc[:, 1].tolist()

        self.x2df = self.y2df


        self.x1 = [round(i.total_seconds()/60, 2) for i in (self.x1df['median_time'] - min([self.x1df['median_time'][0], self.x2df['median_time'][0]])).tolist()]
        self.x2 = [round(i.total_seconds()/60, 2) for i in (self.x2df['median_time'] - min([self.x1df['median_time'][0], self.x2df['median_time'][0]])).tolist()]


    def get_data(self):
        # get frequency and timestamps
        # ___________________________________________
        # get frequency data for isotope 1
        self.y1df = pd.read_csv(self.y1folder)
        self.y1names = self.y1df.iloc[:, 0].tolist()
        self.y1 = self.y1df.iloc[:, 1].tolist()
        self.y1err = self.y1df.iloc[:, 2].tolist()
        self.y1df = self.y1df.sort_values('File')
        self.y1df = self.y1df.reset_index(drop=True)

        # ___________________________________________
        # get frequency data for isotope 2
        self.y2df = pd.read_csv(self.y2folder)
        self.y2names = self.y2df.iloc[:, 0].tolist()
        self.y2 = self.y2df.iloc[:, 1].tolist()
        self.y2err = self.y2df.iloc[:, 2].tolist()
        self.y2df = self.y2df.sort_values('File')
        self.y2df = self.y2df.reset_index(drop=True)

        # ___________________________________________
        # get time stamps for each file for isotope 1
        self.x1df = pd.read_csv(self.x1folder,
                                names=['names', 'date_begin', 'time_begin', 'date_end', 'time_end'],
                                parse_dates = [[1,2], [3,4]])
        self.x1df = self.x1df.sort_values('names')
        self.x1df = self.x1df.reset_index(drop=True)
        self.x1df['timedelta'] = self.x1df['date_end_time_end'] - self.x1df['date_begin_time_begin']
        self.x1names = self.x1df['names'].tolist()
        # calculate the median time stamp for each file also considering the count rate
        self.x1df['median_counts'] = self.get_timestamp_position(self.x1helpfolder, self.x1names)
        self.x1df['median_time'] = (self.x1df['date_begin_time_begin'] +
                                    self.x1df['timedelta'] * self.x1df['median_counts'])


        # ___________________________________________
        # get time stamps for each file for isotope 2
        self.x2df = pd.read_csv(self.x2folder,
                                names=['names', 'date_begin', 'time_begin', 'date_end', 'time_end'],
                                parse_dates = [[1,2], [3,4]])
        self.x2df = self.x2df.sort_values('names')
        self.x2df = self.x2df.reset_index(drop=True)

        self.x2df['timedelta'] = self.x2df['date_end_time_end'] - self.x2df['date_begin_time_begin']
        self.x2names = self.x2df['names'].tolist()
        # calculate the median time stamp for each file also considering the count rate
        self.x2df['median_counts'] = self.get_timestamp_position(self.x2helpfolder, self.x2names)
        self.x2df['median_time'] = (self.x2df['date_begin_time_begin'] +
                                    self.x2df['timedelta'] * self.x2df['median_counts'])
        # the x values for the frequency information is converted into minutes after(/before) the first
        # x1 frequency point started. This had to be changed since the fit wasn't very stable for x values
        # having the correct time converted to seconds since the values were very large (e.g. 1498777404).
        # The following statement reduces the x value to minutes in the 2-3 digit range. Important here is,
        # that x1 and x2 are normalized to the same starting time of MIN(x1[0], x2[0])
        self.x1 = [round(i.total_seconds()/60, 2) for i in (self.x1df['median_time'] - min([self.x1df['median_time'][0], self.x2df['median_time'][0]])).tolist()]
        self.x2 = [round(i.total_seconds()/60, 2) for i in (self.x2df['median_time'] - min([self.x1df['median_time'][0], self.x2df['median_time'][0]])).tolist()]


    def get_timestamp_position(self, folder, x1_or_x2_names):
        # calculates the mean position of ions to arrive during the measurement time
        # --> if beam gate was changed during the run, the average timestamp is not in the middle of the measurement
        median_counts = []
        for i in x1_or_x2_names:
            tmp = folder+i+'_p1_spot_positions.csv'
            self.x1helpdf = pd.read_csv(tmp, names=['event', 'x', 'y', 'time_stamp_us'])
            # we chose the median over the mean to be less sensitive to outliers in case of a drop in count rate
            median_counts.append(float(self.x1helpdf.iloc[:, 0].median() / self.x1helpdf.iloc[:, 0].max()))
        return(median_counts)


    def fit(self):

        if self.lin_extrapolation == False:
            self.fill_tree()

            if self.test:
                out = lmfit.minimize(self.simultaneous_residual, self.param_start,
                                     args=(np.array(self.x1),
                                           np.array(self.x2),
                                           np.array(self.y1),
                                           np.array(self.y2),
                                           np.array(self.y1err),
                                           np.array(self.y2err)))
                print out.params.pretty_print()
                print out.status.pretty_print()
                # print out.covar
                print out.redchi#.pretty_print()
                # print 'aic'
                # print out.aic.pretty_print()

                sys.exit('test')
            else:
                # fit ; timestampt value is reduced since resulting numbers were too large
                iter1, covar, infodict, mesg, ier = scipy.optimize.leastsq(self.simultaneous_residual, self.param_start,
                                                                args=(np.array(self.x1),
                                                                      np.array(self.x2),
                                                                      np.array(self.y1),
                                                                      np.array(self.y2),
                                                                      np.array(self.y1err),
                                                                      np.array(self.y2err)), full_output=True)#, maxfev=1000000
                iter2, covar, infodict, mesg, ier = scipy.optimize.leastsq(self.simultaneous_residual, iter1,
                                                                args=(np.array(self.x1),
                                                                      np.array(self.x2),
                                                                      np.array(self.y1),
                                                                      np.array(self.y2),
                                                                      np.array(self.y1err),
                                                                      np.array(self.y2err)), full_output=True)
                iter3, covar, infodict, mesg, ier = scipy.optimize.leastsq(self.simultaneous_residual, iter2,
                                                                args=(np.array(self.x1),
                                                                      np.array(self.x2),
                                                                      np.array(self.y1),
                                                                      np.array(self.y2),
                                                                      np.array(self.y1err),
                                                                      np.array(self.y2err)), full_output=True)
                iter4, covar, infodict, mesg, ier = scipy.optimize.leastsq(self.simultaneous_residual, iter3,
                                                                args=(np.array(self.x1),
                                                                      np.array(self.x2),
                                                                      np.array(self.y1),
                                                                      np.array(self.y2),
                                                                      np.array(self.y1err),
                                                                      np.array(self.y2err)), full_output=True)
                self.param_result, self.covar, infodict, mesg, ier = scipy.optimize.leastsq(self.simultaneous_residual, iter4,
                                                                args=(np.array(self.x1),
                                                                      np.array(self.x2),
                                                                      np.array(self.y1),
                                                                      np.array(self.y2),
                                                                      np.array(self.y1err),
                                                                      np.array(self.y2err)), full_output=True)
                if ier not in [1,2,3,4]:
                    sys.exit('Fit did not converge')

                # print self.param_result
                # print self.covar

        else:
            init_vals = [((max(self.y1)-min(self.y1))
                          /(max(self.x1)-min(self.x1))),self.y1[0]]
            fit_1, pcov = scipy.optimize.curve_fit(self.linear,
                                       self.x1, self.y1, sigma=self.y1err,
                                       absolute_sigma=True,
                                       p0=init_vals)
            fit_2, pcov = scipy.optimize.curve_fit(self.linear,
                                       self.x1, self.y1, sigma=self.y1err,
                                       absolute_sigma=True,
                                       p0=fit_1)
            fit_3, pcov = scipy.optimize.curve_fit(self.linear,
                                       self.x1, self.y1, sigma=self.y1err,
                                       absolute_sigma=True,
                                       p0=fit_2)
            fit_4, pcov = scipy.optimize.curve_fit(self.linear,
                                       self.x1, self.y1, sigma=self.y1err,
                                       absolute_sigma=True,
                                       p0=fit_3)
            self.param_result, self.covar = scipy.optimize.curve_fit(self.linear,
                                       self.x1, self.y1, sigma=self.y1err,
                                       absolute_sigma=True,
                                       p0=fit_4)


    def fit_function(self, x, param):
        # defines the polinomial function
        return(np.polyval(list(param), x))  # returns a polynomial function of degree len(param)-1 evaluated at x


    def residual(self, p, x, y, yerr):
        # returns residual, which square is minimized later
        return (self.fit_function(x, p) - y )/yerr


    def simultaneous_residual(self, param_start, x1, x2, y1, y2, y1err, y2err):
        # calculate residuals for both datasets
        # print len(param_start),'\n\n\n', [len(param_start) - 1 - j for j in range(len(param_start))]
        if self.test:
            p_s = [param_start['param{}'.format(i)] for i in [len(param_start) - 1 - j for j in range(len(param_start))]]
        else:
            p_s = param_start
        p1 = p_s[:-1]  # Poly1 = R*Poly2
        p2 = [float(i)/p_s[-1] for i in p_s[:-1]]  # Poly1 = R*Poly2
        res1 = self.residual(p1, x1, y1, y1err)
        res2 = self.residual(p2, x2, y2, y2err)
        return np.concatenate((res1, res2))


    def final_residual_1(self, param_start, x1, y1, y1err):
        if self.test:
            p_s = [param_start['param{}'.format(i)] for i in [len(param_start) - 1 - j for j in range(len(param_start))]]
        else:
            p_s = param_start
        p1 = p_s[:-1]  # Poly1 = R*Poly2
        res1 = self.residual(p1, x1, y1, y1err)
        return(res1)


    def final_residual_2(self, param_start, x2, y2, y2err):
        if self.test:
            p_s = [param_start['param{}'.format(i)] for i in [len(param_start) - 1 - j for j in range(len(param_start))]]
        else:
            p_s = param_start
        p2 = [float(i)/p_s[-1] for i in p_s[:-1]]  # Poly1 = R*Poly2
        res2 = self.residual(p2, x2, y2, y2err)
        return(res2)


    def red_chi_sq_calc(self):
        sum_res_sq_1 = 0
        sum_res_sq_2 = 0
        for i in range(len(self.x1)):
            sum_res_sq_1 += (self.final_residual_1(self.param_result, self.x1[i], self.y1[i], self.y1err[i]))**2
        for i in range(len(self.x2)):
            sum_res_sq_2 += (self.final_residual_2(self.param_result, self.x2[i], self.y2[i], self.y2err[i]))**2
        self.red_chi_sq = (sum_res_sq_1+sum_res_sq_2)/(len(self.x1)+len(self.x2)-self.degree-1)


    def ratio_calc(self):
        # calculates the ratio in the middle of the range of interest
        if self.lin_extrapolation == False:
            self.ratio = self.param_result[-1]

            self.ratio_unc = self.param_result_err[-1]

            print 'Degree:              ', self.degree
        else:
            fit_at_y2 = self.linear(self.x2[0], self.param_result[0], self.param_result[1])
            self.ratio = fit_at_y2 / self.y2[0]

            unc_lin = np.sqrt( ( self.param_result_err[0] * self.x2[0] )**2
                               + self.param_result_err[1] ** 2 )
            self.ratio_unc = np.sqrt( ( unc_lin / self.y2[0] )**2
                                      + ( -( fit_at_y2 * self.y2err[0] )
                                          / (self.y2[0])**2 )**2 )
            self.degree = 'linear'
            print 'Degree:              ', self.degree

        print 'Ratio fit parameter: ', self.ratio, '+/-', self.ratio_unc

    def plotting(self):
        # to plot with the real time step instead of the time after the first x1 frequency:
        x1 = self.x1df['median_time'].dt.to_pydatetime().tolist()
        x2 = self.x2df['median_time'].dt.to_pydatetime().tolist()

        mpl.rc('font', family='serif', serif='Utopia')      # Utopia LaTeX font!!
        mpl.rc('text', usetex=False)

        f, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(9,6))
        plt.xlabel('time', fontsize=22, fontweight='bold')
        plt.ylabel('     ', fontsize=22, fontweight='bold')        # makes space for the actual label (next line...)
        f.text(0.03,0.55, 'frequency (Hz)', ha='center', va='center', rotation=90, fontsize=22, fontweight='bold')

        l1 = ax.errorbar(x1, self.y1, yerr=self.y1err, fmt='o', label='{}'.format(self.isotopes[0]))
        l2 = ax2.errorbar(x2, self.y2, yerr=self.y2err, fmt='g o', label='{}'.format(self.isotopes[1]))
        # save data to file
        # data1 = [[self.y1[i], self.y1err[i], x1[i].strftime('%m/%d/%y %H:%M')] for i in range(len(x1))]
        # data2 = [[self.y2[i], self.y2err[i], x2[i].strftime('%m/%d/%y %H:%M')] for i in range(len(x2))]
        # with open('{}{}.csv'.format(self.run_folder, self.isotopes[0]), 'wb') as f:
        #     writer = csv.writer(f)
        #     writer.writerows(data1)
        # with open('{}{}.csv'.format(self.run_folder, self.isotopes[1]), 'wb') as f:
        #     writer = csv.writer(f)
        #     writer.writerows(data2)
        # set range with same delta for both
        big_range = max([max(self.y1)-min(self.y1), max(self.y2)-min(self.y2)])
        max_unc_1 = max(self.y1err)
        max_unc_2 = max(self.y2err)
        big_range_extended = big_range + 0.2*big_range + 2* max([max_unc_1, max_unc_2])
        middle_y1 = min(self.y1) + (max(self.y1) - min(self.y1))/2
        middle_y2 = min(self.y2) + (max(self.y2) - min(self.y2))/2
        ax.set_ylim(middle_y1 - big_range_extended/2, middle_y1 + big_range_extended/2)  # outliers only
        ax2.set_ylim(middle_y2 - big_range_extended/2, middle_y2 + big_range_extended/2)  # most of the data
        # hide the spines between ax and ax2
        ax.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax.xaxis.tick_top()
        ax.tick_params(labeltop='off')  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()
        # put diagonal lines on axis
        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

        # fit
        x_range = max([max(self.x1), max(self.x2)]) - min([min(self.x1), min(self.x2)])
        x = np.linspace(min([min(self.x1), min(self.x2)]) - 0.05*x_range,
                        max([max(self.x1), max(self.x2)]) + 0.05*x_range, 200)
        time_x = [(pd.Timedelta(i, unit='m')+min([self.x1df['median_time'][0], self.x2df['median_time'][0]])).to_pydatetime(warn=False) for i in x]

        if self.lin_extrapolation == False:
            parameter1 = [float(i) for i in list(self.param_result[:-1])]
            parameter2 = [float(i)/float(self.param_result[-1]) for i in list(self.param_result[:-1])]
            l3, = ax.plot(time_x, self.fit_function(np.array(x), parameter1), 'r-', label='fit')   # the star is to unpack the parameter array
            ax2.plot(time_x, self.fit_function(np.array(x), parameter2), 'r-', label='fit')   # the star is to unpack the parameter array
            # ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d, %H'))
            ax.yaxis.set_minor_locator(MultipleLocator(0.1))
            ax.yaxis.set_major_locator(MultipleLocator(0.5))
            ax.yaxis.offsetText.set_fontsize(15)
            ax.tick_params(labelsize=15, top=False, right=False)


            ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
            ax2.yaxis.set_major_locator(MultipleLocator(0.5))
            ax2.yaxis.offsetText.set_fontsize(15)
            ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
            ax2.xaxis.set_minor_locator(mpl.dates.HourLocator())
            ax2.tick_params(labelsize=15, top=False, right=False)
        else:
            l3, = ax.plot([(pd.Timedelta(i, unit='m')
                           +self.x1df['median_time'][0]).to_pydatetime(warn=False) for i in x]
                           , self.linear(np.array(x), self.param_result[0], self.param_result[1]), 'r-', label='fit')   # the star is to unpack the parameter array

        # sys.exit()
        plt.tick_params(labelsize=15, right=False)

        if self.isotopes == ['23Na', '23Mg']:
            self.isotopes = ['$^{23}$Na', '$^{23}$Mg']


        f.legend((l1, l2, l3),('{}'.format(self.isotopes[0]),
                                    '{}'.format(self.isotopes[1]), 'fit')
                      ,fontsize=14,title='{}'.format(self.dirname), bbox_to_anchor=(0.255,0.33))

        plt.tight_layout()


        if self.isotopes == ['$^{23}$Na', '$^{23}$Mg1']:
            self.isotopes = ['23Na', '23Mg']
            anchored_text = mpl.offsetbox.AnchoredText('Ratio: {:1.10f} $\\pm$ {:1.10f}'.format(self.ratio, self.ratio_unc, self.red_chi_sq, self.degree), loc=9)#, fontsize=14)
        else:
            anchored_text = mpl.offsetbox.AnchoredText('Ratio: {:1.10f} $\\pm$ {:1.10f} ($\\chi^2_{{red}}$ = {:1.2f}, polynomial degree = {})'.format(self.ratio, self.ratio_unc, self.red_chi_sq, self.degree), loc='upper right')
        ax.add_artist(anchored_text)

        plt.savefig('{}ratio_deg{}_{}{}_{}{}.pdf'.format(self.run_folder,self.degree, [str(s) for s in re.findall(r'-?\d+\.?\d*', self.isotopes[0])][0], ''.join(x for x in self.isotopes[0] if x.isalpha()), [str(s) for s in re.findall(r'-?\d+\.?\d*', self.isotopes[1])][0], ''.join(x for x in self.isotopes[1] if x.isalpha())))

        # plt.show()
        plt.close()

    def weighted_avg_calc(self, input_list):
        '''
        Takes a list with [[header, header_unc], [x, x_unc.],...] and calculates the weighted mean, inner and outer unc, total unc and the birge ratio.

        The ''internal'' uncertainty includes only measurement uncertainty, that is, uncertainty on N, in computing the uncertainty on the exposure age t. The internal uncertainty could be sometimes called the ''measurement uncertainty,'' although that is less clear and could mean several different things. The ''external'' uncertainty includes all input uncertainties - uncertainties on both N and P in the example above - and is also often called the ''total uncertainty.''

        taken from: https://cosmognosis.wordpress.com/2013/10/25/internal-vs-external-uncertainties-again/
        '''
        weighted_avg = (float(np.sum([float(x[0])/(float(x[1]) ** 2) for x in input_list[1:]]))
                             / float(np.sum([1/(float(x[1]) ** 2) for x in input_list[1:]])))
        internal_unc = float(np.sqrt(float(1/np.sum([1/(float(x[1]) ** 2) for x in input_list[1:]]))))
        external_unc = float(np.sqrt(np.sum([ (float(x[0]) - weighted_avg) ** 2
                                            / (float(x[1]) ** 2) for x in input_list[1:]])
                                            / (np.sum([1/(float(x[1]) ** 2) * (len(input_list)-2)
                                                      for x in input_list[1:]]))))
        total_unc = max([internal_unc, external_unc])
        birge_ratio = external_unc / internal_unc
        return(weighted_avg, total_unc, birge_ratio)


    def save_update(self, z_classes):
        # updates an existing json file with the new calculation
        filename = '{}ratio_dict_{}_{}.json'.format(self.run_folder,
                                                    self.isotopes[0],
                                                    self.isotopes[1])
        if not os.path.exists(filename):
            with open(filename, 'w') as f:
                json.dump(self.ratio_dict, f, sort_keys=True, indent=4)
        else:
            with open(filename) as json_data:
                tmp_dict = json.load(json_data)

            if z_classes not in tmp_dict.keys():
                tmp_dict[z_classes] = {}
            tmp_dict[z_classes]['{}'.format(self.dirname)] = self.ratio_dict[z_classes]['{}'.format(self.dirname)]
            with open(filename, 'w') as f:
                json.dump(tmp_dict, f, sort_keys=True, indent=4)

    def linear(self, x, m, n):
        return(m*x+n)

    def extrapolation(self, merged_dict_file):
        # function to perform the extrapolation to Z=0
        with open(merged_dict_file) as json_data:
            tmp_dict = json.load(json_data)#, parse_float=Decimal)
        z_classes = [str(i) for i in tmp_dict.keys() if str(i) != 'Z0-0']
        if len(z_classes) == 2:
            runs1 = [str(i) for i in tmp_dict[z_classes[0]].keys()]
            runs2 = [str(i) for i in tmp_dict[z_classes[1]].keys()]
            runs = list(set(runs1 + runs2))
        elif len(z_classes) == 3:
            runs1 = [str(i) for i in tmp_dict[z_classes[0]].keys()]
            runs2 = [str(i) for i in tmp_dict[z_classes[1]].keys()]
            runs3 = [str(i) for i in tmp_dict[z_classes[2]].keys()]
            runs = list(set(runs1 + runs2 + runs3))
        elif len(z_classes) == 4:
            runs1 = [str(i) for i in tmp_dict[z_classes[0]].keys()]
            runs2 = [str(i) for i in tmp_dict[z_classes[1]].keys()]
            runs3 = [str(i) for i in tmp_dict[z_classes[2]].keys()]
            runs4 = [str(i) for i in tmp_dict[z_classes[3]].keys()]
            runs = list(set(runs1 + runs2 + runs3 + runs4))
        isotopes = [str(i) for i in tmp_dict[z_classes[0]][runs[0]]['isotopes']]
        tmp_dict['Z0-0'] = {}
        tmp_dict['Z0-0']['extrapolated_list'] = [['run#', 'extr_ratio', 'unc']]
        folder = os.path.dirname(merged_dict_file)+'/Z0-0/'
        if not os.path.exists(folder):
            os.makedirs(folder)
        for run in runs:
            tmp_dict['Z0-0'][run] = {}
            tmp_ratio_list = [['# ions', '$\\nu_c$-ratio', 'unc']]
            for z in z_classes:
                tmp_ratio_list.append([int(z[-1]), tmp_dict[z][run]['ratio'], tmp_dict[z][run]['ratio_unc']])

# fit
            init_vals = [((max([r for z,r,u in tmp_ratio_list[1:]])
                           -min([r for z,r,u in tmp_ratio_list[1:]]))
                          /(max([z for z,r,u in tmp_ratio_list[1:]])
                            -min([z for z,r,u in tmp_ratio_list[1:]]))), tmp_ratio_list[1][1]]
            fit_1, pcov = scipy.optimize.curve_fit(self.linear,
                                       [z for z,r,u in tmp_ratio_list[1:]],
                                       [r for z,r,u in tmp_ratio_list[1:]],
                                       sigma=[u for z,r,u in tmp_ratio_list[1:]],
                                       absolute_sigma=True,
                                       p0=init_vals)
            fit_2, pcov = scipy.optimize.curve_fit(self.linear,
                                       [z for z,r,u in tmp_ratio_list[1:]],
                                       [r for z,r,u in tmp_ratio_list[1:]],
                                       sigma=[u for z,r,u in tmp_ratio_list[1:]],
                                       absolute_sigma=True,
                                       p0=fit_1)
            fit_3, pcov = scipy.optimize.curve_fit(self.linear,
                                       [z for z,r,u in tmp_ratio_list[1:]],
                                       [r for z,r,u in tmp_ratio_list[1:]],
                                       sigma=[u for z,r,u in tmp_ratio_list[1:]],
                                       absolute_sigma=True,
                                       p0=fit_2)
            fit_4, pcov = scipy.optimize.curve_fit(self.linear,
                                       [z for z,r,u in tmp_ratio_list[1:]],
                                       [r for z,r,u in tmp_ratio_list[1:]],
                                       sigma=[u for z,r,u in tmp_ratio_list[1:]],
                                       absolute_sigma=True,
                                       p0=fit_3)
            fit_res, pcov = scipy.optimize.curve_fit(self.linear,
                                       [z for z,r,u in tmp_ratio_list[1:]],
                                       [r for z,r,u in tmp_ratio_list[1:]],
                                       sigma=[u for z,r,u in tmp_ratio_list[1:]],
                                       absolute_sigma=True,
                                       p0=fit_4)
            perr = np.sqrt(np.diag(pcov))
            #save extrapolation


            tmp_dict['Z0-0'][run]['fit'] = tmp_ratio_list
            tmp_dict['Z0-0'][run]['ratio'] = self.linear(0, fit_res[0], fit_res[1])
            if perr[-1] != np.inf:
                if perr[-1] > 0.000001:
                    tmp_dict['Z0-0'][run]['ratio_unc'] = np.sqrt(np.sum([i[2]**2 for i in tmp_ratio_list[1:] if i[0] != 0]))
                else:
                    tmp_dict['Z0-0'][run]['ratio_unc'] = perr[-1]
            else:
                tmp_dict['Z0-0'][run]['ratio_unc'] = np.sqrt(np.sum([i[2]**2 for i in tmp_ratio_list[1:] if i[2] != np.inf]))
            tmp_dict['Z0-0']['extrapolated_list'].append([int(run[-3:]),
                                                          self.linear(0, fit_res[0], fit_res[1]),
                                                          perr[-1]])



            mpl.rc('font', family='serif', serif='Utopia')      # Utopia LaTeX font!!
            mpl.rc('text', usetex=False)

            f, ax = plt.subplots(figsize=(9,6))
            plt.errorbar([z for z,r,u in tmp_ratio_list[1:]],
                         [r for z,r,u in tmp_ratio_list[1:]],
                         yerr=[u for z,r,u in tmp_ratio_list[1:]], fmt='o', label='$\\nu_c$-ratio')
            plt.errorbar([0], [self.linear(0, fit_res[0], fit_res[1])], yerr=tmp_dict['Z0-0'][run]['ratio_unc'], fmt='g *', label='extrapolation', markersize='10')
            x = np.linspace(-0.2,3.2, 101)
            plt.plot(x, self.linear(x, fit_res[0], fit_res[1]), label='fit', color='r')
            tmp_ratio_list.append([0, self.linear(0, fit_res[0], fit_res[1]), tmp_dict['Z0-0'][run]['ratio_unc']])
            plt.xticks([0,1,2,3], [0,1,2,3], fontsize=16)
            plt.yticks(fontsize=16)
            plt.legend(title=run, fontsize=16)
            plt.xlabel(tmp_ratio_list[0][0], fontsize=22, fontweight='bold')
            plt.ylabel(tmp_ratio_list[0][1], fontsize=22, fontweight='bold')
            plt.tight_layout()
            anchored_text = mpl.offsetbox.AnchoredText('R(Z=0) =  {} $\\pm$ {}'.format(self.linear(0, fit_res[0], fit_res[1]), tmp_dict['Z0-0'][run]['ratio_unc']), loc='lower center')
            ax.add_artist(anchored_text)
            plt.savefig(os.path.dirname(merged_dict_file)+'/Z0-0/'+run+'_'+merged_dict_file[-17:-5]+'.pdf')
            plt.close()
#
        tmp_dict['Z0-0']['weighted_avg'], tmp_dict['Z0-0']['weighted_avg_unc'], tmp_dict['Z0-0']['birge_ratio'] = self.weighted_avg_calc([[r, u] for run,r,u in tmp_dict['Z0-0']['extrapolated_list']])

        file = os.path.basename(merged_dict_file).replace('merged', 'extrapolated')
        filename = folder+file
        with open(filename, 'w') as f:
            json.dump(tmp_dict, f, sort_keys=True, indent=4)

        x = []; y = []; yerr = [];
        for run in [i for i in tmp_dict['Z0-0'].keys() if 'run' in i]:
            x.append(int(run[-3:]))
            y.append(tmp_dict['Z0-0'][run]['ratio'])
            yerr.append(tmp_dict['Z0-0'][run]['ratio_unc'])


        f2, ax2 = plt.subplots(figsize=(9,6))
        plt.errorbar(x, y, yerr=yerr, fmt='o', label='extr. $\\nu_c$-ratio')
        plt.xticks(x, x, fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel('run #', fontsize=22, fontweight='bold')
        plt.ylabel('$\\nu_c$ ratio', fontsize=22, fontweight='bold')
        plt.axhline(y=tmp_dict['Z0-0']['weighted_avg'], linewidth=2, color = 'r', label='weighted mean')
        xlim = ax2.get_xlim()
        plt.fill_between(xlim,
                         y1=tmp_dict['Z0-0']['weighted_avg']+tmp_dict['Z0-0']['weighted_avg_unc'],
                         y2=tmp_dict['Z0-0']['weighted_avg']-tmp_dict['Z0-0']['weighted_avg_unc'],
                         color='r', alpha=0.1)
        plt.fill_between(xlim,
                         y1=tmp_dict['Z0-0']['weighted_avg']+tmp_dict['Z0-0']['weighted_avg_unc']*np.sqrt(tmp_dict['Z0-0']['birge_ratio']),
                         y2=tmp_dict['Z0-0']['weighted_avg']-tmp_dict['Z0-0']['weighted_avg_unc']*np.sqrt(tmp_dict['Z0-0']['birge_ratio']),
                         color='r', alpha=0.1)

        ax2.set_xlim(xlim)
        plt.legend(title=merged_dict_file[-17:-5],fontsize=16)
        plt.tight_layout()
        anchored_text = mpl.offsetbox.AnchoredText('weighted mean =  {} $\\pm$ {} ($\\chi^2_{{red}}$ = {:3.2f})'.format(tmp_dict['Z0-0']['weighted_avg'], tmp_dict['Z0-0']['weighted_avg_unc'], tmp_dict['Z0-0']['birge_ratio']), loc='lower center')
        ax2.add_artist(anchored_text)

        if self.mode == 'single_easy':
            plt.savefig('{}{}-{}-poly-fit-ratio.pdf'.format(self.run_folder, self.isotopes[0], self.isotopes[1]))
        else:
            plt.savefig(os.path.dirname(merged_dict_file)+'/Z0-0/extrapolated_ratio_{}.pdf'.format(merged_dict_file[-17:-5]))
        plt.close()

if __name__ == '__main__':
    freq_ratio = Freq_ratio()
    mode = 'single_easy'
    # z_classes = ['Z1-1', 'Z2-2', 'Z3-3']
    z_classes = 'Z1-1'

    if mode == 'single':
        # run_folder = '/Volumes/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/PI-ICR/129Cd/{}/129Cd_run-013/'.format(z_classes)
        # run_folder = '/Volumes/ISOLTRAP/USER/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/PI-ICR/129Cd/{}/129Cd_run-006/'.format(z_classes)
        # isotopes = ['129Cs', '129gCd']
        run_folder = '/Users/jonaskarthein/cernbox/Analysis/131Cs/analyzed/manual/133Cs-131Cs-run02/'
        # run_folder = '/Volumes/ISOLTRAP/USER/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/PI-ICR/127Cd/{}/127Cd_run-002/'.format(z_classes)
        isotopes = ['133Cs', '131Cs']
        freq_ratio.main(run_folder=run_folder, isotopes=isotopes, degree=2, single_or_batch=mode, z_classes=z_classes)
    elif mode == 'batch':
        upper_run_folder = '/Volumes/ISOLTRAP/USER/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/PI-ICR/127Cd/'
        # isotopes = ['129Cs', '129mCd']
        isotopes = ['133Cs', '127gCd']
        freq_ratio.batch(upper_run_folder=upper_run_folder, isotopes=isotopes, degree=2, z_classes=z_classes)
    elif mode == 'merge_single':
        upper_z_folder = '/Volumes/ISOLTRAP/USER/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/PI-ICR/127Cd/'
        isotopes = ['133Cs', '127gCd']
        freq_ratio.merge_single(upper_z_folder, isotopes)
    elif mode == 'extrapolation':
        json_folder = '/Volumes/ISOLTRAP/USER/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/PI-ICR/127Cd/merged_ratio_dict_133Cs_127mCd.json'
        # json_folder = '/Volumes/ISOLTRAP/USER/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/PI-ICR/129Cd/merged_ratio_dict_129Cs_129gCd.json'
        freq_ratio.extrapolation(json_folder)
    elif mode == 'single_easy':
        data_collection_needed = False
        tof_icr = False

        # run_folder = '/Users/jonaskarthein/cernbox/Analysis/131Cs/analyzed/manual/85Rb-133Cs-run03/'
        # isotopes = ['133Cs', '85Rb']

        # run_folder = '/Users/jonaskarthein/cernbox/Analysis/131Cs/analyzed/manual/133Cs-131Cs-run02/'
        # isotopes = ['133Cs', '131Cs']

        run_folder = '/Users/jonaskarthein/cernbox/Analysis/131Cs/analyzed/manual/133Cs-131Cs-comb/'
        isotopes = ['133Cs', '131Cs']

        # run_folder = '/Users/jonaskarthein/cernbox/Analysis/MirrorNuclei/internal_sandwiches/23-part1/'
        # run_folder = '/Users/jonaskarthein/cernbox/Analysis/MirrorNuclei/poly_fit_partial/21-part2/'
        # run_folder = '/Users/jonaskarthein/cernbox/Analysis/MirrorNuclei/poly_fit_partial/23-part2/'
        # isotopes = ['101In_m', '101SrF']
        # freq_ratio.main(run_folder=run_folder, isotopes=isotopes, degree=3, single_or_batch=mode, z_classes=z_classes)
        # isotopes = ['21Na', '21Na-2']
        # isotopes = ['21Ne', '21Ne-2']
        # if data_collection_needed:
        #     if tof_icr:
        #         freq_ratio.get_data_tof_icr(run_folder, isotopes)
        #     else:
        #         freq_ratio.set_folders(run_folder, isotopes)
        #         freq_ratio.get_data()
        #         sys.exit()

        batch = True
        if batch:
            for iso in [isotopes]:
                for deg in [2,3,4,5,6,7,8,9,10,11,12,13,14,15]:
                    freq_ratio.main(run_folder=run_folder, isotopes=iso, degree=deg, single_or_batch=mode, z_classes=z_classes)
        else:
            # isotopes = ['21Ne', '21Na']
            degree = 8
            freq_ratio.main(run_folder=run_folder, isotopes=isotopes, degree=degree, single_or_batch=mode, z_classes=z_classes)



