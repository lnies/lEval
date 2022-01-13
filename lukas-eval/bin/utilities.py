# -*- coding: utf-8 -*-
"""
Created on Mon 11 January 2022
Modified by Lukas.Nies@cern.ch on 21/10/2021
@author: Lukas Nies
@contact: Lukas.Nies@cern.ch
@license: /
"""

import ROOT
from ROOT import RooFit
from ROOT import RooRealVar, RooArgSet, RooArgList, RooDataHist
from ROOT import RooGenericPdf, RooUniform, RooGaussian, RooGaussModel, RooDecay, RooFormulaVar
from ROOT import RooAddPdf, RooMCStudy
from ROOT import Math as RootMath 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy as sc
from scipy.signal import find_peaks
from scipy.optimize import least_squares
from scipy import stats
import math 
import time
import sys
import re


class AME():
    '''
    Base handling AME related stuff
    Params:
        path_to_ame: path to ame mass file
        ame_version: version of AME. Defaults to 2020
    '''
    def __init__(self, path_to_ame, ame_version = 'ame20'):
        self.path_to_ame = path_to_ame
        self.ame_version = ame_version

        # Init masses dataframe
        if self.ame_version == 'ame20':
            self.ame = pd.read_fwf(self.path_to_ame, #usecols=(2,3,4,6,9,10,11,12,18,21,20),
                                      names=['1', 'N-Z', 'N', 'Z', 'A', 'Unnamed: 5', 'element', 'O', 'Unnamed: 8',
                                               'mass_excess', 'mass_excess_err', 'ebinding', 'nan1', 'ebinding_err', 'nan2', 'ET',
                                               'beta_decay_energy', 'beta_decay_energy_err', 'nan18', 'atomic_mass_raw', 'nan20',
                                               'atomic_mass_comma', 'atomic_mass_err'],
                                      widths=(1,3,5,5,5,1,3,4,1,14,12,13,1,10,1,2,13,11,1,3,1,13,12),
                                      header=28,
                                      index_col=False)
        else:
            print(f"(error in AME__init__): Wrong version parsed. Only 'ame20' available.")
        
    def get_value(self, isotope="1H", value="mass", error=False):
        '''
        Returns value for given isotope
        Params:
            isotope: isotope
            value: value to fetch
            error: returns error on value
        '''
        # Split isotope string into mass number and element string
        A = int(re.split("(\\d+)",isotope)[1])
        X = re.split("(\\d+)",isotope)[2]
        if value == 'mass':
            # Get dataframe index of isotope 
            idx = -1
            idx_list = self.ame.index[(self.ame["A"]==A) & (self.ame["element"]==X)].tolist()
            if len(idx_list) < 1:
                print(f"(Error in get_value for A={A}, X={X}): Can't find isotope with given values.")
                return
            elif len(idx_list) > 1:
                print(f"(Error in get_value for A={A}, X={X}): Found multiple isotopes with given values.")
                return
            else:
                idx = idx_list[0]
            #
            if not error:
                # First convert the feteched dataframe entries to str, remove the hashtag, convert to float 
                try:
                    raw = float(str(self.ame.iloc[idx]['atomic_mass_raw']).strip('#'))
                    comma = float(str(self.ame.iloc[idx]['atomic_mass_comma']).strip('#'))/1e6
                except(Exception, TypeError) as err:
                    print(f"(TypeError in get_value for A={A}, X={X}): {err}")
                    return
                return (raw+comma)
            else:
                try:
                    data = float(str(self.ame.iloc[idx]['atomic_mass_err']).strip('#'))
                except(Exception, TypeError) as err:
                    print(f"(TypeError in get_value for A={A}, X={X}): {err}")
                    return
                return data
        else:
            print(f"(Error in get_value: value={value} unknown.")

class MRToFMS(AME):
    '''
    Base class for performing mass extraction from MR-ToF MS data using the C_{ToF} approach.
    Params:
        path_to_ame: path to ame mass file
        ame_version: version of AME. Defaults to 2020
    Inheritance:
        AME: Inherits functionalities from AME 
    '''
    def __init__(self, path_to_ame, ame_version = 'ame20'):
        # Init base class
        AME.__init__(self, path_to_ame = path_to_ame, ame_version = ame_version)
        #
        #   Store init parameters
        #
        # Store fundamental constants
        self.c = 299792458 # m/s
        self.e = 1.602176634e-19 # C
        self.u = 931494.10242 # keV/c**2
        self.u_err = 000.00000028 # MeV/c**2
        # Store some prominent masses
        self.m_p = self.get_value('1H', 'mass')
        self.m_n = self.get_value('1n', 'mass')
        self.m_39K = self.get_value('39K', 'mass')
        self.m_39K_err = self.get_value('39K', 'mass', error=True)
        self.m_85Rb = self.get_value('85Rb', 'mass')
        self.m_85Rb_err = self.get_value('85Rb', 'mass', error=True)
        self.m_133Cs = self.get_value('133Cs', 'mass')
        self.m_133Cs_err = self.get_value('133Cs', 'mass', error=True)

    def calc_weighted_average(self, x, s):
        ###
        ###     Takes array of values (x) plus array of errors (s) and returns weighted average
        ###
        if len(x) != len(s):
            print("Arrays must have same length")
            return 0
        # 
        sum_avg = 0
        sum_w = 0
        for i in range(len(x)):
            w = s[i]**(-2)
            sum_w += w
            sum_avg += w * x[i]
        #
        return sum_avg/sum_w

    def calc_weighted_averge_err(self, s):
        ###
        ###     Takes array s of individual errors
        ###
        sum_w = 0
        for i in range(len(s)):
            sum_w += s[i]**(-2)
        #
        return math.sqrt(1/sum_w)

    def calc_red_chi_square(self, x, s):
        ###
        ###     Calculates reduced chi square for array of values (x) and array of errors (s)
        ###
        if len(x) != len(s):
            print("Arrays must have same length")
            return 0
        # 
        weighted_average = self.calc_weighted_average(x, s)
        #
        chi_square = 0
        for i in range(len(x)):
            chi_square += (x[i]-weighted_average)**2 / s[i]**2
        #
        return chi_square / (len(x)-1)
        
    def calc_C_ToF(self,t,t1,t2):
        ###
        ### Calculation of C_ToF based on Franks 2013 Nature paper https://www.nature.com/articles/nature12226
        ###
        return (2*t-t1-t2)/(2*(t1-t2))

    def calc_C_ToF_err(self, t, t_err, t1, t1_err, t2, t2_err):
        ###
        ### Calculation of C_ToF error based on Franks 2013 Nature paper https://www.nature.com/articles/nature12226
        ###
        del_t = 1 / (t1-t2)
        del_t1 = (-t+t2)/(t1-t2)**2
        del_t2 = (t-t1)/(t1-t2)**2
        #
        return math.sqrt( (del_t*t_err)**2 + (del_t1*t1_err)**2 + (del_t2*t2_err)**2 )
        
    def calc_sqrt_m(self, C_tof, m1, m2):
        ###
        ###     Calculation of the sqrt of the mass of the species of interest
        ###
        Sigma_ref = math.sqrt(m1) + math.sqrt(m2)
        Delta_ref = math.sqrt(m1) - math.sqrt(m2)
        #
        return C_tof * Delta_ref + Sigma_ref/2 

    def calc_sqrt_m_err(self, C_tof, C_tof_err, m1, m1_err, m2, m2_err):
        ###
        ###      Calculation of the err on the sqrt of the mass
        ###
        del_C_tof = math.sqrt(m1) - math.sqrt(m2)
        del_m1 = C_tof + 1/2
        del_m2 = - C_tof + 1/2
        #
        sqrt_m1_err = 1/2 * m1**(-1/2) * m1_err
        sqrt_m2_err = 1/2 * m2**(-1/2) * m2_err
        #
        return math.sqrt( (del_C_tof * C_tof_err)**2 + ( del_m1 * sqrt_m1_err )**2 + ( del_m2 * sqrt_m2_err )**2 )

    def calc_m_err_alternative(self, sqrt_m, sqrt_m_err):
        ### 
        ###     Calculation of the mass error using the error on the sqrt of the mass
        ###
        return 2 * sqrt_m * sqrt_m_err

    def calc_m_err(self, C_tof, C_tof_err, m1, m1_err, m2, m2_err):
        ###
        ###      Calculation of the mass error using error propagation
        ###
        delta_ref = math.sqrt(m1) - math.sqrt(m2)
        sigma_ref = math.sqrt(m1) + math.sqrt(m2)
        #
        del_C_tof = 2 * C_tof * delta_ref**2 + delta_ref * sigma_ref
        del_m1 = C_tof**2 * (1-m1**(-1/2)) + C_tof + 1/4 * (1+m1**(-1/2))
        del_m2 = C_tof**2 * (1-m2**(-1/2)) + C_tof + 1/4 * (1+m2**(-1/2))
        #
        return math.sqrt( (del_C_tof*C_tof_err)**2 + (del_m1 * m1_err)**2 + (del_m2 * m2_err)**2 )

class Peaks:
    """ 
    Wrapper class for finding peaks in an MR-ToF MS spectrum
        dataframe containing the converted .lst content
    """
    # df_file:  
    def __init__(self, df_file):
        self.file = df_file
        self.n_peaks = 0
        self.bins = 500
        self.peak_threshold = 0.002
        self.peak_min_distance = 5
        self.peak_min_height = 10
        self.peak_width_inbins = (3,100)
        self.peak_prominence = None
        self.peak_wlen = None

    def get_binning(self, bins=10):
        """
        Adapts binning to multiples of 0.8ns, assuming that 0.8ns was used to take data (to common case)
        """
        # Get min and max tof from data frame
        minn = self.file.tof.min()
        maxx = self.file.tof.max()
        return round((maxx-minn)/0.8/bins)
    
    def find_peaks(self, bins=10, peak_threshold = 0.002, peak_min_distance = 5, peak_min_height = 10, peak_width_inbins = (3,100), 
                   peak_prominence = None, peak_wlen = None):
        """  """
        #
        self.pos = []
        self.std = []
        self.left_bases = []
        self.right_bases = []
        self.bins = bins
        # faster binning for projections than histograms -> necessary in order to automatically find peaks
        x_proj_for_pfind = self.file.tof.value_counts(bins=self.get_binning(self.bins)).sort_index()
        y_proj_for_pfind = self.file.sweep.value_counts(bins=self.get_binning(self.bins)).sort_index()
        # Do the peak finding
        self.x_proj_peaks, self.peaks_info = sc.signal.find_peaks(x_proj_for_pfind, threshold=peak_threshold,
                                             distance=peak_min_distance,
                                             height=peak_min_height,
                                             width=peak_width_inbins,
                                             prominence=peak_prominence,
                                             wlen=peak_wlen)
        # Calculate some additional meta data for the found peaks
        self.highest_peak = self.peaks_info['peak_heights'].argmax()
        self.n_peaks = len(self.x_proj_peaks)        
        # variables to store earliest left_base and lastest right_base for constraining plot range
        self.earliest_left_base = 1e15
        self.earliest_peak_idx = 1e15
        self.latest_right_base = 0
        self.latest_peak_idx = 0
        for i in range(self.n_peaks):
            # get the peak bases ranges from the peak finder
            left = x_proj_for_pfind.index.mid[self.peaks_info['left_bases'][i]]
            right = x_proj_for_pfind.index.mid[self.peaks_info['right_bases'][i]]
            self.left_bases.append(left)
            self.right_bases.append(right)
            #calculate the median (more accurate due to asym. tails) of the data in the peak ranges
            #better value for the actual peak center than the simple highest point in peak
            peak_pos = self.file.tof[(self.file.tof < right) &
                                  (self.file.tof > left)].median()
            peak_std = self.file.tof[(self.file.tof < right) &
                                  (self.file.tof > left)].std()
            # estimate the mean and sigma from a Gaussian with exponential tail (accounting for asym. peaks)
            try: 
                peak_fit = stats.exponnorm.fit(file_df.tof[(file_df.tof < peak_pos+peak_std) &
                                                           (file_df.tof > peak_pos-peak_std)],
                                               loc=peak_pos)
                # update peak position if a fit was possible
                peak_pos = peak_fit[1]
            except:
                pass
            self.pos.append(peak_pos)
            self.std.append(peak_std)
            # assign earliest and latest bases
            if (left < self.earliest_left_base and not math.isnan(left) and not math.isnan(peak_std)): 
                self.earliest_left_base = left 
                self.earliest_peak_idx = i 
            if (right > self.latest_right_base and not math.isnan(right) and not math.isnan(peak_std)):
                self.latest_right_base = right
                self.latest_peak_idx = i 
                
    def plot(self, bins = 10, focus=False, log=False, silent = False, save = False, path_to_file = "peaks"):
        #
        plt.rcParams["figure.figsize"] = (10,6)
        #
        if self.n_peaks == 0:
            print("Not peaks, no plots :)")
            return 0 
        #            
        xdata = self.file.tof
        n, xe = np.histogram(xdata, bins=self.get_binning(bins))
        cx = 0.5 * (xe[1:] + xe[:-1])
        dx = np.diff(xe)
        plt.errorbar(cx, n, n ** 0.5, fmt="ok", zorder=1)
        plt.plot(xdata, np.zeros_like(xdata)-5, "|", alpha=0.1, label = "ToF Data", zorder = 3)
        if log:
            plt.yscale('log')
        for i in range(self.n_peaks):
            plt.axvline(self.pos[i], c='r', linewidth=1, zorder=3)
        
        xm = np.linspace(xe[0], xe[-1], num=1000)
        plt.legend();
        # plt.xlim(peaks.pos[0]-300, peaks.pos[0]+300)
        # Zoom in on found peaks
        if focus==False:
            plt.xlim(self.earliest_left_base-200, self.latest_right_base+200)
        else:
            plt.xlim(self.pos[focus]-300, self.pos[focus]+300)
        #
        if not silent: 
            plt.show()
            # plt.clf()

        if save:
            plt.savefig(path_to_file+".pdf", dpi=300)
            plt.clf()
        
    def plot2d(self, bins=500, focus=-1, log=False):
        """
        
        """
        plt.rcParams["figure.figsize"] = (5,15)
        tof = self.file.tof
        sweep = self.file.sweep
        plt.plot(tof, sweep, 'o', alpha=0.05, ms=2, label='unbinned data')
        for i in range(self.n_peaks):
            plt.axvline(self.pos[i], c='r', linewidth=1, zorder=3)
        if focus==False:
            plt.xlim(self.earliest_left_base-200, self.latest_right_base+200)
        else:
            plt.xlim(self.pos[focus]-300, self.pos[focus]+300)
        plt.xlabel("Time-of-Flight [ns]")
        plt.ylabel("Rolling sweep number")
        plt.show()

class softCool(Peaks):
    """
    Class for performing software cooling on 2D MR-ToF MS Data
        df_file: dataframe containing the converted .lst content 
        Inherits functionality from the peak finder 
    """
    def __init__(self, df_file, peaks = False):
        """
        Initialize by finding peaks in the ToF
        """ 
        self.file = df_file
        self.coolfile = df_file.copy(deep=True) # copy for storing the cooled spectrum
        if peaks == False:
            self.peaks=Peaks(self.file)
            self.peaks.find_peaks(bins=200,
                            peak_threshold=2,
                            peak_min_distance=5,
                            peak_min_height=100,
                            peak_width_inbins=(0.5, 100),
                            peak_prominence=1,
                            peak_wlen=50)
        else:
            self.peaks = peaks
        self.selected_peak = -1
        self.corr_factors = []
        self.chunk_size = 10
        self.post_cool = False
        
    def select_peak(self, peak_nb):
        """
        
        """
        self.selected_peak = peak_nb
    
    def calc_corr_factors(self, df, tof_cut, chunk_size=10, method="mean"):
        """
        Function for calculating correction factors
        """
        df_cut = df[(df.tof > tof_cut[0]) & (df.tof < tof_cut[1])]
        self.chunk_size = chunk_size
        #
        if method=="mean":
            self.corr_factors = [
                 np.mean(sublist.tof) 
                 for sublist  
                 in 
                    [
                         df_cut[(df_cut.sweep >= i) & (df_cut.sweep < i+self.chunk_size)]
                         for i 
                         in range(0,df.sweep.iloc[-1], self.chunk_size)
                    ]
            ]
        if method=="median":
            self.corr_factors = [
                 np.median(sublist.tof) 
                 for sublist  
                 in 
                    [
                         df_cut[(df_cut.sweep >= i) & (df_cut.sweep < i+self.chunk_size)]
                         for i 
                         in range(0,df.sweep.iloc[-1], self.chunk_size)
                    ]
            ]
        if method=="average":
            self.corr_factors = [
                 np.average(sublist.tof) 
                 for sublist  
                 in 
                    [
                         df_cut[(df_cut.sweep >= i) & (df_cut.sweep < i+self.chunk_size)]
                         for i 
                         in range(0,df.sweep.iloc[-1], self.chunk_size)
                    ]
            ]

    def cool(self, method="mean", tof_cut_left=300, tof_cut_right=300, chunk_size=10, post_cool = False):
        """
        Routine for performing the cooling
        """
        # df to be cooled
        self.post_cool = post_cool
        if not self.post_cool:
            df_to_cool = self.file
        else:
            df_to_cool = self.coolfile
        #
        self.chunk_size = chunk_size
        tof_cut = [self.peaks.pos[self.selected_peak]-tof_cut_left, self.peaks.pos[self.selected_peak]+tof_cut_right]
        self.calc_corr_factors(df_to_cool, tof_cut, self.chunk_size)
        # print(f"Length correction factors: {len(self.corr_factors)}")
        # print(f"Length chunk sizes: {len(range(0,df_to_cool.sweep.iloc[-1], self.chunk_size))}")
        #
        cooled_tofs = [
            [
                # print(cooled.corr_factors[j[1]], j[1])
                row - self.corr_factors[j[1]] + self.corr_factors[0]
                for row 
                in j[0]
            ]
            for j 
            in
            [
                 [df_to_cool[(df_to_cool.sweep >= i) & (df_to_cool.sweep < i+self.chunk_size)].tof, int(i/self.chunk_size)]
                 for i 
                 in range(0,df_to_cool.sweep.iloc[-1], self.chunk_size)
             ]
        ]
        # print(self.coolfile)
        # print(self.coolfile.sweep.tolist())
        # print(f"Length cooled tofs: {[len(tofs) for tofs in cooled_tofs]}, sum: {sum([len(tofs) for tofs in cooled_tofs])}")
        # Flatten 2d array from previous list-comprehension
        self.coolfile.tof = [
             item 
             for sublist 
             in 
             cooled_tofs for item in sublist
        ]
        # self.coolfile.tof = np.hstack(cooled_tofs)

        # Drop empty sweeps
        for idx in self.coolfile[self.coolfile.tof.isnull()].index:
            self.coolfile = self.coolfile.drop(idx)
        # Export the cooled df
        return self.coolfile
        
    def plot2d(self, bins=500, focus=False, log=False):
        """
        Plot the 2d histograms to compare the corrected and not corrected values
        """
        fig, (ax0, ax1) = plt.subplots(1,2,sharey='row', figsize=(7,7))
        tof = self.file.tof
        sweep = self.file.sweep
        # Plot unbinned and un-corrected data
        ax0.plot(tof, sweep, 'o', alpha=0.05, ms=2, label='unbinned data')
        # Plot correction factors
        y_corr = range(0,self.file.sweep.iloc[-1], self.chunk_size)
        x_corr = self.corr_factors
        ax0.plot(x_corr, y_corr, c='r', linewidth=1, zorder=3)
        # Plot corrected data
        tof = self.coolfile.tof
        sweep = self.coolfile.sweep
        ax1.plot(tof, sweep, 'o', alpha=0.05, ms=2, label='unbinned data')
        if self.post_cool:
            y_corr = range(0,self.coolfile.sweep.iloc[-1], self.chunk_size)
            x_corr = self.corr_factors
            ax1.plot(x_corr, y_corr, c='r', linewidth=1, zorder=3)
        #
        # for i in range(self.n_peaks):
        #     plt.axvline(self.pos[i], c='r', linewidth=1, zorder=3)
        if focus==False:
            ax0.set_xlim(self.peaks.earliest_left_base-300, self.peaks.latest_right_base+300)
            ax1.set_xlim(self.peaks.earliest_left_base-300, self.peaks.latest_right_base+300)
        else:
            ax0.set_xlim(self.peaks.pos[focus]-300, self.peaks.pos[focus]+300)
            ax1.set_xlim(self.peaks.pos[focus]-300, self.peaks.pos[focus]+300)

        plt.xlabel("Time-of-Flight [ns]")
        plt.ylabel("Rolling sweep number")

        plt.tight_layout()
        plt.show()

