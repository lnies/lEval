# -*- coding: utf-8 -*-
"""
Created on Mon 11 January 2022
Modified by Lukas.Nies@cern.ch on 21/10/2021
@author: Lukas Nies
@contact: Lukas.Nies@cern.ch
@license: /
"""

import math
import pandas as pd
import sys
import re

# sys.path.insert(0, '/mnt/c/Users/Lukas/cernbox/Data/Analysis_Code/bin/')

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
        #
        if value == 'mass':
            if not error:
                try:
                    raw = float(self.ame["atomic_mass_raw"][(self.ame["A"]==A) & (self.ame["element"]==X)])
                    comma = float(self.ame["atomic_mass_comma"][(self.ame["A"]==A) & (self.ame["element"]==X)])/1e6
                except(Exception, TypeError) as err:
                    print(f"(TypeError in get_value for A={A}, X={X}): {err}")
                    return
                return (raw+comma)
            else:
                try:
                    data = float(self.ame["atomic_mass_err"][(self.ame["A"]==A) & (self.ame["element"]==X)])
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