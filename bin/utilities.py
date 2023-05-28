# -*- coding: utf-8 -*-
"""
Created on Mon 11 January 2022
Modified by Lukas.Nies@cern.ch on 04/02/2021
@author: Lukas Nies
@contact: Lukas.Nies@cern.ch
@license: MIT
"""

import ROOT
from ROOT import RooFit
from ROOT import RooRealVar, RooArgSet, RooArgList, RooDataHist
from ROOT import RooGenericPdf, RooUniform, RooGaussian, RooGaussModel, RooDecay, RooFormulaVar
from ROOT import RooAddPdf, RooMCStudy
from ROOT import Math as RootMath 
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy as sc
from scipy.signal import find_peaks
from scipy.optimize import least_squares
from scipy import stats
import os
import math 
import time
import sys
import re

from process import ProcessorBase, MPANTMpa, MCS6Lst

FILE_LOCATION = os.path.dirname(__file__)+"/"

def get_time_of_measurement(file, as_datetime = False):
    """
    Fetches time of measurement from data file 
    Parameters:
        - file: path to file (only mpa or lst supported)
        - as_datetime: return time as datetime object, otherwise readable string
    """
    # Determine if file is mpa or lst file
    filetype = file.split(".")[1]
    if filetype == 'mpa':
        # Get time of measurement
        raw_file = MPANTMpa()
        pars, data, df = raw_file.process(files=[file])
    elif filetype == 'lst':
        raw_file = MCS6Lst()
        with open(file,'rb') as listfile:
            raw_file.get_time_patch_and_binary(listfile, verbose=False)
        pars = {}
        pars[re.split("/|\.", file)[-2]] = raw_file.pars
    else:
        print(f"Error: filetype {filetype} not supported.")
        return 0
    #
    file_base = re.split("/|\.", file)[-2]
    for key in pars[file_base].keys():
        # Get start time
        if key == 'cmline0':
            #
            start_date = pars[file_base][key].split(" ")[0].split(".")[0]
            start_time = pars[file_base][key].split(" ")[1].split(".")[0]
        # get end time 
        elif 'report-file' in key:
            #
            end_date = pars[file_base][key].split("written")[1].split(" ")[1]
            end_time = pars[file_base][key].split("written")[1].split(" ")[2]
            #
        #
    if as_datetime:
        return datetime.datetime.strptime(start_date+" "+start_time, '%m/%d/%Y %H:%M:%S'), datetime.datetime.strptime(end_date+" "+end_time, '%m/%d/%Y %H:%M:%S')
    else:
        return(start_date+" "+start_time, end_date+" "+end_time)

def get_ref_values(offline_refs_d, parameter = 'mu0', species = '85Rb', n_revs = '2000', 
                   fetch_file_nb = True, skip_file_nb=[], only_file_nb=[],
                  ):
    """
    Browse through all files in passed folder and get values from fit files and meta-data for 
    passed species and number of revs.
    Parameters:
        - offline_refs_d: Root folder of files to browse through
        - parameter: parameter to fetch
        - species: species filter
        - n_revs: number of revs filter
        - fetch_file_nb : Fetches file number instead of time of file
        - skip_file_nb: array of strings of file numbers to skip
        - only_file_nb: array of strings of file numbers to fetch
    """
    # Get number of fit parameters
    x_val = []
    y_val = []
    y_val_err = []
    for file in os.listdir(offline_refs_d):
        path_to_file = offline_refs_d+file
        try:
            file_number = re.split('In_run_|_|\.|revs',file)[1]
            file_species = re.split('In_run_|_|\.|revs',file)[2]
            file_revs = re.split('In_run_|_|\.|revs',file)[3]
        except:
            continue
        #
        if file_species != species or file_revs != n_revs:
            continue
        if file_number in skip_file_nb:
            # print(f"Skipped file {file_number}")
            continue
        if len(only_file_nb) != 0 and file_number not in only_file_nb:
            # print(f"Skipped file {file_number}")
            continue
        #
        file_base = re.split('\.',file)[0]
        file_extension = re.split('\.',file)[1]
        #
        if file_extension == "mpa":
            raw_file = MPANTMpa()
            pars, data, df = raw_file.process(files=[offline_refs_d+file_base+".mpa"])
            for key in pars[file_base].keys():
                if 'report-file' in key:
                    date_array = pars[file_base][key].split("written")[1].split(" ")
                    time = datetime.datetime.strptime(date_array[1]+" "+date_array[2], '%m/%d/%Y %H:%M:%S')
                    #
        if file_extension != 'txt': 
            continue
        #
        fit_file = FitToDict(path_to_file)
        
        if fetch_file_nb:
            x_val.append(file_number)
        if not fetch_file_nb:
            x_val.append(time)
        y_val.append(float(fit_file.get_val(parameter, 'value')))
        y_val_err.append(float(fit_file.get_val(parameter, 'error')))
        #
    return x_val, y_val, y_val_err

class NUBASE():
    '''
    Base handling NUBASE related stuff
    Params:
        path_to_ame: path to ame mass file
        ame_version: version of AME. Defaults to 2020
        nubase_version: version of NUBASE. Defaults to 2020
    '''
    def __init__(self, path_to_ame=FILE_LOCATION+'mass20.txt', path_to_nubase = FILE_LOCATION+'nubase_3.mas20.txt', ame_version = 'ame20', nubase_version = 'nubase20'):
        self.path_to_ame = path_to_ame
        self.ame_version = ame_version
        self.path_to_nubase = path_to_nubase
        self.nubase_version = nubase_version
        # Check for databases (.txt files)
        if not os.path.exists(path_to_ame):
            print(f"(error in NUBASE.__init__): Could not find AME file under {self.path_to_ame}.")
            self.ame = {}
        if not os.path.exists(path_to_nubase):
            print(f"(error in NUBASE.__init__): Could not find NUBASE file under {self.path_to_nubase}.")
            self.nubase = {}            
            
        # Init AME dataframe
        try:
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
                print(f"(error in NUBASE.__init__): Wrong version parsed. Only 'ame20' available.")
            # First use the "#" that indicate extrapolated values to fetch the information whether the mass is extrapolated or not
            self.ame['extrapolated'] = [
                False 
                if len(value.split("#")) < 2 # if this is the case, it's not extrapolated
                else # Otherwise it is extrapolated
                True
                for value 
                in self.ame['mass_excess']
            ]

            # Calculate binding energies (same as ebinding*A) and format mass excess column
            self.ame['mass_excess'] = [
                float(str(value).strip("#"))
                for value 
                in self.ame['mass_excess']
            ]
            self.ame['mass_excess_err'] = [
                float(str(value).strip("#"))
                for value 
                in self.ame['mass_excess_err']
            ]
            self.ame['binding_energy'] = self.ame['mass_excess'] - self.ame['Z'] * self.get_value('1H', 'mass_excess') - self.ame['N'] * self.get_value('1n', 'mass_excess')
            self.ame['binding_energy_err'] = np.sqrt( (self.ame['mass_excess_err'])**2 + (self.ame['Z']*self.get_value('1H', 'mass_excess', error=True))**2 + (self.ame['N']*self.get_value('1n', 'mass_excess', error=True))**2)

        except(Exception) as err: 
            print(f"(error in NUBASE.__init__): Could not load AME: {err}.")

        # Apply mass filters to AME
        # self.apply_mass_filters()
        
        # Init NUBASE dataframe
        try:
            if self.nubase_version == 'nubase20':
                self.nubase = pd.read_fwf(path_to_nubase, #usecols=(2,3,4,6,9,10,11,12,18,21,20),
                                      names=[
                                          "A", "ZZZ", "element", "s", "mass_excess", "mass_excess_err", "excitation_energy", "excitation_energy_err",
                                          "origin", "isom_unc", "isom_inv", "half_life", "half_life_unit", "half_life_err", "Jpi",
                                          "ensdf_year", "discovery", "branch_ratio"
                                      ],
                                      widths=(4,7,5,1,13,11,12,11,2,1,1,9,3,7,14,12,5,90),
                                      header=25,
                                      index_col=False)          
                # Remove mass number from element column
                self.nubase["element"] = [re.split("(\\d+)",row["element"])[2] for idx, row in self.nubase.iterrows()]
            else:
                print(f"(error in NUBASE.__init__): Wrong version parsed. Only 'ame20' available.")
        except(Exception) as err: 
            print(f"(error in NUBASE.__init__): Could not load NUBASE: {err}.")
   


    def apply_mass_filters(self):
        """
        Function to apply mass filters to dataframe from self.ame
        """

        #s1n
        df = self.ame.groupby('Z').\
            apply(lambda df: df.binding_energy.shift(1)-df.binding_energy).reset_index(level=[0,1]).\
            set_index('level_1').rename(columns={'binding_energy':'s1n'})
        self.ame['s1n'] = df['s1n']

        df = self.ame.groupby('Z').\
            apply(lambda df: np.sqrt(df.binding_energy_err.shift(1)**2+df.binding_energy_err**2)).reset_index(level=[0,1]).\
            set_index('level_1').rename(columns={'binding_energy_err':'s1n_err'})
        self.ame['s1n_err'] = df['s1n_err']

        #s2n - This is effectively B(Z,N-2) - B(Z, N)
        df = self.ame.groupby('Z').\
            apply(lambda df: df.binding_energy.shift(2)-df.binding_energy).reset_index(level=[0,1]).\
            set_index('level_1').rename(columns={'binding_energy':'s2n'})
        self.ame['s2n'] = df['s2n']

        df = self.ame.groupby('Z').\
            apply(lambda df: np.sqrt(df.binding_energy_err.shift(2)**2+df.binding_energy_err**2)).reset_index(level=[0,1]).\
            set_index('level_1').rename(columns={'binding_energy_err':'s2n_err'})
        self.ame['s2n_err'] = df['s2n_err']

        #d3n
        df = self.ame.groupby('Z').\
            apply(lambda df: 0.5*((-1)**(df.N))*(-2*df.binding_energy+df.binding_energy.shift(1)+df.binding_energy.shift(-1))).reset_index(level=[0,1]).\
            set_index('level_1')
        df.columns = ['Z','d3n']
        self.ame['d3n'] = df['d3n']

        df = self.ame.groupby('Z').\
            apply(lambda df: np.sqrt((4*df.binding_energy_err**2+df.binding_energy_err.shift(1)**2+df.binding_energy_err.shift(-1)**2)/4)).\
            reset_index(level=[0,1]).\
            set_index('level_1')
        df.columns = ['Z','d3n_err']
        self.ame['d3n_err'] = df['d3n_err']

        #d2n
        df = self.ame.groupby('Z').\
            apply(lambda df: df.s2n-df.s2n.shift(-2)).reset_index(level=[0,1]).\
            set_index('level_1')
        df.columns = ['Z','d2n']
        self.ame['d2n'] = df['d2n']

        df = self.ame.groupby('Z').\
            apply(lambda df: np.sqrt(df.s2n_err**2+df.s2n_err.shift(-2)**2)).\
            reset_index(level=[0,1]).\
            set_index('level_1')
        df.columns = ['Z','d2n_err']
        self.ame['d2n_err'] = df['d2n_err']

        #d2n*
        df = self.ame.groupby('Z').\
            apply(lambda df: df.s2n.shift(2)-df.s2n).reset_index(level=[0,1]).\
            set_index('level_1')
        df.columns = ['Z','d2n_s']
        self.ame['d2n_s'] = df['d2n_s']

        df = self.ame.groupby('Z').\
            apply(lambda df: np.sqrt(df.s2n_err**2+df.s2n_err.shift(2)**2)).\
            reset_index(level=[0,1]).\
            set_index('level_1')
        df.columns = ['Z','d2n_s_err']
        self.ame['d2n_s_err'] = df['d2n_s_err']
        
        #d2n#
        df = self.ame.groupby('Z').\
            apply(lambda df: df.s2n.shift(-2)-df.s2n.shift(-4)).reset_index(level=[0,1]).\
            set_index('level_1')
        df.columns = ['Z','d2n_d']
        self.ame['d2n_d'] = df['d2n_d']

        df = self.ame.groupby('Z').\
            apply(lambda df: np.sqrt(df.s2n_err.shift(-2)**2+df.s2n_err.shift(-4)**2)).\
            reset_index(level=[0,1]).\
            set_index('level_1')
        df.columns = ['Z','d2n_d_err']
        self.ame['d2n_d_err'] = df['d2n_d_err']

        #s1p
        df = self.ame.groupby('N').\
            apply(lambda df: df.binding_energy.shift(1)-df.binding_energy).reset_index(level=[0,1]).\
            set_index('level_1').rename(columns={'binding_energy':'s1p'})
        self.ame['s1p'] = df['s1p']

        df = self.ame.groupby('N').\
            apply(lambda df: np.sqrt(df.binding_energy_err.shift(1)**2+df.binding_energy_err**2)).reset_index(level=[0,1]).\
            set_index('level_1').rename(columns={'binding_energy_err':'s1p_err'})
        self.ame['s1p_err'] = df['s1p_err']

        #s2p
        df = self.ame.groupby('N').\
            apply(lambda df: df.binding_energy.shift(2)-df.binding_energy).reset_index(level=[0,1]).\
            set_index('level_1').rename(columns={'binding_energy':'s2p'})
        self.ame['s2p'] = df['s2p']

        df = self.ame.groupby('N').\
            apply(lambda df: np.sqrt(df.binding_energy_err.shift(2)**2+df.binding_energy_err**2)).reset_index(level=[0,1]).\
            set_index('level_1').rename(columns={'binding_energy_err':'s2p_err'})
        self.ame['s2p_err'] = df['s2p_err']

        #d2p
        df = self.ame.groupby('N').\
            apply(lambda df: df.s2p-df.s2p.shift(-2)).reset_index(level=[0,1]).\
            set_index('level_1')
        df.columns = ['Z','d2p']
        self.ame['d2p'] = df['d2p']

        df = self.ame.groupby('N').\
            apply(lambda df: np.sqrt(df.s2p_err**2+df.s2p_err.shift(-2)**2)).\
            reset_index(level=[0,1]).\
            set_index('level_1')
        df.columns = ['Z','d2p_err']
        self.ame['d2p_err'] = df['d2p_err']

        #d2p*
        df = self.ame.groupby('N').\
            apply(lambda df: df.s2p.shift(2)-df.s2p).reset_index(level=[0,1]).\
            set_index('level_1')
        df.columns = ['Z','d2p_s']
        self.ame['d2p_s'] = df['d2p_s']

        df = self.ame.groupby('N').\
            apply(lambda df: np.sqrt(df.s2p_err**2+df.s2p_err.shift(2)**2)).\
            reset_index(level=[0,1]).\
            set_index('level_1')
        df.columns = ['Z','d2p_s_err']
        self.ame['d2p_s_err'] = df['d2p_s_err']
        
                #d3p
        df = self.ame.groupby('N').\
            apply(lambda df: 0.5*((-1)**(df.Z))*(-2*df.binding_energy+df.binding_energy.shift(1)+df.binding_energy.shift(-1))).reset_index(level=[0,1]).\
            set_index('level_1')
        df.columns = ['Z','d3p']
        self.ame['d3p'] = df['d3p']
        
        df = self.ame.groupby('N').\
            apply(lambda df: np.sqrt((4*df.binding_energy_err**2+df.binding_energy_err.shift(1)**2+df.binding_energy_err.shift(-1)**2)/4)).\
            reset_index(level=[0,1]).\
            set_index('level_1')
        df.columns = ['Z','d3p_err']
        self.ame['d3p_err'] = df['d3p_err']

    def get_value(self, isotope="1H", value="mass", error=False, state=''):
        '''
        Returns value for given isotope
        Params:
            isotope: isotope
            value: value to fetch ("mass", "mass_excess", "excitation_energy")
            error: returns error on value
        '''
        # Split isotope string into mass number and element string
        split_string = re.split("(\\d+)",isotope)
        fetched_value = 0 
        # Iterate through every isotope in the chain of isotopes passed (to handle molecules)
        for i in range(1, len(split_string)-1, 2):
            A = int(split_string[i])
            X = split_string[i+1]
            if value in ['mass', 'mass_excess', 'binding_energy', 'half_life']:
                # Get dataframe index of isotope 
                idx = -1
                idx_list = self.ame.index[(self.ame["A"]==A) & (self.ame["element"]==X)].tolist()
                if len(idx_list) < 1:
                    print(f"(Error in get_value for A={A}, X={X}): Can't find isotope with given values.")
                    return -1
                elif len(idx_list) > 1:
                    print(f"(Error in get_value for A={A}, X={X}): Found multiple isotopes with given values.")
                    return -1
                else:
                    idx = idx_list[0]
                #
                if not error:
                    # First convert the feteched dataframe entries to str, remove the hashtag, convert to float 
                    if value == 'mass':
                        try:
                            raw = float(str(self.ame.iloc[idx]['atomic_mass_raw']).strip('#'))
                            comma = float(str(self.ame.iloc[idx]['atomic_mass_comma']).strip('#'))/1e6
                        except(Exception, TypeError) as err:
                            print(f"(TypeError in get_value for A={A}, X={X}): {err}")
                            return -1
                        fetched_value += raw+comma
                    else:
                        try:
                            fetched_value = float(str(self.ame.iloc[idx][value]).strip('#'))
                        except(Exception, TypeError) as err:
                            print(f"(TypeError in get_value for A={A}, X={X}): {err}")
                            return -1
                else:
                    data = 0
                    try:
                        if value == 'mass': 
                            data = float(str(self.ame.iloc[idx]['atomic_mass_err']).strip('#'))
                        else:
                            data = float(str(self.ame.iloc[idx][value+"_err"]).strip('#'))
                    except(Exception, TypeError) as err:
                        print(f"(TypeError in get_value for A={A}, X={X}): {err}")
                        return -1
                    fetched_value += data**2
                        #
            elif value == 'excitation_energy' and state in ['i', 'j', 'k', 'l', 'm', 'n']:
                # Get dataframe index of isotope 
                idx = -1
                idx_list = self.nubase.index[(self.nubase["A"]==A) & (self.nubase["element"]==X) & (self.nubase["s"]==state)].tolist()
                if len(idx_list) < 1:
                    print(f"(Error in get_value for A={A}, X={X}): Can't find isotope in NUBASE with given values.")
                    return -1
                elif len(idx_list) > 1:
                    print(f"(Error in get_value for A={A}, X={X}): Found multiple isotopes in NUBASE with given values.")
                    return -1
                else:
                    idx = idx_list[0]
                #
                if not error:
                    # First convert the feteched dataframe entries to str, remove the hashtag, convert to float 
                    if value == 'excitation_energy':
                        try:
                            fetched_value = float(str(self.nubase.iloc[idx]['excitation_energy']).strip('#'))
                        except(Exception, TypeError) as err:
                            print(f"(TypeError in get_value for A={A}, X={X}): {err}")
                            return -1
                else:
                    data = 0
                    try:
                        if value == 'excitation_energy': 
                            data = float(str(self.nubase.iloc[idx]['excitation_energy_err']).strip('#'))
                    except(Exception, TypeError) as err:
                        print(f"(TypeError in get_value for A={A}, X={X}): {err}")
                        return -1
                    fetched_value += data**2
            #
            elif value == 'half_life' and state in ['i', 'j', 'k', 'l', 'm', 'n']:
                # Get dataframe index of isotope 
                idx = -1
                idx_list = self.nubase.index[(self.nubase["A"]==A) & (self.nubase["element"]==X) & (self.nubase["s"]==state)].tolist()
                if len(idx_list) < 1:
                    print(f"(Error in get_value for A={A}, X={X}): Can't find isotope in NUBASE with given values.")
                    return -1
                elif len(idx_list) > 1:
                    print(f"(Error in get_value for A={A}, X={X}): Found multiple isotopes in NUBASE with given values.")
                    return -1
                else:
                    idx = idx_list[0]
                #
                if not error:
                    # First convert the feteched dataframe entries to str, remove the hashtag, convert to float 
                    if value == 'half_life':
                        try:
                            fetched_value = str(self.nubase.iloc[idx]['half_life'])
                            # fetched_value = float(str(self.nubase.iloc[idx]['half_life']).strip('#'))
                        except(Exception, TypeError) as err:
                            print(f"(TypeError in get_value for A={A}, X={X}): {err}")
                            return -1
                else:
                    data = 0
                    try:
                        if value == 'half_life': 
                            data = str(self.nubase.iloc[idx]['half_life_err']).strip('#')
                            # data = float(str(self.nubase.iloc[idx]['half_life_err']).strip('#'))
                    except(Exception, TypeError) as err:
                        print(f"(TypeError in get_value for A={A}, X={X}): {err}")
                        return -1
                    fetched_value += data**2
            #
            else:
                print(f"(Error in get_value: value='{value}' and/or state='{state}' unknown.")
                return -1
        #
        if not error:
            return fetched_value
        else: 
            return math.sqrt(fetched_value)

class FitToDict:
    '''
    Class for reading fit files created by the fit functions. Stores fit in .fit dictionary
    Initialization:
    Parameters
        - file_path: path to fit file
    '''
    def __init__(self, file_path, verbose = 0):
        self.fit = {}
        self.line_nb = 0
        self.res_table_line = 0
        self.fit_val_line = 0
        self.file_path = file_path
        # Read file
        self.fit = self.__read(verbose)
        if self.fit == 0:
            print(f"(FitToDict.__init__): File at '{file_path}' does not exist.")
            return 0

    def __initialize(self, verbose = 0):
        '''
        PRIVATE: Init file, read the meta data into dict and save where the results table and fit values table start
        '''
        # open the file
        if os.path.isfile(self.file_path):
            file = open(self.file_path)
            for line in file.readlines():
                # Increment line counter
                self.line_nb += 1
                # get rid of the newline
                line = line[:-1]
                try:
                    # this will break if you have whitespace on the "blank" lines
                    if line:
                        # skip comment lines
                        if line[0] == '#': next
                        # this assumes everything starts on the first column
                        if line[0] == '[':
                            # strip the brackets
                            section = line[1:-1]
                            # create a new section if it doesn't already exist
                            if not section in self.fit:
                                self.fit[section] = {}
                                # Save where which section is
                                if section == 'RESULTS-TABLE':
                                    self.res_table_line = self.line_nb
                                if section == 'FIT-VALUES':
                                    self.fit_val_line = self.line_nb
                        else:
                            # split on first the equal sign
                            (key, val) = line.split('=', 1)
                            # create the attribute as a list if it doesn't
                            # exist under the current section, this will
                            # break if there's no section set yet
                            if not key in self.fit[section]:
                                self.fit[section][key] = val
                            # append the new value to the list
                            #sections[section][key].append(val)
                    #
                except Exception as e:
                    if verbose > 0: 
                        print(str(e) + "line:" +line)
            #
            return 1
            #
        else:
            return 0 

    def __read_tables(self, verbose = 0):
        '''
        Use pandas to read the tables 
        '''
        #
        if 'RESULTS-TABLE' in self.fit:
            if 'FIT-VALUES' in self.fit:
                n_footer = self.line_nb - self.fit_val_line + 1    
            else:
                n_footer = 0 
            if verbose > 0: print(f"res_table_line: {self.res_table_line}\nn_footer: {n_footer}")  
            self.fit['RESULTS-TABLE'] = pd.read_csv(self.file_path, header=self.res_table_line, delimiter=' ', 
                                                    skipfooter=n_footer, engine='python')
        if 'FIT-VALUES' in self.fit: 
            if verbose > 0: print(f"fit_val_line: {self.fit_val_line}")  
            self.fit['FIT-VALUES'] = pd.read_csv(self.file_path, header=self.fit_val_line, delimiter=' ')
    
    def __read(self, verbose = 0):
        '''
        Function for reading a fit into a dict
        Parameters:
            - file_path: Path to file
        Return:
            - Dictionary with meta data, fit results, and fit values for plotting
        '''
        # 
        init = self.__initialize()
        if init != 1:
            return 0
        # read results table
        self.__read_tables(verbose)
        #
        return self.fit

    def get_val(self, key, value = None):
        '''
        Returns value either from meta-data dictionary or from the data frames
        Parameters:
            - key: name of value to be fetched
        '''
        #
        if key in self.fit['META-DATA']:
            return self.fit['META-DATA'][key]
        #
        elif key in self.fit['RESULTS-TABLE']['var'].to_numpy():
            if not value:
                print("Key {key} specified, but no value (e.g. in df for key 'mu0', value could be 'value'")
                return 
            else:
                if value not in self.fit['RESULTS-TABLE'].columns:
                    print("Value {value} not in dataframe")
                    return 
                return(float( self.fit['RESULTS-TABLE'][value][self.fit['RESULTS-TABLE']['var']==key] ))
        #
        else:
            print(f"Key {key} does not exist in the fit dictionary.")
            return

class MRToFUtils(NUBASE):
    '''
    Utility class for performing mass extraction from MR-ToF MS data using the C_{ToF} approach.
    Params:
        path_to_ame: path to ame mass file
        ame_version: version of AME. Defaults to 2020
    Inheritance:
        AME: Inherits functionalities from AME 
    '''
    def __init__(self, path_to_ame=FILE_LOCATION+'mass20.txt', path_to_nubase = FILE_LOCATION+'nubase_3.mas20.txt', ame_version = 'ame20', nubase_version = 'nubase20'):
        # Init base class
        NUBASE.__init__(self, path_to_ame = path_to_ame, ame_version = ame_version, path_to_nubase=path_to_nubase, nubase_version=nubase_version)
        #
        #   Store init parameters
        #
        # Store fundamental constants
        self.c = 299792458 # m/s
        self.e = 1.602176634e-19 # C
        self.e_in_u = 0.00054857990943 # electron mass in atomic mass units
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
        #
        self.tof_calib_loaded = False

    # Functionality

    def load_isoltrap(self, data):
        """
        Load separately ISOLTRAP data, fmt Z, N, mass excess (keV), unc. (keV)
        @author: Vladimir Manea
        """
        isoltrap_import = np.genfromtxt(data,
                                   skip_header=1,
                                   dtype='int,int,float,float')
        isoltrap_table = [list(elem) for elem in isoltrap_import.tolist()]
        df = pd.DataFrame(isoltrap_table, columns=['Z', 'N', 'mass_excess',
                                                        'mass_excess_err'])
        # Add mass number
        df['A'] = df['N'] + df['Z']
        # Calculate binding energy
        df['binding_energy'] = df['mass_excess']-df['Z'] * self.get_value('1H', value="mass_excess") - df['N'] * self.get_value('1n', value="mass_excess")
        df['binding_energy_err'] = (df['mass_excess_err']**2+df['Z']**2*self.get_value('1H', value="mass_excess", error=True)**2+df['N']**2*self.get_value('1n', value="mass_excess", error=True)**2)**0.5
        df = df.sort_values(by=['Z'])
        return df

    def replace_2_in_1(self, df_1, df_2, value = 'mass_excess', self_apply = False, self_apply_df = 'AME'):
        """
        Overwrites in df_1 the data from df_2 (useful, for example, for combining AME and ISOLTRAP data).
        Parameters:
            - value: key of value to replace in df_1
            - self_apply is True, applies to internal dataframe specified by self_apply_df
            - self_apply_df: either 'AME' to apply to self.ame, or 'NUBASE' to apply to self.nubase 
        @author: Vladimir Manea, Maxime Mougeot
        @amended: Lukas Nies
        """
        #
        key = value
        key_err = value+"_err"
        #
        df = df_1.copy(deep=True)
        imxs = df.columns.get_loc(key)
        imxs_err = df.columns.get_loc(key_err)
        for index, row in df_2.iterrows():
            Z = df_2.iat[index,df_2.columns.get_loc('Z')]
            N = df_2.iat[index,df_2.columns.get_loc('N')]
            mxs_new=df_2.iat[index,df_2.columns.get_loc(key)]
            mxs_new_err=df_2.iat[index,df_2.columns.get_loc(key_err)]
            idx = df.index[(df['Z']==Z) & (df['N']==N)].tolist()[0]
            df.iat[idx,imxs] = mxs_new
            df.iat[idx,imxs_err] = mxs_new_err
        # Recalculate binding energies
        if key == 'mass_excess':
            df['binding_energy'] = df[key]-df['Z'] * self.get_value('1H', value=key) - df['N'] * self.get_value('1n', value=key)
            df['binding_energy_err'] = (df[key_err]**2+df['Z']**2*self.get_value('1H', value=key, error=True)**2+df['N']**2*self.get_value('1n', value=key, error=True)**2)**0.5
            #MM added line to avoid extrapolation from being taken into account (see load AME function)
            df['binding_energy'] = df.apply(lambda row: np.nan if np.isnan(row['binding_energy_err']) else row['binding_energy'],axis=1)
        if self_apply and self_apply_df == 'AME':
            self.ame = df
        if self_apply and self_apply_df == 'NUBASE':
            self.nubase = df
        else:
            return df

    # Calculation functions

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

    # Peak identification

    def __calc_tof_params(self, m0, m1, tof0, tof1):
        """
        calculates the ToF parameters based on tofs and masses of two calibrants
        """
        a = (tof0-tof1)/(math.sqrt(m0)-math.sqrt(m1))
        b = tof0 - a * math.sqrt(m0)
        print((tof0-tof1))
        print(math.sqrt(m0))
        return a, b

    def load_tof_calib(self, file):
        """
        Loads ToF calibration parameters from an ISOLTRAP ToF calibration excel sheet. Needs always
        to be in the same format, otherwise the absolute references in this function won't work.
        Tested with a file from 2015 and 2022, both worked.
        """
        data = pd.read_excel(file, 'new tof calibration', header=0)#, index_col=None, usecols = "Q", header = 0, nrows=0)
        self.m0 = float(data.columns[16])
        self.m1 = float(data[data.columns[16]][0])
        self.tofm0_0 = float(data["Unnamed: 28"][0])
        self.tofm1_0 = float(data["Unnamed: 28"][1])
        self.tofm0_1 = float(data.columns[18])
        self.tofm1_1 = float(data[data.columns[18]][0])
        # print(self.m0, self.m1, self.tofm0_0, self.tofm1_0, self.tofm0_1, self.tofm1_1)
        self.a0, self.b0 = self.__calc_tof_params(self.m0, self.m1, self.tofm0_0, self.tofm1_0)
        self.a1, self.b1 = self.__calc_tof_params(self.m0, self.m1, self.tofm0_1, self.tofm1_1)
        self.revN2 = int(data["revolutions at calibration"][0]) # number of revs for calibration
        #
        self.tof_calib_loaded = True

    def __calc_tof_outside_device(self, m, det_loc='EMP2h'):
        """
        Calculates the relevant ToFs outside the MR-ToF MS
        """
        F1 = (self.a0 * np.sqrt(m - self.e_in_u) + self.b0) * 3/4  # Flight time into center of isep cavity, aka pulse down delay. Scaling factor 3/4 depends on location of detector from which total flight time outside of device is determined.
        MCP = 1/3*F1 # Flight time from center of isep cavity to detector, changes between EMP2h and EMP3h
        return F1, MCP

    def calc_ToF(self, m, nrevs = 1000, a0=None, b0=None, a1=None, b1=None, revsN2=None):
        """
        Calculates the ToF for a mass m at nrevs for given calibration parameters. Those can either 
        be passed directly (to be implemented) or loaded through load_tof_calib(self)
        inputs:
            m: scalar or array like, mass in atomic mass units
            nrevs: int number of revs
        """
        #
        if not self.tof_calib_loaded:
            print("(MRToFUtils:calc_ToF)ToF Calibration file not loaded!")
            return False
        #
        F1, MCP = self.__calc_tof_outside_device(m)
        TG1 = ((self.a1 * np.sqrt(m) + self.b1) - F1 - MCP ) / self.revN2 * int(nrevs) 
        tof = TG1 + F1 + MCP 
        return tof

    def recalibrate_tof(self, m1, tof1, nrevs):
        """
        Receives ToF of a known mass at nrevs and recalculates the a1 parameter
        """
        # F1, MCP = self.__calc_tof_outside_device(m)
        # TG1 = tof - F1 - MCP
        # self.a1 = ((TG1 * int(self.revN2) / int(nrevs)) + F1 + MCP - self.b1) / np.sqrt(m) / 1e3 # calculate a1 from absolute ToF and convert to micro-seconds 
        if nrevs != 0:
            m0 = self.m0
            tof0 = self.calc_ToF(m0, nrevs)
            self.revN2 = nrevs
            self.a1, self.b1 = self.__calc_tof_params(m0, m1, tof0, tof1)
        else:
            m0 = self.m0
            tof0 = self.calc_ToF(m0, nrevs)
            self.revN2 = nrevs
            self.a0, self.b0 = self.__calc_tof_params(m0, m1, tof0, tof1)

class MRToFIsotope(MRToFUtils):
    '''
    Class for handling and storing mass data 
    Params:
        - isotope, ref1, ref2: strings of isotopes to be used
        - n_revs: number of revs
        - path_to_ame: path to ame mass file
        - ame_version: version of AME. Defaults to 2020
    Inheritance:
        AME: Inherits functionalities from AME 
        MRToFUtils: Inherits functionalities for calculating masses 
    '''
    def __init__(self, isotope, ref1, ref2, n_revs, state = 'gs', path_to_ame=FILE_LOCATION+'mass20.txt', 
                 path_to_nubase = FILE_LOCATION+'nubase_3.mas20.txt', ame_version = 'ame20', nubase_version = 'nubase20'):
        # Init base class
        MRToFUtils.__init__(self, path_to_ame = path_to_ame, ame_version = ame_version, path_to_nubase=path_to_nubase,
                                  nubase_version=nubase_version)
        # 
        self.isotope = isotope
        self.A = int(re.split("(\\d+)",isotope)[1])
        self.ref1 = ref1
        self.ref2 = ref2 
        self.n_revs = n_revs
        self.state = state
        self.states = []
        self.Var_dict = {}
        #
        self.m_ref1 = self.get_value(self.ref1)
        self.m_ref1_err = self.get_value(self.ref1,error=True)
        self.m_ref2 = self.get_value(self.ref2)
        self.m_ref2_err = self.get_value(self.ref2,error=True)
        #
        self.m_isotope_AME = self.get_value(isotope)
        self.m_isotope_AME_err = self.get_value(isotope,error=True)
        # if isomere is to be calculated
        if self.state != 'gs':
            self.exc_energy_NUBASE = self.get_value(isotope, "excitation_energy", state=self.state) 
            if self.exc_energy_NUBASE == -1:
                print(f"(MRToFIsotope.__init__): Can't find isotope '{self.isotope}-{self.state}' in NUBASE.")
                # return -1   
            self.exc_energy_NUBASE_err = self.get_value(isotope, "excitation_energy", state=state, error=True)
            #
            if self.state == 'm':
                self.isomere_key = 'mu1'
                self.doublet_key = 'E0'
            elif self.state == 'n':
                self.isomere_key = 'mu2'
                self.doublet_key = 'E1'
            else:
                print(f"(MRToFIsotope.__init__): Unknown isomeric state '{self.state}' (only 'm', 'n').")
                # return -1
        #
        self.custom_gs = ''
        self.custom_gs_err = ''
        self.C_tof = ''
        self.C_tof_err = ''
        self.m_isomere = ''
        self.m_isomere_err = ''

    def __infer_states(self, fit_df):
        """
        Infer states from number of peaks passed
        """
        for idx,row in fit_df.iterrows():
            self.Var_dict[row['var']] = row['value']
        #
        dic = self.Var_dict
        for key in dic:
            d = {}
            if(re.match("^mu", key)):
                if len(key.split("-")) == 1:
                    d = {
                        'peak': key.split("-")[0],
                        'state': 'gs'
                    }
                if len(key.split("-")) == 2:
                    d = {
                        'peak': key.split("-")[0],
                        'state': key.split("-")[1]
                    }
                #
                self.states.append(d)

    def __store_tofs(self, file_isotope='', file_ref1='', file_ref2='',
                        t_isotope='', t_ref1='', t_ref2='',
                        t_isotope_err='', t_ref1_err='', t_ref2_err='',
                        centroid = 'mu0', online_ref = '', online_ref2 = '', tweak_tofs = [0,0,0],
                        is_doublet = False, dt = '', dt_err = '',
                        ):
        """
        Stores information passed to calc_mass or calc_exc_energy and loads fit files if passed
        """    
        #
        self.centroid = centroid
        #
        self.file_isotope = file_isotope
        # Store isotope of interest fit file
        if file_isotope != '' and t_isotope == '':
            self.isotope_fit = FitToDict(file_isotope)
            self.__infer_states(self.isotope_fit.fit['RESULTS-TABLE'])
            self.isotope_gs_t = float(self.isotope_fit.get_val(centroid, 'value'))# + tweak_tofs[0]
            self.isotope_gs_t_err = float(self.isotope_fit.get_val(centroid, 'error'))
            # If isomer is passed, store isomere ToF or ToF difference as well
            if self.state != 'gs':
                # If doublet is passed as fit file
                if is_doublet and dt == '':
                    self.isotope_m_dt = float(self.isotope_fit.get_val(f"{centroid}-{self.doublet_key}", 'value')) + tweak_tofs[0]
                    self.isotope_m_dt_err = float(self.isotope_fit.get_val(f"{centroid}-{self.doublet_key}", 'error'))
                else:
                    self.isotope_m_t = float(self.isotope_fit.get_val(self.isomere_key, 'value')) + tweak_tofs[0]
                    self.isotope_m_t_err = float(self.isotope_fit.get_val(self.isomere_key, 'error'))
        # If raw ToFs and not fit files are passed, store ToFs directly
        elif file_isotope == '' and t_isotope != '' and t_isotope_err != '':
            self.isotope_gs_t = t_isotope + tweak_tofs[0]
            self.isotope_gs_t_err = t_isotope_err
            # If isomere ToF is passed, set gs to -1 and store isomere ToF instead 
            if self.state != 'gs':
                # self.isotope_gs_t = -1 
                # self.isotope_gs_t_err = -1 
                # If passed as doublet
                if is_doublet:
                    self.isotope_m_dt =  t_ref1 - self.isotope_gs_t
                    self.isotope_m_dt_err = math.sqrt(t_ref1_err**2 + self.isotope_gs_t_err**2)
                else:
                    self.isotope_m_t = t_isotope + tweak_tofs[0]
                    self.isotope_m_t_err = t_isotope_err


        else:
            print(f"Error input isotope")
            return -1
        #
        self.file_ref1 = file_ref1
        if file_ref1 != '' and t_ref1 == '' and online_ref == '':
            self.ref1_fit = FitToDict(file_ref1)
            self.ref1_t = float(self.ref1_fit.get_val(centroid, 'value')) + tweak_tofs[1]
            self.ref1_t_err = float(self.ref1_fit.get_val(centroid, 'error'))
        elif file_isotope != '' and file_ref1 != '' and t_ref1 == '' and online_ref != '':
            self.ref1_t = float(self.isotope_fit.get_val(online_ref, 'value')) + tweak_tofs[1]
            self.ref1_t_err = float(self.isotope_fit.get_val(online_ref, 'error'))
        elif file_ref1 == '' and t_ref1 != '' and t_ref1_err != '':
            self.ref1_t = t_ref1 + tweak_tofs[1]
            self.ref1_t_err = t_ref1_err
        else:
            if not is_doublet:
                print(f"Error input ref1")
                return
        #
        self.file_ref2 = file_ref2
        if file_ref2 != '' and t_ref2 == '':
            self.ref2_fit = FitToDict(file_ref2)
            self.ref2_t = float(self.ref2_fit.get_val('mu0', 'value')) + tweak_tofs[2]
            self.ref2_t_err = float(self.ref2_fit.get_val('mu0', 'error'))
        elif file_isotope != '' and file_ref2 != '' and t_ref2 == '' and online_ref2 != '':
            self.ref2_t = float(self.isotope_fit.get_val(online_ref2, 'value')) + tweak_tofs[2]
            self.ref2_t_err = float(self.isotope_fit.get_val(online_ref2, 'error'))
        elif file_ref2 == '' and t_ref2 != '' and t_ref2_err != '':
            self.ref2_t = t_ref2 + tweak_tofs[2]
            self.ref2_t_err = t_ref2_err
        else:
            if not is_doublet:
                print(f"Error input ref2")
                return

    def __print_results(self):
        '''
        
        '''

    def calc_exc_energy_from_doublet(self, t0, dt, m0):
        '''
        Calculation of excitation energy based on measured isomer/ground state doublet
        Parameters:
            - t0: ToF of ground state
            - dt: ToF difference between ground state and isomer
            - m0: Ground state mass in keV
        Returns excitation energy in keV
        ''' 
        # Simplified form assuming dt<<t_0
        return 2 * dt/t0 * m0

    def calc_exc_energy_err_from_doublet(self, t0, t0_err, dt, dt_err, m0, m0_err):
        '''
        Calculation of excitation energy error based on measured isomer/ground state doublet
        Parameters:
            - t0, t0_err: ToF of ground state
            - dt, dt_err: ToF difference between ground state and isomer
            - m0, m0_err: Ground state mass in keV
        Returns error on excitation energy in keV
        ''' 
        # Full form
        part1 = ( ( 2*dt/t0**2 + 2/t0 ) * m0 * dt_err ) ** 2
        part2 = ( ( -2*dt**2/t0**3 - 2*dt/t0**2 ) * m0 * t0_err ) ** 2
        part3 = ( ( dt**2/t0**2 + 2*dt/t0 ) * m0_err ) ** 2
        return math.sqrt(part1+part2+part3)
        # Simplified form assuming dt<<t_0 and m0_err<<m0
        # return math.sqrt(dt_err**2 + ((dt * m0_err)/m0)**2) * 2 * m0/t0
    
    def calc_mass(self, file_isotope='', file_ref1='', file_ref2='',
                        t_isotope='', t_ref1='', t_ref2='',
                        t_isotope_err='', t_ref1_err='', t_ref2_err='',
                        centroid = 'mu0', online_ref = '', online_ref2 = '',
                        tweak_tofs = [0,0,0],
                        print_results = False):
        '''
        Calculates mass and mass error from either fit files in form of FitToDict objects passed to method
        or in form of time-of-flights
            - file_isotope, file_ref1, file_ref2: path to fit files to be used
            - t_isotope, t_ref1, t_ref2: time-of-flights to be used. Overwrites the ToFs fetched from fit files
            - centroid: time-of-flight centroid to be used to calculate mass ['mu0', 'numerical_mean']
            - tweak_tofs: array [tof_isotope, tof_ref1, tof_ref2] that add tof to the extracted values from the fit files to tweak the mass and see influence of tof drifts
            - online_ref and online_ref2: centroid of online references 1 and 2
            - print_results: prints short formatted results of calculation
        '''
        # Store ToFs
        self.__store_tofs(file_isotope, file_ref1, file_ref2,t_isotope, t_ref1, t_ref2, t_isotope_err, 
                            t_ref1_err, t_ref2_err, centroid, online_ref, online_ref2, tweak_tofs)
        #
        self.C_tof = self.calc_C_ToF(self.isotope_gs_t, self.ref1_t, self.ref2_t)
        self.C_tof_err = self.calc_C_ToF_err(t=self.isotope_gs_t, t_err=self.isotope_gs_t_err,
                                                       t1=self.ref1_t, t1_err=self.ref1_t_err,
                                                       t2=self.ref2_t, t2_err=self.ref2_t_err)
        #
        self.m_isotope = self.calc_sqrt_m(self.C_tof, self.m_ref1, self.m_ref2)**2
        self.m_isotope_err = self.calc_m_err(self.C_tof, self.C_tof_err, 
                                               self.m_ref1, self.m_ref1_err/self.u ,
                                               self.m_ref2, self.m_ref2_err/self.u)

        self.me_isotope = (self.m_isotope-self.A) * self.u # [keV]
        self.me_isotope_err = self.m_isotope_err * self.u # [keV]

        # 
        if print_results:
            print(f"######################\n\
# Result for {self.isotope}:\n\
######################\n\
# - C_ToF: {self.C_tof:.8f}({self.C_tof_err:.8f})\n\
# - Mass Excess ISOLTRAP: {self.me_isotope:.1f}({self.me_isotope_err:.1f})keV\n\
# - Mass Excess {self.ame_version}: {(self.m_isotope_AME-self.A)*self.u:.1f}({self.m_isotope_AME_err:.1f})keV\n\
# - Mass Difference ISOLTRAP-{self.ame_version}: {abs(self.me_isotope)-abs((self.m_isotope_AME-self.A)*self.u):.1f}keV\n\
######################")
    
    def calc_exc_energy(self, file_isotope='', file_ref1='', file_ref2='',
                        t_isotope='', t_ref1='', t_ref2='',
                        t_isotope_err='', t_ref1_err='', t_ref2_err='',
                        centroid = 'mu0', online_ref = '', online_ref2 = '', tweak_tofs = [0,0,0],
                        custom_gs = '', custom_gs_err = '',
                        is_doublet = False,
                        print_results = False):
        '''
        Calculates the excitation energy for an isomeric state from either fit files in form of FitToDict objects passed to method
        or in form of time-of-flights
            - file_isotope, file_ref1, file_ref2: path to fit files to be used
            - t_isotope, t_ref1, t_ref2: time-of-flights to be used. Overwrites the ToFs fetched from fit files
            - centroid: time-of-flight centroid to be used to calculate mass ['mu0', 'numerical_mean']
            - tweak_tofs: array [tof_isotope, tof_ref1, tof_ref2] that add tof to the extracted values from the fit files to tweak the mass and see influence of tof drifts
            - custom_gs: use ground state mass passsed in [u] to calculate excitation energy. Defaults to AME value
            - custom_gs_err: use ground state mass err passsed in [u] to calculate excitation energy. Defaults to AME value
            - is_doublet: Whether doublet formula is to be used. Does not require reference measurements. Takes file_isotope to calculate doublet. 
            - print_results: prints short formatted results of calculation
        '''
        # Store ToFs
        self.__store_tofs(file_isotope, file_ref1, file_ref2,t_isotope, t_ref1, t_ref2, t_isotope_err, 
                            t_ref1_err, t_ref2_err, centroid, online_ref, online_ref2, tweak_tofs, is_doublet)
        #
        # Calculation of excitation energy via C_tof differences e.g. mass differences
        if not is_doublet:
            self.C_tof = self.calc_C_ToF(self.isotope_m_t, self.ref1_t, self.ref2_t)
            self.C_tof_err = self.calc_C_ToF_err(t=self.isotope_m_t, t_err=self.isotope_m_t_err,
                                                           t1=self.ref1_t, t1_err=self.ref1_t_err,
                                                           t2=self.ref2_t, t2_err=self.ref2_t_err)
            #
            self.m_isomere = self.calc_sqrt_m(self.C_tof, self.m_ref1, self.m_ref2)**2
            self.m_isomere_err = self.calc_m_err(self.C_tof, self.C_tof_err, 
                                                   self.m_ref1, self.m_ref1_err/self.u ,
                                                   self.m_ref2, self.m_ref2_err/self.u)
            
            self.custom_gs = custom_gs
            self.custom_gs_err = custom_gs_err
            if self.custom_gs == '' and self.custom_gs_err == '':
                self.exc_energy = (self.m_isomere-self.m_isotope_AME) * self.u # [keV]
                self.exc_energy_err = math.sqrt(self.m_isomere_err**2 + self.m_isotope_AME_err**2)
            else:
                self.exc_energy = (self.m_isomere-self.custom_gs) * self.u # [keV]
                self.exc_energy_err = math.sqrt(self.custom_gs_err**2 + self.m_isomere_err**2)* self.u # [keV]
        # Calculation via mass-doublet
        else: 
            self.exc_energy = self.calc_exc_energy_from_doublet(t0 = self.isotope_gs_t, dt = self.isotope_m_dt, m0 = self.m_isotope_AME*self.u)
            self.exc_energy_err = self.calc_exc_energy_err_from_doublet(t0 = self.isotope_gs_t, dt = self.isotope_m_dt, m0 = self.m_isotope_AME*self.u,
                                                                    t0_err = self.isotope_gs_t_err, dt_err = self.isotope_m_dt_err, m0_err = self.m_isotope_AME_err,
                                                                    )
        # 
        if print_results:
            print(f"######################\n\
# Result for {self.isotope}-{self.state}:\n\
######################\n\
# - Excitation energy ISOLTRAP: {self.exc_energy:.1f}({self.exc_energy_err:.1f})keV\n\
# - Excitation energy NUBASE: {self.exc_energy_NUBASE:.1f}({self.exc_energy_NUBASE_err:.1f})keV\n\
# - Energy Difference ISOLTRAP-{self.nubase_version}: {abs(self.exc_energy)-abs(self.exc_energy_NUBASE):.1f}keV\n\
######################")
    
    def store_result(self, results_file, overwrite = False, tags=""):
        '''
        Appends results from calc_mass in a results file. Creates new file if file does not exist 
        Parameters:
            - results_file: .csv file to store results in
            - overwrite: checks if entry with isotope, n_revs, and tag exists, and overwrites that entry
            - tags: custom string to add to entry in results_file
        '''
        self.tags = tags
        # Create row to append
        d = {
            'A' : [self.A],
            'isotope': [self.isotope],
            'state': [self.state],
            'n_revs': [self.n_revs],
            'tags': [self.tags],
            'ref1': [self.ref1],
            'ref2': [self.ref2],
            'C_tof': [self.C_tof],
            'C_tof_err': [self.C_tof_err],
            'm_ref1': [self.m_ref1],
            'm_ref1_err': [self.m_ref1_err],
            'm_ref2': [self.m_ref2],
            'm_ref2_err': [self.m_ref2_err],
            'm_isotope_AME': [self.m_isotope_AME],
            'm_isotope_AME_err': [self.m_isotope_AME_err],
            'file_isotope': [self.file_isotope],
            'file_ref1': [self.file_ref1],
            'file_ref2': [self.file_ref2],
        }
        if self.state == "gs":
            append = {
                'm_isotope': [self.m_isotope],
                'm_isotope_err': [self.m_isotope_err],
                'me_isotope': [self.me_isotope],
                'me_isotope_err': [self.me_isotope_err],
            }
        elif self.state == "m" or self.state == "n":
            append = {
                'm_isomere': [self.m_isomere],
                'm_isomere_err': [self.m_isomere_err],
                'exc_energy': [self.exc_energy],
                'exc_energy_err': [self.exc_energy_err],
                'exc_energy_NUBASE': [self.exc_energy_NUBASE],
                'exc_energy_NUBASE_err': [self.exc_energy_NUBASE_err],
            }
        #
        d.update(append)
        #
        if self.custom_gs != '' and self.custom_gs_err != '':
            append = {
                'custom_gs': [self.custom_gs],
                'custom_gs_err': [self.custom_gs_err],
            }
            d.update(append)
        #

        # Load file if exists or create new file
        if not os.path.isfile(results_file):
            print(f"'{results_file}' does not exist, will be created...\n")
            df = pd.DataFrame.from_dict(data=d)
        else:
            df = pd.read_csv(results_file)
            # Check if entry in file exists
            line_exists = False
            idx_list = df.index[(df["isotope"]==self.isotope) & (df["n_revs"]==self.n_revs) & \
                                (df["tags"]==self.tags) & (df["state"]==self.state)].tolist()
            if len(idx_list) != 0:
                line_exists = True
            #
            if line_exists and not overwrite:
                print(f"Entry with tags '{self.tags}' and isotope '{self.isotope}' already exists. To overwrite, set overwrite flag.\n")
                return
            elif line_exists and overwrite:
                print(f"Entry with tags '{self.tags}' and isotope '{self.isotope}' already exists. Will be overwritten.\n")
                for idx in idx_list:
                    for key in d:
                        df.loc[idx,key] = d[key]
            else:
                print(f"Appending to '{results_file}'...\n")
                df2 = pd.DataFrame.from_dict(data=d) 
                df = df.append(df2, ignore_index=True)
        #
        df.to_csv(results_file, index=False)

class TOFPlot():
    """
    Handling all MR-ToF MS related ToF Histograms. Takes data and matplotlib figure/axis.
    :file:
    :fig:
    :ax:
    """
    def __init__(self, file):
        self.file = file # dataframe with tof and sweeps
        #
        # if not external plotting
        # plt.rcParams["figure.figsize"] = figsize
        self.fig, self.ax = plt.subplots(nrows=1, ncols=1, figsize=(8.6,6))
        # plt.close() # close window if created, not needed yet

    def load_utils(self, utils):
        """
        Load an external instance of MRToFUtils
        """
        self.utils = utils

    def get_binning(self, bins=1):
        """
        Adapts binning to multiples of 0.8ns, assuming that 0.8ns was used to take data (to common case)
        """
        # Get min and max tof from data frame
        minn = self.file.tof.min()
        maxx = self.file.tof.max()
        # Avoid having empty binning when maxx is equal to minn
        if minn == maxx:
            maxx += 1
        #
        return round((maxx-minn)/0.8/bins)

    def create_hist1d(self, bins = 10, log=False,
            fs_labels = 25, fs_ticks = 20, figsize = (8.6,6), ylim = (),
            style = 'errorbar', add_vlines = [], legend=False,
        ):
        '''
        Takes figure and axis and fills with 1D Histogram
        Parameters:
        - bins: Number of bins to be rebinned. Default=10
        - focus: if True, sets xlimits to first and last found peak
        - log: if Ture, sets logscale on y-axis
        - silent: if True, shows plot on canvas
        - save: if True, uses path_to_file to save plot as .pdf
        - path_to_file: path to save .pdf in
        - add_vlines: array of tof values where vlines should be added
        '''
        #
        xdata = self.file.tof
        #
        n, xe = np.histogram(xdata, bins=self.get_binning(bins))
        cx = 0.5 * (xe[1:] + xe[:-1])
        dx = np.diff(xe)

        # Plot data
        if style == 'errorbar':
            # use sqrt(n) as error, if n==1 use smaller error to avoid having inifite long error bars in log-scale
            self.ax.errorbar(cx, n, [val ** 0.5 if val != 1 else 0.75 for val in n] ,
                    ecolor='black', elinewidth=1,  
                    fmt="ok", zorder=1, label=f"Data (bins={bins})")
            # self.ax.plot(xdata, np.zeros_like(xdata)-5, "|", alpha=0.1, label = "ToF Data", zorder = 3)
            
        elif style == 'hist':
            self.ax.hist((xdata), bins=self.get_binning(bins=bins), color='grey', edgecolor='black', linewidth=0.1, label=f"Data (bins={bins})")


        # plt.errorbar(cx, n, n ** 0.5, fmt="ok", zorder=1)
        #
        if log:
            self.ax.set_yscale('log')
        #
        
        xm = np.linspace(xe[0], xe[-1], num=1000)

        if legend:
            plt.legend()

        # add vlines
        for vline in add_vlines:
            self.ax.axvline(vline, c='b', linewidth=1, zorder=3, ls = '--')

        # 
        if len(ylim) != 0:
            self.ax.set_ylim(ylim[0], ylim[1])
        
        # Add axis labels
        self.ax.set_xlabel(f'Time-of-Flight (ns)', fontsize=fs_labels)
        self.ax.set_ylabel(f'Counts per bin', fontsize=fs_labels)

        # Set ticks size 
        self.ax.tick_params(axis='both', which='major', labelsize=fs_ticks)

        # return self.fig, self.ax

    def __add_isobar_line(self, vline, text):
        """
        Add vline to axis
        """
        self.ax.axvline(vline, c='b', linewidth=1, zorder=3, ls = '--')
        self.ax.text(vline, 0.9, text, rotation=90, bbox=dict(facecolor='white', alpha=1), 
            transform =self.ax.get_xaxis_transform(),
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=14,
        )

    def add_isobars(self, nrevs, A=None, iso_list=None):
        """
        Add vlines with calculated isobars to the plot
        """
        if iso_list is not None:
            for isobar in iso_list:
                vline = self.utils.calc_ToF(self.utils.get_value(isobar, value='mass', state='gs'), nrevs)*1e3
                self.__add_isobar_line(vline, isobar)
            return

        if A is not None:
            for idx,row in self.utils.ame[self.utils.ame.A==A].iterrows():
                vline = self.utils.calc_ToF(self.utils.get_value(f'{A}{row["element"]}', value='mass', state='gs'), nrevs)*1e3
                self.__add_isobar_line(vline, f'{A}{row["element"]}')
            return

    def add_clusters(self, isotope, nrange, nrevs):
        """
        Add vlines for calculated cluster sizes
        """
        for n in nrange:
            #
            mass = n * self.utils.get_value(isotope, value='mass', state='gs')
            vline = self.utils.calc_ToF(mass, nrevs)*1e3
            self.__add_isobar_line(vline, f'n={n}')
        return


class Peaks(TOFPlot):
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
        #
        TOFPlot.__init__(self, df_file)
    
    def find_peaks(self, bins=10, peak_threshold = None, peak_min_distance = None, peak_min_height = None, peak_width_inbins = None, 
                   peak_prominence = None, peak_wlen = None):
        """  
        Arguments:
            - bins: Rebinning for faster peak finding
            - peak_threshold:
            - peak_min_distance:
            - ...
        """
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
        self.x_proj_peaks, self.peaks_info = sc.signal.find_peaks(x_proj_for_pfind, 
                                             threshold=peak_threshold, # Required vertical distance to its direct neighbouring samples, pretty useless
                                             distance=peak_min_distance, # dinstance in samples, not in value! changes when rebinned! 
                                             height=peak_min_height,
                                             width=peak_width_inbins,
                                             prominence=peak_prominence,
                                             wlen=peak_wlen)
        # Calculate some additional meta data for the found peaks
        self.n_peaks = len(self.x_proj_peaks)
        if self.n_peaks == 0:
            print(f"(Peaks.find_peaks): No peaks with current settings found, try different ones.")
            return
        self.highest_peak = self.peaks_info['peak_heights'].argmax()
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
                
    def plot(self, bins = 10, lines = True, focus=False, log=False, silent = False, 
            fs_labels = 25, fs_ticks = 20, figsize = (8.6,6), ylim = (),
            save = False, path_to_file = "peaks", style = 'errorbar', add_vlines = [],
            external = False, fig = None, ax = None):
        '''
        Plot 1D Histogram with found peaks.
        Parameters:
        - bins: Number of bins to be rebinned. Default=10
        - lines: Draws lines where peaks are found. Default=True
        - focus: if True, sets xlimits to first and last found peak
        - log: if Ture, sets logscale on y-axis
        - silent: if True, shows plot on canvas
        - save: if True, uses path_to_file to save plot as .pdf
        - path_to_file: path to save .pdf in
        - add_vlines: array of tof values where vlines should be added
        '''
        #
        # if not external plotting
        if not external:
            plt.rcParams["figure.figsize"] = figsize
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize)
        #
        if self.n_peaks == 0:
            print("Not peaks, no plots :)")
            return 0 
        #            
        xdata = self.file.tof
        n, xe = np.histogram(xdata, bins=self.get_binning(bins))
        cx = 0.5 * (xe[1:] + xe[:-1])
        dx = np.diff(xe)

        # Plot data
        if style == 'errorbar':
            # use sqrt(n) as error, if n==1 use smaller error to avoid having inifite long error bars in log-scale
            ax.errorbar(cx, n, [val ** 0.5 if val != 1 else 0.75 for val in n] ,
                    ecolor='black', elinewidth=1,  
                    fmt="ok", zorder=1, label=f"Data (bins={bins})")
            # ax.plot(xdata, np.zeros_like(xdata)-5, "|", alpha=0.1, label = "ToF Data", zorder = 3)
            
        elif style == 'hist':
            ax.hist((xdata), bins=self.get_binning(bins=bins), color='grey', edgecolor='black', linewidth=0.1, label=f"Data (bins={bins})")


        # plt.errorbar(cx, n, n ** 0.5, fmt="ok", zorder=1)
        #
        if log:
            ax.set_yscale('log')
        #
        if lines:
            for i in range(self.n_peaks):
                ax.axvline(self.pos[i], c='r', linewidth=1, zorder=3)
        
        xm = np.linspace(xe[0], xe[-1], num=1000)
        plt.legend();
        # plt.xlim(peaks.pos[0]-300, peaks.pos[0]+300)

        # add vlines
        for vline in add_vlines:
            ax.axvline(vline, c='b', linewidth=1, zorder=3, ls = '--')

        
        # Zoom in on found peaks
        if focus:
            ax.set_xlim(self.earliest_left_base-200, self.latest_right_base+200)

        # 
        if len(ylim) != 0:
            ax.set_ylim(ylim[0], ylim[1])
        
        # Add axis labels
        ax.set_xlabel(f'Time-of-Flight [ns]', fontsize=fs_labels)
        ax.set_ylabel(f'Counts per bin', fontsize=fs_labels)

        # Set ticks size 
        ax.tick_params(axis='both', which='major', labelsize=fs_ticks)

        if not external:
            plt.tight_layout()

        if not silent: 
            plt.show()
            # plt.clf()

        # return axis if external plotting is uesd
        if external:
            return fig, ax
        #
        if save:
            plt.savefig(path_to_file+".pdf", dpi=300)
            plt.clf()

    def plottest(self, bins = 10, lines = True, focus=False, log=False, silent = False, 
            fs_labels = 25, fs_ticks = 20, figsize = (8.6,6), ylim = (), legend = False,
            save = False, path_to_file = "peaks", style = 'errorbar', add_vlines = [],
            external = False, fig = None, ax = None):
        '''
        Plot 1D Histogram with found peaks.
        Parameters:
        - bins: Number of bins to be rebinned. Default=10
        - lines: Draws lines where peaks are found. Default=True
        - focus: if True, sets xlimits to first and last found peak
        - log: if Ture, sets logscale on y-axis
        - silent: if True, shows plot on canvas
        - save: if True, uses path_to_file to save plot as .pdf
        - path_to_file: path to save .pdf in
        - add_vlines: array of tof values where vlines should be added
        '''
        #
        if self.n_peaks == 0:
            print("Not peaks, no plots :)")
            return 0 
        # Create plot
        self.create_hist1d(style=style,bins=bins, log=log,
                            fs_labels = fs_labels, fs_ticks = fs_ticks, figsize = figsize, ylim = ylim,
                            # add_vlines = add_vlines,
        )
        if lines:
            for i in range(self.n_peaks):
                self.ax.axvline(self.pos[i], c='r', linewidth=1, zorder=3)
                
        # Zoom in on found peaks
        if focus:
            self.ax.set_xlim(self.earliest_left_base-200, self.latest_right_base+200)

        # 
        if len(ylim) != 0:
            self.ax.set_ylim(ylim[0], ylim[1])
        
        # Rewrite ticks based on first ToF peak found
        originial_xticks = self.ax.get_xticks()
        xticks = [tick+self.pos[0] for tick in originial_xticks-self.pos[0] + abs(originial_xticks-self.pos[0]).min()]
        self.ax.set_xticks(xticks)
        self.ax.set_xticklabels(labels=[f'{tick:.1f}' for tick in xticks-self.pos[0]])

        # Add axis labels
        self.ax.set_xlabel(f'Time-of-Flight (ns) + {self.pos[0]:.1f}ns', fontsize=fs_labels)
        self.ax.set_ylabel(f'Counts per bin', fontsize=fs_labels)

        # Set ticks size 
        self.ax.tick_params(axis='both', which='major', labelsize=fs_ticks)

        self.fig.set_size_inches(figsize)

        if not external:
            plt.tight_layout()

        if not silent: 
            plt.show()
            # plt.clf()

        # return axis if external plotting is uesd
        if external:
            return self.fig, self.ax
        #
        if save:
            plt.savefig(path_to_file, dpi=300)
            # plt.clf()
   
    def plot2d(self, bins=500, y_bins=1, focus=-1, log=False):
        """
        Plot 2D Histogram with found peaks.
        """
        # plt.rcParams["figure.figsize"] = (10,4)
        tof = self.file.tof
        sweep = self.file.sweep

        # Create plot canvas
        fig, ((ax_x, ax_y_hist),(ax_0, ax_y)) = plt.subplots(2,2,sharex='col',sharey='row', figsize=(9,9),
                                                 gridspec_kw={'height_ratios':[1,4],
                                                             'width_ratios':[4,1],
                                                             'hspace': 0.05,
                                                             'wspace':0.05})

        # faster binning for projections than histograms -> necessary in order to automatically find peaks
        x_proj = self.file.tof.value_counts(bins=self.get_binning(bins)).sort_index()
        y_proj = self.file.sweep.value_counts().sort_index()
        # y_proj = self.file.sweep.value_counts(bins=500).sort_index()

        # main plotting
        self.file.plot(x='tof', y='sweep', style='o', alpha=0.15, ms=2, ax=ax_0, label='unbinned data')
        ax_x.semilogy(x_proj.index.mid.to_numpy(), x_proj.to_numpy())
        # ax_y.plot(y_proj.to_numpy(), y_proj.index.mid.to_numpy())
        ax_y.plot(y_proj.to_numpy(), y_proj.index.to_numpy())

        # plt.plot(tof, sweep, 'o', alpha=0.15, ms=2, label='unbinned data')
        for i in range(self.n_peaks):
            ax_0.axvline(self.pos[i], c='r', linewidth=1, zorder=3)
            ax_x.axvline(self.pos[i], c='r', linewidth=1, zorder=3)
        if focus != -1:
            plt.xlim(self.pos[focus]-300, self.pos[focus]+300)

        # Counts histogram
        bar_x = self.file.sweep.value_counts().value_counts().sort_index().index
        bar_y = self.file.sweep.value_counts().value_counts().sort_index()
        ax_y_hist.bar(bar_x, bar_y)
        # ax_y_hist.set_yscale("linear")

        #
        ax_0.set_xlabel(f'Time-of-Flight [ns]', fontsize=20)
        ax_0.set_ylabel(f'Rolling sweep number', fontsize=20)
        ax_x.set_ylabel('# / 0.8 ns', fontsize=20)
        ax_y.set_xlabel('# / 10 sw.', fontsize=20)
        ax_y.xaxis.set_ticks_position('top')
        ax_y.xaxis.set_label_position('top')
        ax_y.set_xlim(0,20)
        #
        plt.show()

class softCool(Peaks, ProcessorBase):
    """
    Class for performing software cooling on 2D MR-ToF MS Data
        df_file: dataframe containing the converted .lst content 
        Inherits functionality from the peak finder 
    """
    def __init__(self, files='', df=''):
        """
        Class for performing software cooling on 2D MR-ToF MS Data
        Parameters:
            file_list: list of .csv-files with dataframe containing the converted .lst content 
        Inherits functionality from the peak finder 
        """ 
        ProcessorBase.__init__(self)
        # Inherits
            # self.files = []
            # self.data = {}
            # self.pars = {}
            # self.df_dict = {}
            # self.step = 0
        #
        self.file_passed = False
        if files != '' and not isinstance(df, pd.DataFrame):
            self.files = files
            self.file_passed = True
            # Read data
            for f in self.files:
                self.df_dict[f] = pd.read_csv(f)

        if files == '' and isinstance(df, pd.DataFrame):
            # Read data
            self.df_dict['df'] = df
        #    
        self.corr_factors = []
        self.chunk_size = 10
        self.post_cool = False
        self.tof = 0
        self.tof_cut_left = 0
        self.tof_cut_right = 0
        self.weighted_average_tof = 0
        self.verbose = 0

    def __prepare_files(self, tof, tof_cut_left=300, tof_cut_right=300, initial_align = True):
        """
        
        """
        #
        self.tof = tof
        self.tof_cut_left = tof_cut_left
        self.tof_cut_right = tof_cut_right
        if initial_align:
            self.__initial_align()
            # 
        # Sum all files
        if self.file_passed:
            self.file = self.add_all(to_csv=False)
        else:
            self.file = self.df_dict['df']
        #
        if not self.post_cool:
            self.coolfile = self.file.copy(deep=True) # copy for storing the cooled spectrum
        # Drop empty sweeps
        for idx in self.coolfile[self.coolfile.tof.isnull()].index:
            self.coolfile = self.coolfile.drop(idx)

    def __initial_align(self):
        """
        Aligns all input files to a weighted average ToF. Onl 
        Parameters:
            - file_list: array of files to be aligned with respect to each other
        """
        #
        weights = []
        averages = []
        weighted_average_tof = 0
        for f in self.df_dict:
            tof_cut = self.df_dict[f][(self.df_dict[f].tof > self.tof-self.tof_cut_left) & (self.df_dict[f].tof < self.tof+self.tof_cut_right)]
            averages.append(np.mean(tof_cut.tof))
            weights.append(len(tof_cut))
            # averages.append(np.mean(self.df_dict[f].tof))
            # weights.append(len(self.df_dict[f]))
        #
        for i in np.arange(len(weights)):
            weighted_average_tof += averages[i] * weights[i]/np.sum(weights)
        #
        if self.verbose > 0:
            print(f"ToF: {self.tof}, Left: {self.tof_cut_left}, Right: {self.tof_cut_right}\n")
            print(f"Weighted average ToF: {weighted_average_tof}.")
            print(f"Averages:\n {averages}.")
            print(f"Weights:\n {weights}.")
        i = 0
        for f in self.df_dict:
            self.df_dict[f].tof += weighted_average_tof - averages[i] 
            if self.verbose > 0:
                print(f"Applied initial ToF correction: {weighted_average_tof - averages[i]} to file {f}.")
            i += 1
    
    def save_cooled_files(self):
        """
        Saves applied correction factors to individual input files and stores them as filename + 'cooled' as file names 
        """
        if not self.file_passed: 
            print(f"(softCool.save_cooled_files): Data not passed as individual files. Can't apply cooling to individual files.")
            return 0
        # if data passed as files, go through master df and save chunks into the individual files
        idx = 0
        for file in self.df_dict:
            # get length of data from original file
            l = len(self.df_dict[file])
            # cut master df into original file length
            sub_df = self.coolfile[idx:idx+l].copy()
            # reset the sweep number
            sub_df["sweep"] = sub_df["sweep"] - sub_df["sweep"].iloc[0]
            # store
            sub_df.to_csv(file.split(".")[0]+"_cooled.csv", index=False)
            idx += l

    def moving_average(self, x, N):
        """
        Moving average filter for input array x with filter length N.
        """ 
        x_averaged = [
            # Normal case away from the edges
            1/(N-1) * np.sum(x[int(i-(N-1)/2):int(i+(N-1)/2)]) 
            if (i-(N-1)/2) >= 0 and (i+(N-1)/2) <= len(x)
            # If around the edges
            else
                # Left edge
                1/(N) * np.sum(x[int(0):int(N)]) 
                if (i-(N-1)/2) < 0
                # Right edge
                else
                1/(N) * np.sum(x[int(len(x)-N):int(len(x))]) 
            # For every element inf the array
            for 
            i in 
            np.arange(0, len(x), 1)
        ]
        #
        return x_averaged

    def calc_corr_factors(self, df, tof_cut, chunk_size=10, method="mean"):
        """
        Function for calculating correction factors
        """
        df_cut = df[(df.tof > tof_cut[0]) & (df.tof < tof_cut[1])]
        self.chunk_size = int(chunk_size)
        #
        if method=="mean" or "moving-average":
            self.corr_factors = [
                # Take mean of slice through sweeps with tof cut 
                np.median(sublist.tof) 
                if len(sublist.tof) != 0
                else 
                0 
                for sublist  
                in 
                [
                    # Apply the sweeps cut
                    df_cut[(df_cut.sweep >= i) & (df_cut.sweep < i+self.chunk_size)]
                    # if the cut is not empty
                    if len(df_cut[(df_cut.sweep >= i) & (df_cut.sweep < i+self.chunk_size)]) != 0
                    # else return a larger slice through the sweeps
                    else 
                        # If away from the edges
                        df_cut[(df_cut.sweep >= i-50) & (df_cut.sweep < i+50+self.chunk_size)]
                        if len(df_cut[(df_cut.sweep >= i-50) & (df_cut.sweep < i+50+self.chunk_size)]) != 0 and i-50 >= 0
                        # Else take slice through the beginning of file
                        else 
                            df_cut[(df_cut.sweep >= i-500) & (df_cut.sweep < i+500+self.chunk_size)]
                            if len(df_cut[(df_cut.sweep >= i-500) & (df_cut.sweep < i+500+self.chunk_size)]) != 0 and i-500 >= 0
                            # Else take slice through the beginning of file
                            else 
                            df_cut[(df_cut.sweep >= 0) & (df_cut.sweep < 500+self.chunk_size)]
                    for i 
                    in range(int(df.sweep.iloc[0]),int(df.sweep.iloc[-1])+1, self.chunk_size)
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
                         in range(0,int(df.sweep.iloc[-1])+1, self.chunk_size)
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
                         in range(0,int(df.sweep.iloc[-1])+1, self.chunk_size)
                    ]
            ]

    def cool(self, tof, tof_cut_left=300, tof_cut_right=300, method="mean", chunk_size=10, 
             post_cool = False, to_csv = False, to_file_csv = False,
             initial_align = False, use_global_mean = False,
             align_with = 0, moving_average_N = 0,
             verbose = 0):
        """
        Routine for performing the cooling
        Parameters:
            - tof: time-of-flight to calculate the correcetion factors around
            - tof_cut_left: left tof cut -> tof-tof_cut_left
            - tof_cut_right: right tof cut -> tof+tof_cut_right
            - method: 'mean', 'median', 'average', 'moving_average'; to calculate correction value for chunk
            - moving_average_N: moving average length is applied to calculating the correction factors,
                              chunk size is automatically set to 1 
            - post_cool: set true for 2nd and more cooling interations
            - initial_align: whether to initially align all files based on a global mean
            - to_csv: if file name givem, saved as csv
            - to_file_csv: if True, save applied cooling to each input file if files instead of df are passed
        Return:
            - Dataframe with cooled ToF
        """
        #    
        self.verbose = verbose
        # Prepare files for cooling
        self.__prepare_files(tof, tof_cut_left=tof_cut_left, tof_cut_right=tof_cut_right, initial_align=initial_align)
        # df to be cooled
        self.post_cool = post_cool
        if self.post_cool:
            df_to_cool = self.coolfile
        else:
            df_to_cool = self.file
        #
        self.chunk_size = chunk_size
        tof_cut = [tof-tof_cut_left, tof+tof_cut_right]
        # If moving average, set chunk size to the size of one sweep to calculate the "signal"
        #   correction factors now represent the signal 
        if method == 'moving-average':
            self.chunk_size = 1
        self.calc_corr_factors(df_to_cool, tof_cut, self.chunk_size)
        # print(self.corr_factors)

        # If moving average, apply moving average filter to "signal"
        if method == 'moving-average':
            self.corr_factors = self.moving_average(self.corr_factors, moving_average_N)
            # print(self.corr_factors)
        # print(f"Length correction factors: {len(self.corr_factors)}")
        # print(f"Length chunk sizes: {len(range(0,int(df_to_cool.sweep.iloc[-1]), int(self.chunk_size)))}")
        
        # Calculate global average
        if use_global_mean:
            mean_tof = np.mean(df_to_cool[(df_to_cool.tof > tof_cut[0]) & (df_to_cool.tof < tof_cut[1])].tof)
        else:
            mean_tof = self.corr_factors[align_with]
        print(f"Len Corr. Factors: {len(self.corr_factors)}")
        
        cooled_tofs = [
            df_to_cool[(df_to_cool.sweep >= i) & (df_to_cool.sweep < i+self.chunk_size)].tof - self.corr_factors[int(i/self.chunk_size)-int(df_to_cool.sweep.iloc[0])] + mean_tof
            # if (self.corr_factors[i] > mean_tof - 10000) 
            # else 
            # df_to_cool[(df_to_cool.sweep >= i) & (df_to_cool.sweep < i+self.chunk_size)].tof
            for i 
            in range(int(df_to_cool.sweep.iloc[0]),int(df_to_cool.sweep.iloc[-1])+1, self.chunk_size)
            # in range(0, len(self.corr_factors), 1)
        ]

        print(f"Len cooled ToFs: {len(cooled_tofs)}")
        print(f"Len hstack cooled ToFs: {len(np.hstack(cooled_tofs))}")
        # print(np.hstack(cooled_tofs))
        # print(self.coolfile)
        # print(self.coolfile.sweep.tolist())
        # print(f"Length cooled tofs: {[len(tofs) for tofs in cooled_tofs]}, sum: {sum([len(tofs) for tofs in cooled_tofs])}")
        # print(f"Length uncooled file: {len(self.file)}")
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

        if to_csv:
            self.coolfile.to_csv(to_csv, index=False)

        if to_file_csv:
            self.save_cooled_files() 

        # Export the cooled df
        return self.coolfile
        
    def plot2d(self, bins=500, focus=False, log=False, alpha=0.25, lw = 2):
        """
        Plot the 2d histograms to compare the corrected and not corrected values
        Parameters:
            - bins: number of bins to rebin
            - focus: peak number to focus on. Defaults to -1 (no focus)
            - log: log scale
            - alpha: alpha for scatter plot
            - linewidth for vertical markers
        """
        fig, (ax0, ax1) = plt.subplots(1,2,sharey='row', sharex=True, figsize=(7,7))
        tof = self.file.tof
        sweep = self.file.sweep
        # Plot unbinned and un-corrected data
        ax0.plot(tof-self.tof, sweep, 'o', alpha=alpha, ms=2, label='unbinned data')
        # Plot correction factors
        y_corr = range(int(self.coolfile.sweep.iloc[0]),int(self.coolfile.sweep.iloc[0])+len(self.corr_factors)*self.chunk_size, self.chunk_size)
        x_corr = self.corr_factors
        ax0.plot(x_corr-self.tof, y_corr, c='r', linewidth=lw, zorder=3)
        # Plot corrected data
        tof = self.coolfile.tof
        sweep = self.coolfile.sweep
        ax1.plot(tof-self.tof, sweep, 'o', alpha=alpha, ms=2, label='unbinned data')
        if self.post_cool:
            y_corr = range(int(self.coolfile.sweep.iloc[0]),int(self.coolfile.sweep.iloc[0])+len(self.corr_factors)*self.chunk_size, self.chunk_size)
            x_corr = self.corr_factors
            ax1.plot(x_corr, y_corr, c='r', linewidth=lw, zorder=3)
        #
        # for i in range(self.n_peaks):
        #     plt.axvline(self.pos[i], c='r', linewidth=lw, zorder=3)

        # Plot file limits
        current_max_sweep = 0
        for key in self.df_dict:
            current_max_sweep += self.df_dict[key].sweep.max()
            ax0.axhline(current_max_sweep, c='grey', linewidth=lw, zorder=1, ls='--')
            ax1.axhline(current_max_sweep, c='grey', linewidth=lw, zorder=1, ls='--')



        # Plot ToF to align around
        ax0.axvline(0, c='b', linewidth=lw, zorder=2, ls='-')
        ax0.axvline( - self.tof_cut_left, c='b', linewidth=lw, zorder=2, ls='--')
        ax0.axvline( + self.tof_cut_right, c='b', linewidth=lw, zorder=2, ls='--')

        ax1.axvline(0, c='b', linewidth=lw, zorder=2, ls='-')

        plt.xlabel(f"Time-of-Flight (ns) - {self.tof:.1f}ns")
        plt.ylabel("Rolling sweep number")

        if focus:
            ax0.set_xlim(-2*self.tof_cut_left, +2*self.tof_cut_right)
            ax1.set_xlim(-2*self.tof_cut_left, +2*self.tof_cut_right)        
        else:
            ax0.set_xlim(-1500, +1500)
            ax1.set_xlim(-1500, +1500)

        plt.tight_layout()
        plt.show()

class ToFExtrapolation:
    """
    Class for simple ToF extrapolation to deal with ToF drifts. Current version implements simple linear ToF 
    drift between measurements. Principle based on "sandwich measurements" similar to ICR measurements in PTMS.
    """
    def __init__(self, isotope_files, ref1_files, ref2_files, isotope, ref1, ref2, n_revs, verbose = 0):
        self.isotope_files = isotope_files
        self.ref1_files = ref1_files
        self.ref2_files = ref2_files
        self.isotope = isotope
        self.ref1 = ref1
        self.ref2 = ref2
        self.n_revs = n_revs

    def pol1d_2points(self, x1, x2, y1, y2):
        '''
        Calculates first degree polynomial through two points.
        Parameters:
            - x1: first x-value
            - x2: second x-value
            - y1: first y-value
            - y2: second y-value
        Returns:
            - m: slope
            - b: y-axis crossing
        '''
        m = (y2-y1)/(x2-x1)
        b = y2 - m*x2
        return m, b

    def pol1d(self, x, m, b):
        '''
        Returns polynomial of first degree.
        Parameters:
            - x: x value
            - m: slope
            - b: y-axis crossing
        '''
        return m*np.array(x)+b

    def calc_lin_fit_params(self, ref_files):
        """
        Calculates simple slope and y-axis crossing between two time/time-of-flight data points
        Parameters:
            - ref_files: array with two entries with paths to reference files (.mpa or .lst)
        Returns:
            - slope
            - y-axis crossing
        """
        if len(ref_files) == 2:
            y_ref = []
            x_ref = []
            fit_file0 = FitToDict(ref_files[0].split(".")[0]+"_fit.txt")
            fit_file1 = FitToDict(ref_files[1].split(".")[0]+"_fit.txt")
            y_ref.append(float(fit_file0.get_val('mu0', 'value')))
            y_ref.append(float(fit_file1.get_val('mu0', 'value')))
            x_ref.append(get_time_of_measurement(ref_files[0].split(".")[0]+".mpa", as_datetime=True))
            x_ref.append(get_time_of_measurement(ref_files[1].split(".")[0]+".mpa", as_datetime=True))
                    
            m, b = self.pol1d_2points(datetime.datetime.timestamp(x_ref[0]), 
                                 datetime.datetime.timestamp(x_ref[1]), 
                                 y_ref[0], y_ref[1])
            
            return m, b
        else:
            print("(ToFExtrapolation.calc_lin_fit_params): WARNING: Only two fit files per references allowed currently.")
            return -1,-1

    def calc_reference_extrapolation(self):
        """ 

        """
        # Get time stamps of measurement files
        x_isotope = []
        for file in self.isotope_files:
            file_base = re.split("/|\.", file)[-2]
            path_to_file = re.split("/In_", file)[0]+"/"
            # Get time of measurement
            x_isotope.append(datetime.datetime.timestamp(
                                get_time_of_measurement(path_to_file+file_base+".mpa", as_datetime=True)
                                )
                            )
        # Get ToF from measurements and referenecs
        self.t_isotope = []
        self.t_ref1 = []
        self.t_ref2 = []
        for file in self.isotope_files:
            fit = FitToDict(file.split(".")[0]+"_fit.txt")
            self.t_isotope.append(fit.get_val("mu0", "value"))
        for file in self.ref1_files:
            fit = FitToDict(file.split(".")[0]+"_fit.txt")
            self.t_ref1.append(fit.get_val("mu0", "value"))
        for file in self.ref2_files:
            fit = FitToDict(file.split(".")[0]+"_fit.txt")
            self.t_ref2.append(fit.get_val("mu0", "value"))

        # Calculate extrapolation parameters
        m_ref1, b_ref1 = self.calc_lin_fit_params(self.ref1_files)
        m_ref2, b_ref2 = self.calc_lin_fit_params(self.ref2_files)

        # Extrapolate
        self.y_ref1_extrapol = self.pol1d(x_isotope, m_ref1, b_ref1)
        self.y_ref1_extrapol_err = [5 for i in x_isotope] # needs to be calculated correctly, currently placeholder
        self.y_ref2_extrapol = self.pol1d(x_isotope, m_ref2, b_ref2)
        self.y_ref2_extrapol_err = [5 for i in x_isotope] # needs to be calculated correctly, currently placeholder

        # Apply extrapolation
        self.y_w_coorection = []
        self.y_wo_coorection_file0 = []
        self.y_wo_coorection_file1 = []
        i = 0
        for file in self.isotope_files:
            # Calculate without correction, earlier file
            self.result_wo_correction_ealier = MRToFIsotope(self.isotope, self.ref1, self.ref2, self.n_revs, state="gs")
            self.result_wo_correction_ealier.calc_mass(file_isotope=file.split(".")[0]+"_fit.txt", 
                           file_ref1 = self.ref1_files[0].split(".")[0]+"_fit.txt",
                           file_ref2 = self.ref2_files[0].split(".")[0]+"_fit.txt",
                            # t_ref1 = y_ref1_extrapol[i], t_ref1_err=y_ref1_extrapol_err[i],
                            # t_ref2 = y_ref2_extrapol[i], t_ref2_err=y_ref2_extrapol_err[i],
                            centroid = 'mu0',
                            tweak_tofs = [0,0,0],
                            print_results = False)
            self.y_wo_coorection_file0.append(
                abs(self.result_wo_correction_ealier.me_isotope)-abs((self.result_wo_correction_ealier.m_isotope_AME-self.result_wo_correction_ealier.A)*self.result_wo_correction_ealier.u)
            )
            # Calculate without correction, later file
            self.result_wo_correction_later = MRToFIsotope(self.isotope, self.ref1, self.ref2, self.n_revs, state="gs")
            self.result_wo_correction_later.calc_mass(file_isotope=file.split(".")[0]+"_fit.txt", 
                           file_ref1 = self.ref1_files[1].split(".")[0]+"_fit.txt",
                           file_ref2 = self.ref2_files[1].split(".")[0]+"_fit.txt",
                            # t_ref1 = y_ref1_extrapol[i], t_ref1_err=y_ref1_extrapol_err[i],
                            # t_ref2 = y_ref2_extrapol[i], t_ref2_err=y_ref2_extrapol_err[i],
                            centroid = 'mu0',
                            tweak_tofs = [0,0,0],
                            print_results = False)
            self.y_wo_coorection_file1.append(
                abs(self.result_wo_correction_later.me_isotope)-abs((self.result_wo_correction_later.m_isotope_AME-self.result_wo_correction_later.A)*self.result_wo_correction_later.u)
            )
            # Calculate with simple correction
            self.result_w_correction = MRToFIsotope(self.isotope, self.ref1, self.ref2, self.n_revs, state="gs")
            self.result_w_correction.calc_mass(file_isotope=file.split(".")[0]+"_fit.txt", 
                            t_ref1 = self.y_ref1_extrapol[i], t_ref1_err=self.y_ref1_extrapol_err[i],
                            t_ref2 = self.y_ref2_extrapol[i], t_ref2_err=self.y_ref2_extrapol_err[i],
                            centroid = 'mu0',
                            tweak_tofs = [0,0,0],
                            print_results = False)
            self.y_w_coorection.append(
                abs(self.result_w_correction.me_isotope)-abs((self.result_w_correction.m_isotope_AME-self.result_w_correction.A)*self.result_w_correction.u)
            )
            #
            i+=1

        # Convert UNIX time back to datetime
        self.x_isotope_time = [datetime.datetime.fromtimestamp(time) for time in x_isotope]

    def plot_extrapolation(self):
        """

        """
        simple_error_plt(y=[self.y_ref1_extrapol], y_err=[self.y_ref1_extrapol_err], x=self.x_isotope_time, label = ["Test"],
                 x_label = "Time", y_label = ["Test"], title = f"{self.n_revs}-files",
                 ref_legend_label='ISOLTRAP21 - AME20 [keV]',
                 x_share = True,
                )

    def plot_comparison(self):
        """

        """
        y = [self.y_w_coorection, self.y_wo_coorection_file0, self.y_wo_coorection_file1]
        y_err = [[5 for i in self.x_isotope_time], [5 for i in self.x_isotope_time], [5 for i in self.x_isotope_time]]
        label = ['With correction', 'Early ref file', 'Late ref file']
        y_label = ['Time-of-Flight [ns]', 'Time-of-Flight [ns]', 'Time-of-Flight [ns]']
        #
        simple_error_plt(y=y, y_err=y_err, x=[self.x_isotope_time for i in range(3)], label = label,
                 x_label = "Time", y_label = y_label, title = "2000revs Files",
                 ref_value = abs(self.result_w_correction.me_isotope)-abs((self.result_w_correction.m_isotope_AME-self.result_w_correction.A)*self.result_w_correction.u),
                 ref_err=self.result_w_correction.me_isotope_err,
                 # ref_err=calc.get_value('105In', 'mass_excess', 'error'),
                 ref_legend_label='ISOLTRAP21 - AME20 [keV]',
                 x_share = False,
                )

#############################