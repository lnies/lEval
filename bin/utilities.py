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

from process import ProcessorBase, MPANTMpa

FILE_LOCATION = os.path.dirname(__file__)+"/"

def get_time_of_measurement(mpa_file, as_datetime = False):
	"""
	Fetches time of measurement from data file 
	Parameters:
		- mpa_file: path to .mpa file
		- as_datetime: return time as datetime object, otherwise readable string
	"""
	# Get time of measurement
	raw_file = MPANTMpa()
	pars, data, df = raw_file.process(files=[mpa_file])
	file_base = re.split("/|\.", mpa_file)[-2]
	for key in pars[file_base].keys():
		if 'report-file' in key:
			date_array = pars[file_base][key].split("written")[1].split(" ")
			if as_datetime:
				time = datetime.datetime.strptime(date_array[1]+" "+date_array[2], '%m/%d/%Y %H:%M:%S')
				return time
			else:
				return(date_array[1]+" "+date_array[2])

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
		except(Exception) as err: 
			print(f"(error in NUBASE.__init__): Could not load AME: {err}.")
		
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
			if value == 'mass' or value == 'mass_excess':
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
					elif value == 'mass_excess':
						try:
							fetched_value = float(str(self.ame.iloc[idx]['mass_excess']).strip('#'))
						except(Exception, TypeError) as err:
							print(f"(TypeError in get_value for A={A}, X={X}): {err}")
							return -1
				else:
					data = 0
					try:
						if value == 'mass': 
							data = float(str(self.ame.iloc[idx]['atomic_mass_err']).strip('#'))
						elif value == 'mass_excess':
							data = float(str(self.ame.iloc[idx]['mass_excess_err']).strip('#'))
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
				return -1   
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
				return -1
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
						centroid = 'mu0', online_ref = '', tweak_tofs = [0,0,0],
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
			self.isotope_gs_t = float(self.isotope_fit.get_val(centroid, 'value')) + tweak_tofs[0]
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
		elif file_isotope != '' and file_ref1 == '' and t_ref1 == '' and online_ref != '':
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
		# Simplified form assuming dt<<t_0
		return math.sqrt(dt_err**2 + ((dt * m0_err)/m0)**2) * 2 * m0/t0
	
	def calc_mass(self, file_isotope='', file_ref1='', file_ref2='',
						t_isotope='', t_ref1='', t_ref2='',
						t_isotope_err='', t_ref1_err='', t_ref2_err='',
						centroid = 'mu0', online_ref = '',
						tweak_tofs = [0,0,0],
						print_results = False):
		'''
		Calculates mass and mass error from either fit files in form of FitToDict objects passed to method
		or in form of time-of-flights
			- file_isotope, file_ref1, file_ref2: path to fit files to be used
			- t_isotope, t_ref1, t_ref2: time-of-flights to be used. Overwrites the ToFs fetched from fit files
			- centroid: time-of-flight centroid to be used to calculate mass ['mu0', 'numerical_mean']
			- tweak_tofs: array [tof_isotope, tof_ref1, tof_ref2] that add tof to the extracted values from the fit files to tweak the mass and see influence of tof drifts
			- print_results: prints short formatted results of calculation
		'''
		# Store ToFs
		self.__store_tofs(file_isotope, file_ref1, file_ref2,t_isotope, t_ref1, t_ref2, t_isotope_err, 
							t_ref1_err, t_ref2_err, centroid, online_ref, tweak_tofs)
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
# - Mass Excess ISOLTRAP: {self.me_isotope:.1f}({self.me_isotope_err:.1f})keV\n\
# - Mass Excess {self.ame_version}: {(self.m_isotope_AME-self.A)*self.u:.1f}({self.m_isotope_AME_err:.1f})keV\n\
# - Mass Difference ISOLTRAP-{self.ame_version}: {abs(self.me_isotope)-abs((self.m_isotope_AME-self.A)*self.u):.1f}keV\n\
######################")
	
	def calc_exc_energy(self, file_isotope='', file_ref1='', file_ref2='',
						t_isotope='', t_ref1='', t_ref2='',
						t_isotope_err='', t_ref1_err='', t_ref2_err='',
						centroid = 'mu0', online_ref = '', tweak_tofs = [0,0,0],
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
							t_ref1_err, t_ref2_err, centroid, online_ref, tweak_tofs, is_doublet)
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
		# Avoid having empty binning when maxx is equal to minn
		if minn == maxx:
			maxx += 1
		#
		return round((maxx-minn)/0.8/bins)
	
	def find_peaks(self, bins=10, peak_threshold = 0.002, peak_min_distance = 5, peak_min_height = 10, peak_width_inbins = (3,100), 
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
		self.x_proj_peaks, self.peaks_info = sc.signal.find_peaks(x_proj_for_pfind, threshold=peak_threshold,
											 distance=peak_min_distance,
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
				
	def plot(self, bins = 10, lines = True, focus=False, log=False, silent = False, save = False, path_to_file = "peaks"):
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
		'''
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
		#
		if log:
			plt.yscale('log')
		#
		if lines:
			for i in range(self.n_peaks):
				plt.axvline(self.pos[i], c='r', linewidth=1, zorder=3)
		
		xm = np.linspace(xe[0], xe[-1], num=1000)
		plt.legend();
		# plt.xlim(peaks.pos[0]-300, peaks.pos[0]+300)
		
		# Zoom in on found peaks
		if focus:
			plt.xlim(self.earliest_left_base-200, self.latest_right_base+200)
		
		# Add axis labels
		plt.xlabel(f'Time-of-Flight [ns]', fontsize=20)
		plt.ylabel(f'Counts per bin', fontsize=20)

		if not silent: 
			plt.show()
			# plt.clf()
		#
		if save:
			plt.savefig(path_to_file+".pdf", dpi=300)
			plt.clf()
		
	def plot2d(self, bins=500, focus=-1, log=False):
		"""
		Plot 2D Histogram with found peaks.
		"""
		# plt.rcParams["figure.figsize"] = (10,4)
		tof = self.file.tof
		sweep = self.file.sweep

		# Create plot canvas
		fig, ((ax_x, blank),(ax_0, ax_y)) = plt.subplots(2,2,sharex='col',sharey='row', figsize=(9,9),
												 gridspec_kw={'height_ratios':[1,4],
															 'width_ratios':[4,1],
															 'hspace': 0.05,
															 'wspace':0.05})

		# faster binning for projections than histograms -> necessary in order to automatically find peaks
		x_proj = self.file.tof.value_counts(bins=self.get_binning(bins)).sort_index()
		y_proj = self.file.sweep.value_counts(bins=500).sort_index()

		# main plotting
		self.file.plot(x='tof', y='sweep', style='o', alpha=0.15, ms=2, ax=ax_0, label='unbinned data')
		ax_x.semilogy(x_proj.index.mid.to_numpy(), x_proj.to_numpy())
		ax_y.plot(y_proj.to_numpy(), y_proj.index.mid.to_numpy())

		# plt.plot(tof, sweep, 'o', alpha=0.15, ms=2, label='unbinned data')
		for i in range(self.n_peaks):
			ax_0.axvline(self.pos[i], c='r', linewidth=1, zorder=3)
			ax_x.axvline(self.pos[i], c='r', linewidth=1, zorder=3)
		if focus != -1:
			plt.xlim(self.pos[focus]-300, self.pos[focus]+300)

		#
		ax_0.set_xlabel(f'Time-of-Flight [ns]', fontsize=20)
		ax_0.set_ylabel(f'Rolling sweep number', fontsize=20)
		ax_x.set_ylabel('# / 0.8 ns', fontsize=20)
		ax_y.set_xlabel('# / 10 sw.', fontsize=20)
		ax_y.xaxis.set_ticks_position('top')
		ax_y.xaxis.set_label_position('top')
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
		if initial_align:
			self.__initial_align(tof, tof_cut_left, tof_cut_right)
			# 
		# Sum all files
		if self.file_passed:
			self.file = self.add_all(to_csv=False)
		else:
			self.file = self.df_dict['df']
		#
		self.coolfile = self.file.copy(deep=True) # copy for storing the cooled spectrum
		# Drop empty sweeps
		for idx in self.coolfile[self.coolfile.tof.isnull()].index:
			self.coolfile = self.coolfile.drop(idx)

	def __initial_align(self, tof, tof_cut_left=300, tof_cut_right=300):
		"""
		Aligns all input files to a weighted average ToF. Onl 
		Parameters:
			- file_list: array of files to be aligned with respect to each other
		"""
		#
		self.tof = tof
		self.tof_cut_left = tof_cut_left
		self.tof_cut_right = tof_cut_right
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
				np.mean(sublist.tof) 
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
		if not self.post_cool:
			df_to_cool = self.file
		else:
			df_to_cool = self.coolfile
		#
		self.chunk_size = chunk_size
		tof_cut = [tof-tof_cut_left, tof+tof_cut_right]
		# If moving average, set chunk size to the size of one sweep to calculate the "signal"
		# 	correction factors now represent the signal 
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
		
	def plot2d(self, bins=500, focus=-1, log=False):
		"""
		Plot the 2d histograms to compare the corrected and not corrected values
		Parameters:
			- bins: number of bins to rebin
			- focus: peak number to focus on. Defaults to -1 (no focus)
			- log: log scale
		"""
		fig, (ax0, ax1) = plt.subplots(1,2,sharey='row', figsize=(7,7))
		tof = self.file.tof
		sweep = self.file.sweep
		# Plot unbinned and un-corrected data
		ax0.plot(tof, sweep, 'o', alpha=0.05, ms=2, label='unbinned data')
		# Plot correction factors
		y_corr = range(int(self.coolfile.sweep.iloc[0]),int(self.coolfile.sweep.iloc[0])+len(self.corr_factors)*self.chunk_size, self.chunk_size)
		x_corr = self.corr_factors
		ax0.plot(x_corr, y_corr, c='r', linewidth=1, zorder=3)
		# Plot corrected data
		tof = self.coolfile.tof
		sweep = self.coolfile.sweep
		ax1.plot(tof, sweep, 'o', alpha=0.05, ms=2, label='unbinned data')
		if self.post_cool:
			y_corr = range(int(self.coolfile.sweep.iloc[0]),int(self.coolfile.sweep.iloc[0])+len(self.corr_factors)*self.chunk_size, self.chunk_size)
			x_corr = self.corr_factors
			ax1.plot(x_corr, y_corr, c='r', linewidth=1, zorder=3)
		#
		# for i in range(self.n_peaks):
		#     plt.axvline(self.pos[i], c='r', linewidth=1, zorder=3)
		if focus!=-1:
			ax0.set_xlim(self.peaks.pos[focus]-300, self.peaks.pos[focus]+300)
			ax1.set_xlim(self.peaks.pos[focus]-300, self.peaks.pos[focus]+300)

		plt.xlabel("Time-of-Flight [ns]")
		plt.ylabel("Rolling sweep number")

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