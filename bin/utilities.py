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

def custom_colors():
	# Color definition from coolors.co
	colors = {
		# https://coolors.co/d1495b-edae49-abe188-1a6fdf-0c0a3e
		'Fall_rgb': {'red':'#D1495B', 'orange':'#EDAE49', 'green':'#ABE188','blue':'#1A6FDF','purple':'#0C0A3E'},
		'Jo-s_favs': 
		{				
						'black': "#000000",
						'red': "#FF2D55", 
						'blue': "#00A2FF",
						'orange': "#FFCC00", 
						'green': "#61D935", 
						'grey': "#C0C0C0", 
						'purple': "#C177DA", 
						'lightblue': "#6FF1E9",
		}
	}
	return colors

def simple_error_plt(y, y_err, x='', x_labels='', \
					 label = ["ISOLTRAP"], x_label='', y_label=[''], title='', \
					 ref_value=None, ref_err=None, ref_legend_label='AME20 Error',
					 x_share = False,
					 ):
	'''
	Simple scatter plot with y error bars.
	Parameters:
	- y: y-values
	- y_err: y-errors
	- x: x-data (exclusive with x_labels)
	- x_labels: array of strings to be used as x-labels
	- x_label: x-axis labeling
	- y_label: y-axis labeling
	- title: plot title
	'''
	colors = custom_colors()
	colors_keys = list(colors['Jo-s_favs'].keys())
	mpl.rc('text', usetex=False)
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.5*1.5, 4.5))
	
	if len(x_labels) != 0 and len(x) == 0:
		x = np.arange(0,len(x_labels),1)
		plt.xticks(x, x_labels, size = 10, rotation = 0)
	elif len(x_labels) == '' and len(x) != '':
		x = x
	else:
		print("Wrong x-data passed. Can only pass either x or x_labels")

	# Loop through list of arrays passed
	twins = []
	for i in np.arange(0,len(y),1):
		#
		if i == 0:
			ax.errorbar(x, y[i], y_err[i],
					   fmt='o', color=colors['Jo-s_favs'][colors_keys[i]], zorder = 2, 
					   label=label[i], 
					   fillstyle='full', mfc="black", linewidth=2, ms =10
			)
			ax.set_ylabel(y_label[i], size=18) #fontweight='bold')
			ax.tick_params(direction='out')
		else:
			if x_share:
				twins.append(ax.twinx())
				twins[i-1].spines['right'].set_position(("axes", 0.8+i*0.2))
				color = colors['Jo-s_favs'][colors_keys[i]]
				twins[i-1].set_ylabel(y_label[i], color=color, fontsize=18)  # we already handled the x-label with ax1
				twins[i-1].tick_params(axis = 'y', direction='in', color=color)
				twins[i-1].yaxis.label.set_color(color)
				ax.spines["right"].set_color(color)
				twins[i-1].spines["right"].set_color(color)
			else:
				twins.append(ax)
				color = colors['Jo-s_favs'][colors_keys[i]]

			twins[i-1].errorbar(x, y[i], y_err[i],
				   fmt='o', color=color, zorder = 2, 
				   label=label[i], 
				   fillstyle='full', mfc=color, linewidth=2, ms =10
			)

			ax.set_zorder(twins[i-1].get_zorder()+1)
			ax.patch.set_visible(False)


	#

	# ax_t = ax.secondary_xaxis('top')
	# ax_t.set_xticks(x)

	ax.set_xlabel(x_label, size=18) #fontweight='bold')
	# ax_r = ax.secondary_yaxis('right')
	# # ax.set_xlabel("Mercury Isotopes", size=14, fontweight='bold')
	# ax_t.tick_params(axis='x', direction='in', labeltop=False)
	# ax_r.tick_params(axis='y', direction='in', labelright=False)


	if ref_value is not None and ref_err is not None:
		# Error band
		ax.fill_between(x, ref_value-ref_err, ref_value+ref_err, facecolor='0.5', alpha=0.5,
				label = ref_legend_label)

	# plt.axhspan(ref_value-ref_err, ref_value+ref_err, facecolor='0.5', alpha=0.5)

	# handles, labels = ax.get_legend_handles_labels()
	# handles = [h[0] for h in handles] # Remove errorbars from legend
	
	plt.legend()#handles, labels, fontsize=12, frameon=False, ncol=1, loc="upper left")
	
	plt.tight_layout()
	plt.title(title, fontsize=20)
	# plt.savefig("./Mercury_Masses_Comparison.pdf", dpi=300)
	plt.show()

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
			elif value == 'excitation_energy' and (state == 'n' or state == 'm'):
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
						if value == 'excitation_energy_err': 
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

	def __initialize(self, verbose = 0):
		'''
		PRIVATE: Init file, read the meta data into dict and save where the results table and fit values table start
		'''
		# open the file
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
			except Exception as e:
				if verbose > 0: 
					print(str(e) + "line:" +line)

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
		self.__initialize()
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
			elif self.state == 'n':
				self.isomere_key = 'mu2'
			else:
				print(f"(MRToFIsotope.__init__): Unknown isomeric state '{self.state}' (only 'm', 'n').")
				return -1
		#
		self.custom_gs = ''
		self.custom_gs_err = ''

	def __store_tofs(self, file_isotope='', file_ref1='', file_ref2='',
						t_isotope='', t_ref1='', t_ref2='',
						t_isotope_err='', t_ref1_err='', t_ref2_err='',
						centroid = 'mu0',
						tweak_tofs = [0,0,0]):
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
			self.isotope_gs_t = float(self.isotope_fit.get_val(centroid, 'value')) + tweak_tofs[0]
			self.isotope_gs_t_err = float(self.isotope_fit.get_val('mu0', 'error'))
			# If isomer is passed, store isomere ToF as well
			if self.state != 'gs':
				self.isotope_m_t = float(self.isotope_fit.get_val(self.isomere_key, 'value')) + tweak_tofs[0]
				self.isotope_m_t_err = float(self.isotope_fit.get_val(self.isomere_key, 'error'))
		# If raw ToFs and not fit files are passed, store ToFs directly
		elif file_isotope == '' and t_isotope != '' and t_isotope_err != '':
			self.isotope_gs_t = t_isotope + tweak_tofs[0]
			self.isotope_gs_t_err = t_isotope_err
			# If isomere ToF is passed, set gs to -1 and store isomere ToF instead 
			if self.state != 'gs':
				self.isotope_gs_t = -1 
				self.isotope_gs_t_err = -1 
				self.isotope_m_t = t_isotope + tweak_tofs[0]
				self.isotope_m_t_err = t_isotope_err
		else:
			print(f"Error input isotope")
			return -1
		#
		self.file_ref1 = file_ref1
		if file_ref1 != '' and t_ref1 == '':
			self.ref1_fit = FitToDict(file_ref1)
			self.ref1_t = float(self.ref1_fit.get_val(centroid, 'value')) + tweak_tofs[1]
			self.ref1_t_err = float(self.ref1_fit.get_val('mu0', 'error'))
		elif file_ref1 == '' and t_ref1 != '' and t_ref1_err != '':
			self.ref1_t = t_ref1 + tweak_tofs[1]
			self.ref1_t_err = t_ref1_err
		else:
			print(f"Error input ref1")
			return
		#
		self.file_ref2 = file_ref2
		if file_ref2 != '' and t_ref2 == '':
			self.ref2_fit = FitToDict(file_ref2)
			self.ref2_t = float(self.ref2_fit.get_val(centroid, 'value')) + tweak_tofs[2]
			self.ref2_t_err = float(self.ref2_fit.get_val('mu0', 'error'))
		elif file_ref2 == '' and t_ref2 != '' and t_ref2_err != '':
			self.ref2_t = t_ref2 + tweak_tofs[2]
			self.ref2_t_err = t_ref2_err
		else:
			print(f"Error input ref2")
			return

	def calc_mass(self, file_isotope='', file_ref1='', file_ref2='',
						t_isotope='', t_ref1='', t_ref2='',
						t_isotope_err='', t_ref1_err='', t_ref2_err='',
						centroid = 'mu0',
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
							t_ref1_err, t_ref2_err, centroid, tweak_tofs)
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
			print(f"Result for {self.isotope}:\n\
	- Mass Excess ISOLTRAP: {self.me_isotope:.1f}({self.me_isotope_err:.1f})keV\n\
	- Mass Excess {self.ame_version}: {(self.m_isotope_AME-self.A)*self.u:.1f}({self.m_isotope_AME_err:.1f})keV\n\
	- Mass Difference ISOLTRAP-{self.ame_version}: {abs(self.me_isotope)-abs((self.m_isotope_AME-self.A)*self.u):.1f}keV\n"
				)
	
	def calc_exc_energy(self, file_isotope='', file_ref1='', file_ref2='',
						t_isotope='', t_ref1='', t_ref2='',
						t_isotope_err='', t_ref1_err='', t_ref2_err='',
						centroid = 'mu0',
						tweak_tofs = [0,0,0],
						custom_gs = '', custom_gs_err = '',
						print_results = False):
		'''
		Calculates the excitation energy for an isomeric state from either fit files in form of FitToDict objects passed to method
		or in form of time-of-flights
			- file_isotope, file_ref1, file_ref2: path to fit files to be used
			- t_isotope, t_ref1, t_ref2: time-of-flights to be used. Overwrites the ToFs fetched from fit files
			- centroid: time-of-flight centroid to be used to calculate mass ['mu0', 'numerical_mean']
			- tweak_tofs: array [tof_isotope, tof_ref1, tof_ref2] that add tof to the extracted values from the fit files to tweak the mass and see influence of tof drifts
			- custom_gs: use ground state mass passsed in [u] to calculate excitation energy. Defaults to AME value
			- print_results: prints short formatted results of calculation
		'''
		# Store ToFs
		self.__store_tofs(file_isotope, file_ref1, file_ref2,t_isotope, t_ref1, t_ref2, t_isotope_err, 
							t_ref1_err, t_ref2_err, centroid, tweak_tofs)
		#
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
			self.exc_energy_err = math.sqrt(self.custom_gs_err**2 + self.m_isotope_AME_err**2)
		# 
		if print_results:
			print(f"Result for {self.isotope}-{self.state}:\n\
	- Excitation energy ISOLTRAP: {self.exc_energy:.1f}({self.exc_energy_err:.1f})keV\n\
	- Excitation energy NUBASE: {self.exc_energy_NUBASE:.1f}({self.exc_energy_NUBASE_err:.1f})keV\n\
	- Energy Difference ISOLTRAP-{self.nubase_version}: {abs(self.exc_energy)-abs(self.exc_energy_NUBASE):.1f}keV\n"
				)
	
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
		x_proj = self.file.tof.value_counts(bins=500).sort_index()
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
	def __init__(self, file_list):
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
		self.files = file_list
		# Read data
		for f in self.files:
			self.df_dict[f] = pd.read_csv(f)
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
		self.file = self.add_all(to_csv=False)
		#
		self.coolfile = self.file.copy(deep=True) # copy for storing the cooled spectrum
		
	def __initial_align(self, tof, tof_cut_left=300, tof_cut_right=300):
		"""
		
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
			tof_cut = self.df_dict[f][(self.df_dict[f].tof > self.tof-self.tof_cut_left) & (self.df_dict[f].tof < self.tof-self.tof_cut_right)]
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
			
	def calc_corr_factors(self, df, tof_cut, chunk_size=10, method="mean"):
		"""
		Function for calculating correction factors
		"""
		df_cut = df[(df.tof > tof_cut[0]) & (df.tof < tof_cut[1])]
		self.chunk_size = int(chunk_size)
		#
		if method=="mean":
			self.corr_factors = [
				 np.mean(sublist.tof) 
				 for sublist  
				 in 
					[
						 df_cut[(df_cut.sweep >= i) & (df_cut.sweep < i+self.chunk_size)]
						 for i 
						 in range(0,int(df.sweep.iloc[-1]), self.chunk_size)
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
						 in range(0,int(df.sweep.iloc[-1]), self.chunk_size)
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
						 in range(0,int(df.sweep.iloc[-1]), self.chunk_size)
					]
			]

	def cool(self, tof, tof_cut_left=300, tof_cut_right=300, method="mean", chunk_size=10, 
			 post_cool = False, to_csv = False, initial_align = False, use_global_mean = False,
			 align_with = 0,
			 verbose = 0):
		"""
		Routine for performing the cooling
		Parameters:
			- tof: time-of-flight to calculate the correcetion factors around
			- tof_cut_left: left tof cut -> tof-tof_cut_left
			- tof_cut_right: right tof cut -> tof+tof_cut_right
			- method: 'mean', 'median', 'average'; to calculate correction value for chunk
			- post_cool: set true for 2nd and more cooling interations
			- initial_align: whether to initially align all files based on a global mean
			- to_csv: if file name givem, saved as csv
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
		self.calc_corr_factors(df_to_cool, tof_cut, self.chunk_size)
		# print(f"Length correction factors: {len(self.corr_factors)}")
		# print(f"Length chunk sizes: {len(range(0,int(df_to_cool.sweep.iloc[-1]), int(self.chunk_size)))}")
		
		# Calculate global average
		if use_global_mean:
			mean_tof = np.mean(df_to_cool[(df_to_cool.tof > tof_cut[0]) & (df_to_cool.tof < tof_cut[1])].tof)
		else:
			mean_tof = self.corr_factors[align_with]
		#
		cooled_tofs = [
			[
				# print(cooled.corr_factors[j[1]], j[1])
				row - self.corr_factors[j[1]] + mean_tof
				for row 
				in j[0]
			]
			for j 
			in
			[
				 [df_to_cool[(df_to_cool.sweep >= i) & (df_to_cool.sweep < i+self.chunk_size)].tof, int(i/self.chunk_size)]
				 for i 
				 in range(0,int(df_to_cool.sweep.iloc[-1]), self.chunk_size)
			 ]
		]
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

		if to_csv != False:
			self.coolfile.to_csv(to_csv, index=False)

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
		y_corr = range(0,int(self.file.sweep.iloc[-1]), self.chunk_size)
		x_corr = self.corr_factors
		ax0.plot(x_corr, y_corr, c='r', linewidth=1, zorder=3)
		# Plot corrected data
		tof = self.coolfile.tof
		sweep = self.coolfile.sweep
		ax1.plot(tof, sweep, 'o', alpha=0.05, ms=2, label='unbinned data')
		if self.post_cool:
			y_corr = range(0,int(self.coolfile.sweep.iloc[-1]), self.chunk_size)
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

