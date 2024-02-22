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
from ROOT import RooGenericPdf, RooUniform, RooPolynomial, RooGaussian, RooGaussModel, RooDecay, RooFormulaVar
from ROOT import RooAddPdf, RooExtendPdf, RooMCStudy
from ROOT import Math as RootMath 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy as sc
from scipy.signal import find_peaks
from scipy.optimize import least_squares
from scipy import stats
import re
import math 
import time
import mmap

# ROOT.gSystem.Load('hyperEMG12_cxx.so')
# ROOT.gSystem.Load('hyperEMG12ex_cxx.so')

from utilities import FitToDict, Peaks, softCool, TOFPlot, custom_colors

class FitMethods(TOFPlot):
	"""
	Generic class with methods inherited by all PDF class models.
	"""
	rootitle = 'generic title'
	roounit = 'geniric unit'
	roodefs = []
	rooCompon = []
	pstate = {}
	rerun = False
	spl_draw = False
	# rx_par = QtCore.QRegExp('par*')
	# rx_unc = QtCore.QRegExp('unc*')
	# rx_name = QtCore.QRegExp('state*')
	par_names = []
	unc_names = []
	fix_names = []

	def __init__(self, lst_file, file_path, name, peaks = 'findit', verbose=False):
		"""
		Init.
		"""
		# Init ToF Plot routine
		TOFPlot.__init__(self, lst_file)
		# Init fit related variables
		self.my_name = name
		self.fit_func_name = ''
		self.this_pdf = RooGaussian()
		self.roi_mean = 0.0
		self.roi_rms = 0.0
		self.roi_counts = 0.0
		self.roi_skewness = 0.0
		#
		self.RooRealVar_dict = {}
		self.RooGenericPdf_dict = {}
		self.Var_dict = {}  
		self.this_pdf = 0
		#
		self.lst_file = lst_file
		self.bins = 100
		self.fitrange = 100
		self.xmin = 0
		self.xmax = 1
		self.results = {}
		self.fits_r = []
		self.limits = {}
		self.meta_data = {}
		self.params = {} 
		# 
		try:
			if 'tof' in lst_file.columns:
				self.xdata_unbinned = lst_file.tof # Only works if data is MR-ToF MS data
			else:
				self.xdata_unbinned = lst_file
		except:
			self.xdata_unbinned = lst_file
		#
		self.dimensions = [0,1] # dimensionality of fit function, e.g. number of shape parameters
		self.n_comps = 1 # multiplicity of components for the fit function
		self.states = [] # array of states passed as limits
		self.numerical_peak = 0
		self.numerical_FWHM_left = 0
		self.numerical_FWHM_right = 0
		self.cps = None
		self.cps_err = None
		#
		self.meta_data_keys = ["range", "swpreset", "cycles", "cmline0", "caloff", "calfact", "time_patch", ]
		#
		#
		# If no peaks are passed, find peaks within the initialization of the function
		if peaks == 'findit':
			self.peaks = Peaks(self.lst_file)
			self.peaks.find_peaks()
			if verbose:
				print("Peaks initialized by class")
				for i in range(self.peaks.n_peaks):
					print(self.peaks.peaks_info)
		elif peaks == 'nopeaks':
			self.peaks = None
		else:
			self.peaks = peaks
		#
		self.base_file_name = file_path

		# Read meta data if lst file path is provided
		# if file_path != False:
		# 	with open(file_path + ".lst",'rb') as listfile:
		# 	    mapped_file = mmap.mmap(listfile.fileno(), 0, access=mmap.ACCESS_READ)
		# 	    for key in self.meta_data_keys:
		# 		    mapped_file.seek(mapped_file.find(f'{key}'.encode('ascii')))
		# 		    self.meta_data[key] = mapped_file.readline().strip('\r\n'.encode('ascii')).decode('ascii')
		
	def minimize(self, data, xmin, xmax, datatype = 'is_th1d', minos = ROOT.kTRUE):
		"""
		Performs the maximum likelihood fit.

		:param data: the histogram to fit
		:param xmin: fit range min
		:param xmax: fit range max
		:return: fitted histogram, fit results
		"""
		x = self.RooRealVar_dict['x']
		# Depending on whether there is a region to be blinded (cut data), define some subregions
		if isinstance(xmin, list):
			ranges=""
			for i in np.arange(len(xmin)):
				lrange=f"win_{self.roodefs[0].GetName()}_range{i}"
				x.setRange(lrange, xmin[i], xmax[i])
				ranges+=lrange
				if i != len(xmin):
					ranges+=","
			# Create full range
			x.setRange(f"win_{self.roodefs[0].GetName()}_rangeFULL", xmin[0], xmax[-1])
			### Fix parameters of NLL to range
			### see tutorial https://root.cern.ch/doc/v624/rf204b__extendedLikelihood__rangedFit_8C.html
			# self.roodefs[0].fixCoefRange(f"win_{self.roodefs[0].GetName()}_rangeFULL")
			# Extend PDF to account for sub ranges
			n_in_int = 0
			for i in np.arange(len(xmin)):
				n_in_int += len(self.lst_file.tof[(self.lst_file.tof > xmin[i]) & (self.lst_file.tof < xmax[i])])

			N = RooRealVar("N", "Extended term", n_in_int, n_in_int, n_in_int)
			# N = RooRealVar("N", "Extended term", n_in_int, n_in_int/2, n_in_int*2)
			# fitmodel = self.roodefs[0]
			fitmodel = RooExtendPdf("extmodel", "Extended model", self.roodefs[0], N, ranges)
			# fitmodel = RooExtendPdf("extmodel", "Extended model", self.roodefs[0], N, f"win_{self.roodefs[0].GetName()}_rangeFULL")


		else:
			ranges = f"win_{self.roodefs[0].GetName()}"
			x.setRange(ranges, xmin, xmax)
			fitmodel = self.roodefs[0]

		# Define data type
		if datatype == 'is_th1d': 
			roohist = RooDataHist('roohist_tof', 'title', RooArgList(x), RooFit.Import(data))
		elif datatype == 'is_roodatahist':
			roohist = data
		else:
			roohist = data

		roodatahist = RooDataHist("data", "data", x, data);


		# Do fit
		result_mlkh = fitmodel.fitTo(roodatahist,
											RooFit.Range(ranges),
											# RooFit.NormRange(f"win_{self.roodefs[0].GetName()}"),
											RooFit.Minos(minos),
											RooFit.PrintEvalErrors(-1),
											# RooFit.SumW2Error(ROOT.kTRUE),
											RooFit.NumCPU(1),
											RooFit.Timer(ROOT.kTRUE),
											RooFit.Save(),
											RooFit.Verbose(ROOT.kFALSE))
		# self.pars_show()
		self.rerun = True
		return roohist, result_mlkh
	
	def get_binning(self, bins=1, pre_bin_size=0.8):
		"""
		Adapts binning to multiples of 0.8ns, assuming that 0.8ns (pre_bin_size) was used to take data (the common case)
		"""
		# Get min and max tof from data frame
		self.bins = bins
		minn = self.xdata_unbinned.min()
		maxx = self.xdata_unbinned.max()
		# Avoid having empty binning when maxx is equal to minn
		if minn == maxx:
			maxx += 1
		#
		return round((maxx-minn)/pre_bin_size/self.bins)

	def data1d2roothist(self, data, bins=1, pre_bin_size=0.8):
		'''
		Converts 1d array of data into ROOT TH1D 
		'''
		# Safe the histogrammed data in arrays
		self.xdata_unbinned = np.array(data)
		# Get min and max tof from data frame
		self.bins = bins
		self.binning = self.get_binning(bins=self.bins, pre_bin_size=pre_bin_size)
		minn = self.xdata_unbinned.min()
		maxx = self.xdata_unbinned.max()
		hist = ROOT.TH1D( 'hist', 'hist converted', self.binning, minn, maxx)
		[hist.Fill(x) for x in self.xdata_unbinned]
		self.y_data, self.x_data = np.histogram(self.xdata_unbinned, bins=self.binning)

		return hist

	def lst2roothist(self, list_file, bins=1, pre_bin_size=0.8):
		"""
		Converts .lst file to root histogram.

		:return: root histogram
		"""
		# Safe the histogrammed data in arrays
		self.xdata_unbinned = self.lst_file.tof
		# Get min and max tof from data frame
		self.bins = bins
		self.binning = self.get_binning(self.bins, pre_bin_size)
		minn = self.lst_file.tof.min()
		maxx = self.lst_file.tof.max()
		hist = ROOT.TH1D( 'hist', 'hist converted', self.binning, minn, maxx)
		[hist.Fill(x) for x in list_file.tof]
		self.y_data, self.x_data = np.histogram(self.xdata_unbinned, bins=self.binning)
		# Calculate counts per shot
		self.cps = np.average(self.lst_file.sweep.value_counts()) 
		self.cps_err = np.std(self.lst_file.sweep.value_counts())
		#
		return hist

	def data2roothist(self, x, y):
		"""
		Converts x and y data into root histogram.
		x data must be evenly spaced
		:return: root histogram
		"""
		# Get min and max tof from data frame
		self.bins = len(x)
		self.binning = len(x)
		self.x_data = x
		self.y_data = y
		minn = x.min()
		maxx = x.max()
		# create empty histogram with right binning
		hist = ROOT.TH1D( 'hist', 'hist converted', self.bins, minn, maxx)
		# fill histogram
		for el in y:
			hist.SetBinContent(i, el)
			i+=1
		return hist

	def find_peak_and_fwhm(self, mu, sigma):
		# Taken form $ROOTSYS/tutorials/fit/langaus.C
		# Seaches for the location (x value) at the maximum of the
		# combined pdf and its full width at half-maximum.
		#
		# The search is probably not very efficient, but it's a first try.
		i = 0;
		MAXCALLS = 10000
		# Search for maximum
		p = mu - 0.1 * sigma
		step = 0.05 * sigma
		lold = -2.0
		l    = -1.0
		while ( (l != lold) and (i < MAXCALLS) ): 
			i += 1
			lold = l
			x = p + step;
			# Evaluate the fit at value x 
			self.RooRealVar_dict['x'].setVal(x)
			l = self.this_pdf.getVal(ROOT.RooArgSet(self.RooRealVar_dict['x']))
			if (l < lold):
				step = -step/10
			p += step
		#
		if (i == MAXCALLS):
			print("fit_Lukas::hyperEmg::find_peak_and_fwhm: could not find peak")
			return -1, -1, -1, -1
		# Maximum found at maxx=x
		maxx = x
		fy = l/2
		# Search for right x location of fy
		p = maxx + sigma
		step = sigma
		lold = -2.0
		l    = -1e300
		i    = 0
		while ( (l != lold) and(i < MAXCALLS) ): 
			i += 1
			lold = l
			x = p + step
			# Evaluate the fit at value x 
			self.RooRealVar_dict['x'].setVal(x)
			l = abs(self.this_pdf.getVal(ROOT.RooArgSet(self.RooRealVar_dict['x'])) - fy)
			if (l > lold):
				step = -step/10
			p += step
		#
		if (i == MAXCALLS):
			print("fit_Lukas::hyperEmg::find_peak_and_fwhm: right FWHM value not found")
			return maxx, -1, -1, -1
		# Right location found 
		fxr = x
		# Search for left x location of fy
		p = maxx - 0.5 *sigma
		step = -sigma
		lold = -2.0
		l    = -1e300
		i    = 0
		while ( (l != lold) and (i < MAXCALLS) ):
			i += 1
			lold = l
			x = p + step
			# Evaluate the fit at value x 
			self.RooRealVar_dict['x'].setVal(x)
			l = abs(self.this_pdf.getVal(ROOT.RooArgSet(self.RooRealVar_dict['x'])) - fy)
			if (l > lold):
				step = -step/10
			p += step
		#
		if (i == MAXCALLS):
			print("fit_Lukas::hyperEmg::find_peak_and_fwhm: left FWHM value not found")
			return maxx, -1, -1, -1
		# 
		fxl = x;
		FWHM = fxr - fxl;
		return maxx, FWHM, fxl, fxr

	def save_fit(self, file_name=False):
		"""
		Saves fit results in plain text file
		"""
		# Write fit file
		if file_name == False: file_name = self.base_file_name+"_fit.txt"
		f = open(file_name, "w")
		f.write(f"RooFit Results\n")
		# If meta data was read from the lst file:
		f.write(f"[META-DATA]\n") 
		f.write(f"file={self.base_file_name}\n")
		for key in self.meta_data:
			f.write(f"{self.meta_data[key]}\n")
		# Fit results
		f.write(f"fit_function={self.fit_func_name}\n")
		f.write(f"binning={self.binning}\n")
		if self.fit_func_name == 'hyperEmg':
			f.write(f"dimensions={self.dimensions}\n")
			f.write(f"n_comps={self.n_comps}\n")
		f.write(f"xmin={self.xmin}\n")
		f.write(f"xmax={self.xmax}\n")
		f.write(f"numerical_peak={self.numerical_peak}\n")
		f.write(f"numerical_FWHM={self.numerical_FWHM}\n")
		f.write(f"numerical_FWHM_left={self.numerical_FWHM_left}\n")
		f.write(f"numerical_FWHM_right={self.numerical_FWHM_right}\n")
		# If counts per second were calculated
		if self.cps != None:
			f.write(f"cps={self.cps}\n")
			f.write(f"cps_err={self.cps_err}\n")
		#
		f.write(f"[RESULTS-TABLE]\n")
		f.write(f"var value error error_lo error_high var_init var_lim_lo var_lim_hi\n")
		for var in self.RooRealVar_dict:
			if var == 'x': continue
			f.write(f"{var} {self.RooRealVar_dict[var].getValV()} {self.RooRealVar_dict[var].getError()} ")
			f.write(f"{self.RooRealVar_dict[var].getErrorLo()} {self.RooRealVar_dict[var].getErrorHi()} ")
			f.write(f"{self.limits[var][0]} {self.limits[var][1]} {self.limits[var][2]}\n")
		#
		f.write(f"[FIT-VALUES]\n")
		f.write(f"tof fit\n")

		# Get fit results for x-axis
		# Get fit ranges
		if isinstance(self.xmin, list):
			plot_xmin = self.xmin[0]
			plot_xmax = self.xmax[-1]
		else:
			plot_xmin = self.xmin
			plot_xmax = self.xmax
		xm = np.linspace(plot_xmin, plot_xmax, num=5000)
		
		y_val = []
		for i in xm:
			self.RooRealVar_dict['x'].setVal(i)
			y_val.append(self.this_pdf.getVal(ROOT.RooArgSet(self.RooRealVar_dict['x'])))

		n, xe = self.y_data, self.x_data

		#
		cx = 0.5 * (xe[1:] + xe[:-1]) 
		dx = np.diff(xe)

		# Normalize values
		integral_cut = sum(y_val) * np.diff(xm)[0]
		left_n_cut = len(xe[xe<plot_xmin])
		right_n_cut = len(xe[xe<plot_xmax])

		n_cut = n[left_n_cut:right_n_cut]        
		y_val = y_val / integral_cut * sum(n_cut) * dx[0]

		it = 0
		for i in xm:
			f.write(f"{i:.3f} {y_val[it]:.9f}\n")
			it+=1
		f.close()

	def load_fit(self, file):
		'''
		Reads fit file that is generated from the FitToDict() class. Stores meta data within class scope.
		Parameters:
			- file: input file
		Return
			- xvalues of fit
			- yvalues of fit 
		'''
		self.fit = FitToDict(file)
		#
		if not 'FIT-VALUES' in self.fit.fit:
			print(f"Fit file {file} has no fit values that were exported.")
			return 0
		#
		xm = np.array(self.fit.fit['FIT-VALUES'].tof)
		y_val = np.array(self.fit.fit['FIT-VALUES'].fit)
		self.xmin = float(self.fit.fit['META-DATA']['xmin'].split(",")[0].replace("[", "").replace("]", ""))
		self.xmax = float(self.fit.fit['META-DATA']['xmax'].split(",")[-1].replace("]", "").replace("[", ""))
		self.numerical_peak = float(self.fit.fit['META-DATA']['numerical_peak'])
		self.numerical_FWHM = float(self.fit.fit['META-DATA']['numerical_FWHM'])
		self.numerical_FWHM_left = float(self.fit.fit['META-DATA']['numerical_FWHM_left'])
		self.numerical_FWHM_right = float(self.fit.fit['META-DATA']['numerical_FWHM_right'])
		self.fit_func_name = self.fit.fit['META-DATA']['fit_function']
		if self.fit_func_name == 'hyperEmg':
			self.n_comps = int(self.fit.fit['META-DATA']['n_comps'])
			self.dimensions = np.fromstring(self.fit.fit['META-DATA']['dimensions'].strip('[]'), sep=',', dtype=int)
		#
		for idx,row in self.fit.fit['RESULTS-TABLE'].iterrows():
			self.Var_dict[row['var']] = row['value']
		# Load states
		self.states = []
		for key in self.Var_dict:
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
		#
		return xm, y_val

	# Plot functions

	def plot(self, bins = 1, log=True, focus=False, xrange_mod = [800,3000], from_file = False, file_out=False, contribs = False,
		silent=False, centroids=False, components=False, carpet=False, legend=True, style='hist', lw_fit = 2,
		fs_legend = 14, fs_xlabel = 20, fs_ylabel = 20, fs_ticks = 15, figsize = (6,4), add_vlines = [], 
		prelim=False, prelimfs = 10, pre_bin_size = 0.8, fitalpha = 0.5, histalpha = 0.75, histlw = 2.0, fitzorder = 2, histzorder = 1,
		tofoffset = True,
		external = False, fig = False, ax = False,
		):
		"""  
		Wrapper for plotting the fit and data
			- bins: number of bins to rebin. Defaults to 1, e.g. no rebinning
			- log: plot y-scale in log. Defaults to false
			- focus: focus xrange on peak number
			- xrange_mod: modify xrange focus. array like with two components
			- from_file: file to read fit from if already fitted earlier .pdf. Defaults to False, fit needs to be executed beforehand 
			- file_out: name to save plot as .pdf. Defaults to False, e.g. no saving
			- contribs: if true, assumes emg components to be stored as contribs (newer fitting method)
			- silent: True shows plot, false does not show plot (for silent printing)
			- centroids: True shows centroids of Gaussian components, as well as location of FWHM of main peak. Defaults to False.
			- components: Plots components to EGH as grey dashed lines
			- carpet: If true, plots carpet 
			- legend: Plot legend if
			- add_vlnies: Array-like, adds vlines at ToF values entered throug array
			- prelim: adds preliminary watermark
			- prelimfs: preliminary watermark font size
			- fitalpha: alpha of fit
			- fitzorder, histzorder: zorders for plot
			- tofoffset: if True, sets ticks to 0 at the first found or fitted peak, writes offset into xlabel
			- external: if rountine is to be used to plot onto an external canvas. requires fig and axis to be passed. axis will be then returned
		"""
		# if not external plotting
		if not external:
			plt.rcParams["figure.figsize"] = figsize
			fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize)
		#
		if len(self.lst_file) == 0:
			print("Fit not excecuted yet or failed.")
			return 0 
		#       
		self.bins = bins

		# If fit is to be read from file
		if from_file:
			# Read fit file
			xm, y_val = self.load_fit(from_file)
			# self.xdata_unbinned = xm

		# Get fit and prepare fit parameters for plotting, from the call_pdf during the same program excecution
		if not from_file:
			# X-Axis for fit
			# Get fit ranges
			if isinstance(self.xmin, list):
				plot_xmin = self.xmin[0]
				plot_xmax = self.xmax[-1]
			else:
				plot_xmin = self.xmin
				plot_xmax = self.xmax
			xm = np.linspace(plot_xmin, plot_xmax, num=5000)
			# Get fit results for x-axis
			y_val = []
			for i in xm:
				self.RooRealVar_dict['x'].setVal(i)
				y_val.append(self.this_pdf.getVal(ROOT.RooArgSet(self.RooRealVar_dict['x'])))
				# print(self.this_pdf.getVal(ROOT.RooArgSet(self.RooRealVar_dict['x'])))	

		# Histogram data
		n, xe = np.histogram(self.xdata_unbinned, bins=self.get_binning(self.bins, pre_bin_size))
		# n = n/len(xdata)
		cx = 0.5 * (xe[1:] + xe[:-1]) 
		dx = np.diff(xe)

		# From infer the ground state / isomeric state character of peaks passed
		#self.__infer_states(self.n_comps, fromm='fit-file')

		# Get fit ranges
		if isinstance(self.xmin, list):
			plot_xmin = self.xmin[0]
			plot_xmax = self.xmax[-1]
		else:
			plot_xmin = self.xmin
			plot_xmax = self.xmax

		# zero the tof if True
		tof_zero = self.numerical_peak if tofoffset else 0


		# Plot data
		if style == 'errorbar':
			# use sqrt(n) as error, if n==1 use smaller error to avoid having inifite long error bars in log-scale
			ax.errorbar(cx - tof_zero, n, [val ** 0.5 if val != 1 else 0.75 for val in n] ,
					ecolor='black', elinewidth=1,  
					fmt="ok", zorder=1, label=f"Data (bins={bins})"
			)
		elif style == 'hist':
			ax.hist((self.xdata_unbinned - tof_zero), bins=self.get_binning(bins=bins, pre_bin_size=pre_bin_size), 
				color='grey', edgecolor='black', linewidth=histlw,  alpha = histalpha, histtype='step', label=f"Data (bins={bins})",
				zorder=histzorder)

		# Normalize values
		integral_cut = sum(y_val) * np.diff(xm)[0]
		left_n_cut = len(xe[xe<plot_xmin])
		right_n_cut = len(xe[xe<plot_xmax])
		n_cut = n[left_n_cut:right_n_cut]        
		y_val = y_val / integral_cut * sum(n_cut) * dx[0]

		# Plot fit	
		ax.plot(xm - tof_zero, 
				y_val, label=f"{self.fit_func_name}({self.dimensions[0]},{self.dimensions[1]})", c='r', linewidth=lw_fit, alpha=fitalpha,
				zorder = fitzorder,
		)
		
		ax.fill_between(xm - tof_zero, y_val, alpha=0.4, color='grey')

		# Plot 'carpet'
		if carpet:
			if log:
				ax.plot(self.xdata_unbinned - tof_zero, 
						 np.zeros_like(self.xdata_unbinned)+0.9, "|", alpha=0.1, zorder = 3)
			ax.plot(self.xdata_unbinned - tof_zero, 
					 np.zeros_like(self.xdata_unbinned)-5, "|", alpha=0.1, zorder = 3)

		# Plot numerical peak position and FWHM
		if centroids:
			ax.axvline(self.numerical_peak - tof_zero, c='r', linewidth=1, zorder=3)
			ax.axvline(self.numerical_FWHM_left - tof_zero, c='blue', linewidth=1, zorder=3)
			ax.axvline(self.numerical_FWHM_right - tof_zero, c='blue', linewidth=1, zorder=3)
			# Plot center of gaussian component
			for comp in range(0,self.n_comps,1):
				if f'mu{comp}' in self.Var_dict.keys():
					ax.axvline(self.Var_dict[f'mu{comp}'] - tof_zero, c='r', linewidth=1, zorder=3)
				if f'mu{comp}-E0' in self.Var_dict.keys():
					ax.axvline(self.Var_dict[f'mu{comp}']+self.Var_dict[f'mu{comp}-E0'] - tof_zero, c='r', linewidth=1, zorder=3)
				if f'mu{comp}-E1' in self.Var_dict.keys():
					ax.axvline(self.Var_dict[f'mu{comp}']+self.Var_dict[f'mu{comp}-E1'] - tof_zero, c='r', linewidth=1, zorder=3)
				else:
					continue
					
		# Plot components
		if components:

			n_ratios = self.n_comps - 1
			n_contribs = np.sum(self.dimensions) - 1 

			# print([self.Var_dict[f'ratio{r}'] for r in np.arange(0,n_ratios,1)])
			# print(sum([self.Var_dict[f'ratio{r}'] for r in np.arange(0,n_ratios,1)]))
			# print(1-sum([self.Var_dict[f'ratio{r}'] for r in np.arange(0,n_ratios,1)]))

			# Loop through components as the ratios are assigned: inner loop must be the different species while outer loop is the egh component
			# Negative egh components
			i_contrib = 0 
			for dim in range(0,self.dimensions[0],1):
				comp = 0
				for state in self.states:
					if state["state"] == 'gs':
						y_val = self.neg_func(x=xm,mu=self.Var_dict[state["peak"]], 
												sigma=self.Var_dict[f'sigma'], 
												ntau=self.Var_dict[f'ntau{dim}'])
					else:
						y_val = self.neg_func(x=xm, mu=self.Var_dict[state["peak"]]+self.Var_dict[state["peak"]+"-"+state["state"]], 
												sigma=self.Var_dict[f'sigma'], 
												ntau=self.Var_dict[f'ntau{dim}'])
					# normalize accounting for fit ratios
					integral_cut = sum(y_val) * np.diff(xm)[0]
					y_val = y_val / integral_cut * sum(n_cut) * dx[0]
					
					if (comp < self.n_comps-1):
						ratio = self.Var_dict[f'ratio{comp}']
					else:
						ratio = 1-sum([self.Var_dict[f'ratio{r}'] for r in np.arange(0,self.n_comps-1,1)])
					comp += 1

					if i_contrib < n_contribs:
						contrib = self.Var_dict[f'contrib{i_contrib}'] 
					else:
						contrib = (1-sum([self.Var_dict[f'contrib{r}'] for r in np.arange(0,n_contribs,1)]))
					#
					y_val *= contrib * ratio

					# print(f"Ratio: {ratio}, Contrib: {contrib}, Dimension {dim}, -, peak" + state["peak"] + state["state"])

					# Plot
					ax.plot(xm - tof_zero, 
							 y_val, c='grey', ls="--", zorder=2, linewidth=2.25,
							 #label=f"Neg. component {comp}:{dim})"
							 )
					# Shade area under component
					ax.fill_between(xm - tof_zero, y_val, alpha=0.4, color='grey')

				#
				i_contrib += 1 
				#
			# Positive egh components
			for dim in range(0,self.dimensions[1],1):
				comp = 0
				for state in self.states:
					if state["state"] == 'gs':
						y_val = self.pos_func(x=xm,mu=self.Var_dict[state["peak"]], 
												sigma=self.Var_dict[f'sigma'], 
												ptau=self.Var_dict[f'ptau{dim}'])
					else:
						y_val = self.pos_func(x=xm, mu=self.Var_dict[state["peak"]]+self.Var_dict[state["peak"]+"-"+state["state"]], 
												sigma=self.Var_dict[f'sigma'], 
												ptau=self.Var_dict[f'ptau{dim}'])
					# normalize accounting for fit ratios
					integral_cut = sum(y_val) * np.diff(xm)[0]
					y_val = y_val / integral_cut * sum(n_cut) * dx[0]
					
					if (comp < self.n_comps-1):
						ratio = self.Var_dict[f'ratio{comp}']
					else:
						ratio = 1-sum([self.Var_dict[f'ratio{r}'] for r in np.arange(0,self.n_comps-1,1)])
					comp += 1

					if i_contrib < n_contribs:
						contrib = self.Var_dict[f'contrib{i_contrib}'] 
					else:
						contrib = (1-sum([self.Var_dict[f'contrib{r}'] for r in np.arange(0,n_contribs,1)]))
					#
					y_val *= contrib * ratio

					# print(f"Ratio: {ratio}, Contrib: {contrib}, Dimension {dim}, +, peak" + state["peak"] + state["state"])
					
					# Plot
					ax.plot(xm - tof_zero, 
							 y_val, c='grey', ls="--", zorder=2, linewidth=2.25,
							 #label=f"Neg. component {comp}:{dim})"
							 )
					# Shade area under component
					ax.fill_between(xm - tof_zero, y_val, alpha=0.4, color='grey')

				#
				i_contrib += 1 
				#

		# add vlines
		for vline in add_vlines:
			ax.axvline(vline-tof_zero, c='blue', linewidth=1, zorder=3, ls = '--')
		
		# Get y axis limits
		ylims = ax.get_ylim()
		if log:
			ax.set_yscale("log")
			ax.set_ylim(0.2,2*ylims[1])

		# Zoom in on found peaks if peaks were searched
		if self.peaks is not None:
			if self.peaks.n_peaks != 0:
				ax.set_xlim(self.peaks.earliest_left_base - xrange_mod[0] - tof_zero, 
						 self.peaks.latest_right_base + xrange_mod[1] - tof_zero)
				if focus:
					ax.set_xlim(self.xmin-tof_zero, self.xmax-tof_zero)
		if self.peaks is None:
			print(self.xmin-tof_zero- xrange_mod[0])
			ax.set_xlim(self.xmin-tof_zero- xrange_mod[0], self.xmax-tof_zero+ xrange_mod[1])

		# Set ticks size 
		ax.tick_params(axis='both', which='major', labelsize=fs_ticks)

		# Add axis labels
		if tofoffset:
			ax.set_xlabel(f'Time-of-Flight (ns) - {tof_zero:.1f} ns', fontsize=fs_xlabel)
		else:
			ax.set_xlabel(f'Time-of-Flight (ns)', fontsize=fs_xlabel)
			ax.get_xaxis().get_major_formatter().set_useOffset(False)
			ax.get_xaxis().get_major_formatter().set_scientific(False)
		ax.set_ylabel(f'Counts per bin', fontsize=fs_ylabel)

        # Check if there are lines passed and plot them
		if len(self.vlines) != 0:
			#
			for vline,text in zip(self.vlines, self.vlines_text):
				self.add_isobar_line(vline-tof_zero, text, external = True, fig = fig, ax = ax, linezorder=99)
			# Rescale y axis
			ylims = ax.get_ylim()
			ax.set_ylim(ylims[0], ylims[1]*10)


		# Format Legend
		if legend:
			plt.legend(fontsize=fs_legend)

		# Print preliminary on top
		if prelim:
			ax.text(0.5, 0.5, 'PRELIMINARY', transform=ax.transAxes,
				fontsize=prelimfs, color='red', alpha=0.3,
				ha='center', va='center', rotation=30)
		if not external:
			plt.tight_layout()

		# Save plot
		if file_out != False:
			print(f"Plot fit save as {file_out}")
			plt.savefig(file_out, dpi=600, transparent=True)


		# Show plot on canvas
		if not silent:
			plt.show()

		# Clear canvas to avoid printing on top of other plot in batch mode
		if silent and not external:
			plt.clf()

		# return axis if external plotting is uesd
		# if external:
		# 	return fig, ax

class Gauss(FitMethods):
	"""
	Class handling the Gaussian model.
	"""
	def __init__(self, lst_file, peaks = 'findit', verbose=False):
		"""
		Initializes the class and number of parameters.
		"""        
		# QtWidgets.QWidget.__init__(self)
		FitMethods.__init__(self, lst_file, 'gauss', 2, peaks=peaks, verbose=verbose)
		# self.setupUi(self)
		# self.update_par_list()

	def call_pdf(self, xmin, xmax, index=0):
		"""
		Setup Variable, Parameters and PDF and performs fitting of the PDF to ROI data.
		:param roi_hist: histogram to fit
		:param xmin: range min for the fit
		:param xmax: range max for the fit
		:param index: peak index for display
		:return list: containing all relevant information about the ROI histogram, PDF, parameters,
					  PDF components and fit results
		"""
		# self.hist = self.lst2roothist(self.lst_file, bins=2000)
		self.hist.GetXaxis().SetRangeUser(xmin, xmax)
		self.compose_title_unit(self.hist.GetXaxis().GetTitle())
		
		# self.set_value(0, self.hist.GetMean(1))
		# self.set_limits(0, xmin, xmax)
		# self.set_value(1, self.hist.GetStdDev(1))
		# self.set_limits(1, 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1))
		#
		self.x = RooRealVar("x", self.rootitle, xmin, xmax, self.roounit)
		self.mean = RooRealVar('mean{}'.format(index), 'mean{}'.format(index),
						  self.hist.GetMean(1), xmin, xmax)
		self.sigma = RooRealVar('sigma{}'.format(index), 'sigma{}'.format(index),
						   self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1))
		self.this_pdf = RooGaussian('gauss_{}'.format(index), 'Gaussian Model {}'.format(index), self.x, self.mean, self.sigma)
		# norm = RooRealVar('events{}'.format(index), "Nb of Events", self.roi_counts)
		# self.ext_this_pdf = RooExtendPdf('ext_gauss{}'.format(index), 'Ext Gaussian', self.this_pdf, norm)
		self.roodefs = [self.this_pdf, self.x, self.mean, self.sigma]
		# self.fix_var(self.roodefs[2:], index)
		self.this_roohist, self.fit_results = self.minimize(self.hist, xmin, xmax)
		return [self.this_roohist, self.this_pdf, self.nb_pars,
				self.x, self.mean, self.sigma, self.pstate,
				self.my_name, self.spl_draw, self.rooCompon, self.fit_results]
	
class LanGauss(FitMethods):
	"""
	Class handling the Landau - Gaussian convolution model.
	Based on: https://root.cern/doc/master/rf208__convolution_8py.html
	"""
	def __init__(self, lst_file, file_path=False, peaks = 'findit', verbose=False):
		"""
		Initializes the class and number of parameters.
		"""        
		# QtWidgets.QWidget.__init__(self)
		FitMethods.__init__(self, lst_file, 'gauss', 2, peaks=peaks, verbose=verbose)
		self.fit_func_name = "Langauss"
		# self.setupUi(self)
		# self.update_par_list()

	def init_limits(self):
		self.limits = {
			'ml':[self.hist.GetMean(1), self.xmin, self.xmax],
			'mg':[self.hist.GetMean(1), self.xmin, self.xmax],
			'sl':[self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1)],
			'sg':[self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1)],
			'ratio':[0.5,0,1],
		}

	def call_pdf(self, xmin, xmax, limits = False, minos = True):
		"""
		Setup Variable, Parameters and PDF and performs fitting of the PDF to ROI data.
		:param roi_hist: histogram to fit
		:param xmin: range min for the fit
		:param xmax: range max for the fit
		:param index: peak index for display
		:return list: containing all relevant information about the ROI histogram, PDF, parameters,
					  PDF components and fit results
		"""
		# self.hist = self.lst2roothist(self.lst_file, bins=1)
		self.xmin=xmin
		self.xmax=xmax
		self.hist.GetXaxis().SetRangeUser(xmin, xmax)
		self.compose_title_unit(self.hist.GetXaxis().GetTitle())

		# Initialize limits
		if not limits:
			self.init_limits()

		# Build variable dicitonary		
		self.RooRealVar_dict['x'] = RooRealVar("x", self.rootitle, xmin, xmax, self.roounit)
		#
		self.RooRealVar_dict['ml'] = RooRealVar('ml', 'mean landau', self.limits['ml'][0], self.limits['ml'][1], self.limits['ml'][2])
		self.RooRealVar_dict['sl'] = RooRealVar('sl', 'sigma landau', self.limits['sl'][0], self.limits['sl'][1], self.limits['sl'][2])
		self.landau = ROOT.RooLandau("landau", "landau", self.RooRealVar_dict['x'], self.RooRealVar_dict['ml'], self.RooRealVar_dict['sl'])
		#
		self.RooRealVar_dict['mg'] = RooRealVar('mg', 'mean gauss', self.limits['mg'][0], self.limits['mg'][1], self.limits['mg'][2])
		self.RooRealVar_dict['sg'] = RooRealVar('sg', 'sigma gauss', self.limits['sg'][0], self.limits['sg'][1], self.limits['sg'][2])
		self.gauss = ROOT.RooGaussian("gauss", "gauss", self.RooRealVar_dict['x'], self.RooRealVar_dict['mg'], self.RooRealVar_dict['sg'])
		#
		self.this_pdf = ROOT.RooFFTConvPdf("lxg", "landau (X) gauss", self.RooRealVar_dict['x'], self.landau, self.gauss)
		
		# norm = RooRealVar('events{}'.format(index), "Nb of Events", self.roi_counts)
		# self.ext_this_pdf = RooExtendPdf('ext_gauss{}'.format(index), 'Ext Gaussian', self.this_pdf, norm)
		self.roodefs = [self.this_pdf, self.RooRealVar_dict]

		# Calculate numerical position of maximum and FWHM
		mu = self.RooRealVar_dict['mg'].getValV()
		sigma = self.RooRealVar_dict['sg'].getValV()
		self.numerical_peak, self.numerical_FWHM, self.numerical_FWHM_left, self.numerical_FWHM_right = self.find_peak_and_fwhm(mu, sigma)
		# Store fit results in dict
		for key in self.RooRealVar_dict:
			self.Var_dict[key] = self.RooRealVar_dict[f'{key}'].getValV()


		# self.fix_var(self.roodefs[2:], index)
		self.this_roohist, self.fit_results = self.minimize(self.hist, xmin, xmax, minos=ROOT.kTRUE)
		#
		return [self.this_roohist, self.this_pdf, self.RooRealVar_dict, self.pstate,
				self.my_name, self.spl_draw, self.rooCompon, self.fit_results]

	def plot(self, figsize=(6,4), silent=False):


		if isinstance(self.xmin, list):
			plot_xmin = self.xmin[0]
			plot_xmax = self.xmax[-1]
		else:
			plot_xmin = self.xmin
			plot_xmax = self.xmax
		xm = np.linspace(plot_xmin, plot_xmax, num=5000)

		y_val = []
		for i in xm:
			self.RooRealVar_dict['x'].setVal(i)
			y_val.append(self.this_pdf.getVal(ROOT.RooArgSet(self.RooRealVar_dict['x'])))

		n, xe = self.y_data, self.x_data

		#
		cx = 0.5 * (xe[1:] + xe[:-1]) 
		dx = np.diff(xe)

		# Normalize values
		integral_cut = sum(y_val) * np.diff(xm)[0]
		left_n_cut = len(xe[xe<plot_xmin])
		right_n_cut = len(xe[xe<plot_xmax])

		n_cut = n[left_n_cut:right_n_cut]        
		y_val = y_val / integral_cut * sum(n_cut) * dx[0]


		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize)
		ax.plot(xm, y_val, color='r')	
		ax.scatter(self.x_data, self.y_data)

		# Show plot on canvas
		if silent == False:
			plt.show()

class Polynomial(FitMethods):
	"""
	Simple polynomial fit
	"""
	def __init__(self, lst_file, file_path=False, peaks = 'findit', verbose=False):
		"""
		Initializes the class and number of parameters.
		"""        
		# QtWidgets.QWidget.__init__(self)
		FitMethods.__init__(self, lst_file, 'polynomial', 2, peaks=peaks, verbose=verbose)
		self.fit_func_name = "polynomial"
		# self.setupUi(self)
		# self.update_par_list()

	def init_limits(self):
		self.limits = {
			'ml':[self.hist.GetMean(1), self.xmin, self.xmax],
			'mg':[self.hist.GetMean(1), self.xmin, self.xmax],
			'sl':[self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1)],
			'sg':[self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1)],
			'ratio':[0.5,0,1],
		}

	def call_pdf(self, xmin, xmax, limits = False, minos = True):
		"""
		Setup Variable, Parameters and PDF and performs fitting of the PDF to ROI data.
		:param roi_hist: histogram to fit
		:param xmin: range min for the fit
		:param xmax: range max for the fit
		:param index: peak index for display
		:return list: containing all relevant information about the ROI histogram, PDF, parameters,
					  PDF components and fit results
		"""
		# self.hist = self.lst2roothist(self.lst_file, bins=1)
		self.xmin=xmin
		self.xmax=xmax
		self.hist.GetXaxis().SetRangeUser(xmin, xmax)
		self.compose_title_unit(self.hist.GetXaxis().GetTitle())

		# Initialize limits
		if not limits:
			self.init_limits()
		else:
			self.limits=limits



		# Build variable dicitonary		
		self.RooRealVar_dict['x'] = RooRealVar("x", self.rootitle, xmin, xmax, self.roounit)
		#
		self.RooRealVar_dict['a0'] = RooRealVar('a0', 'a0', self.limits['a0'][0], self.limits['a0'][1], self.limits['a0'][2])
		self.RooRealVar_dict['a1'] = RooRealVar('a1', 'a1', self.limits['a1'][0], self.limits['a1'][1], self.limits['a1'][2])
		# self.RooRealVar_dict['A'] = RooRealVar('A', 'A', self.limits['a1'][0], self.limits['a1'][1], self.limits['a1'][2])
		#
		self.this_pdf = ROOT.RooPolynomial("pol", "pol", self.RooRealVar_dict['x'], RooArgList(self.RooRealVar_dict['a0'],self.RooRealVar_dict['a1']), lowestOrder = 0)
		
		# norm = RooRealVar('events{}'.format(index), "Nb of Events", self.roi_counts)
		# self.ext_this_pdf = RooExtendPdf('ext_gauss{}'.format(index), 'Ext Gaussian', self.this_pdf, norm)
		self.roodefs = [self.this_pdf, self.RooRealVar_dict]

		# Calculate numerical position of maximum and FWHM
		mu = self.xmin
		sigma = 1
		self.numerical_peak, self.numerical_FWHM, self.numerical_FWHM_left, self.numerical_FWHM_right = self.find_peak_and_fwhm(mu, sigma)
		# Store fit results in dict
		for key in self.RooRealVar_dict:
			self.Var_dict[key] = self.RooRealVar_dict[f'{key}'].getValV()


		# self.fix_var(self.roodefs[2:], index)
		self.this_roohist, self.fit_results = self.minimize(self.hist, xmin, xmax, minos=ROOT.kTRUE)
		#
		return [self.this_roohist, self.this_pdf, self.RooRealVar_dict, self.pstate,
				self.my_name, self.spl_draw, self.rooCompon, self.fit_results]

class multiGauss(FitMethods):
	"""
	Class handling the Gaussian model.
	"""
	def __init__(self, lst_file, file_path = False, peaks = 'findit', verbose=False):
		"""
		Initializes the class and number of parameters.
		"""        
		FitMethods.__init__(self, lst_file, file_path=file_path, name='multiGauss', peaks=peaks, verbose=verbose)
		# self.setupUi(self)
		# self.update_par_list()
		self.fit_func_name = "multiGauss"
		
		if file_path == False: 
			self.base_file_name = 'path_to_file'
		else:
			self.base_file_name = file_path

	def init_limits(self):
		self.limits = {
			'mean0':[self.hist.GetMean(1), self.xmin, self.xmax],
			'mean1':[self.hist.GetMean(1), self.xmin, self.xmax],
			'sigma':[self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1)],
			'ratio':[0.5,0,1],
		}

	def call_pdf(self, xmin, xmax, bins=1, dimensions = [0,1], n_comps = 2, limits=False,
				 minos = ROOT.kFALSE,
		):
		"""
		Setup Variable, Parameters and PDF and performs fitting of the PDF to ROI data.
		:param roi_hist: histogram to fit
		:param xmin: range min for the fit
		:param xmax: range max for the fit
		:return list: containing all relevant information about the ROI histogram, PDF, parameters,
					  PDF components and fit results
		"""
		# Convert histogram
		self.bins = bins
		#
		self.dimensions = dimensions 
		self.n_comps = n_comps
		self.xmin = xmin
		self.xmax = xmax
		self.compose_title_unit(self.hist.GetXaxis().GetTitle())
		# 
		if not limits:
			self.init_limits()
		else:
			self.limits = limits
		#
		
		# self.set_value(0, self.hist.GetMean(1))
		# self.set_limits(0, self.xmin, self.xmax)
		# self.set_value(1, self.hist.GetStdDev(1))
		# self.set_limits(1, 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1))
		#

		# Fill variable dictionary:
		if isinstance(xmin, list):
			self.RooRealVar_dict = {
				'x': RooRealVar("x", self.rootitle, self.xmin[0], self.xmax[-1], self.roounit)
			}
		else:
			self.RooRealVar_dict = {
				'x': RooRealVar("x", self.rootitle, self.xmin, self.xmax, self.roounit)
			}

		self.RooRealVar_dict['sigma'] = RooRealVar('sigma', 'sigma',
						   self.limits["sigma"][0], self.limits["sigma"][1], self.limits["sigma"][2])
		#
		for n in np.arange(0,self.n_comps,1):
			self.RooRealVar_dict[f'mu{n}'] = RooRealVar(f'mu{n}', f'mu{n}',
							  self.limits[f'mu{n}'][0], self.limits[f'mu{n}'][1], self.limits[f'mu{n}'][2])
			self.RooGenericPdf_dict[f'sig{n}'] = RooGaussian(f'gauss{n}', f'Gaussian Model {n}', 	
																							self.RooRealVar_dict['x'], 
																							self.RooRealVar_dict[f'mu{n}'], 
																							self.RooRealVar_dict['sigma'])
		# Build component ratios
		for j in np.arange(0,self.n_comps-1,1): 
			var_name = f"ratio{j}"
			self.RooRealVar_dict[var_name] = RooRealVar(var_name, var_name, self.limits[var_name][0], self.limits[var_name][1], self.limits[var_name][2])


		# put all pdfs together
		all_pdfs = RooArgList()
		for pdf in self.RooGenericPdf_dict:
			all_pdfs.add(self.RooGenericPdf_dict[pdf])
		# put all ratios together
		all_ratios = RooArgList()
		for var in self.RooRealVar_dict:
			m = re.search('ratio*', var)
			if m:
				all_ratios.add(self.RooRealVar_dict[var])

		# Make master pdf
		self.this_pdf = RooAddPdf(f'multiGaussian', f'Multiple Gaussians Superposed', all_pdfs, all_ratios, recursiveFraction = ROOT.kFALSE)

		# norm = RooRealVar('events{}', "Nb of Events", self.roi_counts)
		# self.ext_this_pdf = RooExtendPdf('ext_gauss{}'.format(index), 'Ext Gaussian', self.this_pdf, norm)
		self.roodefs = [self.this_pdf, self.RooRealVar_dict, self.RooGenericPdf_dict]
		# self.fix_var(self.roodefs[2:], index)

		if minos: 
			self.this_roohist, self.fit_results = self.minimize(self.hist, self.xmin, self.xmax, datatype='is_th1d')
		else:
			self.this_roohist, self.fit_results = self.minimize(self.hist, self.xmin, self.xmax, datatype='is_th1d', minos = ROOT.kFALSE)

		# Store fit results in dict
		for key in self.RooRealVar_dict:
			self.Var_dict[key] = self.RooRealVar_dict[f'{key}'].getValV()

		# Calculate numerical position of maximum and FWHM
		mu = self.RooRealVar_dict['mu0'].getValV()
		sigma = self.RooRealVar_dict['sigma'].getValV()
		self.numerical_peak, self.numerical_FWHM, self.numerical_FWHM_left, self.numerical_FWHM_right = self.find_peak_and_fwhm(mu, sigma)

		return [self.this_roohist, self.roodefs, self.this_pdf, self.my_name, self.spl_draw, self.rooCompon, self.fit_results]	
	
class hyperEmg(FitMethods):
	"""
	Class handling the exponentially modified Gaussian as defined in hyper-EMG paper from JLU-GSI group
	https://www.sciencedirect.com/science/article/abs/pii/S1387380616302913
	https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
	"""
	def __init__(self, lst_file, file_path=False, peaks = 'findit', verbose=False):
		"""
		Initializes the class and number of parameters.
		"""    
		#    
		FitMethods.__init__(self, lst_file, file_path=file_path, name='hyperEmg', peaks=peaks, verbose=verbose)
		#
		self.fit_func_name = "hyperEmg"
		# JLU-GSI group definition
		# ground state 
		self.pos_funct = '1/(2*@3)*exp((@2/(TMath::Sqrt(2)*@3))^2-(@0-@1)/@3)*TMath::Erfc(@2/(TMath::Sqrt(2)*@3)-(@0-@1)/(TMath::Sqrt(2)*@2))'
		self.neg_funct = '1/(2*@3)*exp((@2/(TMath::Sqrt(2)*@3))^2+(@0-@1)/@3)*TMath::Erfc(@2/(TMath::Sqrt(2)*@3)+(@0-@1)/(TMath::Sqrt(2)*@2))'
		# isomeric states
		self.pos_funct_on_gs = '1/(2*@3)*exp((@2/(TMath::Sqrt(2)*@3))^2-(@0-(@1+@4))/@3)*TMath::Erfc(@2/(TMath::Sqrt(2)*@3)-(@0-(@1+@4))/(TMath::Sqrt(2)*@2))'
		self.neg_funct_on_gs = '1/(2*@3)*exp((@2/(TMath::Sqrt(2)*@3))^2+(@0-(@1+@4))/@3)*TMath::Erfc(@2/(TMath::Sqrt(2)*@3)+(@0-(@1+@4))/(TMath::Sqrt(2)*@2))'
		self.params = {}
		self.simultaneous = False
		# Wikipedia definition
		# self.pos_funct = '1/(2*@3)*exp(1/(2*@3)*(2*@1+(1/@3)*@2^2-2*@0))*TMath::Erfc((@1+(1/@3)*@2^2-@0)/(TMath::Sqrt(2)*@2))'
		# self.neg_funct = '1/(2*@3)*exp(1/(2*@3)*(2*@1+(1/@3)*@2^2-2*@0))*TMath::Erfc((@1+(1/@3)*@2^2-@0)/(TMath::Sqrt(2)*@2))'

	def pos_func(self, x, mu, sigma, ptau):
		"""
		Analytical function for positive Emg components
		Parameters:
			-x, mu, sigma, ntau: components of Emg function
		"""
		# JLU-GSI definition
		return 1/(2*ptau)*np.exp((sigma/(math.sqrt(2)*ptau))**2-(x-mu)/ptau)*sc.special.erfc(sigma/(math.sqrt(2)*ptau)-(x-mu)/(math.sqrt(2)*sigma))
		# Wikipedia definition
		# return 1/(2*ptau)*np.exp((1/(2*ptau))(2*mu+(sigma**8)/(ptau)-2*x))*sc.special.erfc((mu+(1/ptau)*sigma**2-x)/(1.4142*sigma))
	
	def pos_func_on_gs(self, x, mu, sigma, ptau, exc_energy):
		"""
		Analytical function for positive Emg components
		Parameters:
			-x, mu, sigma, ntau: components of Emg function
			-exc_energy: routine returns function with location of isomer build on ground state: mu = mu0 + E
		"""
		# JLU-GSI definition
		return 1/(2*ptau)*np.exp((sigma/(math.sqrt(2)*ptau))**2-(x-(mu+exc_energy))/ptau)*sc.special.erfc(sigma/(math.sqrt(2)*ptau)-(x-(mu+exc_energy))/(math.sqrt(2)*sigma))
		# Wikipedia definition
		# return 1/(2*ptau)*np.exp((1/(2*ptau))(2*mu+(sigma**8)/(ptau)-2*x))*sc.special.erfc((mu+(1/ptau)*sigma**2-x)/(1.4142*sigma))

	def neg_func(self, x, mu, sigma, ntau):
		"""
		Analytical function for negative Emg components
		Parameters:
			-x, mu, sigma, ntau: components of Emg function
		"""
		# JLU-GSI definition
		return 1/(2*ntau)*np.exp((sigma/(math.sqrt(2)*ntau))**2+(x-mu)/ntau)*sc.special.erfc(sigma/(math.sqrt(2)*ntau)+(x-mu)/(math.sqrt(2)*sigma))
		# Wikipedia definition
		# return 1/(2*ntau)*np.exp((1/(2*ntau))(2*mu+(sigma**8)/(ntau)-2*x))*sc.special.erfc((mu+(1/ntau)*sigma**2-x)/(1.4142*sigma))

	def neg_func_on_gs(self, x, mu, sigma, ntau, exc_energy):
		"""
		Analytical function for negative Emg components
		Parameters:
			-x, mu, sigma, ntau: components of Emg function
			-exc_energy: routine returns function with location of isomer build on ground state: mu = mu0 + E
		"""
		# JLU-GSI definition
		return 1/(2*ntau)*np.exp((sigma/(math.sqrt(2)*ntau))**2+(x-(mu+exc_energy))/ntau)*sc.special.erfc(sigma/(math.sqrt(2)*ntau)+(x-(mu+exc_energy))/(math.sqrt(2)*sigma))
		# Wikipedia definition
		# return 1/(2*ntau)*np.exp((1/(2*ntau))(2*mu+(sigma**8)/(ntau)-2*x))*sc.special.erfc((mu+(1/ntau)*sigma**2-x)/(1.4142*sigma))

	def function_string(self, c0='@0', c1='@1', c2='@2', c3='@3', R='@4',pol="+"):
		"""
		Manual definition of a RooFit like arbitrary string for the hyperEmg fit
		"""
		#
		return f'{R}*1/(2*{c3})*exp(({c2}/(1.4142*{c3}))^2{pol}({c0}-{c1})/{c3})*TMath::Erfc({c2}/(1.4142*{c3}){pol}({c0}-{c1})/(1.4142*{c2}))'

	def build_function_string(self, dimensions = [0,1], params = {}, state = {}, position_mu = '1'):
		#
		funct = ''
		comp = 0
		#
		sigma = '@2'

		# Negative Functions
		for i in range(dimensions[0]):
			
			# Define contribution
			if (comp == dimensions[0]+dimensions[1]-1 and dimensions[0]+dimensions[1] > 1):
				contrib = '(1-'
				#
				for j in range(dimensions[0]+dimensions[1]-1):
					contrib+=f"@{4+2*j}-"
				#
				contrib = contrib[:-1] + ")"
			elif (dimensions[0]+dimensions[1]==1): 
				contrib = f"1"
			else:
				contrib = f"@{4+2*comp}"    
			
			# Define negative tau
			ntau = f"@{3+2*comp}"

			# Define mu and/or excitation energy E
			if state['state'] == 'gs':
				mu = '@'+position_mu

			# If to be fitted as excitation energy
			elif state["state"] != 'gs':
				position_of_exc_parameter = int(2*(dimensions[0]+dimensions[1])-1 + 3)
				mu = f"(@{position_mu}+@{position_of_exc_parameter})"

			else:
				mu = '@'+position_mu
			# Build function
			funct += self.function_string(c1 = mu, c2 = sigma, c3=ntau, R=f"{contrib}", pol="+")+ "+"
			comp+=1
		#
		for i in range(dimensions[1]):
			if (comp == dimensions[0]+dimensions[1]-1 and dimensions[0]+dimensions[1] > 1):
				contrib = '(1-'
				#
				for j in range(dimensions[0]+dimensions[1]-1):
					contrib+=f"@{4+2*j}-"
				#
				contrib = contrib[:-1] + ")"
			elif (dimensions[0]+dimensions[1]==1): 
				contrib = f"1"
			else:
				contrib = f"@{4+2*comp}"    
			
			# Define positive tau
			ptau = f"@{3+2*comp}"
			
			# Define mu and/or excitation energy E
			if state['state'] == 'gs':
				mu = '@'+position_mu

			# If to be fitted as excitation energy
			elif state["state"] != 'gs':
				position_of_exc_parameter = int(2*(dimensions[0]+dimensions[1])-1 + 3)
				mu = f"(@{position_mu}+@{position_of_exc_parameter})"
			
			else:
				mu = f'@{position_mu}'
			# Builf function
			funct += self.function_string(c1 = mu, c2 = sigma, c3=ptau, R=f"{contrib}", pol="-")+ "+"
			comp+=1
		# Add constant background
		# funct = funct[:-1] + f"+@{4+2*comp+1}" # cut away last character which is a "+" and add const. background
		#
		return funct[:-1] # cut away last +
	
	def build_function_string_SAVE(self, dimensions = [0,1], params = {}, state = {}):
		#
		funct = ''
		comp = 0
		#
		if f"sigma" not in params:
				sigma = '@2'
		else:
			sigma = params['sigma']
		# Negative Functions
		for i in range(dimensions[0]):
			
			# Define contribution
			if (comp == dimensions[0]+dimensions[1]-1 and dimensions[0]+dimensions[1] > 1):
				contrib = '(1-'
				#
				if f"contrib{comp-1}" not in params:
					for j in range(dimensions[0]+dimensions[1]-1):
						contrib+=f"@{4+2*j}-"
				else: 
					for j in range(dimensions[0]+dimensions[1]-1):
						temp = f"contrib{j}"
						contrib += f"{params[temp]}-"
				#
				contrib = contrib[:-1] + ")"
			elif (dimensions[0]+dimensions[1]==1): 
				if f"contrib{comp}" not in params:
					contrib = f"1"
				else: 
					temp = f"contrib{comp}"
					contrib = f"{params[temp]}"
			else:
				if f"contrib{comp}" not in params:
					contrib = f"@{4+2*comp}"    
				else: 
					temp = f"contrib{comp}"
					contrib = f"{params[temp]}"
			
			# Define negative tau
			if f"ntau{i}" not in params:
				ntau = f"@{3+2*comp}"
			else:
				temp = f"ntau{i}"
				ntau = f"{params[temp]}"

			# Define mu and/or excitation energy E
			if state['state'] == 'gs' and state['peak'] not in params:
				mu = '@1'
			elif state['state'] == 'gs' and state['peak'] in params:
				temp = state['peak']
				mu = f"{params[temp]}"
			# If to be fitted as excitation energy
			elif state["state"] != 'gs':
				position_of_exc_parameter = int(2*(dimensions[0]+dimensions[1])-1 + 3)
				mu = f"(@1+@{position_of_exc_parameter})"
			# If excitation energies are passed
			elif {state["peak"]}+f'-E0' and state["state"] != 'gs':
				temp = {state["peak"]}+'-E0'
				mu = f"(@1+{params[temp]})"
			elif {state["peak"]}+f'-E1' and state["state"] != 'gs':
				temp = {state["peak"]}+'-E1'
				mu = f"(@1+{params[temp]})"
			else:
				mu = '@1'
			# Builf function
			funct += self.function_string(c1 = mu, c2 = sigma, c3=ntau, R=f"{contrib}", pol="+")+ "+"
			comp+=1
		#
		for i in range(dimensions[1]):
			if (comp == dimensions[0]+dimensions[1]-1 and dimensions[0]+dimensions[1] > 1):
				contrib = '(1-'
				#
				if f"contrib{comp-1}" not in params:
					for j in range(dimensions[0]+dimensions[1]-1):
						contrib+=f"@{4+2*j}-"
				else: 
					for j in range(dimensions[0]+dimensions[1]-1):
						temp = f"contrib{j}"
						contrib += f"{params[temp]}-"
				#
				contrib = contrib[:-1] + ")"
			elif (dimensions[0]+dimensions[1]==1): 
				if f"contrib{comp}" not in params:
					contrib = f"1"
				else: 
					temp = f"contrib{comp}"
					contrib = f"{params[temp]}"
			else:
				if f"contrib{comp}" not in params:
					contrib = f"@{4+2*comp}"    
				else: 
					temp = f"contrib{comp}"
					contrib = f"{params[temp]}"
			
			# Define positive tau
			if f"ptau{i}" not in params:
				ptau = f"@{3+2*comp}"
			else:
				temp = f"ptau{i}"
				ptau = f"{params[temp]}"
			
			# Define mu and/or excitation energy E
			if state['state'] == 'gs' and state['peak'] not in params:
				mu = '@1'
			elif state['state'] == 'gs' and state['peak'] in params:
				temp = state['peak']
				mu = f"{params[temp]}"
			# If to be fitted as excitation energy
			elif state["state"] != 'gs':
				position_of_exc_parameter = int(2*(dimensions[0]+dimensions[1])-1 + 3)
				mu = f"(@1+@{position_of_exc_parameter})"
			# If excitation energies are passed
			elif {state["peak"]}+f'-E0' and state["state"] != 'gs':
				temp = {state["peak"]}+'-E0'
				mu = f"(@1+{params[temp]})"
			elif {state["peak"]}+f'-E1' and state["state"] != 'gs':
				temp = {state["peak"]}+'-E1'
				mu = f"(@1+{params[temp]})"
			else:
				mu = '@1'
			# Builf function
			funct += self.function_string(c1 = mu, c2 = sigma, c3=ptau, R=f"{contrib}", pol="-")+ "+"
			comp+=1
		return funct[:-1] # cut away last character which is a "+"

	def __add_EMG_component(self, j, i, pol='-'):
		# component name
		if pol == "+": # positive sign is negative component
			pdf_name = f"state_{j}-{i}n"
		else:
			pdf_name = f"state_{j}-{i}p"
		print(f'--> Component: {pdf_name}')
		
		# build function
		funct = self.function_string(c0='@0', c1='@1', c2='@2', c3='@3', R='1',pol=pol)
		print(f'--> Function String: {funct}')
		
		# Build RooArgList iteratively depending on which parameters are fixed
		listofRooArgs = RooArgList()
		listofRooArgs.add(self.RooRealVar_dict['x'])
		print("--> Add RooArg: 'x'")
		# Add gaussian center
		mu = self.states[j]["peak"]
		listofRooArgs.add(self.RooRealVar_dict[mu])
		print(f"--> Add RooArg: '{mu}'")
		#
		listofRooArgs.add(self.RooRealVar_dict['sigma'])
		print(f"--> Add RooArg: 'sigma'")
		if pol == '+': # negative component has +  
			tau = f"ntau{i}"
			listofRooArgs.add(self.RooRealVar_dict[f"ntau{i}"])
		else: # positive component has - 
			tau = f"ptau{i}"
			listofRooArgs.add(self.RooRealVar_dict[f"ptau{i}"])
		print(f"--> Add RooArg: '{tau}'")

		state = 'mu0'
		if self.states[j]["state"] != 'gs':
			E = self.states[j]["peak"]+"-"+self.states[j]["state"]
			listofRooArgs.add(self.RooRealVar_dict[E])
			print(f"--> Add RooArg: '{E}'")
		#
			state = E

		self.EMGs_dict[pdf_name] = RooGenericPdf(pdf_name,pdf_name,funct,listofRooArgs)
		# return RooGenericPdf(pdf_name,pdf_name,funct,listofRooArgs)

	def __build_params(self, limits):
		"""
		Takes limits dictionary and puts all fixed fit variables into another dict
		"""
		for limit in limits:
			if limits[limit][0] == limits[limit][1] and limits[limit][0] == limits[limit][2]:
				self.params[limit] = limits[limit][0]

	def __infer_states(self, n_comps, fromm='limits'):
		"""
		Infer states from number of peaks passed
		"""
		if fromm == 'limits':
			dic = self.limits 
		elif fromm == 'fit-file':
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
		if len(dic) != n_comps:
			print("(fit.hyperEmg.__infer_states): Warning: number of infered states does not correspond number of peaks to be fitted!")

	def init_limits(self):
		self.limits = {
			'mean0':[self.hist.GetMean(1), self.xmin, self.xmax],
			'mean1':[self.hist.GetMean(1), self.xmin, self.xmax],
			'sigma':[self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1)],
			'ratio':[0.5,0,1],
		}

	def constraints_from_file(self, file, params_to_fix, for_template_fit=True):
		'''
		Loads and sets constraints for fit from fit file
		Parameters:
			- file: fit file to be loaded
			- params_to_fix: array of strings of parameters to fix
			- for_template_fit: also stores parameters and dict for template fitting
		'''
		#
		self.fit = FitToDict(file)
		#
		if not 'FIT-VALUES' in self.fit.fit:
			print(f"Fit file {file} has no fit values that were exported.")
			return 0
		#
		for idx,row in self.fit.fit['RESULTS-TABLE'].iterrows():
			if row['var'] in params_to_fix:
				print(f"fixed {row['var']} from file")
				self.limits[row['var']] = [row['value'], row['value'], row['value']]
				self.params[row['var']] = row['value']

	def build_RooGenericPdf_dict(self):

		for d in range(self.dimensions[0]+self.dimensions[1]-1): 
			var_name = f"contrib{d}"
			print(f"--> Add RooArg: 'contrib{d}'")
			self.RooRealVar_dict[var_name] = RooRealVar(var_name, var_name, self.limits[var_name][0], self.limits[var_name][1], self.limits[var_name][2])

		# Build one PDF per component, consisting of contributions of different exponetially modified gaussian distributions
		for j in range(len(self.states)):
			pdf_name = f"comp{j}"
			# Start by building the PDF function string
			funct = self.build_function_string(dimensions=self.dimensions, params=self.params, state = self.states[j])
			print(f'--> Function String: {funct}')
			# Build RooArgList iteratively depending on which parameters are fixed
			listofRooArgs = RooArgList()
			listofRooArgs.add(self.RooRealVar_dict['x'])
			print("--> Add RooArg: 'x'")
			# Add gaussian center
			mu = self.states[j]["peak"]
			listofRooArgs.add(self.RooRealVar_dict[mu])
			print(f"--> Add RooArg: '{mu}'")
			#
			listofRooArgs.add(self.RooRealVar_dict['sigma'])
			print(f"--> Add RooArg: 'sigma'")
			# 
			contrib_idx = 0
			for i in range(self.dimensions[0]):
				tau = f"ntau{i}"
				listofRooArgs.add(self.RooRealVar_dict[f"ntau{i}"])
				print(f"--> Add RooArg: '{tau}'")
				if self.dimensions[0]+self.dimensions[1] != 1 and contrib_idx < self.dimensions[0]+self.dimensions[1] - 1:
					listofRooArgs.add(self.RooRealVar_dict[f"contrib{contrib_idx}"])
					print(f"--> Add RooArg: 'contrib{contrib_idx}'")
					contrib_idx += 1
			for i in range(self.dimensions[1]):
				tau = f"ptau{i}"
				listofRooArgs.add(self.RooRealVar_dict[f"ptau{i}"])
				print(f"--> Add RooArg: '{tau}'")
				if self.dimensions[0]+self.dimensions[1] != 1 and contrib_idx < self.dimensions[0]+self.dimensions[1] - 1:
					listofRooArgs.add(self.RooRealVar_dict[f"contrib{contrib_idx}"])
					print(f"--> Add RooArg: 'contrib{contrib_idx}'")
					contrib_idx += 1
			# Reset contrib_idx counter of contribution ratios are to be kept the same for multiple peaks
			if self.simultaneous:
				contrib_idx = 0

			state = 'mu0'
			if self.states[j]["state"] != 'gs':
				E = self.states[j]["peak"]+"-"+self.states[j]["state"]
				listofRooArgs.add(self.RooRealVar_dict[E])
				print(f"--> Add RooArg: '{E}'")
			#
				state = E

			self.RooGenericPdf_dict[pdf_name] = RooGenericPdf(pdf_name,pdf_name,funct,listofRooArgs)
			
		# Build component ratios
			if j < (self.n_comps-1): 
				var_name = f"ratio{j}"
				print(f"--> Add RooArg: 'ratio{j}'")
				self.RooRealVar_dict[var_name] = RooRealVar(var_name, var_name, self.limits[var_name][0], self.limits[var_name][1], self.limits[var_name][2])
			#
			j += 1

	def build_final_hyperEMGpdfs(self):
		# Build one PDF per component, consisting of contributions of different exponetially modified gaussian distributions
		# WARNING:
		# The RooArgList have usually to be passed as attributes of the class, I encountered seg faults when I tried to pass them
		# through the return value of the functions....
		
		for d in range(self.dimensions[0]+self.dimensions[1]-1): 
			var_name = f"contrib{d}"
			self.RooRealVar_dict[var_name] = RooRealVar(var_name, var_name, self.limits[var_name][0], self.limits[var_name][1], self.limits[var_name][2])

		for j in range(len(self.states)):

			EMGs = RooArgList()
			EMG_ratios = RooArgList()
			self.EMGs_dict = {}

			contrib_idx = 0
			for i in range(self.dimensions[0]):
				self.__add_EMG_component(j, i, "+")# negative tail has + sign

			for i in range(self.dimensions[1]):
				self.__add_EMG_component(j, i, "-") # positive tail has - sign

			for d in range(self.dimensions[0]+self.dimensions[1]-1):
				key = f'contrib{d}'
				print(f"--> Add RooArg: 'contrib{d}'")
				EMG_ratios.add(self.RooRealVar_dict[key])

			for emgs in self.EMGs_dict:
				EMGs.add(self.EMGs_dict[emgs])

			self.RooGenericPdf_dict[j] = RooAddPdf(f'hyperEmg-{j}', f'hyperEmg-{j}', EMGs, EMG_ratios, recursiveFraction = ROOT.kFALSE)
		
		#
		# j += 1

	def build_master_string(self, dimensions = [0,1], params = {}, state = {}):
		""" 
		Method that builds a string from n components and all dimensions simultaneously 
		"""
		# add roo vars to dict



		#
		if self.n_comps == 1:
			return self.build_function_string(dimensions=self.dimensions, params=self.params, state = self.states[0], position_mu = '1')
		#
		n_emg = 3 + 2 * (dimensions[0]+dimensions[1]-1)
		n_vars = 3 + 2 * (dimensions[0]+dimensions[1]-1) + 2 * (self.n_comps-1) # 3 variables for base EMG plus two per every extra EMG plus two for every extra hyperEMG
		#

		string = f""

		print(self.states)

		for i in np.arange(1,len(self.states)+1,1):

			position_mu = '1' if i==1 else str(n_emg+2*(i-1)-1)

			if self.states[i-1]['state'] != 'gs': position_mu = '1'  

			hyperEMG_string = self.build_function_string(dimensions=self.dimensions, params=self.params, state = self.states[i-1], position_mu = position_mu)

			if i < self.n_comps:
				contrib = f"@{n_emg+2*i}"
			if i == self.n_comps:
				contrib = '(1-'
				#
				for j in np.arange(1,len(self.states),1):
					contrib+=f"@{n_emg+2*j}-"
				#
				contrib = contrib[:-1] + ")"
				
			string += f"{contrib}*({hyperEMG_string})+"

		print(string[:-1])

		return string[:-1]

	def call_pdf(self, xmin, xmax, dimensions = [1,2], n_comps = 1, simultaneous = False, 
				bins = 1, limits=False, minos = ROOT.kFALSE):
		"""
		Setup Variable, Parameters and PDF and performs fitting of the PDF to ROI data.
		:param roi_hist: histogram to fit
		:param xmin: range min for the fit
		:param xmax: range max for the fit
		:dimensions: array with number of negative and positive components to the hyper-EMG
		:n_comps: number of species withing fit range
		:bins: binning for the histogram to be fitted
		:return list: containing all relevant information about the ROI histogram, PDF, parameters,
					  PDF components and fit results
		"""
		# Test dimension of limits
		self.xmin = xmin
		self.xmax = xmax

		# Convert histogram
		self.bins = bins
		#
		self.binning = self.get_binning(bins=self.bins)
		self.n_comps = n_comps
		self.simultaneous = simultaneous
		self.hist = self.lst2roothist(self.lst_file, bins=1)

		# Initialize limits
		if not limits:
			self.init_limits()
		#  
		else:
			# If limits are passed, save them in dict
			# Careful! This overwrite the constraints set by constraints_from_file
			for limit in limits:
				self.limits[limit] = limits[limit]

		# From limits infer the ground state / isomeric state character of peaks passed
		self.__infer_states(n_comps, fromm='limits')

		# Fill variable dictionary:
		if isinstance(xmin, list):
			self.RooRealVar_dict = {
				'x': RooRealVar("x", self.rootitle, self.xmin[0], self.xmax[-1], self.roounit)
			}
		else:
			self.RooRealVar_dict = {
				'x': RooRealVar("x", self.rootitle, self.xmin, self.xmax, self.roounit)
			}
		# Constant background
		# self.RooRealVar_dict["const_bck"] = RooRealVar("const_bck", "const_bck", self.limits["const_bck"][0], self.limits["const_bck"][1], self.limits["const_bck"][2])
		# self.RooGenericPdf_dict["const_bck"] = RooPolynomial("const_bck","const_bck", self.RooRealVar_dict['x'], RooArgList())
		
		# Define Gaussian components for the fit
		# Shared by all components
		self.RooRealVar_dict["sigma"] = RooRealVar("sigma", "sigma", self.limits["sigma"][0], self.limits["sigma"][1], self.limits["sigma"][2])
		# for ground states and isomers or other peaks
		for j in range(len(self.states)):
			if self.states[j]["state"] == 'gs':
				var_name = self.states[j]["peak"]
			else:
				var_name = self.states[j]["peak"] + "-" + self.states[j]["state"]
			if var_name not in self.RooRealVar_dict:
				self.RooRealVar_dict[var_name] = RooRealVar(var_name, var_name, self.limits[var_name][0], self.limits[var_name][1], self.limits[var_name][2])
		
		# Dimensions for Emg
		self.dimensions = dimensions
		# Definition of negative exponential components
		for i in range(0,dimensions[0],1):
			var_name = f"ntau{i}"
			self.RooRealVar_dict[var_name] = RooRealVar(var_name, var_name, self.limits[var_name][0], self.limits[var_name][1], self.limits[var_name][2])
		# Definition of positive exponential components
		for i in range(0,dimensions[1],1):
			var_name = f"ptau{i}"
			self.RooRealVar_dict[var_name] = RooRealVar(var_name, var_name, self.limits[var_name][0], self.limits[var_name][1], self.limits[var_name][2])
		# Definition of contrib factors: 
		# if not simultaneous: n-1 factors for n components. (+1 contrib for constant background)
		# else: only k-1 variables if k is the sum of both positive and negative dimensions
		if simultaneous:
			n_components = (dimensions[0] + dimensions[1]) - 1
		else:
			n_components = n_comps * (dimensions[0] + dimensions[1]) - 1
		for i in range(0,n_components,1):
			var_name = f"contrib{i}"
			# self.RooRealVar_dict[var_name] = RooRealVar(var_name, var_name, self.limits[var_name][0], self.limits[var_name][1], self.limits[var_name][2])
		
		# Check for fixed parameters from limits passed
		self.__build_params(self.limits)

		# Build all component pdfs where the contribution mixing is done within the generic PDF, not in the RooAddPdf
		# self.build_RooGenericPdf_dict()

		# Alternative: build all components where the contributions are mixed on the RooAddPdf level
		# self.build_final_hyperEMGpdfs()

		# Alternative 2: build very long function string and not use any RooAddPdf

		funct = self.build_master_string(dimensions=dimensions)

		# # #Put all pdfs together
		all_pdfs = RooArgList()
		for pdf in self.RooGenericPdf_dict:
			all_pdfs.add(self.RooGenericPdf_dict[pdf])

		# put all ratio's together
		all_ratios = RooArgList()
		for var in self.RooRealVar_dict:
			m = re.search('ratio*', var)
			if m:
				all_ratios.add(self.RooRealVar_dict[var])

		# Definition hyper-EMG

		# self.this_pdf = RooAddPdf('hyperEmg', 'hyperEmg', all_pdfs, all_ratios, recursiveFraction = ROOT.kFALSE)
		


		for d in range(self.dimensions[0]+self.dimensions[1]-1): 
			var_name = f"contrib{d}"
			self.RooRealVar_dict[var_name] = RooRealVar(var_name, var_name, self.limits[var_name][0], self.limits[var_name][1], self.limits[var_name][2])

		for j in range(self.n_comps-1): 
			var_name = f"ratio{j}"
			self.RooRealVar_dict[var_name] = RooRealVar(var_name, var_name, self.limits[var_name][0], self.limits[var_name][1], self.limits[var_name][2])


		# Build RooArgList iteratively depending on which parameters are fixed
		listofRooArgs = RooArgList()
		listofRooArgs.add(self.RooRealVar_dict['x'])
		print("--> Add RooArg: 'x'")
		# Add gaussian center
		listofRooArgs.add(self.RooRealVar_dict['mu0'])
		print(f"--> Add RooArg: 'mu0'")
		#
		listofRooArgs.add(self.RooRealVar_dict['sigma'])
		print(f"--> Add RooArg: 'sigma'")
		# 
		contrib_idx = 0
		for i in range(self.dimensions[0]):
			tau = f"ntau{i}"
			listofRooArgs.add(self.RooRealVar_dict[f"ntau{i}"])
			print(f"--> Add RooArg: '{tau}'")
			if self.dimensions[0]+self.dimensions[1] != 1 and contrib_idx < self.dimensions[0]+self.dimensions[1] - 1:
				listofRooArgs.add(self.RooRealVar_dict[f"contrib{contrib_idx}"])
				print(f"--> Add RooArg: 'contrib{contrib_idx}'")
				contrib_idx += 1
		for i in range(self.dimensions[1]):
			tau = f"ptau{i}"
			listofRooArgs.add(self.RooRealVar_dict[f"ptau{i}"])
			print(f"--> Add RooArg: '{tau}'")
			if self.dimensions[0]+self.dimensions[1] != 1 and contrib_idx < self.dimensions[0]+self.dimensions[1] - 1:
				listofRooArgs.add(self.RooRealVar_dict[f"contrib{contrib_idx}"])
				print(f"--> Add RooArg: 'contrib{contrib_idx}'")
				contrib_idx += 1

		i = 0
		for state in self.states:

			if i == 0: 
				i+= 1
				continue

			if state['state'] == 'gs':
				mu = self.states[i]["peak"]
				listofRooArgs.add(self.RooRealVar_dict[mu])
				print(f"--> Add RooArg: '{mu}'")

			if state['state'] != 'gs':
				E = self.states[i]["peak"]+"-"+self.states[i]["state"]
				listofRooArgs.add(self.RooRealVar_dict[E])
				print(f"--> Add RooArg: '{E}'")
				# state = E


			ratiostring = f'ratio{i-1}'
			listofRooArgs.add(self.RooRealVar_dict[ratiostring])
			print(f"--> Add RooArg: '{ratiostring}'")
			
			i+= 1

		self.this_pdf = RooGenericPdf('hyperEmg', 'hyperEmg', funct, listofRooArgs)

		#
		self.roodefs = [self.this_pdf, self.RooRealVar_dict, self.RooGenericPdf_dict]
		if minos: 
			self.this_roohist, self.fit_results = self.minimize(self.hist, self.xmin, self.xmax, datatype='is_th1d')
		else:
			self.this_roohist, self.fit_results = self.minimize(self.hist, self.xmin, self.xmax, datatype='is_th1d', minos = ROOT.kFALSE)

		# THIS LINE HAS TO BE HERE FOR SOME UNKNOW REASON TO ME TO AVOID A SEG FAULT IN THE PLOT FUNCTION 
		print(self.this_pdf.getVal(ROOT.RooArgSet(self.RooRealVar_dict['x'])))

		# Calculate numerical position of maximum and FWHM
		mu = self.RooRealVar_dict['mu0'].getValV()
		sigma = self.RooRealVar_dict['sigma'].getValV()
		self.numerical_peak, self.numerical_FWHM, self.numerical_FWHM_left, self.numerical_FWHM_right = self.find_peak_and_fwhm(mu, sigma)
		# Store fit results in dict
		for key in self.RooRealVar_dict:
			self.Var_dict[key] = self.RooRealVar_dict[f'{key}'].getValV()

		return [self.this_roohist, self.roodefs, self.this_pdf, self.my_name, self.spl_draw, self.rooCompon, self.fit_results]	