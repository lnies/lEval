import ROOT
from ROOT import RooFit
from ROOT import RooRealVar, RooArgSet, RooArgList, RooDataHist
from ROOT import RooGenericPdf, RooUniform, RooPolynomial, RooGaussian, RooGaussModel, RooDecay, RooFormulaVar
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
import re
import math 
import time
import mmap

from utilities import FitToDict, Peaks, softCool

class FitMethods():
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
		# QtCore.QObject.__init__(self)
		self.my_name = name
		self.this_pdf = RooGaussian()
		self.roi_mean = 0.0
		self.roi_rms = 0.0
		self.roi_counts = 0.0
		self.roi_skewness = 0.0
		#
		self.RooRealVar_dict = {}
		self.Var_dict = {}  
		#
		self.lst_file = lst_file
		self.fitrange = 100
		self.results = {}
		self.fits_r = []
		self.limits = {}
		self.meta_data = {}
		self.meta_data_keys = ["range", "swpreset", "cycles", "cmline0", "caloff", "calfact", "time_patch", ]
		#
		# If no peaks are passed, find peaks within the initialization of the function
		if peaks == 'findit':
			self.peaks = Peaks(self.lst_file)
			self.peaks.find_peaks()
			if verbose:
				print("Peaks initialized by class")
				for i in range(self.peaks.n_peaks):
					print(self.peaks.peaks_info)
		else:
			self.peaks = peaks
		# Read meta data if lst file path is provided
		# if file_path != False:
		# 	with open(file_path + ".lst",'rb') as listfile:
		# 	    mapped_file = mmap.mmap(listfile.fileno(), 0, access=mmap.ACCESS_READ)
		# 	    for key in self.meta_data_keys:
		# 		    mapped_file.seek(mapped_file.find(f'{key}'.encode('ascii')))
		# 		    self.meta_data[key] = mapped_file.readline().strip('\r\n'.encode('ascii')).decode('ascii')

	def compose_title_unit(self, xtitle):
		"""
		Compose label and unit of RooRealVar by using the X or Y Axis Label of the ROI histogram.
		:param xtitle:
		:return:
		"""
		if '/' in xtitle:
			logging.info('xtitle = {}'.format(xtitle.split('/')[0]))
			self.rootitle = xtitle.split('/')[0]
			self.roounit = xtitle.split('/')[-1]
		else:
			self.rootitle = xtitle
			self.roounit = 'a.u.'
		
	def minimize(self, th1d, xmin, xmax):
		"""
		Performs the maximum likelihood fit.

		:param th1d: the histogram to fit
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
		else:
			ranges = f"win_{self.roodefs[0].GetName()}"
			x.setRange(ranges, xmin, xmax)
		roohist = RooDataHist('roohist_tof', 'title', RooArgList(x), RooFit.Import(th1d))
		print("Test-MINIMIZE")
		result_mlkh = self.roodefs[0].fitTo(roohist,
											RooFit.Range(ranges),
											# RooFit.NormRange(f"win_{self.roodefs[0].GetName()}"),
											RooFit.Minos(ROOT.kTRUE),
											RooFit.PrintEvalErrors(-1),
											# RooFit.SumW2Error(ROOT.kTRUE),
											RooFit.NumCPU(1),
											RooFit.Timer(ROOT.kTRUE),
											RooFit.Save(),
											RooFit.Verbose(ROOT.kFALSE))
		# self.pars_show()
		self.rerun = True
		return roohist, result_mlkh
	
	def get_binning(self, bins=1):
		"""
		Adapts binning to multiples of 0.8ns, assuming that 0.8ns was used to take data (to common case)
		"""
		# Get min and max tof from data frame
		self.bins = bins
		minn = self.lst_file.tof.min()
		maxx = self.lst_file.tof.max()
		return round((maxx-minn)/0.8/self.bins)

	def lst2roothist(self, list_file, bins=1):
		"""
		Converts .lst file to root histogram.

		:return: root histogram
		"""
		# Get min and max tof from data frame
		self.bins = bins
		minn = self.lst_file.tof.min()
		maxx = self.lst_file.tof.max()
		hist = ROOT.TH1D( 'hist', 'hist converted', self.get_binning(bins), minn, maxx)
		[hist.Fill(x) for x in list_file.tof]
		return hist

# def read_fit(path_to_file):
# 	meta_data_keys = ["range", "swpreset", "cycles", "cmline0", "caloff", "calfact", "time_patch",
# 					  "fit_function", "binning", "xmin", "xmax", "numerical_peak", "numerical_FWHM", 
# 					  "numerical_FWHM_left", "numerical_FWHM_right"]

# 	file_number = re.split('In_run_|_',path_to_file)[3]
# 	# Extract fit parameters and save them in dict 
# 	sub_d = {
# 		'file_number': file_number,
# 		'species': re.split('In_run_|revs|_',path_to_file)[4],
# 		'n_revs': re.split('In_run_|revs|_',path_to_file)[5]
# 	}
# 	#

# 	with open(path_to_file,'rb') as fit_file:
# 		mapped_file = mmap.mmap(fit_file.fileno(), 0, access=mmap.ACCESS_READ)
# 		for key in meta_data_keys:
# 			try:
# 				mapped_file.seek(mapped_file.find(f'{key}'.encode('ascii')))
# 				meta_data = mapped_file.readline().strip('\r\n'.encode('ascii')).decode('ascii')
# 				sub_d[key] = re.split('=',meta_data)[1]
# 			except: 
# 				continue
# 		# Find start of fit data 
# 		for num, line in enumerate(fit_file, 1):
# 			if '[RESULTS-TABLE]'.encode('ascii') in line:
# 				header_row =  num
# 		# read fit 
# 		df = pd.read_csv(path_to_file, delimiter=" ", header=header_row)
# 		sub_d['fit'] = df

# 	return sub_d

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
		self.hist = self.lst2roothist(self.lst_file, bins=2000)
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
	
	def plot(self, log=False):
		"""  """
		if len(self.lst_file) == 0:
			print("Fit not excecuted yet or failed.")
			return 0 
		#            
		xdata = self.lst_file.tof
		n, xe = np.histogram(xdata, bins=bins)
		cx = 0.5 * (xe[1:] + xe[:-1])
		dx = np.diff(xe)
		plt.errorbar(cx, n, n ** 0.5, fmt="ok", zorder=1)
		xm = np.linspace(xe[0], xe[-1], num=1000)
		# Plot 'carpet'
		if log:
			plt.plot(xdata, np.zeros_like(xdata)+0.9, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		plt.plot(xdata, np.zeros_like(xdata)-5, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		#
		# Plot fits, loop through all fits
		if log:
			xm = np.linspace(self.mean.getValV()-100, self.mean.getValV()+100, num=1000)
		#
		y_val = [RootMath.gaussian_pdf(x=y, x0=self.mean.getValV(), sigma=self.sigma.getValV()) * len(xdata) * dx[0] for y in xm]
		plt.plot(xm, y_val, ls=":", label=f"init", zorder=2)
		plt.plot(xm, y_val, label=f"fit", c='r', zorder=3)
		
		plt.legend();
		if log:
			plt.yscale("log")
		# Zoom in on found peaks
		if self.peaks.n_peaks != 0:
			plt.xlim(self.peaks.earliest_left_base-200, self.peaks.latest_right_base+200)

class doubleGauss(FitMethods):
	"""
	Class handling the Gaussian model.
	"""
	def __init__(self, lst_file, peaks = 'findit', verbose=False):
		"""
		Initializes the class and number of parameters.
		"""        
		# QtWidgets.QWidget.__init__(self)
		FitMethods.__init__(self, lst_file, 'double_gauss', 2, peaks=peaks, verbose=verbose)
		# self.setupUi(self)
		# self.update_par_list()
		self.hist = self.lst2roothist(self.lst_file, bins=2000)

	def init_limits(self):
		self.limits = {
			'mean0':[self.hist.GetMean(1), self.xmin, self.xmax],
			'mean1':[self.hist.GetMean(1), self.xmin, self.xmax],
			'sigma':[self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1)],
			'ratio':[0.5,0,1],
		}

	def call_pdf(self, xmin, xmax, limits=False, index=0):
		"""
		Setup Variable, Parameters and PDF and performs fitting of the PDF to ROI data.
		:param roi_hist: histogram to fit
		:param xmin: range min for the fit
		:param xmax: range max for the fit
		:param index: peak index for display
		:return list: containing all relevant information about the ROI histogram, PDF, parameters,
					  PDF components and fit results
		"""
		self.xmin = xmin
		self.xmax = xmax
		self.hist.GetXaxis().SetRangeUser(self.xmin, self.xmax)
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
		self.x = RooRealVar("x", self.rootitle, self.xmin, self.xmax, self.roounit)
		self.mean0 = RooRealVar('mean{}0'.format(index), 'mean{}0'.format(index),
						  self.limits["mean0"][0], self.limits["mean0"][1], self.limits["mean0"][2])
		self.mean1 = RooRealVar('mean{}1'.format(index), 'mean{}1'.format(index),
						  self.limits["mean1"][0], self.limits["mean1"][1], self.limits["mean1"][2])
		self.sigma = RooRealVar('sigma{}'.format(index), 'sigma{}'.format(index),
						   self.limits["sigma"][0], self.limits["sigma"][1], self.limits["sigma"][2])
		self.ratio = RooRealVar(f'ratio{index}0', "ratio0", 
							self.limits["ratio"][0], self.limits["ratio"][1], self.limits["ratio"][2])
		sig0 = RooGaussian('gauss_{}0'.format(index), 'Gaussian Model {}0'.format(index), self.x, self.mean0, self.sigma)
		sig1 = RooGaussian('gauss_{}1'.format(index), 'Gaussian Model {}1'.format(index), self.x, self.mean1, self.sigma)

		self.this_pdf = RooAddPdf(f'double_Gaussian{index}', f'Two Gaussians added {index}',
								  RooArgList(sig0, sig1), RooArgList(self.ratio))
		# norm = RooRealVar('events{}'.format(index), "Nb of Events", self.roi_counts)
		# self.ext_this_pdf = RooExtendPdf('ext_gauss{}'.format(index), 'Ext Gaussian', self.this_pdf, norm)
		self.roodefs = [self.this_pdf, self.x, self.mean0, self.mean1, self.sigma]
		# self.fix_var(self.roodefs[2:], index)
		self.this_roohist, self.fit_results = self.minimize(self.hist, self.xmin, self.xmax)
		return [self.this_roohist, self.this_pdf, self.nb_pars,
				self.x, self.mean0, self.mean1, self.sigma, self.pstate,
				self.my_name, self.spl_draw, self.rooCompon, self.fit_results]
	
	def plot(self, bins=2000, log=False):
		"""  """
		if len(self.lst_file) == 0:
			print("Fit not excecuted yet or failed.")
			return 0 
		#            
		xdata = self.lst_file.tof
		n, xe = np.histogram(xdata, bins=bins)
		cx = 0.5 * (xe[1:] + xe[:-1])
		dx = np.diff(xe)
		plt.errorbar(cx, n, n ** 0.5, fmt="ok", zorder=1)
		xm = np.linspace(xe[0], xe[-1], num=1000)
		# Plot 'carpet'
		if log:
			plt.plot(xdata, np.zeros_like(xdata)+0.9, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		plt.plot(xdata, np.zeros_like(xdata)-5, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		#
		# Plot fits, loop through all fits
		if log:
			xm = np.linspace(self.mean0.getValV()-100, self.mean1.getValV()+100, num=1000)
		# Calculate values of PDFs for plotting
		y_val = [(self.ratio.getValV() * RootMath.gaussian_pdf(x=y, x0=self.mean0.getValV(), sigma=self.sigma.getValV()) +
				 (1-self.ratio.getValV()) * RootMath.gaussian_pdf(x=y, x0=self.mean1.getValV(), sigma=self.sigma.getValV()) )
				 for y in xm]
		# Normalize values
		integral_cut = sum(y_val) * np.diff(xm)[0]
		left_n_cut = len(xe[xe<self.xmin])
		right_n_cut = len(xe[xe<self.xmax])
		n_cut = n[left_n_cut:right_n_cut]        
		y_val = y_val / integral_cut * sum(n_cut) * dx[0]
		# Plot
		plt.plot(xm, y_val, ls=":", label=f"init", zorder=2)
		plt.plot(xm, y_val, label=f"fit", c='r', zorder=3)
		
		plt.legend();
		if log:
			plt.yscale("log")
		# Zoom in on found peaks
		if self.peaks.n_peaks != 0:
			plt.xlim(self.peaks.earliest_left_base-200, self.peaks.latest_right_base+200)

		plt.show()

class Egh(FitMethods):
	"""
	Class handling the Exponent-Gaussian hybrid model.
	See article -->
	"""
	def __init__(self, lst_file, peaks = 'findit', verbose=False):
		"""
		Initializes the class and number of parameters.
		"""        
		# QtWidgets.QWidget.__init__(self)
		FitMethods.__init__(self, lst_file, 'Egh', 2, peaks=peaks, verbose=verbose)
		# self.setupUi(self)
		# self.update_par_list()
		self.hist = self.lst2roothist(self.lst_file, bins=2000)
	
	def fit_function(self,x,mean,sigma,tau):
		""" """
		return (np.exp( - ( x-mean )**2 / ( 2*sigma**2 + tau*(x-mean) ) )  )
	
	def calc_pole(self,pmean,psigma,ptau):
		""" """
		return ( -(2*psigma**2)/ptau + pmean )

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
		self.xmin = xmin
		self.xmax = xmax
		self.hist.GetXaxis().SetRangeUser(self.xmin, self.xmax)
		self.compose_title_unit(self.hist.GetXaxis().GetTitle())
		# self.set_value(0, hist.GetMean(1))
		# self.set_limits(0, xmin, xmax)
		# self.set_value(1, hist.GetStdDev(1))
		# self.set_limits(1, 0.01*hist.GetStdDev(1), 10*hist.GetStdDev(1))
		tau_estimate = self.hist.GetStdDev(1)*(self.hist.GetSkewness(1) / 2)**(1/3)
		# self.set_value(2, tau_estimate)
		# self.set_limits(2, -10*tau_estimate, 10*tau_estimate)
		#
		self.x = RooRealVar("x", self.rootitle, self.xmin, self.xmax, self.roounit)
		self.mean = RooRealVar("mean{}".format(index), "mean{}".format(index),
						  self.hist.GetMean(1), self.xmin, self.xmax)
		self.sigma = RooRealVar("sigma{}".format(index), "sigma{}".format(index),
						   self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 100*self.hist.GetStdDev(1))
		self.tau = RooRealVar("tau{}".format(index), "tau{}".format(index),
						 tau_estimate, -10*tau_estimate, 10*tau_estimate)
		self.this_pdf = RooGenericPdf('egh_{}'.format(index),
									  'EGH model {}'.format(index),
									  'exp(((-(@0-@1)^2)/(2*@2^2+@3*(@0-@1))))',
									  RooArgList(self.x, self.mean, self.sigma, self.tau))
		self.roodefs = [self.this_pdf, self.x, self.mean, self.sigma, self.tau]
		# self.fix_var(self.roodefs[2:], index)
		self.this_roohist, self.fit_results = self.minimize(self.hist, self.xmin, self.xmax)
		return [self.this_roohist, self.this_pdf, self.nb_pars,
				self.x, self.mean, self.sigma, self.tau, self.pstate,
				self.my_name, self.spl_draw, self.rooCompon, self.fit_results]
	
	def plot(self, bins = 2000, log=False):
		"""  """
		if len(self.lst_file) == 0:
			print("Fit not excecuted yet or failed.")
			return 0 
		#            
		xdata = self.lst_file.tof
		n, xe = np.histogram(xdata, bins=bins)
		# n = n/len(xdata)
		cx = 0.5 * (xe[1:] + xe[:-1])
		dx = np.diff(xe)
		# Plot data
		# plt.scatter(cx, n, zorder=1)
		plt.errorbar(cx, n, n ** 0.5, fmt="ok", zorder=1)
		# Plot 'carpet'
		if log:
			plt.plot(xdata, np.zeros_like(xdata)+0.9, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		plt.plot(xdata, np.zeros_like(xdata)-5, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		# This function has a pole -> limit plotting 
		pole = self.calc_pole(pmean=self.mean.getValV(), psigma=self.sigma.getValV(), ptau=self.tau.getValV())
		# Check if pole is right or left of pole and limit xrange accordingly
		if pole < self.mean.getValV():
			if log:
				x_range_cut = np.linspace(pole+0.01, self.xmax, num=1000)
			else:
				x_range_cut = xe[xe>pole+0.01]
				# data_cut = n[len(n)-len(x_range_cut):]
		else:
			if log:
				x_range_cut = np.linspace(self.xmin, pole-0.01, num=1000)
			else:
				x_range_cut = xe[xe<pole-0.01]
		data_cut = n[:len(n)-len(x_range_cut)]
		#
		# Plot fits, loop through all fits
		#
		# Get the strength of the fitted function by doing a simple chi2 fit        
		# def fit_func(A, x):
		#     """ """
		#     return A * self.fit_function(x=x, mean=self.mean.getValV(), sigma=self.sigma.getValV(), tau=self.tau.getValV())
		# def fun(A):
		#     return fit_func(A, x_range_cut) - data_cut
		# A0 = 1
		# res = least_squares(fun, A0)
		# Calculate fit 
		# y_val = [self.fit_function(x=y, mean=self.mean.getValV(), sigma=self.sigma.getValV(), tau=self.tau.getValV()) * res.x[0]
		#          for y in x_range_cut]
		
		# Normalize by integration
		y_val = [self.fit_function(x=y, mean=self.mean.getValV(), sigma=self.sigma.getValV(), tau=self.tau.getValV())
								  # * len(xdata) * dx[0] 
								  for y in x_range_cut]
		dx_cut = np.diff(x_range_cut)
		integral_cut = sum(y_val) * dx_cut[0]
		left_n_cut = len(xe[xe<self.xmin])
		right_n_cut = len(xe[xe<self.xmax])
		n_cut = n[left_n_cut:right_n_cut]        
		y_val = y_val / integral_cut * sum(n_cut) * dx[0]

		
		# y_val = [self.fit_function(x=y, mean=self.mean.getValV(), sigma=self.sigma.getValV(), tau=self.tau.getValV())
		#                           * len(xdata) * dx[0] 
		#                           for y in x_range_cut]
		plt.plot(x_range_cut, y_val, ls=":", label=f"init", zorder=2)
		plt.plot(x_range_cut, y_val, label=f"fit", c='r', zorder=3, linewidth=3)
		
		# Get y axis limits
		ylims = plt.ylim()
		
		plt.legend();
		if log:
			plt.yscale("log")
			plt.ylim(0.9,2*ylims[1])
		# Zoom in on found peaks
		if self.peaks.n_peaks != 0:
			plt.xlim(self.peaks.earliest_left_base-200, self.peaks.latest_right_base+200)
		else:
			plt.xlim(self.mean.getValV()-200, self.mean.getValV()+200)

		plt.show()

class Emg(FitMethods):
	"""
	Class handling the exponentially modified Gaussian as defined in 
	https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
	"""
	def __init__(self, lst_file, peaks = 'findit', verbose=False):
		"""
		Initializes the class and number of parameters.
		"""        
		# QtWidgets.QWidget.__init__(self)
		FitMethods.__init__(self, lst_file, 'Egh', 2, peaks=peaks, verbose=verbose)
		# self.setupUi(self)
		# self.update_par_list()
		self.hist = self.lst2roothist(self.lst_file, bins=2000)

	def init_limits(self):
		tau_estimate = self.hist.GetStdDev(1)*(self.hist.GetSkewness(1) / 2)**(1/3)
		llambda_estimate = 1/tau_estimate
		print(llambda_estimate)
		self.limits = {
			'mu':[self.hist.GetMean(1), self.xmin, self.xmax],
			'llambda':[llambda_estimate, 0.001, 10*llambda_estimate],
			'sigma':[self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1)],
		}

	
	def fit_function(self,x,mu,sigma,llambda):
		""" """
		return ( (llambda/2) * np.exp( (llambda/2) * ( 2*mu + llambda*np.power(sigma, 2)-2*x ) ) *
				 (1 - sc.special.erf( (mu+llambda*np.power(sigma, 2)-x)/(math.sqrt(2)*sigma) ) )
			   )
	
	def call_pdf(self, xmin, xmax, limits=False, index=0):
		"""
		Setup Variable, Parameters and PDF and performs fitting of the PDF to ROI data.
		:param roi_hist: histogram to fit
		:param xmin: range min for the fit
		:param xmax: range max for the fit
		:param index: peak index for display
		:return list: containing all relevant information about the ROI histogram, PDF, parameters,
					  PDF components and fit results
		"""
		# Test dimension of limits
		self.xmin = xmin
		self.xmax = xmax
		if isinstance(xmin, list):
			print("(FIT_LUKAS::Emg::call_pdf): Multiple fit-regions not supported in this fit function")
			return 0
		#
		self.compose_title_unit(self.hist.GetXaxis().GetTitle())
		#
		if not limits:
			self.init_limits()
		else:
			self.limits = limits
		#
		# self.set_value(0, hist.GetMean(1))
		# self.set_limits(0, xmin, xmax)
		# self.set_value(1, hist.GetStdDev(1))
		# self.set_limits(1, 0.01*hist.GetStdDev(1), 10*hist.GetStdDev(1))
		# self.set_value(2, tau_estimate)
		# self.set_limits(2, -10*tau_estimate, 10*tau_estimate)
		#
		self.x = RooRealVar("x", self.rootitle, self.xmin, self.xmax, self.roounit)
		self.mu = RooRealVar("mu{}".format(index), "mu{}".format(index),
						  self.limits["mu"][0], self.limits["mu"][1], self.limits["mu"][2])
		self.sigma = RooRealVar("sigma{}".format(index), "sigma{}".format(index),
						  self.limits["sigma"][0], self.limits["sigma"][1], self.limits["sigma"][2])
		self.llambda = RooRealVar("llambda{}".format(index), "llambda{}".format(index),
						  self.limits["llambda"][0], self.limits["llambda"][1], self.limits["llambda"][2])
		self.this_pdf = RooGenericPdf('emg_{}'.format(index),
									  'Emg model {}'.format(index),
									  '@3/2*exp(@3/2*(2*@1+@3*@2^2-2*x))*(1-TMath::Erf((@1+@3*@2^2-x)/(sqrt(2)*@2)))',
									  RooArgList(self.x, self.mu, self.sigma, self.llambda))
		self.roodefs = [self.this_pdf, self.x, self.mu, self.sigma, self.llambda]
		# self.fix_var(self.roodefs[2:], index)
		self.this_roohist, self.fit_results = self.minimize(self.hist, self.xmin, self.xmax)
		return [self.this_roohist, self.this_pdf, self.nb_pars,
				self.x, self.mu, self.sigma, self.llambda, self.pstate,
				self.my_name, self.spl_draw, self.rooCompon, self.fit_results]
	
	def plot(self, bins = 2000, log=False):
		"""  """
		if len(self.lst_file) == 0:
			print("Fit not excecuted yet or failed.")
			return 0 
		#            
		xdata = self.lst_file.tof
		n, xe = np.histogram(xdata, bins=bins)
		# n = n/len(xdata)
		cx = 0.5 * (xe[1:] + xe[:-1])
		dx = np.diff(xe)
		# Plot data
		# plt.scatter(cx, n, zorder=1)
		plt.errorbar(cx, n, n ** 0.5, fmt="ok", zorder=1)
		xm = np.linspace(xe[0], xe[-1], num=5000)
		# Plot 'carpet'
		if log:
			plt.plot(xdata, np.zeros_like(xdata)+0.9, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		plt.plot(xdata, np.zeros_like(xdata)-5, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		# if log:
		# 	xm = np.linspace(self.mu.getValV()-100, self.mu.getValV()+100, num=1000)
		#
		# Plot fits, loop through all fits
		#
		# Get the strength of the fitted function by doing a simple chi2 fit        
		# def fit_func(A, x):
		#     """ """
		#     return A * self.fit_function(x=x, mean=self.mean.getValV(), sigma=self.sigma.getValV(), tau=self.tau.getValV())
		# def fun(A):
		#     return fit_func(A, x_range_cut) - data_cut
		# A0 = 1
		# res = least_squares(fun, A0)
		# Calculate fit 
		# y_val = [self.fit_function(x=y, mean=self.mean.getValV(), sigma=self.sigma.getValV(), tau=self.tau.getValV()) * res.x[0]
		#          for y in x_range_cut]
		
		# Normalize by integration
		y_val = [self.fit_function(x=y, mu=self.mu.getValV(), sigma=self.sigma.getValV(), llambda=self.llambda.getValV())
								  # * len(xdata) * dx[0] 
								  for y in xm]

		# Normalize values
		integral_cut = sum(y_val) * np.diff(xm)[0]
		left_n_cut = len(xe[xe<self.xmin])
		right_n_cut = len(xe[xe<self.xmax])
		n_cut = n[left_n_cut:right_n_cut]        
		y_val = y_val / integral_cut * sum(n_cut) * dx[0]

		
		# y_val = [self.fit_function(x=y, mean=self.mean.getValV(), sigma=self.sigma.getValV(), tau=self.tau.getValV())
		#                           * len(xdata) * dx[0] 
		#                           for y in x_range_cut]
		plt.plot(xm, y_val, ls=":", label=f"init", zorder=2)
		plt.plot(xm, y_val, label=f"fit", c='r', zorder=3, linewidth=3)
		
		# Get y axis limits
		ylims = plt.ylim()
		
		plt.legend();
		if log:
			plt.yscale("log")
			plt.ylim(0.9,2*ylims[1])
		# Zoom in on found peaks
		if self.peaks.n_peaks != 0:
			plt.xlim(self.peaks.earliest_left_base-200, self.peaks.latest_right_base+200)
		else:
			plt.xlim(self.mu.getValV()-200, self.mu.getValV()+200)

		plt.show()

class hyperEmg(FitMethods):
	"""
	Class handling the exponentially modified Gaussian as defined in hyper-EMG paper from JLU-GSI group
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
		self.pos_funct = '1/(2*@3)*exp((@2/(1.4142*@3))^2-(@0-@1)/@3)*TMath::Erfc(@2/(1.4142*@3)-(@0-@1)/(1.4142*@2))'
		self.neg_funct = '1/(2*@3)*exp((@2/(1.4142*@3))^2+(@0-@1)/@3)*TMath::Erfc(@2/(1.4142*@3)+(@0-@1)/(1.4142*@2))'
		self.RooRealVar_dict = {}
		self.RooGenericPdf_dict = {}
		if file_path == False: 
			self.base_file_name = 'path_to_file'
		else:
			self.base_file_name = file_path

	def pos_func(self, x, mu, sigma, ptau):
		"""
		Analytical function for positive Emg components
		"""
		return 1/(2*ptau)*np.exp((sigma/(1.4142*ptau))**2-(x-mu)/ptau)*sc.special.erfc(sigma/(1.4142*ptau)-(x-mu)/(1.4142*sigma))
	
	def neg_func(self, x, mu, sigma, ntau):
		"""
		Analytical function for negative Emg components
		"""
		return 1/(2*ntau)*np.exp((sigma/(1.4142*ntau))**2+(x-mu)/ntau)*sc.special.erfc(sigma/(1.4142*ntau)+(x-mu)/(1.4142*sigma))

	def find_peak_and_fwhm(self):
		# Taken form $ROOTSYS/tutorials/fit/langaus.C
		# Seaches for the location (x value) at the maximum of the
		# combined pdf and its full width at half-maximum.
		#
		# The search is probably not very efficient, but it's a first try.
		i = 0;
		MAXCALLS = 10000
		# Search for maximum
		p = self.RooRealVar_dict['mu0'].getValV() - 0.1 * self.RooRealVar_dict['sigma'].getValV()
		step = 0.05 * self.RooRealVar_dict['sigma'].getValV()
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
		p = maxx + self.RooRealVar_dict['sigma'].getValV()
		step = self.RooRealVar_dict['sigma'].getValV()
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
		p = maxx - 0.5 *self.RooRealVar_dict['sigma'].getValV()
		step = -self.RooRealVar_dict['sigma'].getValV()
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
		f.write(f"dimensions={self.dimensions}\n")
		f.write(f"n_comps={self.n_comps}\n")
		f.write(f"xmin={self.xmin}\n")
		f.write(f"xmax={self.xmax}\n")
		f.write(f"numerical_peak={self.numerical_peak}\n")
		f.write(f"numerical_FWHM={self.numerical_FWHM}\n")
		f.write(f"numerical_FWHM_left={self.numerical_FWHM_left}\n")
		f.write(f"numerical_FWHM_right={self.numerical_FWHM_right}\n")
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
		xm = np.linspace(self.xmin, self.xmax, num=5000)
		# Get fit results for x-axis
		y_val = []
		for i in xm:
			self.RooRealVar_dict['x'].setVal(i)
			f.write(f"{i:.3f} {self.this_pdf.getVal(ROOT.RooArgSet(self.RooRealVar_dict['x'])):.6f}\n")
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
		fitfromfile = FitToDict(file)
		#
		if not 'FIT-VALUES' in fitfromfile.fit:
			print(f"Fit file {file} has no fit values that were exported.")
			return 0
		#
		xm = fitfromfile.fit['FIT-VALUES'].tof
		y_val = fitfromfile.fit['FIT-VALUES'].fit
		self.xmin = float(fitfromfile.fit['META-DATA']['xmin'])
		self.xmax = float(fitfromfile.fit['META-DATA']['xmax'])
		self.numerical_peak = float(fitfromfile.fit['META-DATA']['numerical_peak'])
		self.numerical_FWHM = float(fitfromfile.fit['META-DATA']['numerical_FWHM'])
		self.numerical_FWHM_left = float(fitfromfile.fit['META-DATA']['numerical_FWHM_left'])
		self.numerical_FWHM_right = float(fitfromfile.fit['META-DATA']['numerical_FWHM_right'])
		self.dimensions = np.fromstring(fitfromfile.fit['META-DATA']['dimensions'].strip('[]'), sep=',', dtype=int)
		self.fit_func_name = fitfromfile.fit['META-DATA']['fit_function']
		self.n_comps = int(fitfromfile.fit['META-DATA']['n_comps'])
		#
		for idx,row in fitfromfile.fit['RESULTS-TABLE'].iterrows():
			self.Var_dict[row['var']] = row['value']
		#
		return xm, y_val

	def constraints_from_file(self, file, params_to_fix):
		'''
		Loads and sets constraints for fit from fit file
		Parameters:
			- file: fit file to be loaded
			- params_to_fix: array of strings of parameters to fix
		'''
		#
		fitfromfile = FitToDict(file)
		#
		if not 'FIT-VALUES' in fitfromfile.fit:
			print(f"Fit file {file} has no fit values that were exported.")
			return 0
		#
		for idx,row in fitfromfile.fit['RESULTS-TABLE'].iterrows():
			if row['var'] in params_to_fix:
				self.limits[row['var']] = [row['value'], row['value'], row['value']]


	def call_pdf(self, xmin, xmax, dimensions = [1,2], n_comps = 1, bins = 1, limits=False):
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
		# Convert histogram
		self.bins = bins
		self.binning = self.get_binning(bins=self.bins)
		self.n_comps = n_comps
		self.hist = self.lst2roothist(self.lst_file, bins=1)
		# Test dimension of limits
		self.xmin = xmin
		self.xmax = xmax
		if isinstance(xmin, list):
			print("(FIT_LUKAS::hyperEmg::call_pdf): Multiple fit-regions not supported in this fit function")
			return 0
		#
		self.compose_title_unit(self.hist.GetXaxis().GetTitle())

		# Initialize limits
		if not limits:
			self.init_limits()
		#  
		else:
			# If limits are passed, save them in dict
			# Careful! This overwrite the constraints set by constraints_from_file
			for limit in limits:
				self.limits[limit] = limits[limit]
#
		# Fill variable dictionary:
		self.RooRealVar_dict = {
			'x': RooRealVar("x", self.rootitle, self.xmin, self.xmax, self.roounit)
		}
		# Constant background
		# self.RooRealVar_dict["const_bck"] = RooRealVar("const_bck", "const_bck", self.limits["const_bck"][0], self.limits["const_bck"][1], self.limits["const_bck"][2])
		# self.RooGenericPdf_dict["const_bck"] = RooPolynomial("const_bck","const_bck", self.RooRealVar_dict['x'], RooArgList())
		# Define Gaussian components for the fit
		for i in range(0,n_comps,1):
			var_name = f"mu{i}"
			self.RooRealVar_dict[var_name] = RooRealVar(var_name, var_name, self.limits[var_name][0], self.limits[var_name][1], self.limits[var_name][2])
		self.RooRealVar_dict["sigma"] = RooRealVar("sigma", "sigma", self.limits["sigma"][0], self.limits["sigma"][1], self.limits["sigma"][2])
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
		# Definition of ratio factors: n-1 factors for n components. (+1 ratio for constant background)
		n_components = n_comps * (dimensions[0] + dimensions[1]) - 1
		for i in range(0,n_components,1):
			var_name = f"ratio{i}"
			self.RooRealVar_dict[var_name] = RooRealVar(var_name, var_name, self.limits[var_name][0], self.limits[var_name][1], self.limits[var_name][2])
		# Definition for the negative EMG components
		for i in range(0,dimensions[0],1):
			for j in range(0,n_comps,1):
				pdf_name = f"nsemg{i}{j}"
				self.RooGenericPdf_dict[pdf_name] = RooGenericPdf(pdf_name,pdf_name,self.neg_funct,
														RooArgList(self.RooRealVar_dict['x'], 
																   self.RooRealVar_dict[f'mu{j}'], 
																   self.RooRealVar_dict['sigma'], 
																   self.RooRealVar_dict[f'ntau{i}']))
		# Definition for the positive EMG components
		for i in range(0,dimensions[1],1):
			for j in range(0,n_comps,1):
				pdf_name = f"psemg{i}{j}"
				self.RooGenericPdf_dict[pdf_name] = RooGenericPdf(pdf_name,pdf_name,self.pos_funct,
														RooArgList(self.RooRealVar_dict['x'], 
																   self.RooRealVar_dict[f'mu{j}'], 
																   self.RooRealVar_dict['sigma'], 
																   self.RooRealVar_dict[f'ptau{i}']))
		# Put all pdfs together
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
		self.this_pdf = RooAddPdf('hyperEmg', 'hyperEmg', all_pdfs, all_ratios)
		# self.this_pdf = RooAddPdf('hyperEmg', 'hyperEmg', RooArgList(self.psEmg, self.nsEmg), RooArgList(self.emgratio))
		# 
		self.roodefs = [self.this_pdf, self.RooRealVar_dict, self.RooGenericPdf_dict]
		self.this_roohist, self.fit_results = self.minimize(self.hist, self.xmin, self.xmax)

		# THIS LINE HAS TO BE HERE FOR SOME UNKNOW REASON TO ME TO AVOID A SEG FAULT IN THE PLOT FUNCTION GOD KNOWS WHY I REALLY DON'T KNOW ANYMORE...
		print(self.this_pdf.getVal(ROOT.RooArgSet(self.RooRealVar_dict['x'])))

		# Calculate numerical position of maximum and FWHM
		self.numerical_peak, self.numerical_FWHM, self.numerical_FWHM_left, self.numerical_FWHM_right = self.find_peak_and_fwhm()

		# Store fit results in dict
		for key in self.RooRealVar_dict:
			self.Var_dict[key] = self.RooRealVar_dict[f'{key}'].getValV()

		print(self.Var_dict)

		return [self.this_roohist, self.roodefs, self.this_pdf, self.my_name, self.spl_draw, self.rooCompon, self.fit_results]
	
	def plot(self, bins = 1, log=False, focus=-1, from_file = False, file_out=False, silent=True, 
			 centroids=False, components=False, carpet=False, legend=True):
		"""  
		Wrapper for plotting the fit and data
			- bins: number of bins to rebin. Defaults to 1, e.g. no rebinning
			- log: plot y-scale in log. Defaults to false
			- from_file: file to read fit from if already fitted earlier .pdf. Defaults to False, fit needs to be executed beforehand 
			- file_out: name to save plot as .pdf. Defaults to False, e.g. no saving
			- silent: True shows plot, false does not show plot (for silent printing)
			- centroids: True shows centroids of Gaussian components, as well as location of FWHM of main peak. Defaults to False.
			- components: Plots components to EGH as grey dashed lines
			- carpet: If true, plots carpet 
			- legend: Plot legend if
		"""
		if len(self.lst_file) == 0:
			print("Fit not excecuted yet or failed.")
			return 0 
		#       
		self.bins = bins
		xdata = self.lst_file.tof
		n, xe = np.histogram(xdata, bins=self.get_binning(self.bins))
		# n = n/len(xdata)
		cx = 0.5 * (xe[1:] + xe[:-1]) 
		dx = np.diff(xe)
		# Plot data
		# plt.scatter(cx, n, zorder=1)
		
		# Get fit and prepare fit parameters for plotting,
		#	- Either from the call_pdf during the same program excecution
		# 	- or from the fit file that saves the (not normalized) fit
		if not from_file:
			# X-Axis for fit
			xm = np.linspace(self.xmin, self.xmax, num=5000)
			# Get fit results for x-axis
			y_val = []
			for i in xm:
				self.RooRealVar_dict['x'].setVal(i)
				y_val.append(self.this_pdf.getVal(ROOT.RooArgSet(self.RooRealVar_dict['x'])))	
		else:
			# Read fit file
			xm, y_val = self.load_fit(from_file)

		# Plot data
		plt.errorbar(cx - self.numerical_peak,
					 n, n ** 0.5, fmt="ok", zorder=1, label=f"Data (bins={bins})")

		# Normalize values
		integral_cut = sum(y_val) * np.diff(xm)[0]
		left_n_cut = len(xe[xe<self.xmin])
		right_n_cut = len(xe[xe<self.xmax])
		n_cut = n[left_n_cut:right_n_cut]        
		y_val = y_val / integral_cut * sum(n_cut) * dx[0]

		# Plot fit	
		plt.plot(xm - self.numerical_peak, 
				 y_val, label=f"{self.fit_func_name}({self.dimensions[0]},{self.dimensions[1]})", c='r', zorder=3, linewidth=3)
		
		# Plot 'carpet'
		if carpet:
			if log:
				plt.plot(xdata - self.numerical_peak, 
						 np.zeros_like(xdata)+0.9, "|", alpha=0.1, zorder = 3)
			plt.plot(xdata - self.numerical_peak, 
					 np.zeros_like(xdata)-5, "|", alpha=0.1, zorder = 3)

		# Plot numerical peak position and FWHM
		if centroids:
			plt.axvline(self.numerical_peak - self.numerical_peak, c='r', linewidth=1, zorder=3)
			plt.axvline(self.numerical_FWHM_left - self.numerical_peak, c='blue', linewidth=1, zorder=3)
			plt.axvline(self.numerical_FWHM_right - self.numerical_peak, c='blue', linewidth=1, zorder=3)
			# Plot center of gaussian component
			for comp in range(0,self.n_comps,1):
				plt.axvline(self.Var_dict[f'mu{comp}'] - self.numerical_peak, c='r', linewidth=1, zorder=3)

		# Plot components
		if components:
			i_ratio = 0 
			n_ratios = self.n_comps * (self.dimensions[0] + self.dimensions[1]) - 1
			# Loop through components as the ratios are assigned: inner loop must be the different species while outer loop is the egh component
			for dim in range(0,self.dimensions[0],1):
				for comp in range(0,self.n_comps,1):
					y_val = self.neg_func(x=xm,mu=self.Var_dict[f'mu{comp}'], 
												sigma=self.Var_dict[f'sigma'], 
												ntau=self.Var_dict[f'ntau{dim}'])
					# normalize accounting for fit ratios
					integral_cut = sum(y_val) * np.diff(xm)[0]
					y_val = y_val / integral_cut * sum(n_cut) * dx[0]
					if (self.n_comps == 0 and (self.dimensions[0]+self.dimensions[1]) == 1):
						ratio = 1
					else:
						if i_ratio != n_ratios:
							y_val *= self.Var_dict[f'ratio{i_ratio}']
						else:
							y_val *= 1-sum([self.Var_dict[f'ratio{r}'] for r in np.arange(0,n_ratios,1)])
					i_ratio += 1
					# Plot
					plt.plot(xm - self.numerical_peak, 
							 y_val, label=f"Neg. component {comp}:{dim})", c='grey', ls="--", zorder=2, linewidth=1.75)

			for dim in range(0,self.dimensions[1],1):
				for comp in range(0,self.n_comps,1):
					y_val = self.pos_func(x=xm,mu=self.Var_dict[f'mu{comp}'], 
												sigma=self.Var_dict[f'sigma'], 
												ptau=self.Var_dict[f'ptau{dim}'])
					# normalize accounting for fit ratios
					integral_cut = sum(y_val) * np.diff(xm)[0]
					y_val = y_val / integral_cut * sum(n_cut) * dx[0]
					if (self.n_comps == 0 and (self.dimensions[0]+self.dimensions[1]) == 1):
						ratio = 1
					else:
						if i_ratio != n_ratios:
							y_val *= self.Var_dict[f'ratio{i_ratio}']
						else:
							y_val *= 1-sum([self.Var_dict[f'ratio{r}'] for r in np.arange(0,n_ratios,1)])
					i_ratio += 1
					# Plot
					plt.plot(xm - self.numerical_peak, 
							 y_val, label=f"Pos. component {comp}:{dim})", c='grey', ls="--", zorder=2, linewidth=1.75)

		# Get y axis limits
		ylims = plt.ylim()
		if log:
			plt.yscale("log")
			plt.ylim(0.1,2*ylims[1])

		# Zoom in on found peaks
		if self.peaks.n_peaks != 0:
			plt.xlim(self.peaks.earliest_left_base - 400 - self.numerical_peak, 
					 self.peaks.latest_right_base + 400 - self.numerical_peak)
			if focus != -1:
				plt.xlim(self.peaks.pos[focus] - 600 - self.numerical_peak, 
						 self.peaks.pos[focus] + 600 - self.numerical_peak)

		# Add axis labels
		plt.xlabel(f'Time-of-Flight [ns] - {self.numerical_peak:.1f}ns', fontsize=20)
		plt.ylabel(f'Counts per bin', fontsize=20)

		# Format Legend
		if legend:
			plt.legend(fontsize=20)

		# Save plot
		if file_out != False:
			print(f"Plot fit save as {file_out}")
			plt.savefig(file_out, dpi=300)

		# Show plot on canvas
		if silent == False:
			plt.show()

		# Clear canvas to avoid printing on top of other plot in batch mode
		if silent:
			plt.clf()

class hyperEmg2(FitMethods):
	"""
	Class handling the exponentially modified Gaussian as defined in 
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
		self.pos_funct = '1/(2*@3)*exp((@2/(1.4142*@3))^2-(@0-@1)/@3)*TMath::Erfc(@2/(1.4142*@3)-(@0-@1)/(1.4142*@2))'
		self.neg_funct = '1/(2*@3)*exp((@2/(1.4142*@3))^2+(@0-@1)/@3)*TMath::Erfc(@2/(1.4142*@3)+(@0-@1)/(1.4142*@2))'
		self.RooRealVar_dict = {}
		self.RooGenericPdf_dict = {}
		if file_path == False: 
			self.base_file_name = 'path_to_file'
		else:
			self.base_file_name = file_path

	def init_limits(self):
		tau_estimate = self.hist.GetStdDev(1)*(self.hist.GetSkewness(1) / 2)**(1/3)
		llambda_estimate = 1/tau_estimate
		print(llambda_estimate)
		self.limits = {
			'mu':[self.hist.GetMean(1), self.xmin, self.xmax],
			'llambda':[llambda_estimate, 0.001, 10*llambda_estimate],
			'sigma':[self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1)],
		}

	def function_string(self, c0='@0', c1='@1', c2='@2', c3='@3', R='@4',pol="+"):
		return f'{R}*1/(2*{c3})*exp(({c2}/(1.4142*{c3}))^2{pol}({c0}-{c1})/{c3})*TMath::Erfc({c2}/(1.4142*{c3}){pol}({c0}-{c1})/(1.4142*{c2}))'

	def build_function_string(self, dimensions = [0,1]):
		#
		funct = ''
		# 
		comp = 0
		#
		for i in range(dimensions[0]):
			if (comp == dimensions[0]+dimensions[1]-1 and dimensions[0]+dimensions[1] > 1):
				ratio = '(1-'
				for j in range(dimensions[0]+dimensions[1]-1):
					ratio+=f"@{4+2*j}-"
				ratio = ratio[:-1] + ")"
			else: 
				ratio = f"@{4+2*comp}"
			funct += self.function_string(c3=f"@{3+2*comp}", R=f"{ratio}", pol="+")+ "+"
			comp+=1
		#
		for i in range(dimensions[1]):
			if (comp == dimensions[0]+dimensions[1]-1 and dimensions[0]+dimensions[1] > 1):
				ratio = '(1-'
				for j in range(dimensions[0]+dimensions[1]-1):
					ratio+=f"@{4+2*j}-"
				ratio = ratio[:-1] + ")"
			else: 
				ratio = f"@{4+2*comp}"
			funct += self.function_string(c3=f"@{3+2*comp}", R=f"{ratio}", pol="-")+ "+"
			comp+=1

		return funct[:-1] # cut away last character which is a "+"

	def pos_func(self, x, mu, sigma, ptau):
		"""
		Analytical function for positive Emg components
		"""
		return 1/(2*ptau)*np.exp((sigma/(1.4142*ptau))**2-(x-mu)/ptau)*sc.special.erfc(sigma/(1.4142*ptau)-(x-mu)/(1.4142*sigma))
	
	def neg_func(self, x, mu, sigma, ntau):
		"""
		Analytical function for negative Emg components
		"""
		return 1/(2*ntau)*np.exp((sigma/(1.4142*ntau))**2+(x-mu)/ntau)*sc.special.erfc(sigma/(1.4142*ntau)+(x-mu)/(1.4142*sigma))

	def find_peak_and_fwhm(self):
		# Taken form $ROOTSYS/tutorials/fit/langaus.C
		# Seaches for the location (x value) at the maximum of the
		# combined pdf and its full width at half-maximum.
		#
		# The search is probably not very efficient, but it's a first try.
		i = 0;
		MAXCALLS = 10000
		# Search for maximum
		p = self.RooRealVar_dict['mu0'].getValV() - 0.1 * self.RooRealVar_dict['sigma'].getValV()
		step = 0.05 * self.RooRealVar_dict['sigma'].getValV()
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
		p = maxx + self.RooRealVar_dict['sigma'].getValV()
		step = self.RooRealVar_dict['sigma'].getValV()
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
		p = maxx - 0.5 *self.RooRealVar_dict['sigma'].getValV()
		step = -self.RooRealVar_dict['sigma'].getValV()
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
		f.write(f"xmin={self.xmin}\n")
		f.write(f"xmax={self.xmax}\n")
		f.write(f"numerical_peak={self.numerical_peak}\n")
		f.write(f"numerical_FWHM={self.numerical_FWHM}\n")
		f.write(f"numerical_FWHM_left={self.numerical_FWHM_left}\n")
		f.write(f"numerical_FWHM_right={self.numerical_FWHM_right}\n")
		#
		f.write(f"[RESULTS-TABLE]\n")
		f.write(f"var value error error_lo error_high var_init var_lim_lo var_lim_hi\n")
		for var in self.RooRealVar_dict:
			if var == 'x': continue
			f.write(f"{var} {self.RooRealVar_dict[var].getValV()} {self.RooRealVar_dict[var].getError()} ")
			f.write(f"{self.RooRealVar_dict[var].getErrorLo()} {self.RooRealVar_dict[var].getErrorHi()} ")
			f.write(f"{self.limits[var][0]} {self.limits[var][1]} {self.limits[var][2]}\n")
		f.close()

	def call_pdf(self, xmin, xmax, dimensions = [1,2], n_comps = 1, bins = 1, limits=False):
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
		# Convert histogram
		self.bins = bins
		self.binning = self.get_binning(bins=self.bins)
		self.n_comps = n_comps
		self.hist = self.lst2roothist(self.lst_file, bins=self.bins)
		# Test dimension of limits
		self.xmin = xmin
		self.xmax = xmax
		if isinstance(xmin, list):
			print("(FIT_LUKAS::hyperEmg::call_pdf): Multiple fit-regions not supported in this fit function")
			return 0
		#
		self.compose_title_unit(self.hist.GetXaxis().GetTitle())
		#
		if not limits:
			self.init_limits()
		else:
			self.limits = limits
		# Build variable dictionary:
		self.RooRealVar_dict = {
			'x': RooRealVar("x", self.rootitle, self.xmin, self.xmax, self.roounit)
		}
		# Define Gaussian components for the fit
		for i in range(0,n_comps,1):
			var_name = f"mu{i}"
			self.RooRealVar_dict[var_name] = RooRealVar(var_name, var_name, self.limits[var_name][0], self.limits[var_name][1], self.limits[var_name][2])
		self.RooRealVar_dict["sigma"] = RooRealVar("sigma", "sigma", self.limits["sigma"][0], self.limits["sigma"][1], self.limits["sigma"][2])
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
		# Definition of ratio factors: n-1 factors for n components. 
		n_components = (dimensions[0] + dimensions[1]) - 1
		for i in range(0,n_components,1):
			var_name = f"ratio{i}"
			self.RooRealVar_dict[var_name] = RooRealVar(var_name, var_name, self.limits[var_name][0], self.limits[var_name][1], self.limits[var_name][2])
		# Build component 0 and 1 PDFs based on the ratios passed 
		# Start by building the PDF function string
		funct = self.build_function_string(dimensions = dimensions)

		for j in range(0,n_comps,1):
			# Build argument lists
			arglist = RooArgList(self.RooRealVar_dict['x'], self.RooRealVar_dict[f'mu{j}'], self.RooRealVar_dict['sigma']) 
			for i in range(0,dimensions[0],1):
				arglist.add(self.RooRealVar_dict[f"ntau{i}"])
				arglist.add(self.RooRealVar_dict[f"ratio{i}"])
			for i in range(0,dimensions[1],1):
				arglist.add(self.RooRealVar_dict[f"ptau{i}"])
				if i != (dimensions[1] - 1): 
					arglist.add(self.RooRealVar_dict[f"ratio{i+dimensions[0]}"])
			# Build pdf
			pdf_name = f"comp{j}"
			self.RooGenericPdf_dict[pdf_name] = RooGenericPdf(pdf_name,pdf_name,funct,arglist)

			# Build component ratios
			if j < (n_comps-1): 
				var_name = f"ratio_comp{j}"
				self.RooRealVar_dict[var_name] = RooRealVar(var_name, var_name, self.limits[var_name][0], self.limits[var_name][1], self.limits[var_name][2])

		# Put all pdfs together
		all_pdfs = RooArgList()
		for pdf in self.RooGenericPdf_dict:
			all_pdfs.add(self.RooGenericPdf_dict[pdf])

		# Build list of argument ratios
		# put all ratio's together
		all_ratios = RooArgList()
		for var in self.RooRealVar_dict:
			m = re.search('ratio_comp*', var)
			if m:
				all_ratios.add(self.RooRealVar_dict[var])

		# Definition hyper-EMG
		self.this_pdf = RooAddPdf('hyperEmg', 'hyperEmg', all_pdfs, all_ratios)
		# self.this_pdf = RooAddPdf('hyperEmg', 'hyperEmg', RooArgList(self.psEmg, self.nsEmg), RooArgList(self.emgratio))
		# 
		self.roodefs = [self.this_pdf, self.RooRealVar_dict, self.RooGenericPdf_dict]
		self.this_roohist, self.fit_results = self.minimize(self.hist, self.xmin, self.xmax)

		# THIS LINE HAS TO BE HERE FOR SOME UNKNOW REASON TO ME TO AVOID A SEG FAULT IN THE PLOT FUNCTION GOD KNOWS WHY I REALLY DON'T KNOW ANYMORE...
		print(self.this_pdf.getVal(ROOT.RooArgSet(self.RooRealVar_dict['x'])))

		# Calculate numerical position of maximum and FWHM
		self.numerical_peak, self.numerical_FWHM, self.numerical_FWHM_left, self.numerical_FWHM_right = self.find_peak_and_fwhm()

		return [self.this_roohist, self.roodefs, self.this_pdf, self.my_name, self.spl_draw, self.rooCompon, self.fit_results]

	def plot(self, bins = 1, log=False, save=False, silent=True):
		"""  """
		if len(self.lst_file) == 0:
			print("Fit not excecuted yet or failed.")
			return 0 
		#       
		# Convert histogram
		self.bins = bins
		xdata = self.lst_file.tof
		n, xe = np.histogram(xdata, bins=self.get_binning(bins=self.bins))
		# n = n/len(xdata)
		cx = 0.5 * (xe[1:] + xe[:-1])
		dx = np.diff(xe)
		# Plot data
		# plt.scatter(cx, n, zorder=1)
		plt.errorbar(cx, n, n ** 0.5, fmt="ok", zorder=1, label=f"Data (bins={bins})")
		xm = np.linspace(self.xmin, self.xmax, num=5000)
		# Plot 'carpet'
		if log:
			plt.plot(xdata, np.zeros_like(xdata)+0.9, "|", alpha=0.1, zorder = 3)
		plt.plot(xdata, np.zeros_like(xdata)-5, "|", alpha=0.1, zorder = 3)
		# if log:
		# 	xm = np.linspace(self.mu.getValV()-100, self.mu.getValV()+100, num=1000)
		#
		y_val = []
		for i in xm:
			self.RooRealVar_dict['x'].setVal(i)
			y_val.append(self.this_pdf.getVal(ROOT.RooArgSet(self.RooRealVar_dict['x'])))			
		# Normalize values
		integral_cut = sum(y_val) * np.diff(xm)[0]
		left_n_cut = len(xe[xe<self.xmin])
		right_n_cut = len(xe[xe<self.xmax])
		n_cut = n[left_n_cut:right_n_cut]        
		y_val = y_val / integral_cut * sum(n_cut) * dx[0]
		# Plot fit	
		plt.plot(xm, y_val, label=f"{self.fit_func_name}({self.dimensions[0]},{self.dimensions[1]})", c='r', zorder=3, linewidth=3)
		
		# Plot numerical peak position and FWHM
		plt.axvline(self.numerical_peak, c='r', linewidth=1, zorder=3)
		plt.axvline(self.numerical_FWHM_left, c='blue', linewidth=1, zorder=3)
		plt.axvline(self.numerical_FWHM_right, c='blue', linewidth=1, zorder=3)
		# Plot center of gaussian component
		for comp in range(0,self.n_comps,1):
			plt.axvline(self.RooRealVar_dict[f'mu{comp}'].getValV(), c='r', linewidth=1, zorder=3)

		# Plot components
		n_ratios = (self.dimensions[0] + self.dimensions[1]) - 1
		for comp in range(0,self.n_comps,1):
			i_ratio = 0 
			for dim in range(0,self.dimensions[0],1):
				y_val = self.neg_func(x=xm,mu=self.RooRealVar_dict[f'mu{comp}'].getValV(), 
											sigma=self.RooRealVar_dict[f'sigma'].getValV(), 
											ntau=self.RooRealVar_dict[f'ntau{dim}'].getValV())
				# normalize accounting for fit ratios
				integral_cut = sum(y_val) * np.diff(xm)[0]
				y_val = y_val / integral_cut * sum(n_cut) * dx[0]
				if (self.n_comps == 0 and (self.dimensions[0]+self.dimensions[1]) == 1):
					ratio = 1
				else:
					if i_ratio != n_ratios:
						y_val *= self.RooRealVar_dict[f'ratio{i_ratio}'].getValV()
					else:
						y_val *= 1-sum([self.RooRealVar_dict[f'ratio{r}'].getValV() for r in np.arange(0,n_ratios,1)])
				#
				if comp < (self.n_comps-1):
					y_val *= self.RooRealVar_dict[f'ratio_comp{comp}'].getValV()
				else:
					y_val *= 1-sum([self.RooRealVar_dict[f'ratio_comp{r}'].getValV() for r in np.arange(0,self.n_comps-1,1)])

				i_ratio += 1
				# Plot
				plt.plot(xm, y_val, label=f"Neg. component {comp}:{dim})", c='grey', ls="--", zorder=2, linewidth=1.75)

			for dim in range(0,self.dimensions[1],1):
				y_val = self.pos_func(x=xm,mu=self.RooRealVar_dict[f'mu{comp}'].getValV(), 
											sigma=self.RooRealVar_dict[f'sigma'].getValV(), 
											ptau=self.RooRealVar_dict[f'ptau{dim}'].getValV())
				# normalize accounting for fit ratios
				integral_cut = sum(y_val) * np.diff(xm)[0]
				y_val = y_val / integral_cut * sum(n_cut) * dx[0]
				if (self.n_comps == 0 and (self.dimensions[0]+self.dimensions[1]) == 1):
					ratio = 1
				else:
					if i_ratio != n_ratios:
						y_val *= self.RooRealVar_dict[f'ratio{i_ratio}'].getValV()
					else:
						y_val *= 1-sum([self.RooRealVar_dict[f'ratio{r}'].getValV() for r in np.arange(0,n_ratios,1)])
				#
				if comp < (self.n_comps-1):
					y_val *= self.RooRealVar_dict[f'ratio_comp{comp}'].getValV()
				else:
					y_val *= 1-sum([self.RooRealVar_dict[f'ratio_comp{r}'].getValV() for r in np.arange(0,self.n_comps-1,1)])
				#
				i_ratio += 1
				# Plot
				plt.plot(xm, y_val, label=f"Pos. component {comp}:{dim})", c='grey', ls="--", zorder=2, linewidth=1.75)

		# Get y axis limits
		ylims = plt.ylim()
		plt.legend()
		if log:
			plt.yscale("log")
			plt.ylim(0.9,2*ylims[1])
		# Zoom in on found peaks
		if self.peaks.n_peaks != 0:
			plt.xlim(self.peaks.earliest_left_base-200, self.peaks.latest_right_base+200)
		else:
			plt.xlim(self.mu.getValV()-200, self.mu.getValV()+200)

		if save:
			print(f"Plot fit save as {self.base_file_name}.pdf")
			plt.savefig(self.base_file_name+"_fit.pdf", dpi=300)

		if silent == False:
			plt.show()

		# Clear canvas to avoid printing on top of other plot in batch mode
		plt.clf()

class conv_Emg(FitMethods):
	"""
	Class handling the exponentially modified Gaussian as defined in 
	https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
	"""
	def __init__(self, lst_file, peaks = 'findit', verbose=False):
		"""
		Initializes the class and number of parameters.
		"""        
		# QtWidgets.QWidget.__init__(self)
		FitMethods.__init__(self, lst_file, 'Egh', 2, peaks=peaks, verbose=verbose)
		# self.setupUi(self)
		# self.update_par_list()
		self.hist = self.lst2roothist(self.lst_file, bins=2000)

	def init_limits(self):
		tau_estimate = self.hist.GetStdDev(1)*(self.hist.GetSkewness(1) / 2)**(1/3)
		llambda_estimate = 1/tau_estimate
		print(llambda_estimate)
		self.limits = {
			'mu':[self.hist.GetMean(1), self.xmin, self.xmax],
			'llambda':[llambda_estimate, 0.001, 10*llambda_estimate],
			'sigma':[self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1)],
		}
	
	def call_pdf(self, xmin, xmax, limits=False, index=0):
		"""
		Setup Variable, Parameters and PDF and performs fitting of the PDF to ROI data.
		:param roi_hist: histogram to fit
		:param xmin: range min for the fit
		:param xmax: range max for the fit
		:param index: peak index for display
		:return list: containing all relevant information about the ROI histogram, PDF, parameters,
					  PDF components and fit results
		"""
		# Test dimension of limits
		self.xmin = xmin
		self.xmax = xmax
		if isinstance(xmin, list):
			print("(FIT_LUKAS::Emg::call_pdf): Multiple fit-regions not supported in this fit function")
			return 0
		#
		self.compose_title_unit(self.hist.GetXaxis().GetTitle())
		#
		if not limits:
			self.init_limits()
		else:
			self.limits = limits
		#
		# self.set_value(0, hist.GetMean(1))
		# self.set_limits(0, xmin, xmax)
		# self.set_value(1, hist.GetStdDev(1))
		# self.set_limits(1, 0.01*hist.GetStdDev(1), 10*hist.GetStdDev(1))
		# self.set_value(2, tau_estimate)
		# self.set_limits(2, -10*tau_estimate, 10*tau_estimate)
		#
		self.x = RooRealVar("x", self.rootitle, self.xmin, self.xmax, self.roounit)
		self.mu = RooRealVar("mu{}".format(index), "mu{}".format(index),
						  self.limits["mu"][0], self.limits["mu"][1], self.limits["mu"][2])
		self.sigma = RooRealVar("sigma{}".format(index), "sigma{}".format(index),
						  self.limits["sigma"][0], self.limits["sigma"][1], self.limits["sigma"][2])
		self.tau0 = RooRealVar("tau{}0".format(index), "tau{}0".format(index),
						  self.limits["tau0"][0], self.limits["tau0"][1], self.limits["tau0"][2])
		self.tau1 = RooRealVar("tau{}1".format(index), "tau{}1".format(index),
						  self.limits["tau1"][0], self.limits["tau1"][1], self.limits["tau1"][2])
		# self.tau2 = RooRealVar("tau{}2".format(index), "tau{}2".format(index),
		# 				  self.limits["tau2"][0], self.limits["tau2"][1], self.limits["tau2"][2])
		self.ratio_exp0 = RooRealVar(f'ratio_exp0{index}', "ratio_exp0", 
							self.limits["ratio_exp0"][0], self.limits["ratio_exp0"][1], self.limits["ratio_exp0"][2])
		# self.ratio_exp1 = RooRealVar(f'ratio_exp1{index}', "ratio_exp1", 
		# 					self.limits["ratio_exp1"][0], self.limits["ratio_exp1"][1], self.limits["ratio_exp1"][2])
		gm = RooGaussian(f'g{index}', f'Gaussian Model {index}', self.x, self.mu, self.sigma)
		basis_func0 = RooGenericPdf("basis0", "basis0","(1/@1)*exp(-(@0/@1))",RooArgList(self.x,self.tau0))
		basis_func1 = RooGenericPdf("basis1", "basis1","(1/@1)*exp(-(@0/@1))",RooArgList(self.x,self.tau1))
		# basis_func2 = RooGenericPdf("basis2", "basis2","(1/@1)*exp(-(@0/@1))",RooArgList(self.x,self.tau2))
		sum_pdf = RooAddPdf(f'exponetial_sum{index}', f'exponetial_sum {index}',
										  RooArgList(basis_func0, basis_func1), RooArgList(self.ratio_exp0))
										  # RooArgList(basis_func0, basis_func1, basis_func2), RooArgList(self.ratio_exp0, self.ratio_exp1))
		self.x.setBins(10000, "cache")


		self.this_pdf = ROOT.RooFFTConvPdf("conv_Emg2", "Exp2 (X) gauss", self.x, gm, sum_pdf)


		self.roodefs = [self.this_pdf, self.x, self.mu, self.sigma, self.tau0, self.tau1, self.ratio_exp0]
		# self.roodefs = [self.this_pdf, self.x, self.mu, self.sigma, self.tau0, self.tau1, self.tau2, self.ratio_exp0, self.ratio_exp1]
		# self.fix_var(self.roodefs[2:], index)
		self.this_roohist, self.fit_results = self.minimize(self.hist, self.xmin, self.xmax)

		# THIS LINE HAS TO BE HERE FOR SOME UNKNOW REASON TO ME TO AVOID A SEG FAULT IN THE PLOT FUNCTION GOD KNOWS WHY I REALLY DON'T KNOW ANYMORE...
		print(self.this_pdf.getVal(ROOT.RooArgSet(self.x, self.mu, self.sigma, self.tau0, self.tau1, self.ratio_exp0)))
		print(self.this_pdf.getVal(ROOT.RooArgSet(self.x, self.mu, self.sigma, self.tau0, self.tau1, self.ratio_exp0)))
		# print(self.this_pdf.getVal(ROOT.RooArgSet(self.x, self.mu, self.sigma, self.tau0, self.tau1, self.tau2, self.ratio_exp0, self.ratio_exp1)))


		# frame = self.x.frame(RooFit.Title("landau (x) gauss convolution"))
		# cloned_pdf = ROOT.TF1(self.this_pdf.asTF(RooArgSet(self.x)))
		# self.hist.Draw()
		# self.this_pdf.plotOn(frame)
		# c = ROOT.TCanvas("rf208_convolution", "rf208_convolution", 600, 600)
		# frame.Draw()
		# c.Draw()

		# time.sleep(30)

		return [self.this_roohist, self.this_pdf, self.nb_pars,
				self.x, self.mu, self.sigma, self.tau0, self.tau1, self.ratio_exp0, self.pstate,
				self.my_name, self.spl_draw, self.rooCompon, self.fit_results]
		# return [self.this_roohist, self.this_pdf, self.nb_pars,
		# 		self.x, self.mu, self.sigma, self.tau0, self.tau1, self.tau2, self.ratio_exp0, self.ratio_exp1, self.pstate,
		# 		self.my_name, self.spl_draw, self.rooCompon, self.fit_results]

class conv_doubleEmg(FitMethods):
	"""
	Class handling the exponentially modified Gaussian as defined in 
	https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
	"""
	def __init__(self, lst_file, peaks = 'findit', verbose=False):
		"""
		Initializes the class and number of parameters.
		"""        
		# QtWidgets.QWidget.__init__(self)
		FitMethods.__init__(self, lst_file, 'Egh', 2, peaks=peaks, verbose=verbose)
		# self.setupUi(self)
		# self.update_par_list()
		self.hist = self.lst2roothist(self.lst_file, bins=2000)

	def init_limits(self):
		tau_estimate = self.hist.GetStdDev(1)*(self.hist.GetSkewness(1) / 2)**(1/3)
		llambda_estimate = 1/tau_estimate
		print(llambda_estimate)
		self.limits = {
			'mu':[self.hist.GetMean(1), self.xmin, self.xmax],
			'llambda':[llambda_estimate, 0.001, 10*llambda_estimate],
			'sigma':[self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1)],
		}
	
	def call_pdf(self, xmin, xmax, limits=False, index=0):
		"""
		Setup Variable, Parameters and PDF and performs fitting of the PDF to ROI data.
		:param roi_hist: histogram to fit
		:param xmin: range min for the fit
		:param xmax: range max for the fit
		:param index: peak index for display
		:return list: containing all relevant information about the ROI histogram, PDF, parameters,
					  PDF components and fit results
		"""
		# Test dimension of limits
		self.xmin = xmin
		self.xmax = xmax
		if isinstance(xmin, list):
			print("(FIT_LUKAS::Emg::call_pdf): Multiple fit-regions not supported in this fit function")
			return 0
		#
		self.compose_title_unit(self.hist.GetXaxis().GetTitle())
		#
		if not limits:
			self.init_limits()
		else:
			self.limits = limits
		#
		# self.set_value(0, hist.GetMean(1))
		# self.set_limits(0, xmin, xmax)
		# self.set_value(1, hist.GetStdDev(1))
		# self.set_limits(1, 0.01*hist.GetStdDev(1), 10*hist.GetStdDev(1))
		# self.set_value(2, tau_estimate)
		# self.set_limits(2, -10*tau_estimate, 10*tau_estimate)
		###
		self.x = RooRealVar("x", self.rootitle, self.xmin, self.xmax, self.roounit)
		###
		self.mu0 = RooRealVar("mu{}0".format(index), "mu{}0".format(index),
						  self.limits["mu0"][0], self.limits["mu0"][1], self.limits["mu0"][2])
		self.mu1 = RooRealVar("mu{}1".format(index), "mu{}1".format(index),
						  self.limits["mu1"][0], self.limits["mu1"][1], self.limits["mu1"][2])
		self.sigma = RooRealVar("sigma{}".format(index), "sigma{}".format(index),
						  self.limits["sigma"][0], self.limits["sigma"][1], self.limits["sigma"][2])
		self.ratio_gauss = RooRealVar(f'ratio{index}_gauss', "ratio_gauss", 
							self.limits["ratio_gauss"][0], self.limits["ratio_gauss"][1], self.limits["ratio_gauss"][2])
		gm0 = RooGaussian(f'g{index}0', f'Gaussian Model {index}0', self.x, self.mu0, self.sigma)
		gm1 = RooGaussian(f'g{index}1', f'Gaussian Model {index}1', self.x, self.mu1, self.sigma)
		gm_sum = RooAddPdf(f'gaussian_sum{index}', f'gaussian_sum {index}',
										  RooArgList(gm0, gm1), RooArgList(self.ratio_gauss))
		###
		self.tau0 = RooRealVar("tau{}0".format(index), "tau{}0".format(index),
						  self.limits["tau0"][0], self.limits["tau0"][1], self.limits["tau0"][2])
		self.tau1 = RooRealVar("tau{}1".format(index), "tau{}1".format(index),
						  self.limits["tau1"][0], self.limits["tau1"][1], self.limits["tau1"][2])
		self.ratio_exp = RooRealVar(f'ratio{index}_exp', "ratio_exp", 
							self.limits["ratio_exp"][0], self.limits["ratio_exp"][1], self.limits["ratio_exp"][2])
		basis_func0 = RooGenericPdf("basis0", "basis0","(1/@1)*exp(-(@0/@1))",RooArgList(self.x,self.tau0))
		basis_func1 = RooGenericPdf("basis1", "basis1","(1/@1)*exp(-(@0/@1))",RooArgList(self.x,self.tau1))
		exp_sum = RooAddPdf(f'exp_sum{index}', f'exp_sum {index}',
										  RooArgList(basis_func0, basis_func1), RooArgList(self.ratio_exp))
		###
		self.x.setBins(10000, "cache")
		self.this_pdf = ROOT.RooFFTConvPdf("conv_doubleEmg2", "Exp2 (X) gauss2", self.x, gm_sum, exp_sum)
		self.roodefs = [self.this_pdf, self.x, self.mu0, self.mu1, self.sigma, self.tau0, self.tau1, self.ratio_exp, self.ratio_gauss]
		# self.fix_var(self.roodefs[2:], index)
		self.this_roohist, self.fit_results = self.minimize(self.hist, self.xmin, self.xmax)

		# THIS LINE HAS TO BE HERE FOR SOME UNKNOW REASON TO ME TO AVOID A SEG FAULT IN THE PLOT FUNCTION GOD KNOWS WHY I REALLY DON'T KNOW ANYMORE...
		print(self.this_pdf.getVal(ROOT.RooArgSet(self.x, self.mu0, self.mu1, self.sigma, self.tau0, self.tau1, self.ratio_exp, self.ratio_gauss)))


		# frame = self.x.frame(RooFit.Title("landau (x) gauss convolution"))
		# cloned_pdf = ROOT.TF1(self.this_pdf.asTF(RooArgSet(self.x)))
		# self.hist.Draw()
		# self.this_pdf.plotOn(frame)
		# c = ROOT.TCanvas("rf208_convolution", "rf208_convolution", 600, 600)
		# frame.Draw()
		# c.Draw()

		# time.sleep(30)

		return [self.this_roohist, self.this_pdf, self.nb_pars,
				self.x, self.mu0, self.mu1, self.sigma, self.tau0, self.tau1, self.ratio_exp, self.ratio_gauss, self.pstate,
				self.my_name, self.spl_draw, self.rooCompon, self.fit_results]

	def plot(self, bins = 2000, log=False):
		"""  """
		if len(self.lst_file) == 0:
			print("Fit not excecuted yet or failed.")
			return 0 
		#            
		xdata = self.lst_file.tof
		n, xe = np.histogram(xdata, bins=bins)
		# n = n/len(xdata)
		cx = 0.5 * (xe[1:] + xe[:-1])
		dx = np.diff(xe)
		# Plot data
		# plt.scatter(cx, n, zorder=1)
		plt.errorbar(cx, n, n ** 0.5, fmt="ok", zorder=1)
		xm = np.linspace(self.xmin, self.xmax, num=5000)
		# Plot 'carpet'
		if log:
			plt.plot(xdata, np.zeros_like(xdata)+0.9, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		plt.plot(xdata, np.zeros_like(xdata)-5, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		# if log:
		# 	xm = np.linspace(self.mu.getValV()-100, self.mu.getValV()+100, num=1000)
		#
		# Clone pdf, otherwise access results in seg fault
		# cloned_pdf = self.this_pdf.clone("Test")
		y_val = []
		for i in xm:
			self.x.setVal(i)
			# x_copy = self.x
			# print(self.this_pdf.getVal(ROOT.RooArgSet(x_copy)))
			y_val.append(self.this_pdf.getVal(ROOT.RooArgSet(self.x, self.mu0, self.mu1, self.sigma, self.tau0, self.tau1, self.ratio_exp, self.ratio_gauss)))

		# Normalize values
		integral_cut = sum(y_val) * np.diff(xm)[0]
		left_n_cut = len(xe[xe<self.xmin])
		right_n_cut = len(xe[xe<self.xmax])
		n_cut = n[left_n_cut:right_n_cut]        
		y_val = y_val / integral_cut * sum(n_cut) * dx[0]

		
		# y_val = [self.fit_function(x=y, mean=self.mean.getValV(), sigma=self.sigma.getValV(), tau=self.tau.getValV())
		#                           * len(xdata) * dx[0] 
		#                           for y in x_range_cut]
		plt.plot(xm, y_val, ls=":", label=f"init", zorder=2)
		plt.plot(xm, y_val, label=f"fit", c='r', zorder=3, linewidth=3)
		
		# Get y axis limits
		ylims = plt.ylim()
		
		plt.legend();
		if log:
			plt.yscale("log")
			plt.ylim(0.9,2*ylims[1])
		# Zoom in on found peaks
		if self.peaks.n_peaks != 0:
			plt.xlim(self.peaks.earliest_left_base-200, self.peaks.latest_right_base+200)
		else:
			plt.xlim(self.mu.getValV()-200, self.mu.getValV()+200)

		plt.show()
# 
class conv_doubleEmg_2comp(FitMethods):
	"""
	Class handling the exponentially modified Gaussian as defined in 
	https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
	"""
	def __init__(self, lst_file, peaks = 'findit', verbose=False):
		"""
		Initializes the class and number of parameters.
		"""        
		# QtWidgets.QWidget.__init__(self)
		FitMethods.__init__(self, lst_file, 'Egh', 2, peaks=peaks, verbose=verbose)
		# self.setupUi(self)
		# self.update_par_list()
		self.hist = self.lst2roothist(self.lst_file, bins=2000)

	def init_limits(self):
		tau_estimate = self.hist.GetStdDev(1)*(self.hist.GetSkewness(1) / 2)**(1/3)
		llambda_estimate = 1/tau_estimate
		print(llambda_estimate)
		self.limits = {
			'mu':[self.hist.GetMean(1), self.xmin, self.xmax],
			'llambda':[llambda_estimate, 0.001, 10*llambda_estimate],
			'sigma':[self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1)],
		}
	
	def call_pdf(self, xmin, xmax, limits=False, index=0):
		"""
		Setup Variable, Parameters and PDF and performs fitting of the PDF to ROI data.
		:param roi_hist: histogram to fit
		:param xmin: range min for the fit
		:param xmax: range max for the fit
		:param index: peak index for display
		:return list: containing all relevant information about the ROI histogram, PDF, parameters,
					  PDF components and fit results
		"""
		# Test dimension of limits
		self.xmin = xmin
		self.xmax = xmax
		if isinstance(xmin, list):
			print("(FIT_LUKAS::Emg::call_pdf): Multiple fit-regions not supported in this fit function")
			return 0
		#
		self.compose_title_unit(self.hist.GetXaxis().GetTitle())
		#
		if not limits:
			self.init_limits()
		else:
			self.limits = limits
		#
		# self.set_value(0, hist.GetMean(1))
		# self.set_limits(0, xmin, xmax)
		# self.set_value(1, hist.GetStdDev(1))
		# self.set_limits(1, 0.01*hist.GetStdDev(1), 10*hist.GetStdDev(1))
		# self.set_value(2, tau_estimate)
		# self.set_limits(2, -10*tau_estimate, 10*tau_estimate)
		###
		self.x = RooRealVar("x", self.rootitle, self.xmin, self.xmax, self.roounit)
		###
		self.mu0 = RooRealVar("mu{}0".format(index), "mu{}0".format(index),
						  self.limits["mu0"][0], self.limits["mu0"][1], self.limits["mu0"][2])
		self.mu1 = RooRealVar("mu{}1".format(index), "mu{}1".format(index),
						  self.limits["mu1"][0], self.limits["mu1"][1], self.limits["mu1"][2])
		self.sigma = RooRealVar("sigma{}".format(index), "sigma{}".format(index),
						  self.limits["sigma"][0], self.limits["sigma"][1], self.limits["sigma"][2])
		gm0 = RooGaussian(f'g{index}0', f'Gaussian Model {index}0', self.x, self.mu0, self.sigma)
		gm1 = RooGaussian(f'g{index}1', f'Gaussian Model {index}1', self.x, self.mu1, self.sigma)
		print("Test-Gauss")
		###
		self.tau0 = RooRealVar("tau{}0".format(index), "tau{}0".format(index),
						  self.limits["tau0"][0], self.limits["tau0"][1], self.limits["tau0"][2])
		self.tau1 = RooRealVar("tau{}1".format(index), "tau{}1".format(index),
						  self.limits["tau1"][0], self.limits["tau1"][1], self.limits["tau1"][2])
		self.ratio_exp = RooRealVar(f'ratio{index}_exp', "ratio_exp", 
							self.limits["ratio_exp"][0], self.limits["ratio_exp"][1], self.limits["ratio_exp"][2])
		basis_func0 = RooGenericPdf("basis0", "basis0","(1/@1)*exp(-(@0/@1))",RooArgList(self.x,self.tau0))
		basis_func1 = RooGenericPdf("basis1", "basis1","(1/@1)*exp(-(@0/@1))",RooArgList(self.x,self.tau1))
		exp_sum= RooAddPdf(f'exp_sum{index}', f'exp_sum {index}',
										  RooArgList(basis_func0, basis_func1), RooArgList(self.ratio_exp))
		print("Test-Exp")
		###
		self.x.setBins(10000, "cache")
		self.comp0 = ROOT.RooFFTConvPdf("conv_doubleEmg_comp0", "Exp2 (X) gauss2 comp0", self.x, gm0, exp_sum)
		self.comp1 = ROOT.RooFFTConvPdf("conv_doubleEmg_comp1", "Exp2 (X) gauss2 comp1", self.x, gm1, exp_sum)
		print("Test-FFT")
		###
		self.ratio_comp = RooRealVar(f'ratio{index}_comp', "ratio_comp", 
							self.limits["ratio_comp"][0], self.limits["ratio_comp"][1], self.limits["ratio_comp"][2])
		self.this_pdf = RooAddPdf(f'exp_sum{index}', f'exp_sum {index}',
										  RooArgList(self.comp0, self.comp1), RooArgList(self.ratio_comp))
		print("Test-PDF")
		###


		###
		self.roodefs = [self.this_pdf, self.x, self.mu0, self.mu1, self.sigma, self.tau0, self.tau1, self.ratio_exp, self.ratio_comp]
		# self.fix_var(self.roodefs[2:], index)
		self.this_roohist, self.fit_results = self.minimize(self.hist, self.xmin, self.xmax)

		# THIS LINE HAS TO BE HERE FOR SOME UNKNOW REASON TO ME TO AVOID A SEG FAULT IN THE PLOT FUNCTION GOD KNOWS WHY I REALLY DON'T KNOW ANYMORE...
		print(self.this_pdf.getVal(ROOT.RooArgSet(self.x, self.mu0, self.mu1, self.sigma, self.tau0, self.tau1, self.ratio_exp, self.ratio_comp)))


		# frame = self.x.frame(RooFit.Title("landau (x) gauss convolution"))
		# cloned_pdf = ROOT.TF1(self.this_pdf.asTF(RooArgSet(self.x)))
		# self.hist.Draw()
		# self.this_pdf.plotOn(frame)
		# c = ROOT.TCanvas("rf208_convolution", "rf208_convolution", 600, 600)
		# frame.Draw()
		# c.Draw()

		# time.sleep(30)

		return [self.this_roohist, self.this_pdf, self.nb_pars,
				self.x, self.mu0, self.mu1, self.sigma, self.tau0, self.tau1, self.ratio_exp, self.ratio_comp, self.pstate,
				self.my_name, self.spl_draw, self.rooCompon, self.fit_results]

	def plot(self, bins = 2000, log=False):
		"""  """
		if len(self.lst_file) == 0:
			print("Fit not excecuted yet or failed.")
			return 0 
		#            
		xdata = self.lst_file.tof
		n, xe = np.histogram(xdata, bins=bins)
		# n = n/len(xdata)
		cx = 0.5 * (xe[1:] + xe[:-1])
		dx = np.diff(xe)
		# Plot data
		# plt.scatter(cx, n, zorder=1)
		plt.errorbar(cx, n, n ** 0.5, fmt="ok", zorder=1)
		xm = np.linspace(self.xmin, self.xmax, num=5000)
		# Plot 'carpet'
		if log:
			plt.plot(xdata, np.zeros_like(xdata)+0.9, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		plt.plot(xdata, np.zeros_like(xdata)-5, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		# if log:
		# 	xm = np.linspace(self.mu.getValV()-100, self.mu.getValV()+100, num=1000)
		#
		# Clone pdf, otherwise access results in seg fault
		# cloned_pdf = self.this_pdf.clone("Test")
		y_val = []
		for i in xm:
			self.x.setVal(i)
			# x_copy = self.x
			# print(self.this_pdf.getVal(ROOT.RooArgSet(x_copy)))
			y_val.append(self.this_pdf.getVal(ROOT.RooArgSet(self.x, self.mu0, self.mu1, self.sigma, self.tau0, self.tau1, self.ratio_exp, self.ratio_comp)))

		# Normalize values
		integral_cut = sum(y_val) * np.diff(xm)[0]
		left_n_cut = len(xe[xe<self.xmin])
		right_n_cut = len(xe[xe<self.xmax])
		n_cut = n[left_n_cut:right_n_cut]        
		y_val = y_val / integral_cut * sum(n_cut) * dx[0]

		
		# y_val = [self.fit_function(x=y, mean=self.mean.getValV(), sigma=self.sigma.getValV(), tau=self.tau.getValV())
		#                           * len(xdata) * dx[0] 
		#                           for y in x_range_cut]
		plt.plot(xm, y_val, ls=":", label=f"init", zorder=2)
		plt.plot(xm, y_val, label=f"fit", c='r', zorder=3, linewidth=3)
		
		# Get y axis limits
		ylims = plt.ylim()
		
		plt.legend();
		if log:
			plt.yscale("log")
			plt.ylim(0.9,2*ylims[1])
		# Zoom in on found peaks
		if self.peaks.n_peaks != 0:
			plt.xlim(self.peaks.earliest_left_base-200, self.peaks.latest_right_base+200)
		else:
			plt.xlim(self.mu.getValV()-200, self.mu.getValV()+200)

		plt.show()
# 
class doubleEmg(FitMethods):
	"""
	Class handling the exponentially modified Gaussian as defined in 
	https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
	"""
	def __init__(self, lst_file, peaks = 'findit', verbose=False):
		"""
		Initializes the class and number of parameters.
		"""        
		# QtWidgets.QWidget.__init__(self)
		FitMethods.__init__(self, lst_file, 'Egh', 2, peaks=peaks, verbose=verbose)
		# self.setupUi(self)
		# self.update_par_list()
		self.hist = self.lst2roothist(self.lst_file, bins=2000)

	def init_limits(self):
		tau_estimate = self.hist.GetStdDev(1)*(self.hist.GetSkewness(1) / 2)**(1/3)
		llambda_estimate = 1/tau_estimate
		self.limits = {
			'mu0':[self.hist.GetMean(1), self.xmin, self.xmax],
			'mu1':[self.hist.GetMean(1), self.xmin, self.xmax],
			'llambda':[llambda_estimate, 0.0001, 10*llambda_estimate],
			'sigma':[self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1)],
			'ratio':[0.5,0,1],
		}

	
	def fit_function(self,x,mu,sigma,llambda):
		""" """
		return ( (llambda/2) * np.exp( (llambda/2) * ( 2*mu + llambda*np.power(sigma, 2)-2*x ) ) *
				 (1 - sc.special.erf( (mu+llambda*np.power(sigma, 2)-x)/(math.sqrt(2)*sigma) ) )
			   )
	
	def call_pdf(self, xmin, xmax, limits=False, index=0):
		"""
		Setup Variable, Parameters and PDF and performs fitting of the PDF to ROI data.
		:param roi_hist: histogram to fit
		:param xmin: range min for the fit
		:param xmax: range max for the fit
		:param index: peak index for display
		:return list: containing all relevant information about the ROI histogram, PDF, parameters,
					  PDF components and fit results
		"""
		self.xmin = xmin
		self.xmax = xmax
		#
		if isinstance(xmin, list):
			_xmin = xmin[0]
			_xmax = xmax[len(xmax)-1]
		else:
			_xmin = xmin
			_xmax = xmax
		self.hist.GetXaxis().SetRangeUser(_xmin, _xmax)
		self.compose_title_unit(self.hist.GetXaxis().GetTitle())
		#
		if not limits:
			self.init_limits()
		else:
			self.limits = limits
		#
		# self.set_value(0, hist.GetMean(1))
		# self.set_limits(0, xmin, xmax)
		# self.set_value(1, hist.GetStdDev(1))
		# self.set_limits(1, 0.01*hist.GetStdDev(1), 10*hist.GetStdDev(1))
		# self.set_value(2, tau_estimate)
		# self.set_limits(2, -10*tau_estimate, 10*tau_estimate)
		#
		self.x = RooRealVar("x", self.rootitle, _xmin, _xmax, self.roounit)
		#
		self.mu0 = RooRealVar("mu{}0".format(index), "mu0{}".format(index),
						  self.limits["mu0"][0], self.limits["mu0"][1], self.limits["mu0"][2])
		self.mu1 = RooRealVar("mu{}1".format(index), "mu{}1".format(index),
						  self.limits["mu1"][0], self.limits["mu1"][1], self.limits["mu1"][2])
		self.sigma = RooRealVar("sigma{}".format(index), "sigma{}".format(index),
						  self.limits["sigma"][0], self.limits["sigma"][1], self.limits["sigma"][2])
		self.llambda = RooRealVar("llambda{}".format(index), "llambda{}".format(index),
						  self.limits["llambda"][0], self.limits["llambda"][1], self.limits["llambda"][2])
		self.ratio = RooRealVar(f'ratio{index}0', "ratio0", 
						  self.limits["ratio"][0], self.limits["ratio"][1], self.limits["ratio"][2])
		sig0 = RooGenericPdf('emg_{}0'.format(index),
									  'Emg model {}0'.format(index),
									  '@3/2*exp(@3/2*(2*@1+@3*@2^2-2*x))*(1-TMath::Erf((@1+@3*@2^2-x)/(sqrt(2)*@2)))',
									  RooArgList(self.x, self.mu0, self.sigma, self.llambda))
		sig1 = RooGenericPdf('emg_{}1'.format(index),
									  'Emg model {}1'.format(index),
									  '@3/2*exp(@3/2*(2*@1+@3*@2^2-2*x))*(1-TMath::Erf((@1+@3*@2^2-x)/(sqrt(2)*@2)))',
									  RooArgList(self.x, self.mu1, self.sigma, self.llambda))
		self.this_pdf = RooAddPdf(f'double_Emg{index}', f'Two EMGs added {index}',
								  RooArgList(sig0, sig1), RooArgList(self.ratio))
		self.roodefs = [self.this_pdf, self.x, self.mu0, self.mu1, self.sigma, self.llambda]
		# self.fix_var(self.roodefs[2:], index)
		self.this_roohist, self.fit_results = self.minimize(self.hist, self.xmin, self.xmax)
		return [self.this_roohist, self.this_pdf, self.nb_pars,
				self.x, self.mu0, self.mu1, self.sigma, self.llambda, self.pstate,
				self.my_name, self.spl_draw, self.rooCompon, self.fit_results]
	
	def plot(self, bins = 2000, log=False):
		"""  """
		if len(self.lst_file) == 0:
			print("Fit not excecuted yet or failed.")
			return 0 
		#            
		xdata = self.lst_file.tof
		n, xe = np.histogram(xdata, bins=bins)
		# n = n/len(xdata)
		cx = 0.5 * (xe[1:] + xe[:-1])
		dx = np.diff(xe)
		# Plot data
		# plt.scatter(cx, n, zorder=1)
		plt.errorbar(cx, n, n ** 0.5, fmt="ok", zorder=1)
		xm = np.linspace(xe[0], xe[-1], num=5000)
		# Plot 'carpet'
		if log:
			plt.plot(xdata, np.zeros_like(xdata)+0.9, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		plt.plot(xdata, np.zeros_like(xdata)-5, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		# if log:
		# 	xm = np.linspace(self.mu.getValV()-100, self.mu.getValV()+100, num=1000)
		#
		# Plot fits, loop through all fits
		#
		# Get the strength of the fitted function by doing a simple chi2 fit        
		# def fit_func(A, x):
		#     """ """
		#     return A * self.fit_function(x=x, mean=self.mean.getValV(), sigma=self.sigma.getValV(), tau=self.tau.getValV())
		# def fun(A):
		#     return fit_func(A, x_range_cut) - data_cut
		# A0 = 1
		# res = least_squares(fun, A0)
		# Calculate fit 
		# y_val = [self.fit_function(x=y, mean=self.mean.getValV(), sigma=self.sigma.getValV(), tau=self.tau.getValV()) * res.x[0]
		#          for y in x_range_cut]
		
		# Normalize by integration
		y_val = [(self.ratio.getValV() * self.fit_function(x=y, mu=self.mu0.getValV(), sigma=self.sigma.getValV(), llambda=self.llambda.getValV())) +
				 ((1 - self.ratio.getValV())*self.fit_function(x=y, mu=self.mu1.getValV(), sigma=self.sigma.getValV(), llambda=self.llambda.getValV()))
								  # * len(xdata) * dx[0] 
								  for y in xm]

		# Normalize values
		integral_cut = sum(y_val) * np.diff(xm)[0]
		# If several regions are defined
		if isinstance(self.xmin, list):
			left_n_cut = len(xe[xe<self.xmin[0]])
			right_n_cut = len(xe[xe<self.xmax[len(self.xmax)-1]])
		else:
			left_n_cut = len(xe[xe<self.xmin])
			right_n_cut = len(xe[xe<self.xmax])
		n_cut = n[left_n_cut:right_n_cut]        
		y_val = y_val / integral_cut * sum(n_cut) * dx[0]

		
		# y_val = [self.fit_function(x=y, mean=self.mean.getValV(), sigma=self.sigma.getValV(), tau=self.tau.getValV())
		#                           * len(xdata) * dx[0] 
		#                           for y in x_range_cut]
		plt.plot(xm, y_val, ls=":", label=f"init", zorder=2)
		plt.plot(xm, y_val, label=f"fit", c='r', zorder=3, linewidth=3)
		
		# Get y axis limits
		ylims = plt.ylim()
		
		plt.legend();
		if log:
			plt.yscale("log")
			plt.ylim(0.9,2*ylims[1])
		# Zoom in on found peaks
		if self.peaks.n_peaks != 0:
			plt.xlim(self.peaks.earliest_left_base-200, self.peaks.latest_right_base+200)
		else:
			plt.xlim(self.mu0.getValV()-200, self.mu1.getValV()+200)

		plt.show()

class doubleEgh(FitMethods):
	"""
	Class handling the Exponent-Gaussian hybrid model.
	See article -->
	"""
	def __init__(self, lst_file, peaks = 'findit', verbose=False):
		"""
		Initializes the class and number of parameters.
		"""        
		# QtWidgets.QWidget.__init__(self)
		FitMethods.__init__(self, lst_file, 'doubleEgh', 2, peaks=peaks, verbose=verbose)
		# self.setupUi(self)
		# self.update_par_list()
		self.hist = self.lst2roothist(self.lst_file, bins=2000)

	def init_limits(self):
		tau_estimate = self.hist.GetStdDev(1)*(self.hist.GetSkewness(1) / 2)**(1/3)
		self.limits = {
			'mean0':[self.hist.GetMean(1), self.xmin, self.xmax],
			'mean1':[self.hist.GetMean(1), self.xmin, self.xmax],
			'tau':[tau_estimate,-10*tau_estimate, 10*tau_estimate],
			'sigma':[self.hist.GetStdDev(1), 0.01*self.hist.GetStdDev(1), 10*self.hist.GetStdDev(1)],
			'ratio':[0.5,0,1],
		}
	
	def fit_function(self,x,mean,sigma,tau):
		""" """
		return (np.exp( - ( x-mean )**2 / ( 2*sigma**2 + tau*(x-mean) ) )  )
	
	def calc_pole(self,pmean,psigma,ptau):
		""" """
		return ( -(2*psigma**2)/ptau + pmean )

	def call_pdf(self, xmin, xmax, limits = False, index=0):
		"""
		Setup Variable, Parameters and PDF and performs fitting of the PDF to ROI data.
		:param roi_hist: histogram to fit
		:param xmin: range min for the fit
		:param xmax: range max for the fit
		:param index: peak index for display
		:return list: containing all relevant information about the ROI histogram, PDF, parameters,
					  PDF components and fit results
		"""
		self.xmin = xmin
		self.xmax = xmax
		self.hist.GetXaxis().SetRangeUser(self.xmin, self.xmax)
		self.compose_title_unit(self.hist.GetXaxis().GetTitle())
		# 
		if not limits:
			self.init_limits()
		else:
			self.limits = limits
		#
		# self.set_value(0, hist.GetMean(1))
		# self.set_limits(0, xmin, xmax)
		# self.set_value(1, hist.GetStdDev(1))
		# self.set_limits(1, 0.01*hist.GetStdDev(1), 10*hist.GetStdDev(1))
		# self.set_value(2, tau_estimate)
		# self.set_limits(2, -10*tau_estimate, 10*tau_estimate)
		#
		# Definitions for the variables
		self.x = RooRealVar("x", self.rootitle, self.xmin, self.xmax, self.roounit)
		self.mean0 = RooRealVar("mean{}0".format(index), "mean{}0".format(index),
						  self.limits["mean0"][0], self.limits["mean0"][1], self.limits["mean0"][2])
		self.mean1 = RooRealVar("mean{}1".format(index), "mean{}1".format(index),
						  self.limits["mean1"][0], self.limits["mean1"][1], self.limits["mean1"][2])
		self.sigma = RooRealVar("sigma{}".format(index), "sigma{}".format(index),
						   self.limits["sigma"][0], self.limits["sigma"][1], self.limits["sigma"][2])
		self.tau = RooRealVar("tau{}".format(index), "tau{}".format(index),
						 self.limits["tau"][0], self.limits["tau"][1], self.limits["tau"][2])
		self.ratio = RooRealVar(f'ratio{index}0', "ratio0", 0.5, 0., 1., " ")
		# Definition for the components
		sig0 = RooGenericPdf('egh_{}0'.format(index),
									  'EGH model {}0'.format(index),
									  'exp(((-(@0-@1)^2)/(2*@2^2+@3*(@0-@1))))',
									  RooArgList(self.x, self.mean0, self.sigma, self.tau))
		sig1 = RooGenericPdf('egh_{}1'.format(index),
									  'EGH model {}1'.format(index),
									  'exp(((-(@0-@1)^2)/(2*@2^2+@3*(@0-@1))))',
									  RooArgList(self.x, self.mean1, self.sigma, self.tau))
		# Add PDFs
		self.this_pdf = RooAddPdf(f'double_Egh{index}', f'Two EGHs added {index}',
								  RooArgList(sig0, sig1), RooArgList(self.ratio))
		self.roodefs = [self.this_pdf, self.x, self.mean0, self.mean1, self.sigma, self.tau]
		# self.fix_var(self.roodefs[2:], index)
		# Do fit
		self.this_roohist, self.fit_results = self.minimize(self.hist, self.xmin, self.xmax)
		# Return results
		return [self.this_roohist, self.this_pdf, self.nb_pars,
				self.x, self.mean0, self.mean1, self.sigma, self.tau, self.pstate,
				self.my_name, self.spl_draw, self.rooCompon, self.fit_results]           

	def plot(self, bins=2000, log=False):
		"""  """
		if len(self.lst_file) == 0:
			print("Fit not excecuted yet or failed.")
			return 0 
		#            
		xdata = self.lst_file.tof
		n, xe = np.histogram(xdata, bins=bins)
		cx = 0.5 * (xe[1:] + xe[:-1])
		dx = np.diff(xe)
		plt.errorbar(cx, n, n ** 0.5, fmt="ok", zorder=1)
		xm = np.linspace(xe[0], xe[-1], num=1000)
		# Plot 'carpet'
		if log:
			plt.plot(xdata, np.zeros_like(xdata)+0.9, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		plt.plot(xdata, np.zeros_like(xdata)-5, "|", alpha=0.1, label = "ToF Data", zorder = 3)
		#
		# Plot fits, loop through all fits
		if log:
			xm = np.linspace(self.mean0.getValV()-100, self.mean1.getValV()+100, num=1000)
		# Calculate values of PDFs for plotting
		y_val = [(self.ratio.getValV() * RootMath.gaussian_pdf(x=y, x0=self.mean0.getValV(), sigma=self.sigma.getValV()) +
				 (1-self.ratio.getValV()) * RootMath.gaussian_pdf(x=y, x0=self.mean1.getValV(), sigma=self.sigma.getValV()) )
				 for y in xm]
		# Normalize values
		integral_cut = sum(y_val) * np.diff(xm)[0]
		left_n_cut = len(xe[xe<self.xmin])
		right_n_cut = len(xe[xe<self.xmax])
		n_cut = n[left_n_cut:right_n_cut]        
		y_val = y_val / integral_cut * sum(n_cut) * dx[0]
		# Plot
		plt.plot(xm, y_val, ls=":", label=f"init", zorder=2)
		plt.plot(xm, y_val, label=f"fit", c='r', zorder=3)
		
		plt.legend();
		if log:
			plt.yscale("log")
		# Zoom in on found peaks
		if self.peaks.n_peaks != 0:
			plt.xlim(self.peaks.earliest_left_base-200, self.peaks.latest_right_base+200)

		plt.show()
