# -*- coding: utf-8 -*-
"""
Created on Mon 08 August 2023
@author: Lukas Nies
@contact: Lukas.Nies@cern.ch
@license: MIT
"""

import ROOT
from ROOT import RooFit
from ROOT import RooRealVar, RooArgList, RooArgSet, RooAddPdf, RooGenericPdf, RooDataHist
from ROOT import Math as RootMath 
import pandas as pd
import numpy as np
import time

def lst2roothist(df, bins=1):
	"""
	Converts .lst file to root histogram.
	:return: root histogram
	"""
	# Safe the histogrammed data in arrays
	xdata_unbinned = df.tof
	# Get min and max tof from data frame
	pre_bin_size = 0.8 # intrinsic binning from TDC in ns
	binning = round((df.tof.max()-df.tof.min())/pre_bin_size/bins)
	# Fill TH1D
	hist = ROOT.TH1D( 'hist', 'hist converted', binning, df.tof.min(), df.tof.max())
	[hist.Fill(x) for x in df.tof]
	#
	return hist


def main():
	
	print('Example hyperEMG fit.')

	# Load data 
	df = pd.read_csv('./data.csv')
	# df = pd.read_csv('./data.csv')
	histTH1D = lst2roothist(df)

	# Initial guesses and fit range
	mu0_init = 55445260 # ns
	mu1_init = 55445260 + 2000 # ns
	xmin = mu0_init - 250 # ns
	xmax = mu1_init + 250 # ns


	# Create variables to be fitted
	x = RooRealVar("x", "x", xmin, xmax, 'ns') # range of x given by fit range
	x.setRange('x-range', xmin, xmax)

	mu0_var = RooRealVar("mu0", "mu0", mu0_init, mu0_init-200, mu0_init+200, 'ns')
	mu1_var = RooRealVar("mu1", "mu1", mu1_init, mu1_init-200, mu1_init+200, 'ns')
	sigma_var = RooRealVar("sigma", "sigma", 40, 10, 100, 'ns')
	ntau0_var = RooRealVar("ntau0", "ntau0", 30, 10, 100, 'ns') 
	ptau0_var = RooRealVar("ptau0", "ptau0", 50, 10, 250, 'ns')
	ptau1_var = RooRealVar("ptau1", "ptau1", 300, 250, 500, 'ns')
	contrib0_var = RooRealVar("contrib0", "contrib0", 0.8, 0.5, 0.99, '%')
	contrib1_var = RooRealVar("contrib1", "contrib1", 0.2, 0.01, 0.5, '%')
	ratio_var = RooRealVar("ratio", "ratio", 0.3, 0.1, 0.5, '%')


	## Method 1: build individual hyperEMGs and add them later
	# This method is very very slow

	# hyperEMG12_string = '@4*1/(2*@3)*exp((@2/(1.4142*@3))^2+(@0-@1)/@3)*TMath::Erfc(@2/(1.4142*@3)+(@0-@1)/(1.4142*@2))+@6*1/(2*@5)*exp((@2/(1.4142*@5))^2-(@0-@1)/@5)*TMath::Erfc(@2/(1.4142*@5)-(@0-@1)/(1.4142*@2))+(1-@4-@6)*1/(2*@7)*exp((@2/(1.4142*@7))^2-(@0-@1)/@7)*TMath::Erfc(@2/(1.4142*@7)-(@0-@1)/(1.4142*@2))'

	# hyperEMG_pdf1 = RooGenericPdf('hyperEMG1','hyperEMG1', hyperEMG12_string, RooArgList(x, mu0_var, sigma_var, ntau0_var, contrib0_var, ptau0_var, contrib1_var, ptau1_var)) 
	# hyperEMG_pdf2 = RooGenericPdf('hyperEMG2','hyperEMG2', hyperEMG12_string, RooArgList(x, mu1_var, sigma_var, ntau0_var, contrib0_var, ptau0_var, contrib1_var, ptau1_var)) 

	# fitmodel = RooAddPdf('model','model', RooArgList(hyperEMG_pdf1,hyperEMG_pdf2), RooArgList(ratio_var))


	## Method 2: build one single master fit model through string
	
	hyperEMG12_string = '@4*1/(2*@3)*exp((@2/(1.4142*@3))^2+(@0-@1)/@3)*TMath::Erfc(@2/(1.4142*@3)+(@0-@1)/(1.4142*@2))+@6*1/(2*@5)*exp((@2/(1.4142*@5))^2-(@0-@1)/@5)*TMath::Erfc(@2/(1.4142*@5)-(@0-@1)/(1.4142*@2))+(1-@4-@6)*1/(2*@7)*exp((@2/(1.4142*@7))^2-(@0-@1)/@7)*TMath::Erfc(@2/(1.4142*@7)-(@0-@1)/(1.4142*@2))'
	hyperEMG12_string2 = '@4*1/(2*@3)*exp((@2/(1.4142*@3))^2+(@0-@8)/@3)*TMath::Erfc(@2/(1.4142*@3)+(@0-@8)/(1.4142*@2))+@6*1/(2*@5)*exp((@2/(1.4142*@5))^2-(@0-@8)/@5)*TMath::Erfc(@2/(1.4142*@5)-(@0-@8)/(1.4142*@2))+(1-@4-@6)*1/(2*@7)*exp((@2/(1.4142*@7))^2-(@0-@8)/@7)*TMath::Erfc(@2/(1.4142*@7)-(@0-@8)/(1.4142*@2))'
	fitmodel_string = '@9*('+hyperEMG12_string+")+(1-@9)*("+hyperEMG12_string2+")"

	print(fitmodel_string)

	longarglist = RooArgList(x, mu0_var, sigma_var, ntau0_var, contrib0_var, ptau0_var, contrib1_var, ptau1_var)
	longarglist.add(mu1_var)
	longarglist.add(ratio_var)

	fitmodel = RooGenericPdf('model','model', fitmodel_string, longarglist)

	# Fit data 

	roodata = RooDataHist("data", "data", x, histTH1D)

	result = fitmodel.fitTo(roodata, RooFit.Range('x-range'),
								RooFit.Minos(ROOT.kTRUE),
								RooFit.PrintEvalErrors(-1),
								RooFit.NumCPU(1),
								RooFit.Timer(ROOT.kTRUE),
								RooFit.Save(),
								RooFit.Verbose(ROOT.kFALSE))


	frame = x.frame(ROOT.RooFit.Name("xframe"), ROOT.RooFit.Title("RooPlot"), ROOT.RooFit.Bins(20))
	
	roodata.plotOn(frame)#, ROOT.RooFit.DataError("SumW2"))
	fitmodel.plotOn(frame)

	canvas = ROOT.TCanvas("canvas")
	canvas.cd()

	frame.Draw()

	time.sleep(10000)

if __name__ == '__main__':
  main()