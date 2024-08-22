# -*- coding: utf-8 -*-
"""
Created on Mon 08 August 2023
@author: Lukas Nies
@contact: Lukas.Nies@cern.ch
@license: MIT
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
import time
import sys
sys.path.insert(0, '../../bin/')
from utilities import Peaks, MRToFUtils
from fit import hyperEmg


def main():
	file_path = './data'
	data_df = pd.read_csv(file_path+".csv")

	bins=20
	peaks=Peaks(data_df)
	peaks.find_peaks(bins=bins, peak_min_distance=1000)

	fit = hyperEmg(data_df, file_path=file_path, peaks=peaks)
	fit.lst2roothist(data_df)

	first_peak = 45140066.1
	second_peak = 45140066.1 + 1000
	third_peak = 45140066.1 + 1450
	fourth_peak = 45140066.1 + 2930
	fifth_peak = 45140066.1 + 4150
	shoulder = 45140066.1 + 300

	limits = {
    'mu0':[first_peak, first_peak-50, first_peak+50],
    'mu1':[shoulder, shoulder -25, shoulder +52],
    'mu2':[second_peak, second_peak -100, second_peak +100],
    'mu3':[third_peak, third_peak -100, third_peak +100],
    'mu4':[fourth_peak, fourth_peak -100, fourth_peak +100],
    'mu5':[fifth_peak, fifth_peak -100, fifth_peak +100],
    #
    'sigma':[20, 5, 60],
    #
    'ntau0':[25 , 5 , 1000],
    'ptau0':[25 ,2.5 , 100],
    'ptau1':[250 ,100 , 750],

    'contrib0':[0.1, 0.001, 0.3], # 1st neg. comp. of early isotope
    'contrib1':[0.8, 0.4, 1], # 1st neg. comp. of late isotope

    'ratio0':[0.65, 0.55, 0.99],
    'ratio1':[0.05, 0.01, 0.2],
    'ratio2':[0.1, 0.01, 0.40],
    'ratio3':[0.025, 0.01, 0.3],
    'ratio4':[0.05, 0.01, 0.3],
	}	
	result = fit.call_pdf(first_peak-200, fifth_peak+300, dimensions = [1,2], n_comps = 6, 
						simultaneous = True, bins = bins, limits=limits, minos=False)

	fit.plot(bins=20, log=True, from_file=False, focus=False, file_out = file_path+"_fit.pdf", histlw = 2, xrange_mod = [1200,3000],
			 silent=False, centroids=False, components=False, carpet=False, legend=True, style='hist', residuals=True, prelim=True)

	fit.save_fit(file_path+"_fit.txt")

# #############
if __name__ == '__main__':
  main()