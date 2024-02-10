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
from process import MCS6Lst

def main():
	
	print('Example plot MR-ToF MS data.')

	print('Processing data, this might take a while.')
	# Process data
	proc = MCS6Lst()
	files = [
		'./data.lst', 
	]
	# proc.process(files, to_csv=True) # Use this if the data has not yet been converted to .csv

	# Load AME into a helper class
	path_to_ame = './mass20.txt' # path to 2020 AME file 
	calc = MRToFUtils(path_to_ame)

	# Load MR-ToF MS data 
	df = pd.read_csv('./data.csv')
	# Load Time-of-Flight calibration
	ToF_Calib = './tof_calib.xlsm'
	calc.load_tof_calib(ToF_Calib)

	# Find peaks in data
	peaks=Peaks(df)
	peaks.find_peaks(bins=20)
	# Plot data with found peaks
	peaks.plot(focus=False, log=True, lines=True, bins=20)
	# Plot data in 2D
	peaks.plot2d(x_bins=20, hist2d_y_bins = 1000, y_bins = 100, focus=False, log=True)
	# Plot data with all A=88 isobars, where data was taken at 2000 revolutions
	peaks.load_utils(calc) # first load time-of-flight calibration into peaks class
	peaks.add_isobars(A=88, nrevs=2000)
	peaks.plot(focus=False, log=True, lines=True, bins=20)
	# Recalibrate time-of-flight calibration based on known species in the spectrum
	calc.recalibrate_tof(calc.get_value('88Sr', value='mass', state='gs'), peaks.pos[0]/1e3, 2000)
	peaks.add_isobars(iso_list=['88Rb', '88Sr'], nrevs=2000) # add only certain isotopes
	peaks.add_isobars(iso_list=['69Br19F'], nrevs=2000) # you can also add molecules
	peaks.plot(focus=False, log=True, lines=True, bins=20) # plot again

if __name__ == '__main__':
  main()