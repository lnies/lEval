# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 17:43:00 2016

@author: mamougeo

modified on Mon Mar 27 2017: Define list_ini_file to read the ini files.
This function also replaces the shell script list_ini_files.sh

modified on Thu Mar 30 2017: add find_peak function 

"""

import sys,os
import numpy as np
import matplotlib.pyplot as plt
from fit_routine import *
import glob as gb 
#from logger import Logger


def read_ini_file(filename_ini):
    """
    Load analysis ini file
    :param filename_ini:
    :return: numpy array containing the parameters to be used by the fit routine
    """
    return np.genfromtxt(filename_ini, dtype='object')

def list_ini_file(folder_path):
	pattern = folder_path + '*.ini'
	return gb.glob(pattern)
	
def main(argv):
    """
    Main function of the script. Performs the fit for the file specified in ini file
    :param argv: Command line arguments : Ini_file target_directory
    :return: None
    """
    list_ini = list_ini_file(argv[0])
    for file in list_ini:
        routine_parameters = read_ini_file(file)
        mrtof_analysis(routine_parameters)
		#print routine_parameters
        #for par in routine_parameters[:]:
        #    if par[4] in ['Gaussian']:
        #        print "start"
        #        mrtof_analysis(par)
                


if __name__ == "__main__":
    main(sys.argv[1:])
