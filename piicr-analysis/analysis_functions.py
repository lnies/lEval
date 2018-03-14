# ---------------------------------------------------------------------------
# Written by Jonas Karthein in 2016/2017. Questions to jonas.karthein@cern.ch
# ---------------------------------------------------------------------------

import math
import matplotlib as mpl
mpl.use('Qt5Agg')

import matplotlib.pyplot as plt
import numpy as np
from piicr_daq import piicr_daq
from analysis_functions import *
import csv
import sys, os
import glob
import datetime
import pandas as pd
from read_write_functions import *
from python_plotter_functions import *
import xml.etree.ElementTree as ET
import shutil
import platform
import subprocess
import time
from cross_checks_functions import weighted_mean


def piicr_add_time_stamps(file_path):
    '''Adds timestamps for raw piicr data from .tdc to .txt'''
    if platform.system() == 'Windows':
        exe_str = 'G:\\Experiments\\ISOLTRAP\\Software\\PI-ICR\\Python-DAQ\\remote_piicr_add_timestamps.exe'
        # file_path = 'E:\\Jonas\\test'     # folder where raw data .txt and .tdc is placed
        process = subprocess.Popen(exe_str+' '+file_path,  stderr=subprocess.PIPE)
        time.sleep(5)   # just a waiting time until all files are manipulated...writing is not so fast it seems :D
        # process.kill()  # should be fine now, it closes automatically (but needs time to write the files, the code is faster than saving)


def folder_creation_analysis(folder_name, isotopes):
    os.chdir(folder_name)
    txt_files = []                      # looks for all files with .txt ending
    for file in glob.glob("*.txt"):
        txt_files.append(file)

    for i in isotopes:
        if platform.system() == 'Windows':
            if os.path.isdir(folder_name+'\\'+i) == False:
                    os.mkdir(folder_name+'\\'+i)
            if os.path.isdir(folder_name+'\\'+i+'\\c') == False:
                os.mkdir(folder_name+'\\'+i+'\\c')
            if os.path.isdir(folder_name+'\\'+i+'\\p1p2') == False:
                os.mkdir(folder_name+'\\'+i+'\\p1p2')
            if os.path.isdir(folder_name+'\\'+i+'\\freq_config') == False:
                os.mkdir(folder_name+'\\'+i+'\\freq_config')

            for j in txt_files:
                if (i+'_c') in j:   # move center files
                    shutil.move(folder_name+'\\'+j, folder_name+'\\'+i+'\\c\\'+j)
                if i in j:          # move p1p2 files
                    if os.path.isfile(folder_name+'\\'+j) == True:
                        shutil.move(folder_name+'\\'+j, folder_name+'\\'+i+'\\p1p2\\'+j)

            out_files = []
            for file in glob.glob("*.out"):
                if str(file[:len(i)]) == str(i):
                    out_files.append(file)
            for out in out_files:
                if os.path.isfile(folder_name+'\\'+out) == True:
                    shutil.move(folder_name+'\\'+out, folder_name+'\\'+i+'\\freq_config\\'+out)

            # if os.path.isfile(folder_name+'\\'+i+'.out') == True:
            #     shutil.move(folder_name+'\\'+i+'.out', folder_name+'\\'+i+'\\freq_config\\'+i+'.out')



        elif platform.system() == 'Darwin': # MacOS
            if os.path.isdir(folder_name+'/'+i) == False:
                os.mkdir(folder_name+'/'+i)
            if os.path.isdir(folder_name+'/'+i+'/c') == False:
                os.mkdir(folder_name+'/'+i+'/c')
            if os.path.isdir(folder_name+'/'+i+'/p1p2') == False:
                os.mkdir(folder_name+'/'+i+'/p1p2')
            if os.path.isdir(folder_name+'/'+i+'/freq_config') == False:
                os.mkdir(folder_name+'/'+i+'/freq_config')

            for j in txt_files:
                if (i+'_c') in j:   # move center files
                    shutil.move(folder_name+'/'+j, folder_name+'/'+i+'/c/'+j)
                if i in j:          # move p1p2 files
                    if os.path.isfile(folder_name+'/'+j) == True:
                        shutil.move(folder_name+'/'+j, folder_name+'/'+i+'/p1p2/'+j)

            out_files = []
            for file in glob.glob("*.out"):
                if str(file[:len(i)]) == str(i):
                    out_files.append(file)
            for out in out_files:
                if os.path.isfile(folder_name+'/'+out) == True:
                    shutil.move(folder_name+'/'+out, folder_name+'/'+i+'/freq_config/'+out)


    if platform.system() == 'Darwin':
        os.chdir(folder_name)
        tdc_files = []                      # looks for all files with .tdc ending
        for file in glob.glob("*.tdc"):
            tdc_files.append(file)
        if os.path.isdir(folder_name+'/bin') == False:
            os.mkdir(folder_name+'/bin')
        for k in tdc_files:
            if os.path.isfile(folder_name+'/'+k) == True:
                shutil.move(folder_name+'/'+k, folder_name+'/bin/'+k)

        os.chdir(folder_name)
        bin_files = []                      # looks for all files with .tdc ending
        for file in glob.glob("*.bin"):
            bin_files.append(file)
        if os.path.isdir(folder_name+'/bin') == False:
            os.mkdir(folder_name+'/bin')
        for l in bin_files:
            if os.path.isfile(folder_name+'/'+l) == True:
                shutil.move(folder_name+'/'+l, folder_name+'/bin/'+l)


    elif platform.system() == 'Windows': # MacOS
        os.chdir(folder_name)
        tdc_files = []                      # looks for all files with .tdc ending
        for file in glob.glob("*.tdc"):
            tdc_files.append(file)
        if os.path.isdir(folder_name+'\\bin') == False:
            os.mkdir(folder_name+'\\bin')
        for k in tdc_files:
            if os.path.isfile(folder_name+'\\'+k) == True:
                shutil.move(folder_name+'\\'+k, folder_name+'\\bin\\'+k)

        os.chdir(folder_name)
        bin_files = []                      # looks for all files with .tdc ending
        for file in glob.glob("*.bin"):
            bin_files.append(file)
        if os.path.isdir(folder_name+'\\bin') == False:
            os.mkdir(folder_name+'\\bin')
        for l in bin_files:
            if os.path.isfile(folder_name+'\\'+l) == True:
                shutil.move(folder_name+'\\'+l, folder_name+'\\bin\\'+l)



def angle_calculator(centerXY, all_peak_positions_p1, all_peak_positions_p2, file_name, upper_folder_path, folder_name_0, latex_name):
    """The function calculates the phase between P1 nd P2."""
    axis_scaling_factor = 0.031746032

    vector_p1 = []
    vector_p2 = []
    for i in range(0, len(all_peak_positions_p1)-1, 1):
        hilf_p1 = [0, 0]
        hilf_p1[0] = all_peak_positions_p1[i+1][1] - centerXY[1][1]     # P1_x - C_x
        hilf_p1[1] = all_peak_positions_p1[i+1][3] - centerXY[1][3]     # P1_y - C_y
        vector_p1.append(hilf_p1)
        hilf_p2 = [0, 0]
        hilf_p2[0] = all_peak_positions_p2[i+1][1] - centerXY[1][1]     # P2_x - C_x
        hilf_p2[1] = all_peak_positions_p2[i+1][3] - centerXY[1][3]     # P2_y - C_y
        vector_p2.append(hilf_p2)

    #   lengths of both vectors
    length_vector_p1 = []
    length_vector_p2 = []
    for i in range(len(vector_p1)):
        hilf_length_p1 = 0
        hilf_length_p1 = math.sqrt(vector_p1[i][0] ** 2 + vector_p1[i][1] ** 2)
        length_vector_p1.append(hilf_length_p1)
        hilf_length_p2 = 0
        hilf_length_p2 = math.sqrt(vector_p2[i][0] ** 2 + vector_p2[i][1] ** 2)
        length_vector_p2.append(hilf_length_p2)

    vector_lengths = [['Inj. #', 'length of vector / ch.', 'length of vector / ch.']]
    for i in range(len(length_vector_p1)):
        vector_lengths.append([i+1, length_vector_p1[i], length_vector_p2[i]])

    python_plot(vector_lengths, 'Lengths of P1 & P2 vectors', '%s_p1p2_vector_lengths' % (file_name), 'no', ['P1', 'P2'], 'no', 'gauss', 'Utopia', 'red', ['blue', 'red'], 1, 10, 5, 3, '2lines-no-error', 50, 1, 'on', 'off', '', 0, 0)
    vector_lengths = [['Inj. #', 'length of P1 vector / ch.', 'length of P2 vector / ch.', 'length of P1 vector / mm', 'length of P2 vector / mm']]
    for i in range(len(length_vector_p1)):
        vector_lengths.append([i+1, length_vector_p1[i], length_vector_p2[i], length_vector_p1[i] * axis_scaling_factor, length_vector_p2[i] * axis_scaling_factor])
    write_csv(upper_folder_path, vector_lengths, 'vector_length_%s' % folder_name_0)

    #   math.atan2(y, x) use for angle calculation. angle between 0 and 2pi. --> phase between p1 and p2  -#-#-#- 11.05.17: added cases to correct angle to around 2pi
    p1p2_angle = []
    for i in range(len(vector_p1)):
        hilf_angle = [all_peak_positions_p1[i+1][0][:-3], 0]
        hilf_angle[1] = 2 * math.pi - math.atan2(vector_p2[i][1], vector_p2[i][0]) + math.atan2(vector_p1[i][1], vector_p1[i][0])
        if hilf_angle[1] < math.pi:     # correct for cases where angle is not +/- 2pi
            hilf_angle[1] += 2 * math.pi
        elif hilf_angle[1] > 3 * math.pi:
            hilf_angle[1] -= 2 * math.pi
        p1p2_angle.append(hilf_angle)

    #   derivative calculation for error; see https://en.wikipedia.org/wiki/Atan2
    d_angle_p1_x = []     # derivative of angle for P1_x
    d_angle_p2_x = []
    d_angle_p1_y = []
    d_angle_p2_y = []
    d_angle_c_x = []
    d_angle_c_y = []
    for i in range(len(vector_p1)):
        hilf_err_angle_p1_x = 0
        hilf_err_angle_p1_x = vector_p1[i][1] / length_vector_p1[i] / length_vector_p1[i]
        d_angle_p1_x.append(hilf_err_angle_p1_x)
        hilf_err_angle_p2_x = 0
        hilf_err_angle_p2_x = - vector_p2[i][1] / length_vector_p2[i] / length_vector_p2[i]
        d_angle_p2_x.append(hilf_err_angle_p2_x)
        hilf_err_angle_p1_y = 0
        hilf_err_angle_p1_y = - vector_p1[i][0] / length_vector_p1[i] / length_vector_p1[i]
        d_angle_p1_y.append(hilf_err_angle_p1_y)
        hilf_err_angle_p2_y = 0
        hilf_err_angle_p2_y = vector_p2[i][0] / length_vector_p2[i] / length_vector_p2[i]
        d_angle_p2_y.append(hilf_err_angle_p2_y)
        hilf_err_angle_c_x = 0
        hilf_err_angle_c_x = vector_p2[i][1] / length_vector_p2[i] / length_vector_p2[i] - vector_p1[i][1] / length_vector_p1[i] / length_vector_p1[i]
        d_angle_c_x.append(hilf_err_angle_c_x)
        hilf_err_angle_c_y = 0
        hilf_err_angle_c_y = vector_p1[i][0] / length_vector_p1[i] / length_vector_p1[i] - vector_p2[i][0] / length_vector_p2[i] / length_vector_p2[i]
        d_angle_c_y.append(hilf_err_angle_c_y)

    #   error calculation using error of position fit (not FWHM/standard deviation!)
    p1p2_angle_error = []
    for i in range(len(vector_p1)):
        hilf_angle_error = [all_peak_positions_p1[i+1][0], 0]
        hilf_angle_error[1] = math.sqrt((((d_angle_p1_x[i]) ** 2) * ((all_peak_positions_p1[i+1][2]) ** 2)) +
                                        (((d_angle_p2_x[i]) ** 2) * ((all_peak_positions_p2[i+1][2]) ** 2)) +
                                        (((d_angle_p1_y[i]) ** 2) * ((all_peak_positions_p1[i+1][4]) ** 2)) +
                                        (((d_angle_p2_y[i]) ** 2) * ((all_peak_positions_p2[i+1][4]) ** 2)) +
                                        (((d_angle_c_x[i]) ** 2) * ((centerXY[1][2]) ** 2)) +
                                        (((d_angle_c_y[i]) ** 2) * ((centerXY[1][4]) ** 2)))
        p1p2_angle_error.append(hilf_angle_error)


    # for a calculation of nu_+ and nu_- we also calculate the individual phases

    p1_angle = []   # corresponds to nu_-
    for i in range(len(vector_p1)):
        hilf_angle_p1 = [all_peak_positions_p1[i+1][0], 0]
        hilf_angle_p1[1] = math.pi + math.atan2(vector_p1[i][1], vector_p1[i][0])
        p1_angle.append(hilf_angle_p1)
    p2_angle = []   # corresponds to nu_+
    for i in range(len(vector_p2)):
        hilf_angle_p2 = [all_peak_positions_p1[i+1][0], 0]
        hilf_angle_p2[1] = math.pi - math.atan2(vector_p2[i][1], vector_p2[i][0])
        p2_angle.append(hilf_angle_p2)
    p1_angle_error = []
    for i in range(len(vector_p1)):
        hilf_angle_error_p1 = [all_peak_positions_p1[i+1][0], 0]
        hilf_angle_error_p1[1] = math.sqrt((((d_angle_p1_x[i]) ** 2) * ((all_peak_positions_p1[i+1][2]) ** 2)) +
                                        (((d_angle_p2_x[i]) ** 2) * ((all_peak_positions_p2[i+1][2]) ** 2)) +
                                        (((d_angle_c_x[i]) ** 2) * ((centerXY[1][2]) ** 2)))
        p1_angle_error.append(hilf_angle_error_p1)
    p2_angle_error = []
    for i in range(len(vector_p2)):
        hilf_angle_error_p2 = [all_peak_positions_p1[i+1][0], 0]
        hilf_angle_error_p2[1] = math.sqrt((((d_angle_p1_y[i]) ** 2) * ((all_peak_positions_p1[i+1][4]) ** 2)) +
                                        (((d_angle_p2_y[i]) ** 2) * ((all_peak_positions_p2[i+1][4]) ** 2)) +
                                        (((d_angle_c_y[i]) ** 2) * ((centerXY[1][4]) ** 2)))
        p2_angle_error.append(hilf_angle_error_p2)

    angle_info = [['Name', 'P1-P2-angle / radian', 'P1-P2-angle-unc.', 'P1-angle', 'P1-angle-unc.', 'P2-angle', 'P2-angle-unc.']]
    angle_plot = [['Run #', 'P1-P2-angle / radian', 'unc.']]
    for i in range(len(p1p2_angle)):
        angle_info.append([p1p2_angle[i][0], p1p2_angle[i][1], p1p2_angle_error[i][1], p1_angle[i][1], p2_angle[i][1], p1_angle_error[i][1], p2_angle_error[i][1]])
        angle_plot.append([int(p1p2_angle[i][0][-3:]), p1p2_angle[i][1], p1p2_angle_error[i][1]])

    # python_plot(angle_plot, 'P1-P2-angle for %s' % (latex_name), '%s_p1_p2_angle' % (folder_name_0), 'no', '%s P1-P2-angle' % (latex_name), 'no', 'linear', 'Utopia', 'red', 'green', 'full', 1, 5, 1, 'scatter', 20, 3, 'on', 'off', '', '', '')
    write_csv(upper_folder_path, angle_info, 'angle_info_%s' % folder_name_0)

    return(p1p2_angle, p1p2_angle_error, p1_angle, p2_angle, p1_angle_error, p2_angle_error)


def angle_calculator_two(centerXY, all_peak_positions_p1, all_peak_positions_p2, file_name, upper_folder_path, folder_name_0, latex_name):
    """The function calculates the phase between P1 nd P2."""
    axis_scaling_factor = 0.031746032

    vector_p1 = []
    vector_p2_dom = []
    vector_p2_rec = []
    for i in range(0, len(all_peak_positions_p1)-1, 1):
        hilf_p1 = [0, 0]
        hilf_p1[0] = all_peak_positions_p1[i+1][1] - centerXY[1][1]     # P1_x - C_x
        hilf_p1[1] = all_peak_positions_p1[i+1][3] - centerXY[1][3]     # P1_y - C_y
        vector_p1.append(hilf_p1)
        hilf_p2_dom = [0, 0]
        hilf_p2_rec = [0, 0]
        hilf_p2_dom[0] = all_peak_positions_p2[i+1][1] - centerXY[1][1]     # P2_x - C_x
        hilf_p2_dom[1] = all_peak_positions_p2[i+1][3] - centerXY[1][3]     # P2_y - C_y
        hilf_p2_rec[0] = all_peak_positions_p2[i+1][5] - centerXY[1][1]     # P2_x - C_x
        hilf_p2_rec[1] = all_peak_positions_p2[i+1][7] - centerXY[1][3]     # P2_y - C_y
        vector_p2_dom.append(hilf_p2_dom)
        vector_p2_rec.append(hilf_p2_rec)

    #   lengths of both vectors
    length_vector_p1 = []
    length_vector_p2_dom = []
    length_vector_p2_rec = []
    for i in range(len(vector_p1)):
        hilf_length_p1 = 0
        hilf_length_p1 = math.sqrt(vector_p1[i][0] ** 2 + vector_p1[i][1] ** 2)
        length_vector_p1.append(hilf_length_p1)
        hilf_length_p2_dom = 0
        hilf_length_p2_dom = math.sqrt(vector_p2_dom[i][0] ** 2 + vector_p2_dom[i][1] ** 2)
        length_vector_p2_dom.append(hilf_length_p2_dom)
        hilf_length_p2_rec = 0
        hilf_length_p2_rec = math.sqrt(vector_p2_rec[i][0] ** 2 + vector_p2_rec[i][1] ** 2)
        length_vector_p2_rec.append(hilf_length_p2_rec)

    vector_lengths_dom = [['Inj. #', 'length of p1 vector / ch.', 'length of p2 (dominant) vector / ch.']]
    vector_lengths_rec = [['Inj. #', 'length of p1 vector / ch.', 'length of p2 (recessive) vector / ch.']]
    for i in range(len(length_vector_p1)):
        vector_lengths_dom.append([i+1, length_vector_p1[i], length_vector_p2_dom[i]])
        vector_lengths_rec.append([i+1, length_vector_p1[i], length_vector_p2_rec[i]])

    python_plot(vector_lengths_dom, 'Lengths of P1 & P2 (dominant) vectors', '%s_p1p2_vector_length_dom' % (file_name), 'no', ['P1', 'P2'], 'no', 'gauss', 'Utopia', 'red', ['blue', 'red'], 1, 10, 5, 3, '2lines-no-error', 50, 1, 'on', 'off', '', 0, 0)
    python_plot(vector_lengths_rec, 'Lengths of P1 & P2 (recessive) vectors', '%s_p1p2_vector_length_rec' % (file_name), 'no', ['P1', 'P2'], 'no', 'gauss', 'Utopia', 'red', ['blue', 'red'], 1, 10, 5, 3, '2lines-no-error', 50, 1, 'on', 'off', '', 0, 0)

    vector_lengths = [['Inj. #', 'length of P1 vector / ch.', 'length of P2 (dominant) vector / ch.', 'length of P2 (recessive) vector / ch.', 'length of P1 vector / mm', 'length of P2 (dominant) vector / mm', 'length of P2 (recessive) vector / mm']]
    for i in range(len(length_vector_p1)):
        vector_lengths.append([i+1, length_vector_p1[i], length_vector_p2_dom[i], length_vector_p2_rec[i], length_vector_p1[i] * axis_scaling_factor, length_vector_p2_dom[i] * axis_scaling_factor, length_vector_p2_rec[i] * axis_scaling_factor])
    write_csv(upper_folder_path, vector_lengths, 'vector_length_%s' % folder_name_0)

    #   math.atan2(y, x) use for angle calculation. angle between 0 and 2pi. --> phase between p1 and p2  -#-#-#- 11.05.17: added cases to correct angle to around 2pi
    p1p2_angle_dom = []         # dominant peak
    for i in range(len(vector_p1)):
        hilf_angle = [all_peak_positions_p1[i+1][0][:-3], 0]
        hilf_angle[1] = 2 * math.pi - math.atan2(vector_p2_dom[i][1], vector_p2_dom[i][0]) + math.atan2(vector_p1[i][1], vector_p1[i][0])
        if hilf_angle[1] < math.pi:     # correct for cases where angle is not +/- 2pi
            hilf_angle[1] += 2 * math.pi
        elif hilf_angle[1] > 3 * math.pi:
            hilf_angle[1] -= 2 * math.pi
        p1p2_angle_dom.append(hilf_angle)

    p1p2_angle_rec = []         # recessive peak
    for i in range(len(vector_p1)):
        hilf_angle = [all_peak_positions_p1[i+1][0][:-3], 0]
        hilf_angle[1] = 2 * math.pi - math.atan2(vector_p2_rec[i][1], vector_p2_rec[i][0]) + math.atan2(vector_p1[i][1], vector_p1[i][0])
        if hilf_angle[1] < math.pi:     # correct for cases where angle is not +/- 2pi
            hilf_angle[1] += 2 * math.pi
        elif hilf_angle[1] > 3 * math.pi:
            hilf_angle[1] -= 2 * math.pi
        p1p2_angle_rec.append(hilf_angle)

    #   derivative calculation for error; see https://en.wikipedia.org/wiki/Atan2
    d_angle_p1_x = []     # derivative of angle for P1_x
    d_angle_p2_dom_x = []
    d_angle_p2_rec_x = []
    d_angle_p1_y = []
    d_angle_p2_dom_y = []
    d_angle_p2_rec_y = []
    d_angle_c_x = []
    d_angle_c_y = []
    for i in range(len(vector_p1)):
        hilf_err_angle_p1_x = 0
        hilf_err_angle_p1_x = vector_p1[i][1] / length_vector_p1[i] / length_vector_p1[i]
        d_angle_p1_x.append(hilf_err_angle_p1_x)
        hilf_err_angle_p2_dom_x = 0
        hilf_err_angle_p2_dom_x = - vector_p2_dom[i][1] / length_vector_p2_dom[i] / length_vector_p2_dom[i]
        d_angle_p2_dom_x.append(hilf_err_angle_p2_dom_x)
        hilf_err_angle_p2_rec_x = 0
        hilf_err_angle_p2_rec_x = - vector_p2_rec[i][1] / length_vector_p2_rec[i] / length_vector_p2_rec[i]
        d_angle_p2_rec_x.append(hilf_err_angle_p2_rec_x)
        hilf_err_angle_p1_y = 0
        hilf_err_angle_p1_y = - vector_p1[i][0] / length_vector_p1[i] / length_vector_p1[i]
        d_angle_p1_y.append(hilf_err_angle_p1_y)
        hilf_err_angle_p2_dom_y = 0
        hilf_err_angle_p2_dom_y = vector_p2_dom[i][0] / length_vector_p2_dom[i] / length_vector_p2_dom[i]
        d_angle_p2_dom_y.append(hilf_err_angle_p2_dom_y)
        hilf_err_angle_p2_rec_y = 0
        hilf_err_angle_p2_rec_y = vector_p2_rec[i][0] / length_vector_p2_rec[i] / length_vector_p2_rec[i]
        d_angle_p2_rec_y.append(hilf_err_angle_p2_rec_y)
        hilf_err_angle_c_x = 0
        hilf_err_angle_c_x = vector_p2_rec[i][1] / length_vector_p2_rec[i] / length_vector_p2_rec[i] - vector_p1[i][1] / length_vector_p1[i] / length_vector_p1[i]
        d_angle_c_x.append(hilf_err_angle_c_x)
        hilf_err_angle_c_y = 0
        hilf_err_angle_c_y = vector_p1[i][0] / length_vector_p1[i] / length_vector_p1[i] - vector_p2_rec[i][0] / length_vector_p2_rec[i] / length_vector_p2_rec[i]
        d_angle_c_y.append(hilf_err_angle_c_y)

    #   error calculation using error of position fit (not FWHM/standard deviation!)
    p1p2_angle_dom_error = []       # dominant
    for i in range(len(vector_p1)):
        hilf_angle_dom_error = [all_peak_positions_p1[i+1][0], 0]
        hilf_angle_dom_error[1] = math.sqrt((((d_angle_p1_x[i]) ** 2) * ((all_peak_positions_p1[i+1][2]) ** 2)) +
                                        (((d_angle_p2_dom_x[i]) ** 2) * ((all_peak_positions_p2[i+1][2]) ** 2)) +
                                        (((d_angle_p1_y[i]) ** 2) * ((all_peak_positions_p1[i+1][4]) ** 2)) +
                                        (((d_angle_p2_dom_y[i]) ** 2) * ((all_peak_positions_p2[i+1][4]) ** 2)) +
                                        (((d_angle_c_x[i]) ** 2) * ((centerXY[1][2]) ** 2)) +
                                        (((d_angle_c_y[i]) ** 2) * ((centerXY[1][4]) ** 2)))
        p1p2_angle_dom_error.append(hilf_angle_dom_error)

    p1p2_angle_rec_error = []       # recessive
    for i in range(len(vector_p1)):
        hilf_angle_error = [all_peak_positions_p1[i+1][0], 0]
        hilf_angle_error[1] = math.sqrt((((d_angle_p1_x[i]) ** 2) * ((all_peak_positions_p1[i+1][2]) ** 2)) +
                                        (((d_angle_p2_rec_x[i]) ** 2) * ((all_peak_positions_p2[i+1][6]) ** 2)) +
                                        (((d_angle_p1_y[i]) ** 2) * ((all_peak_positions_p1[i+1][4]) ** 2)) +
                                        (((d_angle_p2_rec_y[i]) ** 2) * ((all_peak_positions_p2[i+1][8]) ** 2)) +
                                        (((d_angle_c_x[i]) ** 2) * ((centerXY[1][2]) ** 2)) +
                                        (((d_angle_c_y[i]) ** 2) * ((centerXY[1][4]) ** 2)))
        p1p2_angle_rec_error.append(hilf_angle_error)

    # for a calculation of nu_+ and nu_- we also calculate the individual phases

    p1_angle = []   # corresponds to nu_-
    for i in range(len(vector_p1)):
        hilf_angle_p1 = [all_peak_positions_p1[i+1][0], 0]
        hilf_angle_p1[1] = math.pi + math.atan2(vector_p1[i][1], vector_p1[i][0])
        p1_angle.append(hilf_angle_p1)
    p2_angle_dom = []   # corresponds to nu_+
    for i in range(len(vector_p2_dom)):
        hilf_angle_p2 = [all_peak_positions_p1[i+1][0], 0]
        hilf_angle_p2[1] = math.pi - math.atan2(vector_p2_dom[i][1], vector_p2_dom[i][0])
        p2_angle_dom.append(hilf_angle_p2)
    p2_angle_rec = []   # corresponds to nu_+
    for i in range(len(vector_p2_rec)):
        hilf_angle_p2 = [all_peak_positions_p1[i+1][0], 0]
        hilf_angle_p2[1] = math.pi - math.atan2(vector_p2_rec[i][1], vector_p2_rec[i][0])
        p2_angle_rec.append(hilf_angle_p2)

    p1_angle_error = []
    for i in range(len(vector_p1)):
        hilf_angle_error_p1 = [all_peak_positions_p1[i+1][0], 0]
        hilf_angle_error_p1[1] = math.sqrt((((d_angle_p1_x[i]) ** 2) * ((all_peak_positions_p1[i+1][2]) ** 2)) +
                                        (((d_angle_p2_dom_x[i]) ** 2) * ((all_peak_positions_p2[i+1][6]) ** 2)) +       # [6] due to higher unc. in rec
                                        (((d_angle_c_x[i]) ** 2) * ((centerXY[1][2]) ** 2)))
        p1_angle_error.append(hilf_angle_error_p1)
    p2_angle_dom_error = []
    for i in range(len(vector_p2_dom)):
        hilf_angle_error_p2 = [all_peak_positions_p1[i+1][0], 0]
        hilf_angle_error_p2[1] = math.sqrt((((d_angle_p1_y[i]) ** 2) * ((all_peak_positions_p1[i+1][4]) ** 2)) +
                                        (((d_angle_p2_dom_y[i]) ** 2) * ((all_peak_positions_p2[i+1][4]) ** 2)) +
                                        (((d_angle_c_y[i]) ** 2) * ((centerXY[1][4]) ** 2)))
        p2_angle_dom_error.append(hilf_angle_error_p2)
    p2_angle_rec_error = []
    for i in range(len(vector_p2_rec)):
        hilf_angle_error_p2 = [all_peak_positions_p1[i+1][0], 0]
        hilf_angle_error_p2[1] = math.sqrt((((d_angle_p1_y[i]) ** 2) * ((all_peak_positions_p1[i+1][4]) ** 2)) +
                                        (((d_angle_p2_rec_y[i]) ** 2) * ((all_peak_positions_p2[i+1][8]) ** 2)) +
                                        (((d_angle_c_y[i]) ** 2) * ((centerXY[1][4]) ** 2)))
        p2_angle_rec_error.append(hilf_angle_error_p2)

    angle_info = [['Name', 'P1-P2-angle (dominant) / radian', 'P1-P2-angle-unc.', 'P1-P2-angle (recessive) / radian', 'P1-P2-angle-unc.', 'P1-angle', 'P1-angle-unc.', 'P2-angle (dominant)', 'P2-angle-unc.', 'P2-angle (recessive)', 'P2-angle-unc.']]
    for i in range(len(p1p2_angle_dom)):
        angle_info.append([p1p2_angle_dom[i][0], p1p2_angle_dom[i][1], p1p2_angle_dom_error[i][1], p1p2_angle_rec[i][1], p1p2_angle_rec_error[i][1], p1_angle[i][1], p1_angle_error[i][1], p2_angle_dom[i][1], p2_angle_dom_error[i][1], p2_angle_rec[i][1], p2_angle_rec_error[i][1]])

    write_csv(upper_folder_path, angle_info, 'angle_info_%s' % folder_name_0)

    return(p1p2_angle_dom, p1p2_angle_dom_error, p1p2_angle_rec, p1p2_angle_rec_error, p1_angle, p1_angle_error, p2_angle_dom, p2_angle_dom_error, p2_angle_rec, p2_angle_rec_error)




def piicr_cyclotron_frequency(p1p2_angle, p1p2_angle_error, p1_angle, p2_angle, p1_angle_error, p2_angle_error, red_cyc_freq, mag_freq, cyc_acc_time, corrected_cyc_acc_time):
    """
    The function calculates PI-ICR cyclotron frequency.

    Times are in micro seconds, frequencies in Hertz.
    """
    n_plus_rounds = []
    n_minus_rounds = []
    if isinstance(red_cyc_freq, (list, )):
        for i in range(len(cyc_acc_time)):
            #   calculated rounds after n+ excitation:
            hilf1 = math.floor(corrected_cyc_acc_time[i] * 0.000001 * red_cyc_freq[i])
            n_plus_rounds.append(hilf1)
            #   calculated rounds after n- excitation:
            hilf2 = math.floor(corrected_cyc_acc_time[i] * 0.000001 * mag_freq[i])
            n_minus_rounds.append(hilf2)
        #   freqency calculation
        piicr_cyc_freq = []
        for i in range(len(p1p2_angle)):
            hilf_cyc_freq = [p1p2_angle[i][0], 0, 0]
            hilf_cyc_freq[1] = (p1p2_angle[i][1] + 2 * math.pi * (n_plus_rounds[i] + n_minus_rounds[i])) / (2 * math.pi * corrected_cyc_acc_time[i] * 0.000001)
            hilf_cyc_freq[2] = p1p2_angle_error[i][1] / (2 * math.pi * cyc_acc_time[i] * 0.000001)
            piicr_cyc_freq.append(hilf_cyc_freq)
        #   nu_- calculation
        piicr_mag_freq = []
        for i in range(len(p1_angle)):
            hilf_mag_freq = [p1_angle[i][0], 0, 0]
            hilf_mag_freq[1] = (p1_angle[i][1] + 2 * math.pi * n_minus_rounds[i]) / (2 * math.pi * corrected_cyc_acc_time[i] * 0.000001)
            hilf_mag_freq[2] = p1_angle_error[i][1] / (2 * math.pi * cyc_acc_time[i] * 0.000001)
            piicr_mag_freq.append(hilf_mag_freq)
        #   nu_+ calculation
        piicr_red_cyc_freq = []
        for i in range(len(p2_angle)):
            hilf_red_cyc_freq = [p2_angle[i][0], 0, 0]
            hilf_red_cyc_freq[1] = (p2_angle[i][1] + 2 * math.pi * n_plus_rounds[i]) / (2 * math.pi * corrected_cyc_acc_time[i] * 0.000001)
            hilf_red_cyc_freq[2] = p2_angle_error[i][1] / (2 * math.pi * cyc_acc_time[i] * 0.000001)
            piicr_red_cyc_freq.append(hilf_red_cyc_freq)

    else:
        #   calculated rounds after n+ excitation:
        n_plus_rounds = math.floor(float(corrected_cyc_acc_time) * 0.000001 * float(red_cyc_freq))
        #   calculated rounds after n- excitation:
        n_minus_rounds = math.floor(float(corrected_cyc_acc_time) * 0.000001 * float(mag_freq))
    #   freqency calculation
        piicr_cyc_freq = []
        for i in range(len(p1p2_angle)):
            hilf_cyc_freq = [p1p2_angle[i][0], 0, 0]
            hilf_cyc_freq[1] = (p1p2_angle[i][1] + 2 * math.pi * (n_plus_rounds + n_minus_rounds)) / (2 * math.pi * float(corrected_cyc_acc_time) * 0.000001)
            hilf_cyc_freq[2] = p1p2_angle_error[i][1] / (2 * math.pi * float(cyc_acc_time) * 0.000001)
            piicr_cyc_freq.append(hilf_cyc_freq)
        #   nu_- calculation
        piicr_mag_freq = [[u'run #', u'Magnetron freq. / Hz', u'err']]
        for i in range(len(p1_angle)):
            hilf_mag_freq = [i, 0, 0] # [p1_angle[i][0], 0, 0]
            hilf_mag_freq[1] = (p1_angle[i][1] + 2 * math.pi * n_minus_rounds) / (2 * math.pi * float(corrected_cyc_acc_time) * 0.000001)
            hilf_mag_freq[2] = p1_angle_error[i][1] / (2 * math.pi * float(cyc_acc_time) * 0.000001)
            piicr_mag_freq.append(hilf_mag_freq)
        #   nu_+ calculation
        piicr_red_cyc_freq = [[u'run #', u'Red. cyclotron freq. / Hz', u'err']]
        for i in range(len(p2_angle)):
            hilf_red_cyc_freq = [i, 0, 0] # [p2_angle[i][0], 0, 0]
            hilf_red_cyc_freq[1] = (p2_angle[i][1] + 2 * math.pi * n_plus_rounds) / (2 * math.pi * float(corrected_cyc_acc_time) * 0.000001)
            hilf_red_cyc_freq[2] = p2_angle_error[i][1] / (2 * math.pi * float(cyc_acc_time) * 0.000001)
            piicr_red_cyc_freq.append(hilf_red_cyc_freq)
    return(piicr_cyc_freq, piicr_mag_freq, piicr_red_cyc_freq)


def piicr_cyclotron_frequency_2(p1p2_angle, p1p2_angle_error, cyc_acc_time, n_acc):

    piicr_cyc_freq = []
    for i in range(len(p1p2_angle)):
        piicr_cyc_freq.append([p1p2_angle[i][0],
                              ((float(p1p2_angle[i][1]) + 2 * math.pi * (float(n_acc) - 1)) / (2 * math.pi * float(cyc_acc_time) * 0.000001)),
                              (float(p1p2_angle_error[i][1]) / (2 * math.pi * float(cyc_acc_time) * 0.000001))])
    return(piicr_cyc_freq)


def piicr_cyclotron_frequency_two(p1p2_angle_1, p1p2_angle_1_error, p1p2_angle_2, p1p2_angle_2_error, cyc_acc_time, n_acc):

    piicr_cyc_freq = []
    for i in range(len(p1p2_angle_1)):
        piicr_cyc_freq.append([p1p2_angle_1[i][0],
                              ((float(p1p2_angle_1[i][1]) + 2 * math.pi * (float(n_acc) - 1)) / (2 * math.pi * float(cyc_acc_time) * 0.000001)),
                              (float(p1p2_angle_1_error[i][1]) / (2 * math.pi * float(cyc_acc_time) * 0.000001)),
                              ((float(p1p2_angle_2[i][1]) + 2 * math.pi * (float(n_acc) - 1)) / (2 * math.pi * float(cyc_acc_time) * 0.000001)),
                              (float(p1p2_angle_2_error[i][1]) / (2 * math.pi * float(cyc_acc_time) * 0.000001))
                              ])
    return(piicr_cyc_freq)



def get_FWHM_parameter(all_peak_pos_c_file_path):
    '''Function to calculate mean FWHM for parameter fixing'''
    if platform.system() == 'Windows':
        all_peak_pos_p1_file_path = os.path.dirname(os.path.dirname(all_peak_pos_c_file_path)) + '\\p1p2\\_all_peak_positions_p1.csv'
        all_peak_pos_p2_file_path = os.path.dirname(os.path.dirname(all_peak_pos_c_file_path)) + '\\p1p2\\_all_peak_positions_p2.csv'
    elif platform.system() == 'Darwin': # MacOS
        all_peak_pos_p1_file_path = os.path.dirname(os.path.dirname(all_peak_pos_c_file_path)) + '/p1p2/_all_peak_positions_p1.csv'
        all_peak_pos_p2_file_path = os.path.dirname(os.path.dirname(all_peak_pos_c_file_path)) + '/p1p2/_all_peak_positions_p2.csv'

    all_peak_pos_file_paths = [all_peak_pos_c_file_path, all_peak_pos_p1_file_path, all_peak_pos_p2_file_path]
    all_peak_pos = {}

    for i in all_peak_pos_file_paths:
        all_peak_pos['{}'.format(os.path.splitext(os.path.basename(i))[0])] = pd.read_csv(i)

    X_FWHM_c = all_peak_pos['_all_peak_positions_c']['X-FWHM'].tolist()
    X_FWHM_unc_c = all_peak_pos['_all_peak_positions_c']['X-FWHM-error'].tolist()
    X_FWHM_p1 = all_peak_pos['_all_peak_positions_p1']['X-FWHM'].tolist()
    X_FWHM_unc_p1 = all_peak_pos['_all_peak_positions_p1']['X-FWHM-error'].tolist()
    X_FWHM_p2 = all_peak_pos['_all_peak_positions_p2']['X-FWHM'].tolist()
    X_FWHM_unc_p2 = all_peak_pos['_all_peak_positions_p2']['X-FWHM-error'].tolist()
    Y_FWHM_c = all_peak_pos['_all_peak_positions_c']['Y-FWHM'].tolist()
    Y_FWHM_unc_c = all_peak_pos['_all_peak_positions_c']['Y-FWHM-error'].tolist()
    Y_FWHM_p1 = all_peak_pos['_all_peak_positions_p1']['Y-FWHM'].tolist()
    Y_FWHM_unc_p1 = all_peak_pos['_all_peak_positions_p1']['Y-FWHM-error'].tolist()
    Y_FWHM_p2 = all_peak_pos['_all_peak_positions_p2']['Y-FWHM'].tolist()
    Y_FWHM_unc_p2 = all_peak_pos['_all_peak_positions_p2']['Y-FWHM-error'].tolist()

    list_fwhm_x_c = [['X-FWHM c', 'X FWHM-error c'], [X_FWHM_c[0], X_FWHM_unc_c[0]]]
    list_fwhm_y_c = [['Y-FWHM c', 'Y FWHM-error c'], [Y_FWHM_c[0], Y_FWHM_unc_c[0]]]

    list_fwhm_x_p1 = [['X-FWHM p1', 'X-FWHM-error p1']]
    list_fwhm_x_p2 = [['X-FWHM p2', 'X-FWHM-error p2']]
    list_fwhm_y_p1 = [['Y-FWHM p1', 'Y-FWHM-error p1']]
    list_fwhm_y_p2 = [['Y-FWHM p2', 'Y-FWHM-error p2']]

    for i in range(len(X_FWHM_p1)):
        list_fwhm_x_p1.append([X_FWHM_p1[i], X_FWHM_unc_p1[i]])
        list_fwhm_x_p2.append([X_FWHM_p2[i], X_FWHM_unc_p2[i]])
        list_fwhm_y_p1.append([Y_FWHM_p1[i], Y_FWHM_unc_p1[i]])
        list_fwhm_y_p2.append([Y_FWHM_p2[i], Y_FWHM_unc_p2[i]])

    sigma_x_c = weighted_mean(list_fwhm_x_c)[0]
    sigma_x_p1 = weighted_mean(list_fwhm_x_p1)[0]
    sigma_x_p2 = weighted_mean(list_fwhm_x_p2)[0]
    sigma_y_c = weighted_mean(list_fwhm_y_c)[0]
    sigma_y_p1 = weighted_mean(list_fwhm_y_p1)[0]
    sigma_y_p2 = weighted_mean(list_fwhm_y_p2)[0]

    FWHM_parameter_fix_dict = {'c': {'fixed': True,
                                     'sigmax': sigma_x_c,
                                     'sigmay': sigma_y_c,
                                     'original data x': list_fwhm_x_c,
                                     'original data y': list_fwhm_y_c},
                               'p1': {'fixed': True,
                                      'sigmax': sigma_x_p1,
                                      'sigmay': sigma_y_p1,
                                      'original data x': list_fwhm_x_p1,
                                      'original data y': list_fwhm_y_p1},
                               'p2': {'fixed': True,
                                      'sigmax': sigma_x_p2,
                                      'sigmay': sigma_y_p2,
                                      'original data x': list_fwhm_x_p2,
                                      'original data y': list_fwhm_y_p2},
                               'isotope': os.path.basename(os.path.dirname(os.path.dirname(all_peak_pos_c_file_path)))}

    return(FWHM_parameter_fix_dict)


def analysis_main(upper_folder_path, fit_range, tof_cut, z_class_analysis, FWHM_parameter_fix, isomer_analysis):
    '''The function analysis PI-ICR raw data and calculates the frequencies.'''

    # ---------------------------------------------------------------------------
    # ADVANCED SETTINGS -- Only change if you know what you do!!!
    # ---------------------------------------------------------------------------

    axis_scaling_factor = 0.031746032
    dont_show = 1           # if =1 the plots are not displayed, for showing =0! file_name = '133Cs_c_091'  # if only one file should be looked at, uncomment! if peak_finder_plot = 1 this number will automatically set to 0 since they don't work together.
    pattern = 1             # chose if you want pattern one (='1'; p1+p2 in file) or pattern two (='2'; p1+p2 in file) should be calculated. If only one pattern is saved in the file (e.g. the center) just put a '3' there
    sigma = 10              # changes the measurement window; eqal the number of sigma of the gaussian distribution calculated for the totaltime
    bins = 60  # 85         # changes the bin number per axis of the histogram plot
    cut_counter = 1         # adds a nametag to certain files when you enable cutting
    trigger_on_off = 1      # when triggered (-> MCP time != 0) then trigger = 1, else set it to zero (e.g. for darkcount tests)
    clear_one = 0           # when there are engough events with only 1 event at each delayline set it to 1, elso to zero (if 1 time window gets calculated)
    peak_finder_plot = 0    # when neither 1 nor 3 peaks for both could be found an error will appear (e.g. in line 485 no parameter found). To better check and correct for that (adjust the peakfinder function in width or in height) set peak_finder_plot to one and you will get a plot with the found peaks. Normally you don't need that, so set it then to 0
    all_peak_positions = [] # list to save all fit positions. Ordering: file_name, X-Fit-Pos, X-Fit-Pos-err, Y-Fit-Pos, Y-Fit-Pos-Err
    hist2d_spines_off = 1   # if you want the axis lines for the 2d Histogram to be switched off set this variable to '1', if not set it to '0'
    hist2d_min_max = 750    # this variable sets the +- axis range for the 2d Histogram. Set 750 for the whole MCP
    simple_plot_on = 0      # if you want an additional plot that dispays the whole detector area set it to '1', otherwise to '0'
    surface3D_plot_on = 0   # if you want a 3D plot set it to '1', otherwise to '0'
    trigger_splitter = 0    # if = '0' it will go through all files in one calculation. if = '1' it will cut the raw data is portions of #trigger_number and calculates for each
    trigger_numbers = 1000  # if trigger_splitter = '1' it will cut the raw data in portions of #trigger_numbers trigger.
    color_map = 2           # 1 for blue, 2 for multicolor
    plot_type = 'scatter'   # 'scatter' or 'histogram'
    file_or_folder = 'folder'               # in case you just want to analyze only one single file e.g. for cutting an isomer out
    file_name_input = '123Cd-p2-127.txt'    # if file_or_folder == 'file' input here the fil name WITH the ending .txt
    modes = ['c', 'p1', 'p2']
    analysis = 'new'        # 'old' ==  without time info in the data and using the mpant and tmp files; 'new' == with time info in the data and using the xml files
    path, folder_name_0 = os.path.split(upper_folder_path)  # saves the folder name
    latex_name = '$^{%s}$%s' % (filter(str.isdigit, folder_name_0), filter(str.isalpha, folder_name_0))

    if isomer_analysis:
        one_or_two_peaks_to_fit = 'two'
    else:
        one_or_two_peaks_to_fit = 'one'
    # ---------------------------------------------------------------------------
    # DON'T CHANGE
    # ---------------------------------------------------------------------------

    if tof_cut[0] == '' and tof_cut[1] == '':
        tof_cut_on = 0  # DON'T CHANGE HERE        # if you want a TOF cut with TOF_min/max, set it to '1', if you want to see the whole data and the whole TOF range set it to '0'
        tof_min = 58    # DON'T CHANGE HERE        # only valid if tof_cut_on == 1! unit is microseconds
        tof_max = 65    # DON'T CHANGE HERE        # only valid if tof_cut_on == 1! unit is microseconds
    else:
        tof_cut_on = 1  # DON'T CHANGE HERE        # if you want a TOF cut with TOF_min/max, set it to '1', if you want to see the whole data and the whole TOF range set it to '0'
        tof_min = float(tof_cut[0])    # DON'T CHANGE HERE        # only valid if tof_cut_on == 1! unit is microseconds
        tof_max = float(tof_cut[1])    # DON'T CHANGE HERE        # only valid if tof_cut_on == 1! unit is microseconds

    if fit_range[0] == '' and fit_range[1] == '' and fit_range[2] == '' and fit_range[3] == '':
        fit_range = [-750, 750, -750, 750]
    else:
        fit_range = [float(fit_range[0]), float(fit_range[1]), float(fit_range[2]), float(fit_range[3])]

    if z_class_analysis[0] != '' and z_class_analysis[1] != '':
        z_class_analysis = ['yes', float(z_class_analysis[0]), float(z_class_analysis[1])]
    # ---------------------------------------------------------------------------
    # Get peak-positions for c,p1,p2
    # ---------------------------------------------------------------------------

    all_peak_positions_c = []
    all_peak_positions_p1 = []
    all_peak_positions_p2 = []

    for mode in modes:          # get the peak positions either from a pre-calculated file or if it doesn't exist, calculate new
        all_peak_positions = []  # reset variable
        if mode == 'c':
            if platform.system() == 'Windows':
                data_folder_path = '%s\\%s' % (upper_folder_path, mode)  # Macintosh /Volumes/dfs/DATA/2017/2017_Feb-PI-ICR/2017-02-cross-checks/P1-P2/  # Windows C://Users/jkarthei/cernbox/Python/test2      #  G://Experiments/ISOLTRAP/DATA/2017/2017_Feb-PI-ICR/2017-02-long-run
            elif platform.system() == 'Darwin': # MacOS
                data_folder_path = '%s/%s' % (upper_folder_path, mode)  # Macintosh /Volumes/dfs/DATA/2017/2017_Feb-PI-ICR/2017-02-cross-checks/P1-P2/  # Windows C://Users/jkarthei/cernbox/Python/test2      #  G://Experiments/ISOLTRAP/DATA/2017/2017_Feb-PI-ICR/2017-02-long-run
        else:
            if platform.system() == 'Windows':
                data_folder_path = '%s\\p1p2' % (upper_folder_path)  # Macintosh /Volumes/dfs/DATA/2017/2017_Feb-PI-ICR/2017-02-cross-checks/P1-P2/  # Windows C://Users/jkarthei/cernbox/Python/test2      #  G://Experiments/ISOLTRAP/DATA/2017/2017_Feb-PI-ICR/2017-02-long-run
            elif platform.system() == 'Darwin': # MacOS
                data_folder_path = '%s/p1p2' % (upper_folder_path)  # Macintosh /Volumes/dfs/DATA/2017/2017_Feb-PI-ICR/2017-02-cross-checks/P1-P2/  # Windows C://Users/jkarthei/cernbox/Python/test2      #  G://Experiments/ISOLTRAP/DATA/2017/2017_Feb-PI-ICR/2017-02-long-run
        path, folder_name = os.path.split(data_folder_path)  # saves the folder name
        os.chdir(data_folder_path)

        if os.path.isfile('_all_peak_positions_%s.csv' % (mode)):
            if mode == 'c':
                all_peak_positions_c = read_csv(data_folder_path, '_all_peak_positions_%s' % (mode))
                all_peak_positions_c[1:] = [c[:1] + map(float, c[1:]) for c in all_peak_positions_c[1:]]  # convert numbers from string into float
                time_info = read_csv(data_folder_path, '_time_info_c')
            elif mode == 'p1':
                all_peak_positions_p1 = read_csv(data_folder_path, '_all_peak_positions_%s' % (mode))
                all_peak_positions_p1[1:] = [c[:1] + map(float, c[1:]) for c in all_peak_positions_p1[1:]]  # convert numbers from string into float
                time_info = read_csv(data_folder_path, '_time_info_p1p2')
            elif mode == 'p2':
                all_peak_positions_p2 = read_csv(data_folder_path, '_all_peak_positions_%s' % (mode))
                all_peak_positions_p2[1:] = [c[:1] + map(float, c[1:]) for c in all_peak_positions_p2[1:]]  # convert numbers from string into float
                time_info = read_csv(data_folder_path, '_time_info_p1p2')
        else:
            if mode == 'c':
                all_peak_positions_c = piicr_daq(data_folder_path, axis_scaling_factor, dont_show, 3, sigma, bins, [-750, 750, -750, 750], cut_counter, tof_cut_on, tof_min, tof_max, trigger_on_off, clear_one, peak_finder_plot, all_peak_positions, hist2d_spines_off, hist2d_min_max, simple_plot_on, surface3D_plot_on, 0, trigger_numbers, file_or_folder, file_name_input, color_map, analysis, z_class_analysis, FWHM_parameter_fix['c'], 'one')
                time_info = read_csv(data_folder_path, '_time_info_c')
            elif mode == 'p1':
                all_peak_positions_p1 = piicr_daq(data_folder_path, axis_scaling_factor, dont_show, 1, sigma, bins, [-750, 750, -750, 750], cut_counter, tof_cut_on, tof_min, tof_max, trigger_on_off, clear_one, peak_finder_plot, all_peak_positions, hist2d_spines_off, hist2d_min_max, simple_plot_on, surface3D_plot_on, trigger_splitter, trigger_numbers, file_or_folder, file_name_input, color_map, analysis, z_class_analysis, FWHM_parameter_fix['p1'], 'one')
                time_info = read_csv(data_folder_path, '_time_info_p1p2')
            elif mode == 'p2':
                all_peak_positions_p2 = piicr_daq(data_folder_path, axis_scaling_factor, dont_show, 2, sigma, bins, fit_range, cut_counter, tof_cut_on, tof_min, tof_max, trigger_on_off, clear_one, peak_finder_plot, all_peak_positions, hist2d_spines_off, hist2d_min_max, simple_plot_on, surface3D_plot_on, trigger_splitter, trigger_numbers, file_or_folder, file_name_input, color_map, analysis, z_class_analysis, FWHM_parameter_fix['p2'], one_or_two_peaks_to_fit)
                time_info = read_csv(data_folder_path, '_time_info_p1p2')
                print '\n\n\nall preak positions ::\n'
                print all_peak_positions_p2
                print '\n\n\n'

    # ---------------------------------------------------------------------------
    # Get excitation frequencies and phase accumulation time
    # ---------------------------------------------------------------------------

    if platform.system() == 'Windows':
        file_path_1 = '%s\\freq_config' % (upper_folder_path)  # Macintosh /Volumes/dfs/DATA/2017/2017_Feb-PI-ICR/2017-02-cross-checks/P1-P2/  # Windows C://Users/jkarthei/cernbox/Python/test2      #  G://Experiments/ISOLTRAP/DATA/2017/2017_Feb-PI-ICR/2017-02-long-run
    elif platform.system() == 'Darwin': # MacOS
        file_path_1 = '%s/freq_config' % (upper_folder_path)  # Macintosh /Volumes/dfs/DATA/2017/2017_Feb-PI-ICR/2017-02-cross-checks/P1-P2/  # Windows C://Users/jkarthei/cernbox/Python/test2      #  G://Experiments/ISOLTRAP/DATA/2017/2017_Feb-PI-ICR/2017-02-long-run

    path, folder_name_1 = os.path.split(data_folder_path)  # saves the folder name
    os.chdir(file_path_1)

    if os.path.isfile('piicr_excitation.csv'):
        piicr_excitation = read_csv(file_path_1, 'piicr_excitation')
    else:
        piicr_excitation = read_freq_config(file_path_1, folder_name_0)

    # ---------------------------------------------------------------------------
    # Calculate angle and piicr frequencies
    # ---------------------------------------------------------------------------

    os.chdir(upper_folder_path)
    if len(all_peak_positions_p2[1]) == 10: # one spot
        p1p2_angle, p1p2_angle_error, p1_angle, p2_angle, p1_angle_error, p2_angle_error = angle_calculator(all_peak_positions_c, all_peak_positions_p1, all_peak_positions_p2, folder_name_0, upper_folder_path, folder_name_0, latex_name)
        # piicr_cyc_freq, piicr_mag_freq, piicr_red_cyc_freq = piicr_cyclotron_frequency(p1p2_angle, p1p2_angle_error, p1_angle, p2_angle, p1_angle_error, p2_angle_error, piicr_excitation[1][1], piicr_excitation[1][2], piicr_excitation[1][3], piicr_excitation[1][3])
        piicr_cyc_freq = piicr_cyclotron_frequency_2(p1p2_angle, p1p2_angle_error, piicr_excitation[1][3], piicr_excitation[1][4])


    elif len(all_peak_positions_p2[0]) == 14:   # two spots, e.g. isomers
        p1p2_angle_1, p1p2_angle_1_error, p1p2_angle_2, p1p2_angle_2_error, p1_angle, p1_angle_error, p2_angle_1, p2_angle_1_error, p2_angle_2, p2_angle_2_error = angle_calculator_two(all_peak_positions_c, all_peak_positions_p1, all_peak_positions_p2, folder_name_0, upper_folder_path, folder_name_0, latex_name)
        piicr_cyc_freq = piicr_cyclotron_frequency_two(p1p2_angle_1, p1p2_angle_1_error, p1p2_angle_2, p1p2_angle_2_error, piicr_excitation[1][3], piicr_excitation[1][4])

    # ---------------------------------------------------------------------------
    # Prepare plotting
    # ---------------------------------------------------------------------------

    if len(piicr_cyc_freq[0]) == 3:
        piicr_cyc_freq_new = [[u'run #', u'cyclotron freq. / Hz', u'err']]
        piicr_cyc_freq_new_t = [[u'time', u'cyclotron freq. / Hz', u'err']]

        temp_time_info = ''
        # print piicr_cyc_freq
        # print time_info
        for i in range(len(piicr_cyc_freq)):
            hilf6 = [i, piicr_cyc_freq[i][1], piicr_cyc_freq[i][2]]
            piicr_cyc_freq_new.append(hilf6)
            # for j in time_info:
            #     print '###############'
            #     print j
            #     if str('_'.join(piicr_cyc_freq[i][0][:9].split('_', 2)[:2])) == j[0]:
            #         temp_time_info = j
            #         print temp_time_info
            # print piicr_cyc_freq[i][1], piicr_cyc_freq[i][2], i
            # hilf17 = [datetime.datetime.strptime(temp_time_info[1]+' '+temp_time_info[2], '%m/%d/%Y %H:%M:%S'), piicr_cyc_freq[i][1], piicr_cyc_freq[i][2]]
            # piicr_cyc_freq_new_t.append(hilf17)
        python_plot(piicr_cyc_freq_new, '$\\nu_{\mathrm{c}}$ analysis for %s' % (latex_name), '%s' % (folder_name_0), 'no', '%s cyc. freq.' % (latex_name), 'no', 'linear', 'Utopia', 'red', 'green', 'full', 1, 5, 1, 'scatter', 20, 3, 'on', 'off', '', '', '')
        # python_plot(piicr_cyc_freq_new_t, '%s analysis' % (latex_name), '%s_time' % (folder_name_0), 'no', '%s cyc. freq.' % (latex_name), 'no', 'linear', 'Utopia', 'red', 'green', 'full', 1, 5, 1, 'scatter', 20, 3, 'on', 'off', '', '', '')
        # python_plot(piicr_mag_freq, '$\\nu_-$ analysis for %s' % (latex_name), '%s_mag_freq' % (folder_name_0), 'no', '%s mag. freq.' % (latex_name), 'no', 'linear', 'Utopia', 'red', 'green', 'full', 1, 5, 1, 'scatter', 20, 3, 'on', 'off', '', '', '')
        # python_plot(piicr_red_cyc_freq, '$\\nu_+$ analysis for %s' % (latex_name), '%s_red_cyc_freq' % (folder_name_0), 'no', '%s red. cyc. freq.' % (latex_name), 'no', 'linear', 'Utopia', 'red', 'green', 'full', 1, 5, 1, 'scatter', 20, 3, 'on', 'off', '', '', '')
        piicr_cyc_freq.insert(0, ['File', 'PI-ICR cyclotron frequency / Hz', 'Err / Hz'])
        write_csv(upper_folder_path, piicr_cyc_freq, 'piicr_cyc_freq_%s' % folder_name_0)
        # write_csv(upper_folder_path, piicr_mag_freq, 'piicr_mag_freq_%s' % folder_name_0)
        # write_csv(upper_folder_path, piicr_red_cyc_freq, 'piicr_red_cyc_freq_%s' % folder_name_0)

    elif len(piicr_cyc_freq[0]) == 5:
        piicr_cyc_freq_dom = [[u'run #', u'cyclotron freq. (dominant) / Hz', u'err']]
        piicr_cyc_freq_rec = [[u'run #', u'cyclotron freq. (recessive) / Hz', u'err']]

        temp_time_info = ''

        for i in range(len(piicr_cyc_freq)):
            hilf6 = [i, piicr_cyc_freq[i][1], piicr_cyc_freq[i][2]]
            piicr_cyc_freq_dom.append(hilf6)
            hilf7 = [i, piicr_cyc_freq[i][3], piicr_cyc_freq[i][4]]
            piicr_cyc_freq_rec.append(hilf7)

        python_plot(piicr_cyc_freq_dom, '$\\nu_{\mathrm{c}}$ analysis for %s (dominant)' % (latex_name), '%s' % (folder_name_0+'_dom'), 'no', '%s cyc. freq.' % (latex_name), 'no', 'linear', 'Utopia', 'red', 'green', 'full', 1, 5, 1, 'scatter', 20, 3, 'on', 'off', '', '', '')

        python_plot(piicr_cyc_freq_rec, '$\\nu_{\mathrm{c}}$ analysis for %s (recessive)' % (latex_name), '%s' % (folder_name_0+'_rec'), 'no', '%s cyc. freq.' % (latex_name), 'no', 'linear', 'Utopia', 'red', 'green', 'full', 1, 5, 1, 'scatter', 20, 3, 'on', 'off', '', '', '')

        piicr_cyc_freq.insert(0, ['File', 'PI-ICR cyclotron frequency (domminant) / Hz', 'Err / Hz', 'PI-ICR cyclotron frequency (recessive) / Hz', 'Err / Hz'])

        write_csv(upper_folder_path, piicr_cyc_freq, 'piicr_cyc_freq_%s' % folder_name_0)

        piicr_cyc_freq_ground = [[i[0], i[1], i[2]] for i in piicr_cyc_freq]
        piicr_cyc_freq_isomer = [[i[0], i[3], i[4]] for i in piicr_cyc_freq]

        write_csv(upper_folder_path, piicr_cyc_freq_ground, 'piicr_cyc_freq_%s_ground' % folder_name_0)
        write_csv(upper_folder_path, piicr_cyc_freq_isomer, 'piicr_cyc_freq_%s_isomer' % folder_name_0)

if __name__ == '__main__':
    # ---------------------------------------------------------------------------
    # Please change before use to the correct folder -- Instructions in REAM-ME-first
    # ---------------------------------------------------------------------------

    upper_folder_path = '/Volumes/Macintosh HD/Users/jonaskarthein/cernbox/Python/8587Rb_over_night_23.04.17_failed/85Rb'
    fit_range = [-750, 750, -750, 750]
    # ---------------------------------------------------------------------------
    # Run skript
    # ---------------------------------------------------------------------------

    # analysis_main(upper_folder_path, fit_range)



    all_peak_pos_file_paths = '/Users/jonaskarthein/cernbox/Python/Doktor/2017-11-02_Cd-analysis/PI-ICR/127Cd/127Cd_run-002/127mCd/c/_all_peak_positions_c.csv'
    get_FWHM_parameter(all_peak_pos_file_paths)
