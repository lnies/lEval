# ---------------------------------------------------------------------------
# Written by Jonas Karthein in 2016/2017. Questions to jonas.karthein@cern.ch
# ---------------------------------------------------------------------------

from read_write_functions import *
from python_plotter_functions import *
import math
import numpy as np
import matplotlib as mpl
mpl.use('Qt5Agg')
# import peakutils
import scipy.optimize
# from peakutils.plot import plot as pplot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import multiprocessing
import time
import ROOT as root
from collections import Counter
from array import array
import sys
import pandas as pd


def load_and_pattern(file_name, pattern, analysis_old_new):
    """
    Loading function.

    The function loads the given text file that is found in the
    folder (finding is not included here) and saves the data depending on the
    way the raw data was saved in the .txt (e.g. p1 and p2 alternating). In the
    example case it would only save p1 OR p2.
    """
    load_file = []
    with open('%s.txt' % file_name, 'rb') as infile:
        load_file = [[str(i) for i in line.strip().split()] for line in infile]

    if analysis_old_new == 'new':
        time_info_hilf = load_file[-1]
        load_file.pop(-1)   # delete time stamp
        load_file.pop(-1)   # when the time stamps are added, 5 empty lines are also added (we don't know why)
        load_file.pop(-1)
        load_file.pop(-1)
        load_file.pop(-1)
        load_file.pop(-1)
        time_info = [i for i in time_info_hilf]
        # print time_info
        if not len(load_file) % 10 == 0:
            load_file.pop(0)   # check if there was a bad triggering at the beginning
            load_file.pop(0)
            load_file.pop(0)
            load_file.pop(0)
            load_file.pop(0)
            print('File info                ::  Bad trigger found, first ejection will be deleted!')
    elif analysis_old_new == 'old':
        time_info = []


    signals = 0
    for index in range(0, len(load_file), 5):
        signals += max([float(load_file[index+1][0]), float(load_file[index+2][0]), float(load_file[index+3][0]), float(load_file[index+4][0])])
    counts_info = [len(load_file)/5, signals]


    load_file = [map(int, x) for x in load_file]
    print 'Length of input list     :: ', len(load_file)
    results = []
    if pattern == 1:
        for i in range(0, len(load_file), 10):
            for j in range(5):
                results.append(load_file[i+j])
        file_name = file_name.replace('p1p2_', '')
        file_name += str('_p1')
    elif pattern == 2:
        for i in range(0, len(load_file), 10):
            for j in range(5):
                results.append(load_file[(i+j+int(5))])
        file_name = file_name.replace('p1p2_', '')
        file_name += str('_p2')
    elif pattern == 3:
        for i in range(len(load_file)):
            # print len(load_file)
            results.append(load_file[i])
    else:
        for i in range(0, len(load_file), 5*pattern[1]):
            if (i+4+((int(pattern[0])-1)*5)) < len(load_file):
                for j in range(5):
                    results.append(load_file[i+j+((int(pattern[0])-1)*5)])
        file_name = file_name.replace('p1p2_', '')
        file_name += str('_p%s' % (pattern[0]))

    return(results, file_name, time_info, counts_info)


def mm_scale(au_scale, axis_scaling_factor):
    """The function converts the scale given in channels into millimeters."""
    return (au_scale*axis_scaling_factor)


def update_ax2(ax1, ax2):
    """The function updates the x axis from channels into millimeters."""
    y1, y2 = ax1.get_xlim()
    ax2.set_xlim(mm_scale(y1), mm_scale(y2))
    ax2.figure.canvas.draw()


def count_distribution(list_of_lists):
    '''returns list with sorted #events with 0, 1, 2, ... counts'''
    temp_df = pd.DataFrame(list_of_lists, columns=['ej', 'x', 'y', 'tof'])
    series = pd.cut(temp_df['ej'], bins=range(0,int(temp_df.max(0)['ej']+1))).value_counts()
    counts_bin = pd.cut(series, bins=range(-1,10+1), labels=[str(i) for i in range(10+1)]).value_counts().sort_index()
    return(counts_bin.tolist())


# ------------------------
# Measurement window + error: takes all 1 hit per event shots and counts (t_x1 + t_x2)-t_MCP and (t_y1 + t_y2)- t_MCP

# comment c1: there are always 4 hits on the mcp when there is only one at the detector on some other events there are also some...
# ...totally odd times like the first two events of MCP so I ignored them. So how I chose between the to other events? I tried my...
# ..."window" routine with both, the third and the fourth mcp count and calculated the the average detector time and the gaussian...
# ...errors: for the third event I get with the test file window_x=1514, window_y=1498, err_window_x=14, err_window_y=24.4 and for...
# ...both 23 counts, which sounds quite reasonable to me. But for the fourth MCP time I get with the test file window_x=538, ...
# ...window_y=522, err_window_x=976.7, err_window_y=976.8 and for both 23 counts. So this number doesn't stay constant in time, ...
# ...but changes quite drastically (look at the big "errors"). Therefore I set the value to the third MCP time. BUT: Further studies needed!

# when no measurement window can be calculated set clear_one to 0!
# ------------------------


def window(clear_one, results, window_x, window_y, window_xy, err_window_x, err_window_y, err_window_xy, data_folder_path, file_name):
    """
    Measurement window calculation.

    The function calculates the measurement window an theirs errors. It
    takes, when there are engough 1 hit events (normally not the case) and
    the variable 'clear_one' is set to '1', all 1 hit per event shots and
    counts (t_x1 + t_x2)-t_MCP and (t_y1 + t_y2)- t_MCP. If this is not the
    case it takes the given start values (round 1500+-25 for x&y and 3000+-30
    for xy) set by experience. In a second 'iteration' it calculates the
    proper average windows for the current file. This is then be used for the
    next one if there is one and so on. Since the acceptance of the window is
    very broad (sigma*err_window = +-5*25), this shouldn't matter so much.
    """
    import numpy
    import copy

    if clear_one == 1:
        window_x = 0
        window_y = 0
        counter = 0
        sum_x = [['SumX', 'Ejection#']]  # save all window_x info of each ejection
        sum_y = [['SumY', 'Ejection#']]  # save all window_y info of each ejection
        sum_x_y = [['SumX', 'SumY', 'Ejection#']]
        for i in range(0, len(results), 5):
            if results[i][0] == 1 and results[i+1][0] == 1 and results[i+2][0] == 1 and results[i+3][0] == 1:
                counter += 1
                window_x = window_x + ((results[i][1] + results[i+1][1]) - 2*results[i+4][1])           # look for comment c1 above
                window_y = window_y + ((results[i+2][1] + results[i+3][1]) - 2*results[i+4][1])         # look for comment c1 above
                sum_x.append([((results[i][1] + results[i+1][1]) - 2*results[i+4][1]), ((i+5)/5)])
                sum_y.append([((results[i+2][1] + results[i+3][1]) - 2*results[i+4][1]), ((i+5)/5)])
                sum_x_y.append([((results[i][1] + results[i+1][1]) - 2*results[i+4][1]), ((results[i+2][1] + results[i+3][1]) - 2*results[i+4][1]), ((i+5)/5)])
        if counter > 0:
            window_x = window_x/counter
            window_y = window_y/counter
        else:
            window_x = 1600 # 1516.0  # 1545
            window_y = 1600 # 1506.0  # 1533

        err_window_x = 0
        err_window_y = 0
        counter = 0

        for i in range(0, len(results), 5):
            if results[i][0] == 1 and results[i+1][0] == 1 and results[i+2][0] == 1 and results[i+3][0] == 1:
                counter += 1
                err_window_x += (((results[i][1] + results[i+1][1]) - 2*results[i+4][1]) - window_x)**2
                err_window_y += (((results[i+2][1] + results[i+3][1]) - 2*results[i+4][1]) - window_y)**2
        if counter > 0:
            err_window_x = numpy.sqrt(err_window_x/counter)
            err_window_y = numpy.sqrt(err_window_y/counter)
        else:
            err_window_x = 30.0
            err_window_y = 30.0

        window_xy = 0
        counter = 0

        for i in range(0, len(results), 5):
            if results[i][0] == 1 and results[i+1][0] == 1 and results[i+2][0] == 1 and results[i+3][0] == 1:
                counter += 1
                window_xy = window_xy + ((results[i][1] + results[i+1][1] + results[i+2][1] + results[i+3][1]) - 4*results[i+4][1])           # look for comment c1 above
        if counter > 0:
            window_xy = window_xy/counter
        else:
            window_xy = 3200.0

        err_window_xy = 0
        counter = 0

        for i in range(0, len(results), 5):
            if results[i][0] == 1 and results[i+1][0] == 1 and results[i+2][0] == 1 and results[i+3][0] == 1:
                counter += 1
                err_window_xy += (((results[i][1] + results[i+1][1] + results[i+2][1] + results[i+3][1]) - 4*results[i+4][1]) - window_xy)**2
        if counter > 0:
            err_window_xy = numpy.sqrt(err_window_xy/counter)
        else:
            err_window_xy = 40.0
        # print sum_x, sum_y, sum_x_y
        write_csv(data_folder_path, sum_x, '%s_sum_x' % (file_name))
        write_csv(data_folder_path, sum_y, '%s_sum_y' % (file_name))
        python_plot(sum_x, 'SumX', '%s_sum_x_hist' % (file_name), 'no', '', 'no', 'gauss', 'Utopia', 'red', 'green', 'full', 0, 5, 1, 'histogram', 60, 5, 'off', 'off', '', '', '')
        python_plot(sum_y, 'SumY', '%s_sum_y_hist' % (file_name), 'no', '', 'no', 'gauss', 'Utopia', 'red', 'green', 'full', 0, 5, 1, 'histogram', 60, 5, 'off', 'off', '', '', '')
        python_plot(sum_x_y, 'SumX-SumY', '%s_sum_x_y_hist' % (file_name), 'no', '', 'no', 'gauss', 'Utopia', 'red', 'green', 'full', 0, 5, 1, '2Dhistogram', 60, 5, 'off', 'off', '', '', '')

        return(window_x, window_y, window_xy, err_window_x, err_window_y, err_window_xy, sum_x, sum_y)

    elif clear_one == 0:
        if window_x == 0:
            window_x = 1600.0  # 1545
            window_y = 1600.0  # 1533
            window_xy = 3200.0  # 3079
            err_window_x = 30.0
            err_window_y = 30.0
            err_window_xy = 40.0
            sum_x = [['SumX', 'Ejection#']]  # save all window_x info of each ejection
            sum_y = [['SumY', 'Ejection#']]  # save all window_y info of each ejection
        else:
            window_x = copy.copy(window_x)
            window_y = copy.copy(window_y)
            window_xy = copy.copy(window_xy)
            err_window_x = copy.copy(err_window_x)
            err_window_y = copy.copy(err_window_y)
            err_window_xy = copy.copy(err_window_xy)
            sum_x = [['SumX', 'Ejection#']]  # save all window_x info of each ejection
            sum_y = [['SumY', 'Ejection#']]  # save all window_y info of each ejection

        return(window_x, window_y, window_xy, err_window_x, err_window_y, err_window_xy, sum_x, sum_y)


def save_fit_info(fit_range_save, pattern, folder_name, all_peak_positions, cut_counter, data_folder_path):
    """Saves all fit information into a file."""
    import csv

    if fit_range_save[0] == -750 and fit_range_save[1] == 750 and fit_range_save[2] == -750 and fit_range_save[3] == 750 and pattern == 1:
        write_csv(data_folder_path, all_peak_positions, '_all_peak_positions_p1')
    elif fit_range_save[0] == -750 and fit_range_save[1] == 750 and fit_range_save[2] == -750 and fit_range_save[3] == 750 and pattern == 2:
        write_csv(data_folder_path, all_peak_positions, '_all_peak_positions_p2')
    elif fit_range_save[0] == -750 and fit_range_save[1] == 750 and fit_range_save[2] == -750 and fit_range_save[3] == 750 and pattern == 3:
        write_csv(data_folder_path, all_peak_positions, '_all_peak_positions_%s' %(folder_name))
    else:
        if cut_counter == 1 and pattern == 1:
            write_csv(data_folder_path, all_peak_positions, '_all_peak_positions_p1_cut1')
        elif cut_counter == 1 and pattern == 2:
            write_csv(data_folder_path, all_peak_positions, '_all_peak_positions_p2_cut1')
        elif cut_counter == 1 and pattern == 3:
            write_csv(data_folder_path, all_peak_positions, '_all_peak_positions_%s_cut1' %(folder_name))
        elif cut_counter == 2 and pattern == 1:
            write_csv(data_folder_path, all_peak_positions, '_all_peak_positions_p1_cut2')
        elif cut_counter == 2 and pattern == 2:
            write_csv(data_folder_path, all_peak_positions, '_all_peak_positions_p2_cut2')
        elif cut_counter == 2 and pattern == 3:
            write_csv(data_folder_path, all_peak_positions, '_all_peak_positions_%s_cut2' %(folder_name))
        elif cut_counter == 3 and pattern == 1:
            write_csv(data_folder_path, all_peak_positions, '_all_peak_positions_p1_cut3')
        elif cut_counter == 3 and pattern == 2:
            write_csv(data_folder_path, all_peak_positions, '_all_peak_positions_p2_cut3')
        elif cut_counter == 3 and pattern == 3:
            write_csv(data_folder_path, all_peak_positions, '_all_peak_positions_%s_cut3' %(folder_name))


# ------------------------
# X and Y position: Detector-length * (t_x1 - t_MCP)/(t_x1 + t_x2 - 2*t_MCP) = X, same for Y. Detector-length...
# ...can be calibrated by darkcount measurement (MCP propertiers are known)
# ------------------------

def zerolistmaker(n):
    '''Creates a list of certain length filled with zeros.'''
    listofzeros = [0] * (n + 1)
    return listofzeros

def convert_cartesian_to_polar(x,y):
    '''Converts a list of x and y coordinates respectively into polar.'''
    r = []
    phi = []

    for i in range(len(x)):
        r.append(math.sqrt((x[i] ** 2 + y[i] ** 2)))
        phi.append(math.atan2(y[i], x[i]))

    return(r, phi)


def rec_0(input_list, q):
# def rec_0(results, trigger_on_off, tof_cut_on, tof_min, tof_max, window_x, window_y, err_window_x, err_window_y, window_xy, err_window_xy, sigma):
    """
    Event reconstruction function.

    The function reconstructs all possible events in one ejection with 1 time
    stamp in all 5 lines --> event only possible if times are within the
    measurement window(= on the detector) and 1 time per X1, X2, Y1, Y2 and MCP
    channel is found. If data is lost in one or more channels the event is not
    reconstructed.
    Input_list consists of the raw_data in the first entry and a dictionary with further values.
    """

    trigger_on_off = input_list[1]['trigger_on_off']
    tof_cut_on = input_list[1]['tof_cut_on']
    tof_min = input_list[1]['tof_min']  #/0.000025
    tof_max = input_list[1]['tof_max']  #/0.000025
    window_x = input_list[1]['window_x']
    window_y = input_list[1]['window_y']
    err_window_x = input_list[1]['err_window_x']
    err_window_y = input_list[1]['err_window_y']
    window_xy = input_list[1]['window_xy']
    err_window_xy = input_list[1]['err_window_xy']
    sigma = input_list[1]['sigma']

    results = input_list[0]

    cpu_core = input_list[2]
    # trigger_on_off = 1
    # tof_cut_on = 1
    # tof_min = 58/0.000025
    # tof_max = 68/0.000025
    # window_x = 1600
    # window_y = 1600
    # err_window_x = 30
    # err_window_y = 30
    # window_xy = 3200
    # err_window_xy = 50
    # sigma = 7


    spot_positions = []         # [ *ersetzen* for i in 10*range(len(results))]        # 1st column = number of ejection, 2nd column x position, 3rd column = y position
    short_save = [0., 0., 0.]   # nur ein Zwischenspeicher (buffer variable)
    short_save_2 = [1., 0.]     # nur ein Zwischenspeicher (buffer variable)
    check_event_taken = []
    double_taken_event_counter = 0
    temp_results = []



    if trigger_on_off == 1:
        for i in range(0, len(results), 5):                     # check each ejection
            check_event_taken.append(zerolistmaker(results[i][0]))      # Matrix with size of results to check if the event was used for a match
            check_event_taken.append(zerolistmaker(results[i+1][0]))
            check_event_taken.append(zerolistmaker(results[i+2][0]))
            check_event_taken.append(zerolistmaker(results[i+3][0]))
            check_event_taken.append(zerolistmaker(results[i+4][0]))
            if results[i][0] > 0 and results[i+1][0] > 0 and results[i+2][0] > 0 and results[i+3][0] > 0:     # make sure that there is at least one hit at the detector, then build every possible combination (even when a event is missing) and store all reasonable values
                for j in range(1, results[i][0]+1):                                      # x1
                    if results[i][j] > 250000:                                                 # excludes trigger-errors: saves time
                        for k in range(1, results[i+1][0]+1):                            # x2
                            if results[i+1][k] > 250000:                                       # excludes trigger-errors: saves time
                                for l in range(1, results[i+2][0]+1):                    # y1
                                    if results[i+2][l] > 250000:                               # excludes trigger-errors: saves time
                                        for m in range(1, results[i+3][0]+1):            # y2
                                            if results[i+3][m] > 250000:                       # excludes trigger-errors: saves time
                                                for n in range(1, results[i+4][0]+1):    # MCP
                                                    if results[i+4][n] > 250000:               # excludes trigger-errors: saves time
                                                        if (results[i][j]+results[i+1][k]-2*results[i+4][n]) < (window_x + sigma*err_window_x) and \
                                                            (results[i][j]+results[i+1][k]-2*results[i+4][n]) > (window_x - sigma*err_window_x) and \
                                                            (results[i+2][l]+results[i+3][m]-2*results[i+4][n]) < (window_y + sigma*err_window_y) and \
                                                            (results[i+2][l]+results[i+3][m]-2*results[i+4][n]) > (window_y - sigma*err_window_y) and \
                                                            (results[i][j] + results[i+1][k] + results[i+2][l] + results[i+3][m] - 4*results[i+4][n]) < (window_xy + sigma*err_window_xy) and \
                                                                (results[i][j] + results[i+1][k] + results[i+2][l] + results[i+3][m] - 4*results[i+4][n]) > (window_xy - sigma*err_window_xy):                 # check for x pairs and y pairs in their x/y-measurement window
                                                            short_save = [0., 0., 0., 0.]
                                                            short_save[0] = (float(i)+len(results)*(cpu_core-1))/5
                                                            short_save[1] = (float(results[i][j])-float(results[i+1][k]))/2
                                                            short_save[2] = (float(results[i+2][l])-float(results[i+3][m]))/2
                                                            short_save[3] = float(results[i+4][n])*0.000025    # ch to micro seconds
                                                            if short_save[1] < window_x and short_save[1] > -window_x and short_save[2] < window_y and short_save[2] > -window_y:              # check if the matched signal is also on the mcp
                                                                # if tof_cut_on == 1:
                                                                #     if short_save[3] > tof_min and short_save[3] < tof_max:
                                                                #         spot_positions.append(short_save)
                                                                #         q.put(short_save)
                                                                #         check_event_taken[i][j] += 1
                                                                #         check_event_taken[i+1][k] += 1
                                                                #         check_event_taken[i+2][l] += 1
                                                                #         check_event_taken[i+3][m] += 1
                                                                #         check_event_taken[i+4][n] += 1
                                                                #         short_save_2 = [1., 0.]                  # next 15 lines save an accepted event into temp_results so that the function window can calculate the time windows again
                                                                #         short_save_2[1] = results[i][j]
                                                                #         temp_results.append(short_save_2)
                                                                #         short_save_2 = [1., 0.]
                                                                #         short_save_2[1] = results[i+1][k]
                                                                #         temp_results.append(short_save_2)
                                                                #         short_save_2 = [1., 0.]
                                                                #         short_save_2[1] = results[i+2][l]
                                                                #         temp_results.append(short_save_2)
                                                                #         short_save_2 = [1., 0.]
                                                                #         short_save_2[1] = results[i+3][m]
                                                                #         temp_results.append(short_save_2)
                                                                #         short_save_2 = [1., 0.]
                                                                #         short_save_2[1] = results[i+4][n]
                                                                #         temp_results.append(short_save_2)
                                                                #         if check_event_taken[i][j] > 1 or check_event_taken[i+1][k] > 1 or check_event_taken[i+2][l] > 1 or check_event_taken[i+3][m] > 1 or check_event_taken[i+4][n] > 1:
                                                                #         # print 'Double counting in ejection %s' % ((i+1)/5)
                                                                #             double_taken_event_counter += 1
                                                                # else:
                                                                spot_positions.append(short_save)
                                                                q.put(short_save)
                                                                check_event_taken[i][j] += 1
                                                                check_event_taken[i+1][k] += 1
                                                                check_event_taken[i+2][l] += 1
                                                                check_event_taken[i+3][m] += 1
                                                                check_event_taken[i+4][n] += 1
                                                                short_save_2 = [1., 0.]                  # next 15 lines save an accepted event into temp_results so that the function window can calculate the time windows again
                                                                short_save_2[1] = results[i][j]
                                                                temp_results.append(short_save_2)
                                                                short_save_2 = [1., 0.]
                                                                short_save_2[1] = results[i+1][k]
                                                                temp_results.append(short_save_2)
                                                                short_save_2 = [1., 0.]
                                                                short_save_2[1] = results[i+2][l]
                                                                temp_results.append(short_save_2)
                                                                short_save_2 = [1., 0.]
                                                                short_save_2[1] = results[i+3][m]
                                                                temp_results.append(short_save_2)
                                                                short_save_2 = [1., 0.]
                                                                short_save_2[1] = results[i+4][n]
                                                                temp_results.append(short_save_2)
                                                                if check_event_taken[i][j] > 1 or check_event_taken[i+1][k] > 1 or check_event_taken[i+2][l] > 1 or check_event_taken[i+3][m] > 1 or check_event_taken[i+4][n] > 1:
                                                                # print 'Double counting in ejection %s' % ((i+1)/5)
                                                                    double_taken_event_counter += 1
    elif trigger_on_off == 0:                                           # darkcount routine
        for i in range(0, len(results), 5):                               # check each ejection
            if results[i][0] > 0 and results[i+1][0] > 0 and results[i+2][0] > 0 and results[i+3][0] > 0:     # make sure that there is at least one hit at the detector, then build every possible combination (even when a event is missing) and store all reasonable values
                for j in range(1, results[i][0]+1):                      # x1
                    for k in range(1, results[i+1][0]+1):                # x2
                        for l in range(1, results[i+2][0]+1):            # y1
                            for m in range(1, results[i+3][0]+1):        # y2
                                if (results[i][j]+results[i+1][k]) < (window_x + sigma*err_window_x) and (results[i][j]+results[i+1][k]) > (window_x - sigma*err_window_x) and (results[i+2][l]+results[i+3][m]) < (window_y + sigma*err_window_y) and (results[i+2][l]+results[i+3][m]) > (window_y - sigma*err_window_y) and (results[i][j] + results[i+1][k] + results[i+2][l] + results[i+3][m]) < (window_xy + sigma*err_window_xy) and (results[i][j] + results[i+1][k] + results[i+2][l] + results[i+3][m]) > (window_xy - sigma*err_window_xy):                 # check for x pairs and y pairs in their x/y-measurement window
                                    short_save = [0., 0., 0.]
                                    short_save[0] = i
                                    short_save[1] = (float(results[i][j])-float(results[i+1][k]))/2
                                    short_save[2] = (float(results[i+2][l])-float(results[i+3][m]))/2
                                    if short_save[1] < window_x and short_save[1] > -window_x and short_save[2] < window_y and short_save[2] > -window_y:              # use this if statement instead the one of the next line to find the spot
                                        spot_positions.append(short_save)
                                        q.put(short_save)
    # return(spot_positions, check_event_taken, double_taken_event_counter, temp_results)
    return(spot_positions)


def devide(list_in, parts):
    '''
    The function divides a list into equal parts of minimum size 10
    (to keep p1 and p2 always together) and returns it as a dict
    of list of lists.
    '''
    length = len(list_in)
    part_length = int(length/10/parts)*10
    dict_parts = {}
    for i in range(parts-1):
        hilf = []
        for j in range(i*part_length, i*part_length+part_length):
            hilf.append(list_in[j])
        dict_parts["results_part{0}".format(i+1)]=hilf
    hilf = []
    for j in range((parts-1)*part_length, length):
        hilf.append(list_in[j])
    dict_parts["results_part{0}".format(parts)]=hilf
    return(dict_parts)


def multiprocess_rec_0(import_dict):

    number_cpu_cores = multiprocessing.cpu_count()
    print '\n### Fit\n\nCPU cores used for fit   :: ', number_cpu_cores

    dict_parts = devide(import_dict['results'], number_cpu_cores)

    m = multiprocessing.Manager()
    q = m.Queue()        # shared queue (= shared memory)

    t1 = time.time()
    if number_cpu_cores == 2:
        p1 = multiprocessing.Process(target=rec_0, args=([dict_parts['results_part1'], import_dict, 1], q))
        p2 = multiprocessing.Process(target=rec_0, args=([dict_parts['results_part2'], import_dict, 2], q))

        p1.start()
        p2.start()

        p1.join()
        p2.join()
    elif number_cpu_cores == 4:
        p1 = multiprocessing.Process(target=rec_0, args=([dict_parts['results_part1'], import_dict, 1], q))
        p2 = multiprocessing.Process(target=rec_0, args=([dict_parts['results_part2'], import_dict, 2], q))
        p3 = multiprocessing.Process(target=rec_0, args=([dict_parts['results_part3'], import_dict, 3], q))
        p4 = multiprocessing.Process(target=rec_0, args=([dict_parts['results_part4'], import_dict, 4], q))

        p1.start()
        p2.start()
        p3.start()
        p4.start()

        p1.join()
        p2.join()
        p3.join()
        p4.join()
    elif number_cpu_cores == 8:
        p1 = multiprocessing.Process(target=rec_0, args=([dict_parts['results_part1'], import_dict, 1], q))
        p2 = multiprocessing.Process(target=rec_0, args=([dict_parts['results_part2'], import_dict, 2], q))
        p3 = multiprocessing.Process(target=rec_0, args=([dict_parts['results_part3'], import_dict, 3], q))
        p4 = multiprocessing.Process(target=rec_0, args=([dict_parts['results_part4'], import_dict, 4], q))
        p5 = multiprocessing.Process(target=rec_0, args=([dict_parts['results_part5'], import_dict, 5], q))
        p6 = multiprocessing.Process(target=rec_0, args=([dict_parts['results_part6'], import_dict, 6], q))
        p7 = multiprocessing.Process(target=rec_0, args=([dict_parts['results_part7'], import_dict, 7], q))
        p8 = multiprocessing.Process(target=rec_0, args=([dict_parts['results_part8'], import_dict, 8], q))

        p1.start()
        p2.start()
        p3.start()
        p4.start()
        p5.start()
        p6.start()
        p7.start()
        p8.start()

        p1.join()
        p2.join()
        p3.join()
        p4.join()
        p5.join()
        p6.join()
        p7.join()
        p8.join()

    spots = []
    while q.empty() is False:
        spots_hilf = q.get()
        spots.append(spots_hilf)

    return(spots)
# -----------------------------------------------
# N2
# -----------------------------------------------


def rec_1(results, trigger_on_off, tof_cut_on, tof_min, tof_max, window_x, window_y, err_window_x, err_window_y, window_xy, err_window_xy, sigma, check_event_taken):
    """
    Event reconstruction function.

    The function reconstructs all possible events in one ejection with only 3
    matching time stamp in the X1, X2, Y1, Y2 lines --> event only possible if
    times are within the measurement window(= on the detector) and 3 time per
    X1, X2, Y1, Y2 are matching to a MCP time. If data is lost in two or more
    channels the event is not reconstructed. The function checks if a point
    was already reconstructed by a lower order reconstruction, so points
    won't be counted twice.
    """
    import copy

    spot_positions_reconstructed = []  # [ *ersetzen* for i in 10*range(len(results))]        # 1st column = number of ejection, 2nd column x position, 3rd column = y position
    short_save_rec = [0., 0., 0.]  # nur ein Zwischenspeicher (buffer variable)
    check_event_taken_rec = copy.copy(check_event_taken)
    double_taken_event_counter_rec1 = 0

    if trigger_on_off == 1:
        for i in range(0, len(results), 5):                     # check each ejection
            if results[i][0] > 0 and results[i+1][0] > 0 and results[i+2][0] > 0:     # make sure that there is at least one hit at the detector, then build every possible combination (even when a event is missing) and store all reasonable values
                for j in range(1, results[i][0]+1):                                      # x1
                    if results[i][j] > 1000000:                                                 # excludes trigger-errors: saves time
                        for k in range(1, results[i+1][0]+1):                            # x2
                            if results[i+1][k] > 1000000:                                       # excludes trigger-errors: saves time
                                for l in range(1, results[i+2][0]+1):                    # y1
                                    if results[i+2][l] > 1000000:                               # excludes trigger-errors: saves time
                                        for n in range(1, results[i+4][0]+1):    # MCP
                                            if results[i+4][n] > 1000000:               # excludes trigger-errors: saves time
                                                if (results[i][j]+results[i+1][k]-2*results[i+4][n]) < (window_x + sigma*err_window_x) and \
                                                    (results[i][j]+results[i+1][k]-2*results[i+4][n]) > (window_x - sigma*err_window_x) and \
                                                    ((results[i+2][l] - results[i+4][n]) > 0 and
                                                        (results[i+2][l] - results[i+4][n]) < (window_y - 0.3*window_y)) and \
                                                    (results[i][j] + results[i+1][k] + results[i+2][l] + (window_y-results[i+2][l]+2*results[i+4][n]) - 4*results[i+4][n]) < (window_xy + sigma*err_window_xy) and \
                                                        (results[i][j] + results[i+1][k] + results[i+2][l] + (window_y-results[i+2][l]+2*results[i+4][n]) - 4*results[i+4][n]) > (window_xy - sigma*err_window_xy):
                                                    short_save_rec = [0., 0., 0.]
                                                    short_save_rec[0] = i
                                                    short_save_rec[1] = (float(results[i][j])-float(results[i+1][k]))/2
                                                    short_save_rec[2] = (float(results[i+2][l])-float((window_y-results[i+2][l]+2*results[i+4][n])))/2
                                                    if short_save_rec[1] < window_x and short_save_rec[1] > -window_x and short_save_rec[2] < window_y and short_save_rec[2] > -window_y:              # check if the matched signal is also on the mcp
                                                        check_event_taken_rec[i][j] += 1
                                                        check_event_taken_rec[i+1][k] += 1
                                                        check_event_taken_rec[i+2][l] += 1
                                                        check_event_taken_rec[i+4][n] += 1
                                                        if check_event_taken_rec[i][j] > 1 or check_event_taken_rec[i+1][k] > 1 or check_event_taken_rec[i+2][l] > 1 or check_event_taken_rec[i+4][n] > 1:
                                                        # print 'Double counting in ejection reconstruction %s'%((i+1)/5)
                                                            double_taken_event_counter_rec1 += 1
                                                        else:
                                                            spot_positions_reconstructed.append(short_save_rec)
    if trigger_on_off == 1:
        for i in range(0, len(results), 5):                     # check each ejection
            if results[i][0] > 0 and results[i+1][0] > 0 and results[i+3][0] > 0:     # make sure that there is at least one hit at the detector, then build every possible combination (even when a event is missing) and store all reasonable values
                for j in range(1, results[i][0]+1):                                      # x1
                    if results[i][j] > 1000000:                                                 # excludes trigger-errors: saves time
                        for k in range(1, results[i+1][0]+1):                            # x2
                            if results[i+1][k] > 1000000:                                       # excludes trigger-errors: saves time
                                for m in range(1, results[i+3][0]+1):            # y2
                                    if results[i+3][m] > 1000000:                       # excludes trigger-errors: saves time
                                        for n in range(1, results[i+4][0]+1):    # MCP
                                            if results[i+4][n] > 1000000:               # excludes trigger-errors: saves time
                                                if (results[i][j]+results[i+1][k]-2*results[i+4][n]) < (window_x + sigma*err_window_x) and \
                                                    (results[i][j]+results[i+1][k]-2*results[i+4][n]) > (window_x - sigma*err_window_x) and \
                                                    (results[i+3][m] - results[i+4][n]) > 0 and \
                                                    (results[i+3][m] - results[i+4][n]) < (window_y - 0.3*window_y) and \
                                                    (results[i][j] + results[i+1][k] + results[i+3][m] + (window_y-results[i+3][m]+2*results[i+4][n]) - 4*results[i+4][n]) < (window_xy + sigma*err_window_xy) and\
                                                        (results[i][j] + results[i+1][k] + results[i+3][m] + (window_y-results[i+3][m]+2*results[i+4][n]) - 4*results[i+4][n]) > (window_xy - sigma*err_window_xy):
                                                    short_save_rec = [0., 0., 0.]
                                                    short_save_rec[0] = i
                                                    short_save_rec[1] = (float(results[i][j])-float(results[i+1][k]))/2
                                                    short_save_rec[2] = (float((window_y-results[i+3][m]+2*results[i+4][n]))-float(results[i+3][m]))/2
                                                    if short_save_rec[1] < window_x and short_save_rec[1] > -window_x and short_save_rec[2] < window_y and short_save_rec[2] > -window_y:              # check if the matched signal is also on the mcp
                                                        check_event_taken_rec[i][j] += 1
                                                        check_event_taken_rec[i+1][k] += 1
                                                        check_event_taken_rec[i+3][m] += 1
                                                        check_event_taken_rec[i+4][n] += 1
                                                        if check_event_taken_rec[i][j] > 1 or check_event_taken_rec[i+1][k] > 1 or check_event_taken_rec[i+3][m] > 1 or check_event_taken_rec[i+4][n] > 1:
                                                        # print 'Double counting in ejection reconstruction %s'%((i+1)/5)
                                                            double_taken_event_counter_rec1 += 1
                                                        else:
                                                            spot_positions_reconstructed.append(short_save_rec)
    if trigger_on_off == 1:
        for i in range(0, len(results), 5):                     # check each ejection
            if results[i][0] > 0 and results[i+2][0] > 0 and results[i+3][0] > 0:     # make sure that there is at least one hit at the detector, then build every possible combination (even when a event is missing) and store all reasonable values
                for j in range(1, results[i][0]+1):                                      # x1
                    if results[i][j] > 1000000:                                                 # excludes trigger-errors: saves time
                        for l in range(1, results[i+2][0]+1):                    # y1
                            if results[i+2][l] > 1000000:                               # excludes trigger-errors: saves time
                                for m in range(1, results[i+3][0]+1):            # y2
                                    if results[i+3][m] > 1000000:                       # excludes trigger-errors: saves time
                                        for n in range(1, results[i+4][0]+1):    # MCP
                                            if results[i+4][n] > 1000000:               # excludes trigger-errors: saves time
                                                if (results[i+2][l]+results[i+3][m]-2*results[i+4][n]) < (window_y + sigma*err_window_y) and \
                                                    (results[i+2][l]+results[i+3][m]-2*results[i+4][n]) > (window_y - sigma*err_window_y) and \
                                                    ((results[i][j]-results[i+4][n]) > 0 and
                                                        (results[i][j]-results[i+4][n]) < (window_x - 0.3*window_x)) and \
                                                    (results[i+2][l]+results[i+3][m] + results[i][j] + (window_x-results[i][j]+2*results[i+4][n]) - 4*results[i+4][n]) < (window_xy + sigma*err_window_xy) and \
                                                        (results[i+2][l]+results[i+3][m] + results[i][j] + (window_x-results[i][j]+2*results[i+4][n]) - 4*results[i+4][n]) > (window_xy - sigma*err_window_xy):
                                                    short_save_rec = [0., 0., 0.]
                                                    short_save_rec[0] = i
                                                    short_save_rec[1] = (float(results[i][j])-float((window_x-results[i][j]+2*results[i+4][n])))/2
                                                    short_save_rec[2] = (float(results[i+2][l])-float(results[i+3][m]))/2
                                                    if short_save_rec[1] < window_x and short_save_rec[1] > -window_x and short_save_rec[2] < window_y and short_save_rec[2] > -window_y:              # check if the matched signal is also on the mcp
                                                        check_event_taken_rec[i][j] += 1
                                                        check_event_taken_rec[i+2][l] += 1
                                                        check_event_taken_rec[i+3][m] += 1
                                                        check_event_taken_rec[i+4][n] += 1
                                                        if check_event_taken_rec[i][j] > 1 or check_event_taken_rec[i+2][l] > 1 or check_event_taken_rec[i+3][m] > 1 or check_event_taken_rec[i+4][n] > 1:
                                                        # print 'Double counting in ejection reconstruction %s'%((i+1)/5)
                                                            double_taken_event_counter_rec1 += 1
                                                        else:
                                                            spot_positions_reconstructed.append(short_save_rec)
    if trigger_on_off == 1:
        for i in range(0, len(results), 5):                     # check each ejection
            if results[i+1][0] > 0 and results[i+2][0] > 0 and results[i+3][0] > 0:     # make sure that there is at least one hit at the detector, then build every possible combination (even when a event is missing) and store all reasonable values
                for k in range(1, results[i+1][0]+1):                            # x2
                    if results[i+1][k] > 1000000:                                       # excludes trigger-errors: saves time
                        for l in range(1, results[i+2][0]+1):                    # y1
                            if results[i+2][l] > 1000000:                               # excludes trigger-errors: saves time
                                for m in range(1, results[i+3][0]+1):            # y2
                                    if results[i+3][m] > 1000000:                       # excludes trigger-errors: saves time
                                        for n in range(1, results[i+4][0]+1):    # MCP
                                            if results[i+4][n] > 1000000:               # excludes trigger-errors: saves time
                                                if (results[i+2][l]+results[i+3][m]-2*results[i+4][n]) < (window_y + 0.5*sigma*err_window_y) and \
                                                    (results[i+2][l]+results[i+3][m]-2*results[i+4][n]) > (window_y - sigma*err_window_y) and \
                                                    ((results[i+1][k]-results[i+4][n]) > 0 and
                                                        (results[i+1][k]-results[i+4][n]) < (window_x - 0.3*window_x)) and \
                                                    (results[i+2][l]+results[i+3][m] + results[i+1][k] + (window_x-results[i+1][k]+2*results[i+4][n]) - 4*results[i+4][n]) < (window_xy + sigma*err_window_xy) and \
                                                        (results[i+2][l]+results[i+3][m] + results[i+1][k] + (window_x-results[i+1][k]+2*results[i+4][n]) - 4*results[i+4][n]) > (window_xy - sigma*err_window_xy):
                                                    short_save_rec = [0., 0., 0.]
                                                    short_save_rec[0] = i
                                                    short_save_rec[1] = (float((window_x-results[i+1][k]+2*results[i+4][n]))-float(results[i+1][k]))/2
                                                    short_save_rec[2] = (float(results[i+2][l])-float(results[i+3][m]))/2
                                                    if short_save_rec[1] < window_x and short_save_rec[1] > -window_x and short_save_rec[2] < window_y and short_save_rec[2] > -window_y:              # check if the matched signal is also on the mcp
                                                        check_event_taken_rec[i+1][k] += 1
                                                        check_event_taken_rec[i+2][l] += 1
                                                        check_event_taken_rec[i+3][m] += 1
                                                        check_event_taken_rec[i+4][n] += 1
                                                        if check_event_taken_rec[i+1][k] > 1 or check_event_taken_rec[i+2][l] > 1 or check_event_taken_rec[i+3][m] > 1 or check_event_taken_rec[i+4][n] > 1:
                                                        # print 'Double counting in ejection reconstruction %s'%((i+1)/5)
                                                            double_taken_event_counter_rec1 += 1
                                                        else:
                                                            spot_positions_reconstructed.append(short_save_rec)

    return(spot_positions_reconstructed, check_event_taken_rec, double_taken_event_counter_rec1)   # , temp_results)


def rec_2(results, trigger_on_off, tof_cut_on, tof_min, tof_max, window_x, window_y, err_window_x, err_window_y, window_xy, err_window_xy, sigma, check_event_taken_rec):
    """
    Event reconstruction function.

    The function reconstructs all possible events in one ejection with only 1
    matching time stamp in the X1, X2 and Y1, Y2 respectively --> event only
    possible if times are within the measurement window(= on the detector) and
    1 X and 1 Y time is matching to a MCP time. If data is lost in two or more
    channels the event is not reconstructed. The function checks if a point
    was already reconstructed by a lower order reconstruction, so points
    won't be counted twice.
    """
    import copy

    spot_positions_reconstructed_2 = []     # [ *ersetzen* for i in 10*range(len(results))]        # 1st column = number of ejection, 2nd column x position, 3rd column = y position
    short_save_rec_2 = [0., 0., 0.]         # nur ein Zwischenspeicher (buffer variable)
    check_event_taken_rec_2 = copy.copy(check_event_taken_rec)
    double_taken_event_counter_rec2 = 0

    if trigger_on_off == 1:
        for i in range(0, len(results), 5):                     # check each ejection
            if results[i][0] > 0 and results[i+1][0] > 0 and results[i+2][0] > 0 and results[i+3][0] > 0:     # make sure that there is at least one hit at the detector, then build every possible combination (even when a event is missing) and store all reasonable values
                for j in range(1, results[i][0]+1):                                      # x1
                    if results[i][j] > 1000000:                                                 # excludes trigger-errors: saves time
                        for k in range(1, results[i+1][0]+1):                            # x2
                            if results[i+1][k] > 1000000:                                       # excludes trigger-errors: saves time
                                for l in range(1, results[i+2][0]+1):                    # y1
                                    if results[i+2][l] > 1000000:                               # excludes trigger-errors: saves time
                                        for m in range(1, results[i+3][0]+1):            # y2
                                            if results[i+3][m] > 1000000:                       # excludes trigger-errors: saves time
                                                for n in range(1, results[i+4][0]+1):    # MCP
                                                    if results[i+4][n] > 1000000:               # excludes trigger-errors: saves time
                                                        if ((results[i][j]-results[i+4][n]) > 0 and
                                                            (results[i][j]-results[i+4][n]) < (window_x)) and \
                                                            ((results[i+2][l]-results[i+4][n]) > 0 and
                                                                (results[i+2][l]-results[i+4][n]) < (window_y)) and \
                                                            ((results[i][j] + (window_y-results[i][j]+2*results[i+4][n]) + results[i+2][l] + (window_y-results[i+2][l]+2*results[i+4][n]) - 4*results[i+4][n]) < (window_xy + sigma*err_window_xy) and
                                                                (results[i][j] + (window_y-results[i][j]+2*results[i+4][n]) + results[i+2][l] + (window_y-results[i+2][l]+2*results[i+4][n]) - 4*results[i+4][n]) > (window_xy - sigma*err_window_xy)):
                                                            short_save_rec_2 = [0., 0., 0.]
                                                            short_save_rec_2[0] = i
                                                            short_save_rec_2[1] = (float(results[i][j])-float((window_x-results[i][j]+2*results[i+4][n])))/2
                                                            short_save_rec_2[2] = (float(results[i+2][l])-float((window_y-results[i+2][l]+2*results[i+4][n])))/2
                                                            if short_save_rec_2[1] < window_x and short_save_rec_2[1] > -window_x and short_save_rec_2[2] < window_y and short_save_rec_2[2] > -window_y:              # check if the matched signal is also on the mcp
                                                                check_event_taken_rec_2[i][j] += 1
                                                                check_event_taken_rec_2[i+2][l] += 1
                                                                check_event_taken_rec_2[i+4][n] += 1
                                                                if check_event_taken_rec_2[i][j] > 1 or check_event_taken_rec_2[i+2][l] > 1 or check_event_taken_rec_2[i+4][n] > 1:
                                                                # print 'Double counting in ejection reconstruction %s'%((i+1)/5)
                                                                    double_taken_event_counter_rec2 += 1
                                                                else:
                                                                    spot_positions_reconstructed_2.append(short_save_rec_2)
                                                        if ((results[i+1][k]-results[i+4][n]) > 0 and
                                                            (results[i+1][k]-results[i+4][n]) < (window_x)) and \
                                                            ((results[i+2][l]-results[i+4][n]) > 0 and
                                                                (results[i+2][l]-results[i+4][n]) < (window_y)) and \
                                                            ((results[i+1][k] + (window_y-results[i+1][k]+2*results[i+4][n]) + results[i+2][l] + (window_y-results[i+2][l]+2*results[i+4][n]) - 4*results[i+4][n]) < (window_xy + sigma*err_window_xy) and
                                                                (results[i+1][k] + (window_y-results[i+1][k]+2*results[i+4][n]) + results[i+2][l] + (window_y-results[i+2][l]+2*results[i+4][n]) - 4*results[i+4][n]) > (window_xy - sigma*err_window_xy)):
                                                            short_save_rec_2 = [0., 0., 0.]
                                                            short_save_rec_2[0] = i
                                                            short_save_rec_2[1] = (float(results[i+1][k])-float((window_x-results[i+1][k]+2*results[i+4][n])))/2
                                                            short_save_rec_2[2] = (float(results[i+2][l])-float((window_y-results[i+2][l]+2*results[i+4][n])))/2
                                                            if short_save_rec_2[1] < window_x and short_save_rec_2[1] > -window_x and short_save_rec_2[2] < window_y and short_save_rec_2[2] > -window_y:              # check if the matched signal is also on the mcp
                                                                check_event_taken_rec_2[i+1][k] += 1
                                                                check_event_taken_rec_2[i+2][l] += 1
                                                                check_event_taken_rec_2[i+4][n] += 1
                                                                if check_event_taken_rec_2[i+1][k] > 1 or check_event_taken_rec_2[i+2][l] > 1 or check_event_taken_rec_2[i+4][n] > 1:
                                                                # print 'Double counting in ejection reconstruction %s'%((i+1)/5)
                                                                    double_taken_event_counter_rec2 += 1
                                                                else:
                                                                    spot_positions_reconstructed_2.append(short_save_rec_2)
                                                        if ((results[i][j]-results[i+4][n]) > 0 and
                                                            (results[i][j]-results[i+4][n]) < (window_x)) and \
                                                            ((results[i+3][m]-results[i+4][n]) > 0 and
                                                                (results[i+3][m]-results[i+4][n]) < (window_y)) and \
                                                            ((results[i][j] + (window_y-results[i][j]+2*results[i+4][n]) + results[i+3][m] + (window_y-results[i+3][m]+2*results[i+4][n]) - 4*results[i+4][n]) < (window_xy + sigma*err_window_xy) and
                                                                (results[i][j] + (window_y-results[i][j]+2*results[i+4][n]) + results[i+3][m] + (window_y-results[i+3][m]+2*results[i+4][n]) - 4*results[i+4][n]) > (window_xy - sigma*err_window_xy)):
                                                            short_save_rec_2 = [0., 0., 0.]
                                                            short_save_rec_2[0] = i
                                                            short_save_rec_2[1] = (float(results[i][j])-float((window_y-results[i][j]+2*results[i+4][n])))/2
                                                            short_save_rec_2[2] = (float(results[i+3][m])-float((window_y-results[i+3][m]+2*results[i+4][n])))/2
                                                            if short_save_rec_2[1] < window_x and short_save_rec_2[1] > -window_x and short_save_rec_2[2] < window_y and short_save_rec_2[2] > -window_y:              # check if the matched signal is also on the mcp
                                                                check_event_taken_rec_2[i][j] += 1
                                                                check_event_taken_rec_2[i+3][m] += 1
                                                                check_event_taken_rec_2[i+4][n] += 1
                                                                if check_event_taken_rec_2[i][j] > 1 or check_event_taken_rec_2[i+3][m] > 1 or check_event_taken_rec_2[i+4][n] > 1:
                                                                # print 'Double counting in ejection reconstruction %s'%((i+1)/5)
                                                                    double_taken_event_counter_rec2 += 1
                                                                else:
                                                                    spot_positions_reconstructed_2.append(short_save_rec_2)
                                                        if ((results[i+1][k]-results[i+4][n]) > 0 and
                                                            (results[i+1][k]-results[i+4][n]) < (window_x)) and \
                                                            ((results[i+3][m]-results[i+4][n]) > 0 and
                                                                (results[i+3][m]-results[i+4][n]) < (window_y)) and \
                                                            ((results[i+1][k] + (window_y-results[i+1][k]+2*results[i+4][n]) + results[i+3][m] + (window_y-results[i+3][m]+2*results[i+4][n]) - 4*results[i+4][n]) < (window_xy + sigma*err_window_xy) and
                                                                (results[i+1][k] + (window_y-results[i+1][k]+2*results[i+4][n]) + results[i+3][m] + (window_y-results[i+3][m]+2*results[i+4][n]) - 4*results[i+4][n]) > (window_xy - sigma*err_window_xy)):
                                                            short_save_rec_2 = [0., 0., 0.]
                                                            short_save_rec_2[0] = i
                                                            short_save_rec_2[1] = (float(results[i+1][k])-float((window_y-results[i+1][k]+2*results[i+4][n])))/2
                                                            short_save_rec_2[2] = (float(results[i+3][m])-float((window_y-results[i+3][m]+2*results[i+4][n])))/2
                                                            if short_save_rec_2[1] < window_x and short_save_rec_2[1] > -window_x and short_save_rec_2[2] < window_y and short_save_rec_2[2] > -window_y:              # check if the matched signal is also on the mcp
                                                                check_event_taken_rec_2[i+1][k] += 1
                                                                check_event_taken_rec_2[i+3][m] += 1
                                                                check_event_taken_rec_2[i+4][n] += 1
                                                                if check_event_taken_rec_2[i+1][k] > 1 or check_event_taken_rec_2[i+3][m] > 1 or check_event_taken_rec_2[i+4][n] > 1:
                                                                # print 'Double counting in ejection reconstruction %s'%((i+1)/5)
                                                                    double_taken_event_counter_rec2 += 1
                                                                else:
                                                                    spot_positions_reconstructed_2.append(short_save_rec_2)


#         histo2D(spot_positions_reconstructed_2, file_name, pattern, 3)
    return(spot_positions_reconstructed_2, check_event_taken_rec_2, double_taken_event_counter_rec2)   # , temp_results)


# ------------------------
# TOF distribution plot
# ------------------------
def tof_plot(spot_positions, bins, file_name, dont_show):
    """Function to plot the time of flight distribution of dataset."""
    mpl.rc('font', family='serif', serif='Utopia')      # Utopia LaTeX font!!
    mpl.rc('text', usetex=False)

    tof = [float(x[3]) for x in spot_positions]   # in microseconds
    fig = plt.figure(17, figsize=(8, 8), facecolor='white')
    ax = fig.add_subplot(111)
    plt.hist(tof, bins, histtype='step', stacked=True, fill=False)  # , range=[43, 58]
    plt.xlabel('time of flight / $\mu$s', fontsize=22)
    plt.ylabel('counts', fontsize=22)
    plt.title('Time of flight distribution', fontsize=26)
    ax.text(ax.get_xticks()[0]*1.01, ax.get_yticks()[-1]*0.95, '%s' % (file_name), fontsize=10)
#            plt.ylim(0, 20)
    plt.savefig('%s_tof.pdf' % file_name)

    if dont_show == 0:
        plt.show()
    elif dont_show == 1:
        plt.close()


def ion_cloud_3d(spot_positions, axis_scaling_factor, dont_show, file_name):
    """
    Nice 3D ioncloud plot.

    The function gives a nice combination between the 2D MCP plot and the
    time of flight information. It therefore reconstructs a picture of the
    3D ioncloud.
    """
    mpl.rc('font', family='serif', serif='Utopia')      # Utopia LaTeX font!!
    mpl.rc('text', usetex=False)

    xss = [x[1]*axis_scaling_factor for x in spot_positions]  # in nanoseconds
    yss = [x[2]*axis_scaling_factor for x in spot_positions]  # in nanoseconds
    zss = [x[3]*0.000025 for x in spot_positions]   # in mikroseconds
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xss, yss, zss, zdir='xss', s=2, c='green')
    ax.set_xlabel('X direction / mm', fontsize=22)
    ax.set_ylabel('Y direction / mm', fontsize=22)
    ax.set_zlabel('TOF / $\mu$s', fontsize=22)
    ax.set_xlim([-20, 20])
    ax.set_ylim([-20, 20])

    ax.view_init(0, 0)
    #     plt.draw()
    #     plt.pause(.001)
    # for angle in range(0, 360):
    #     ax.view_init(28, angle)
    #     plt.draw()
    #     plt.pause(.001)


    plt.savefig('%s_ioncloud.pdf' % file_name)

    if dont_show == 0:
        plt.show()
    elif dont_show == 1:
        plt.close()


def spot_fit_root(file_name, spot_positions, all_peak_positions, mean_x, mean_y, fit_range, nll_plot, plot_fit, iteration, sigma_import):
    '''
    Function to fit position to 2D unbinned csv file.

    :param fit_range: vector containing [0] = cut_x_min, [1] = cut_x_max, [2] = cut_y_min, [3] = cut_y_max
    '''
    root.gErrorIgnoreLevel = root.kInfo
    root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)


    # import data into TTree
    tree = root.TTree('tree', 'tree')
    tree.ReadFile('%s_spot_positions.csv' % (file_name), 'inj/D:x/D:y/D:t/D', ',')

    # set vector properties
    inj = root.RooRealVar('inj', 'inj', 0, 0.0, 5000.0)
    x = root.RooRealVar('x', 'x', 0, fit_range[0], fit_range[1])
    y = root.RooRealVar('y', 'y', 0, fit_range[2], fit_range[3])
    t = root.RooRealVar('t', 't', 1000000.0, 1000000., 3000000.0)

    # set dataset
    ds = root.RooDataSet('x-y', 'x-y-data', tree, root.RooArgSet(x, y))

    # set fitting model
    meanx = root.RooRealVar('meanx', 'meanx', mean_x, mean_x-50.0, mean_x+50.0)
    meany = root.RooRealVar('meany', 'meany', mean_y, mean_y-50.0, mean_y+50.0)

    if sigma_import['fixed'] == True:
        sigmax = root.RooRealVar('sigmax', 'sigmax', float(sigma_import['sigmax']))
        sigmay = root.RooRealVar('sigmay', 'sigmay', float(sigma_import['sigmay']))
        sigmax.setConstant(True)
        sigmay.setConstant(True)
    else:
        sigmax = root.RooRealVar('sigmax', 'sigmax', 50.0, 10.0, 100.0)
        sigmay = root.RooRealVar('sigmay', 'sigmay', 50.0, 10.0, 100.0)

    gaussx = root.RooGaussian('gaussx', 'Gaussian distribution', x, meanx, sigmax)
    gaussy = root.RooGaussian('gaussy', 'Gaussian distribution', y, meany, sigmay)
    gaussxy = root.RooProdPdf('gxy', 'gx*gy', root.RooArgList(gaussx,gaussy))

    # actual fit:
    result = gaussxy.fitTo(ds, root.RooFit.PrintLevel(-1), root.RooFit.Verbose(False)) # , root.RooFit.NumCPU(4))

    if nll_plot == 'yes':
        # plot projections and fits
        c = root.TCanvas('c', 'results', 1000, 700)
        c.Divide(2, 2)

        c.cd(1)
        framex = x.frame()
        ds.plotOn(framex)
        gaussxy.plotOn(framex)
        gaussx.paramOn(framex)
        framex.Draw()
        framex.SetTitle('X-Proj. & unbinned 2D-max.likelihood fit')
        framex.SetXTitle('channels')

        c.cd(2)
        framey = y.frame()
        ds.plotOn(framey)
        gaussxy.plotOn(framey)
        gaussy.paramOn(framey)
        framey.Draw()
        framey.SetTitle('Y-Proj. & unbinned 2D-max.likelihood fit')
        framey.SetXTitle('channels')

        # plot nll's
        c.cd(3)
        nllx = gaussxy.createNLL(ds) # , root.RooFit.NumCPU(4))
        profile_llmeanx = nllx.createProfile(root.RooArgSet(meanx))
        pllframex = meanx.frame()
        nllx.plotOn(pllframex, root.RooFit.ShiftToZero())
        profile_llmeanx.plotOn(pllframex, root.RooFit.LineColor(root.kRed))
        pllframex.SetTitle('X-NLL analysis')
        pllframex.SetXTitle('channels')
        # pllframex.SetMinimum(0)
        # pllframex.SetMaximum(3)
        pllframex.Draw()


        c.cd(4)
        nlly = gaussxy.createNLL(ds) # , root.RooFit.NumCPU(4))
        profile_llmeany = nlly.createProfile(root.RooArgSet(meany))
        pllframey = meany.frame()
        nlly.plotOn(pllframey, root.RooFit.ShiftToZero())
        profile_llmeany.plotOn(pllframey, root.RooFit.LineColor(root.kRed))
        pllframey.SetTitle('Y-NLL analysis')
        pllframey.SetXTitle('channels')
        # pllframey.SetMinimum(0)
        # pllframey.SetMaximum(3)
        pllframey.Draw()
        c.Update()

        pdffile = file_name + '_x-y-projection-fits.pdf'
        c.SaveAs(pdffile)

    # fit results
    x_pos = meanx.getValV()
    x_pos_err = meanx.getError()
    y_pos = meany.getValV()
    y_pos_err = meany.getError()
    x_sigma = sigmax.getValV()
    x_sigma_err = sigmax.getError()
    y_sigma = sigmay.getValV()
    y_sigma_err = sigmay.getError()

    print 'Fit results in iter. {}   ::  x = {}, y = {}'.format(iteration, x_pos, y_pos)

    if plot_fit == 'yes':
        only_spot_positions = [['X position / ch.', 'Y position / ch.']]
        for i in range(len(spot_positions)):
            only_spot_positions.append([float(spot_positions[i][1]), float(spot_positions[i][2])])
        python_plot(only_spot_positions, '', '%s_mcp_iter%s' % (file_name, iteration), 'no', '', 'no', 'gauss', 'Utopia', 'red', 'black', 'full', 0, 8, 6, 'mcp', 60, 1, 'on', 'off', 'mcp', [x_sigma, x_sigma_err, y_sigma, y_sigma_err], [x_pos, x_pos_err, y_pos, y_pos_err])
    # python_plot(only_spot_positions, 'MCP-PS', '%s_mcp_2' % file_name, 'no', '', 'no', 'gauss', 'Utopia', 'red', 'black', 'full', 0, 8, 6, 'mcp', 60, 9, 'on', 'off', 'mcp', [x_sigma, x_sigma_err, y_sigma, y_sigma_err], [x_pos, x_pos_err, y_pos, y_pos_err])

    # python_plot(only_spot_positions, 'MCP-PS', '%s_mcp_2' % file_name, 'no', '', 'no', 'gauss', 'Utopia', 'red', 'black', 'full', 0, 8, 6, '2dhistogram-mcp', 60, 1, 'on', 'off', 'mcp', [x_sigma, x_sigma_err, y_sigma, y_sigma_err], [x_pos, x_pos_err, y_pos, y_pos_err])

    #root.gApplication.Run()

    # output
    all_peak_positions = [file_name, x_pos, x_pos_err, y_pos, y_pos_err, x_sigma, x_sigma_err, y_sigma, y_sigma_err, ds.numEntries()]

    return(all_peak_positions)


def z_class_reduction(spot_positions, z_class_analysis):
    '''
    :param: spot_positions      = csv list with entries of ejection number, x-pos and y-pos of reconstructed spot as well as its ToF
    :param: z_class_analysis    = list with three entries: string z-class-analysis 'yes'/'no', min number of ions and max number of ions

    Function which reduces the spot positions list to a list which contains only spots with given number of ions in the trap.
    '''
    ejections = []
    for entry in spot_positions:
        ejections.append(entry[0])

    ejections_dict = dict(Counter(ejections))   # count how many ions per ejection were found

    reduced_spot_positions = []
    for key in ejections_dict:
        if ejections_dict[key] >= z_class_analysis[1] and ejections_dict[key] <= z_class_analysis[2]:
            for entry in spot_positions:
                if key == entry[0]:
                    reduced_spot_positions.append(entry)

    return(reduced_spot_positions)


def root_1D_unbinned_gauss_fit(liste, file_name):
    '''Function takes 1D list with unbinned data and returns a list with gaussian fit parameters: x, x_unc., sigma, sigma_unc.'''
    root.gErrorIgnoreLevel = root.kInfo
    root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)



    import_list = [float(x) for x in liste]

    # get the x where most of the ions arrive
    hist, hilf_edges = np.histogram(import_list, bins=2000, range=(0,200))

    mpl.rc('font', family='serif', serif='Utopia')      # Utopia LaTeX font!!
    mpl.rc('text', usetex=False)
    fig = plt.figure(1110002, figsize=(9,6))

    plt.plot(hilf_edges[1:], hist)

    plt.xlabel('time of flight / $\mu$s', fontsize=22)
    plt.ylabel('counts', fontsize=22)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    x_edges = []
    for i in range(len(hilf_edges)-1):
        x_edges.append(np.mean([hilf_edges[i], hilf_edges[i+1]]))
    x_max = x_edges[index_max_value(hist)]

    plt.plot([x_max], [max(hist)], 'r *')
    plt.tight_layout()
    plt.savefig('{}_tof_total.pdf'.format(file_name))
    plt.close()
    # calculate mean of limited range of x interval
    temp_list = []
    for i in import_list:
        if i < x_max + 5 and i > x_max - 5:
            temp_list.append(i)

    mean_x = np.mean(temp_list)
    std_x = np.std(temp_list)


    # import list to ROOT
    # f = root.TFile( 'test.root', 'recreate' )
    tree = root.TTree( 'tree', 'tree' )

    x = array('d', [ 0. ])
    tree.Branch('x', x, 'x/D')
    for i in range(len(import_list)):
        x[0] = import_list[i]
        tree.Fill()
    # unbinned max. likelihood fit
    x = root.RooRealVar('x', 'x', 0, 200, 'us')
    ds = root.RooDataSet('x-data', 'x-data', root.RooArgSet(x), root.RooFit.Import(tree))
    meant = root.RooRealVar('meant', 'meant', mean_x, mean_x-5, mean_x+5)#, 'us')
    sigmax = root.RooRealVar('sigmax', 'sigmax', std_x, 0.01, 5)
    gaussx = root.RooGaussian('gaussx', 'Gaussian distribution', x, meant, sigmax)
    x.setRange('range', mean_x-10, mean_x+10)
    result = gaussx.fitTo(ds, root.RooFit.Range('range'), root.RooFit.PrintLevel(-1), root.RooFit.Verbose(False)) # , root.RooFit.NumCPU(4))
    # fit results
    x_pos = meant.getValV()
    x_pos_err = meant.getError()
    x_sigma = sigmax.getValV()
    x_sigma_err = sigmax.getError()
    # plot
    # c = root.TCanvas('c', 'x', 1000, 700)
    # # tree.Draw('x')
    # tframe = x.frame()
    # ds.plotOn(tframe, root.RooFit.Binning(2000))
    # gaussx.plotOn(tframe,root.RooFit.LineColor(root.kRed))
    # tframe.Draw()
    # c.Update()
    # pdffile = 'file_name' + '_hist.pdf'
    # c.SaveAs(pdffile)
    # root.gApplication.Run()
    # f.Write()
    # f.Close()
    return([x_pos, x_pos_err, x_sigma, x_sigma_err])


def index_max_value(liste):
########## find index of maximum in list
    index_max = 0
    max_value = -99999999999.
    for i in range(len(liste)):
        if liste[i] > max_value:
            max_value = liste[i]
            index_max = i
    return(index_max)


def two_spot_fit_root(file_name, all_peak_positions, mean_x_1, mean_y_1, mean_x_2, mean_y_2, plot_fit, iteration, sigma_dominant_peak):
    '''
    Function to fit position to 2D unbinned csv file.

    :param fit_range: vector containing [0] = cut_x_min, [1] = cut_x_max, [2] = cut_y_min, [3] = cut_y_max
    '''
    root.gErrorIgnoreLevel = root.kWarning
    root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)

    axis_scaling_factor = 0.031746032
    fit_status_meaning = {0: 'converged',
                          1: 'Covariance was mad  epos defined',
                          2: 'Hesse is invalid',
                          3: 'Edm is above max',
                          4: 'Reached call limit',
                          5: 'Any other failure'}

    # import data into TTree
    tree_1 = root.TTree('tree_1', 'tree_1')
    tree_1.ReadFile('%s_spot_positions_dominant.csv' % (file_name), 'inj/D:x/D:y/D:t/D', ',')
    tree_2 = root.TTree('tree_2', 'tree_2')
    tree_2.ReadFile('%s_spot_positions_recessive.csv' % (file_name), 'inj/D:x/D:y/D:t/D', ',')

    # set vector properties
    inj = root.RooRealVar('inj', 'inj', 0, 0.0, 5000.0)
    x = root.RooRealVar('x', 'x', 0, -750, 750)
    y = root.RooRealVar('y', 'y', 0, -750, 750)
    t = root.RooRealVar('t', 't', 1000000.0, 1000000., 3000000.0)

    # set dataset
    ds = root.RooDataSet('ds', 'x-y-data', tree_1, root.RooArgSet(x, y))
    cat = root.RooCategory('cat', 'Category')
    cat.defineType('data_dominant')
    cat.defineType('data_recessive')
    cat.setLabel('data_dominant')
    ds.addColumn(cat)

    ds_2 = root.RooDataSet('ds_2', 'x-y-data', tree_2, root.RooArgSet(x, y))
    cat.setLabel('data_recessive')
    ds_2.addColumn(cat)
    ds.append(ds_2)


    # set fitting model
    meanx1 = root.RooRealVar('meanx1', 'meanx1', mean_x_1, mean_x_1-150.0, mean_x_1+150.0)
    meany1 = root.RooRealVar('meany1', 'meany1', mean_y_1, mean_y_1-150.0, mean_y_1+150.0)

    meanx2 = root.RooRealVar('meanx2', 'meanx2', mean_x_2, mean_x_2-150.0, mean_x_2+150.0)
    meany2 = root.RooRealVar('meany2', 'meany2', mean_y_2, mean_y_2-150.0, mean_y_2+150.0)

    if sigma_dominant_peak['fixed'] == True:
        sigmax = root.RooRealVar('sigmax', 'sigmax', float(sigma_dominant_peak['sigmax']))
        sigmay = root.RooRealVar('sigmay', 'sigmay', float(sigma_dominant_peak['sigmay']))
        sigmax.setConstant(True)
        sigmay.setConstant(True)
    else:
        sigmax = root.RooRealVar('sigmax', 'sigmax', 50.0, 10.0, 100.0)
        sigmay = root.RooRealVar('sigmay', 'sigmay', 50.0, 10.0, 100.0)

    gaussx1 = root.RooGaussian('gaussx1', 'Gaussian distribution', x, meanx1, sigmax)
    gaussy1 = root.RooGaussian('gaussy1', 'Gaussian distribution', y, meany1, sigmay)
    g1 = root.RooProdPdf('g1', 'gx*gy1', root.RooArgList(gaussx1,gaussy1))
    gaussx2 = root.RooGaussian('gaussx1', 'Gaussian distribution', x, meanx2, sigmax)
    gaussy2 = root.RooGaussian('gaussy1', 'Gaussian distribution', y, meany2, sigmay)
    g2 = root.RooProdPdf('g2', 'gx*gy2', root.RooArgList(gaussx2,gaussy2))

    simPdf = root.RooSimultaneous('simPdf', 'Simultaneous PDF', cat)
    simPdf.addPdf(g1, 'data_dominant')
    simPdf.addPdf(g2, 'data_recessive')

    # two_gauss = root.RooAddPdf('2g', 'g1+g2', root.RooArgList(gaussxy1,gaussxy2))
    # two_gauss.Print()


    # actual fit:
    result = simPdf.fitTo(ds, root.RooFit.Minos(1), root.RooFit.Save(), root.RooFit.PrintLevel(-1), root.RooFit.Verbose(False)) # , root.RooFit.NumCPU(4))


    # fit results
    x_pos_1 = meanx1.getValV()
    x_pos_err_1 = meanx1.getError()
    y_pos_1 = meany1.getValV()
    y_pos_err_1 = meany1.getError()
    x_pos_2 = meanx2.getValV()
    x_pos_err_2 = meanx2.getError()
    y_pos_2 = meany2.getValV()
    y_pos_err_2 = meany2.getError()
    x_sigma = sigmax.getValV()
    x_sigma_err = sigmax.getError()
    y_sigma = sigmay.getValV()
    y_sigma_err = sigmay.getError()

    print 'Fit results in iter. {}   ::  status = {}, x(dom) = {:.2f} mm, y(dom) = {:.2f} mm, x(rec) = {:.2f} mm, y(rec) = {:.2f} mm'.format(iteration, fit_status_meaning[result.status()], x_pos_1*axis_scaling_factor, y_pos_1*axis_scaling_factor, x_pos_2*axis_scaling_factor, y_pos_2*axis_scaling_factor)

    pll_plot = 'yes'
    spot_pos_dom = pd.read_csv('%s_spot_positions_dominant.csv' % (file_name)).values.tolist()
    spot_pos_rec = pd.read_csv('%s_spot_positions_recessive.csv' % (file_name)).values.tolist()


    if pll_plot == 'yes' and plot_fit == 'yes':
        c4 = root.TCanvas('c4', 'results', 1000, 700)
        c4.Divide(2, 2)


        # plot nll's
        c4.cd(1)
        nllx1 = simPdf.createNLL(ds) # , root.RooFit.NumCPU(4))
        profile_llmeanx1 = nllx1.createProfile(root.RooArgSet(meanx1))
        pllframex1 = meanx1.frame()#root.RooFit.Range(x_pos_1-3*x_sigma, x_pos_1+3*x_sigma))
        nllx1.plotOn(pllframex1, root.RooFit.ShiftToZero())
        profile_llmeanx1.plotOn(pllframex1, root.RooFit.LineColor(root.kRed))
        pllframex1.SetTitle('LL and profileLL in x position (dominant)')
        pllframex1.SetXTitle('channels')
        pllframex1.SetYTitle('projection of -log(likelihood)')
        pllframex1.Draw()
        c4.Modified()
        c4.Update()


        c4.cd(2)
        nlly1 = simPdf.createNLL(ds) # , root.RooFit.NumCPU(4))
        profile_llmeany1 = nlly1.createProfile(root.RooArgSet(meany1))
        pllframey1 = meany1.frame()#root.RooFit.Range(y_pos_1-3*y_sigma, y_pos_1+3*y_sigma))
        nlly1.plotOn(pllframey1, root.RooFit.ShiftToZero())
        profile_llmeany1.plotOn(pllframey1, root.RooFit.LineColor(root.kRed))
        pllframey1.SetTitle('LL and profileLL in y position (dominant)')
        pllframey1.SetXTitle('channels')
        pllframey1.SetYTitle('projection of -log(likelihood)')
        pllframey1.Draw()
        c4.Modified()
        c4.Update()


        c4.cd(3)
        nllx2 = simPdf.createNLL(ds) # , root.RooFit.NumCPU(4))
        profile_llmeanx2 = nllx2.createProfile(root.RooArgSet(meanx2))
        pllframex2 = meanx2.frame()#root.RooFit.Range(x_pos_2-3*x_sigma, x_pos_2+3*x_sigma))
        nllx2.plotOn(pllframex2, root.RooFit.ShiftToZero())
        profile_llmeanx2.plotOn(pllframex2, root.RooFit.LineColor(root.kRed))
        pllframex2.SetTitle('LL and profileLL in x position (recessive)')
        pllframex2.SetXTitle('channels')
        pllframex2.SetYTitle('projection of -log(likelihood)')
        pllframex2.Draw()
        c4.Modified()
        c4.Update()


        c4.cd(4)
        nlly2 = simPdf.createNLL(ds) # , root.RooFit.NumCPU(4))
        profile_llmeany2 = nlly2.createProfile(root.RooArgSet(meany2))
        pllframey2 = meany2.frame()#root.RooFit.Range(y_pos_2-3*y_sigma, y_pos_2+3*y_sigma))
        nlly2.plotOn(pllframey2, root.RooFit.ShiftToZero())
        profile_llmeany2.plotOn(pllframey2, root.RooFit.LineColor(root.kRed))
        pllframey2.SetTitle('LL and profileLL in y position (recessive)')
        pllframey2.SetXTitle('channels')
        pllframey2.SetYTitle('projection of -log(likelihood)')
        pllframey2.Draw()
        c4.Modified()
        c4.Update()

        c4.Paint()

        pdffile = file_name + '_PLL.png'        # pdf doesn't work, switched to png
        c4.SaveAs(pdffile)


    if plot_fit == 'yes':
        spots = spot_pos_dom + spot_pos_rec
        only_spot_positions = [['X position / ch.', 'Y position / ch.']]
        for i in range(len(spots)):
            only_spot_positions.append([float(spots[i][1]), float(spots[i][2])])

        python_plot(only_spot_positions, '', '%s_mcp_iter%s' % (file_name, iteration), 'no', '', 'no', 'gauss', 'Utopia', 'red', 'black', 'full', 0, 8, 6, 'mcp', 60, 1, 'on', 'off', 'mcp', [x_sigma, x_sigma_err, y_sigma, y_sigma_err], [x_pos_1, x_pos_err_1, y_pos_1, y_pos_err_1, x_pos_2, x_pos_err_2, y_pos_2, y_pos_err_2, fit_status_meaning[result.status()]])

    all_peak_positions = [file_name, x_pos_1, x_pos_err_1, y_pos_1, y_pos_err_1, x_pos_2, x_pos_err_2, y_pos_2, y_pos_err_2, x_sigma, x_sigma_err, y_sigma, y_sigma_err, ds.numEntries()]

    return(all_peak_positions)


if __name__ == '__main__':

    test = [9.406500000000001, 9.24965, 9.765, 9.674900000000001, 9.182500000000001, 8.67375, 9.524475, 8.8002, 9.616325, 9.702425, 9.52105, 53.2659, 8.822925, 9.29995, 10.3392, 9.138225, 8.67225, 9.375475, 9.2508, 9.5458, 8.991100000000001, 9.0881, 9.085975000000001, 9.040775, 9.76345, 9.656600000000001, 9.15975, 9.11665, 9.840250000000001, 9.3325, 9.464725, 9.391, 8.801300000000001, 9.52935, 9.256825000000001, 9.1248, 9.446875, 8.914825, 9.27045, 9.935725, 9.018075, 8.894400000000001, 9.27265, 8.848700000000001, 9.8879, 9.222025, 9.0929, 9.084475000000001, 8.69215, 9.252225000000001, 10.070975, 8.908900000000001, 9.052675, 9.016475, 9.298025, 9.360325, 9.075525, 9.552325, 9.479725, 9.1927, 9.582, 9.16505, 9.407175, 9.35, 9.567, 9.33755, 9.364475, 9.465275, 9.57255, 9.40795, 9.239625, 9.442675000000001, 9.34535, 10.46465, 9.53635, 8.616100000000001, 9.440975, 9.190475000000001, 9.145, 9.321675, 9.462375, 9.748700000000001, 9.29315, 8.767800000000001, 9.552475000000001, 9.6072, 9.375350000000001, 9.3581, 9.3147, 9.75065, 9.54175, 9.7756, 9.285525, 9.470425, 9.48715, 9.311250000000001, 27.996075, 9.64215, 9.708625, 9.357750000000001, 9.864600000000001, 9.189825, 9.300475, 9.2505, 9.479125, 9.579375, 9.756425, 9.3438, 9.460475, 9.38965, 10.184675, 9.16885, 9.392800000000001, 9.214025000000001, 9.369075, 26.2204, 9.55025, 9.1483, 9.914875, 8.721300000000001, 9.339925000000001, 9.501825, 10.506525, 9.2507, 9.092425, 9.525625, 9.161225, 9.57905, 9.389475000000001, 9.2775, 8.81535, 9.668725, 9.329025, 9.174975, 9.88685, 9.539225, 9.4773, 8.887175000000001, 10.57765, 9.27255, 8.7935, 9.541775000000001, 9.595775, 9.24735, 9.1987, 9.006300000000001, 9.350425000000001, 9.3507, 9.052100000000001, 9.038350000000001, 9.50985, 9.558575000000001, 9.60145, 9.105500000000001, 9.484075, 9.369950000000001, 9.775825000000001, 8.9411, 9.029875, 9.182450000000001, 9.742650000000001, 9.077825, 9.343200000000001, 26.920650000000002, 9.31485, 9.21175, 9.44945, 9.220375, 9.2276, 9.315975, 9.92295, 7.6025, 9.444975000000001, 41.330525, 9.512025000000001, 9.335825, 8.669275, 9.809225, 9.285175, 9.16215, 9.510925, 9.617575, 9.27675, 9.3831, 9.307925000000001, 9.307925000000001, 9.307925000000001, 9.307925000000001, 10.359375, 9.06615, 8.5494, 9.343325, 8.606125, 9.696900000000001, 9.494525000000001, 8.813925000000001, 9.106250000000001, 9.6441, 9.925500000000001, 9.713275000000001, 9.193975, 11.272725000000001, 9.61035, 9.057825000000001, 9.485175, 9.346225, 9.2848, 9.214425, 9.578875, 9.2712, 9.7359, 9.212150000000001, 8.998750000000001, 9.3125, 9.724475, 9.54955, 9.423925, 9.370800000000001, 9.33285, 9.221300000000001, 9.1896, 9.490450000000001, 9.002225000000001, 10.9859, 9.298175, 9.157575, 9.458525, 9.49285, 8.963625, 9.343775, 9.470600000000001, 9.401425, 9.512825000000001, 9.300600000000001, 12.150275, 9.585550000000001, 9.527475, 9.518650000000001, 9.5, 9.755725, 9.55725, 9.015425, 9.223025, 9.384025000000001, 9.308300000000001, 9.353325, 9.52285, 9.85725, 9.373275, 9.459325, 9.32595, 9.440575, 9.2598, 9.69505, 9.2644, 9.32085, 9.868350000000001, 9.0683, 9.787925000000001, 9.408100000000001, 9.505375, 9.153975, 9.742525, 10.746625, 9.73135, 9.662550000000001, 9.25375, 9.58455, 9.693800000000001, 9.311125, 9.846, 9.18995, 9.40375, 9.338225000000001, 9.299075, 184.748825, 9.214825000000001, 9.202625000000001, 9.723775, 9.2193, 9.270050000000001, 9.922075000000001, 9.50915, 8.849625, 9.7711, 9.20175, 9.251925, 9.34685, 9.187875, 9.63245, 10.069675, 9.098025, 9.0092, 9.697825, 9.142925, 9.042950000000001, 45.65305, 9.513, 9.444, 9.769350000000001, 9.52075, 9.38765, 64.344925, 9.680925, 9.264975, 9.241725, 9.301275, 9.939475, 9.436175, 9.11755, 8.324125, 9.06685, 9.5724, 9.5408, 8.202675000000001, 9.274475, 9.564275, 9.661475000000001, 9.507900000000001, 9.15305, 10.1395, 9.108600000000001, 9.0799, 9.2363, 9.18945, 9.2123, 9.511625, 9.516475, 10.164425, 9.54285, 9.2646, 9.4276, 9.6883, 8.569075, 9.8285, 9.4608, 194.24785, 9.16435, 9.1042, 9.192450000000001, 9.2927, 9.25995, 9.573625, 11.18525, 9.22525, 9.380700000000001, 9.462875, 9.6942, 9.789575000000001, 9.351600000000001, 9.44975, 8.9145, 10.9574, 9.35335, 9.481875, 8.963650000000001, 9.1516, 9.22565, 9.260875, 9.43125, 9.6074, 27.778475, 9.500975, 8.5469, 9.09765, 9.200275000000001, 9.765975000000001, 9.516275, 9.2448, 9.5907, 9.297600000000001, 9.118450000000001, 9.380625, 61.652225, 8.941075, 9.496875000000001, 9.134725000000001, 9.68265, 9.21325, 9.180075, 9.403450000000001, 9.509025000000001, 9.5082, 9.991375, 9.031, 9.408900000000001, 9.364975000000001, 44.16655, 9.04995, 9.514850000000001, 9.084525000000001, 10.19645, 9.43835, 32.254425000000005, 9.1907, 9.481475, 9.69705]
    root_1D_unbinned_gauss_fit(test)
