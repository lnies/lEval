
# ---------------------------------------------------------------------------
# Written by Jonas Karthein in 2016/2017/2018. Questions to jonas.karthein@cern.ch
# ---------------------------------------------------------------------------

import sys, os, time, json, math, glob, csv, multiprocessing, scipy.optimize, platform
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Qt5Agg')

import matplotlib.pyplot as plt
import array as ARRAY
import copy as COPY
import ROOT as root
from mpl_toolkits.mplot3d import axes3d
from collections import Counter

from read_write_functions import *
from python_plotter_functions import *
from ellipse_selector import *


def piicr_daq(data_folder_path, axis_scaling_factor, dont_show, pattern, sigma, bins, fit_range, cut_counter, tof_cut_on, tof_min, tof_max, trigger_on_off, clear_one, peak_finder_plot, all_peak_positions, hist2d_spines_off, hist2d_min_max, simple_plot_on, surface3D_plot_on, trigger_splitter, trigger_numbers, file_or_folder, file_name_input, color_map, analysis_old_new, z_class_analysis, FWHM_parameter_fix, one_or_two_peaks_to_fit, ellipse_selection):
    axis_scaling_factor = 0.031746032       # convert a.u. to mm
    origin_cartesian = [0, -31.5]       # in case the origin for a ring shape analysis is not at 0/0 change it here
    analysis_dict = {}
    counts_info_dict = {}

    if FWHM_parameter_fix == {}:
        FWHM_parameter_fix_dict = {'fixed': False}     # dict which controls the fixing of the FWHM in the spot_fit_root
    else:
        FWHM_parameter_fix_dict = FWHM_parameter_fix

    nll_plot = 'no'
    plot_fit = 'no'

    input_list = {}

    fit_range_save = fit_range

    all_peak_positions.append(['File name', 'X-Pos', 'X-Pos-error', 'Y-Pos', 'Y-Pos-error', 'X-FWHM', 'X-FWHM-error', 'Y-FWHM', 'Y-FWHM-error', 'Counts'])

    abspath = os.path.abspath(__file__)     # saves the folder of the script
    dname = os.path.dirname(abspath)


    if peak_finder_plot == 1:           # they don't work nicely together
        dont_show = 0

    exception_counter = 0               # continues if error occures, but counts them

    txt_files = []                      # looks for all files with .txt ending

    if file_or_folder == 'file':
        os.chdir(data_folder_path)
        txt_files.append(file_name_input)
    elif file_or_folder == 'folder':
        os.chdir(data_folder_path)
        for file in glob.glob("*.txt"):
            txt_files.append(file)
    print 'Files in the batch       :: ', txt_files
    x_fit_no_success = 0
    y_fit_no_success = 0

    path, folder_name = os.path.split(data_folder_path)  # saves the folder name

    time_info = []

    for i in txt_files:
        try:
            file_name = i[:-4]
            print('\n### Loading\n\nFile name                ::  %s\nPattern                  ::  %s') % (file_name, pattern)
            results, file_name, time_info_hilf, counts_info = load_and_pattern(file_name, pattern, analysis_old_new)

            counts_info_dict['{}'.format(file_name)] = {}
            counts_info_dict['{}'.format(file_name)]['#ej'] = counts_info[0]           # total number of ejections in file --> p1+p2
            counts_info_dict['{}'.format(file_name)]['#signals'] = counts_info[1]      # maximal number of signals arriving at the delay lines

            hilf_153053 = [file_name[:-3]]
            hilf_153053.extend(time_info_hilf)
            time_info.append(hilf_153053)

            try:                      # when a window with proper data is alread calculated
                window_x
            except NameError:
                window_x = 0
                window_y = 0
                window_xy = 0
                err_window_x = 0
                err_window_y = 0
                err_window_xy = 0
    # ------------------------
    # measurement window guess (only from first file)
    # ------------------------
            window_x, window_y, window_xy, err_window_x, err_window_y, err_window_xy, sum_x, sum_y = window(clear_one, results, window_x, window_y, window_xy, err_window_x, err_window_y, err_window_xy, data_folder_path, file_name)
            # print 'Guess:          X-Window = %s (%s) ; Y-Window = %s (%s) ; X-Y-Window = %s (%s)' % (window_x, err_window_x, window_y, err_window_y, window_xy, err_window_xy)

    # ------------------------
    # spot calculation without reconstruction
    # ------------------------
            if trigger_splitter == 0:
                hilf_file_name = file_name

                input_list['trigger_on_off'] = trigger_on_off
                input_list['tof_cut_on'] = tof_cut_on
                input_list['tof_min'] = tof_min
                input_list['tof_max'] = tof_max
                input_list['window_x'] = window_x
                input_list['window_y'] = window_y
                input_list['err_window_x'] = err_window_x
                input_list['err_window_y'] = err_window_y
                input_list['window_xy'] = window_xy
                input_list['err_window_xy'] = err_window_xy
                input_list['sigma'] = sigma

                input_list['results'] = results

                spot_positions_rec = multiprocess_rec_0(input_list)
                counts_info_dict['{}'.format(file_name)]['#rec-spots'] = len(spot_positions_rec)      # number of reconstructed spots
                counts_info_dict['{}'.format(file_name)]['#rec-spots-distribution'] = count_distribution(spot_positions_rec)

                # spot_positions, check_event_taken, double_taken_event_counter, temp_results = rec_0(results, trigger_on_off, tof_cut_on, tof_min, tof_max, 1600, 1600, 30, 30, 3200, 50, 4)     # weird noise --> fixed window! 27.04.2017
                # spot_positions, check_event_taken, double_taken_event_counter, temp_results = rec_0(results, trigger_on_off, tof_cut_on, tof_min, tof_max, window_x, window_y, err_window_x, err_window_y, window_xy, err_window_xy, sigma)
                # window_x, window_y, window_xy, err_window_x, err_window_y, err_window_xy, sum_x, sum_y = window(1, temp_results, window_x, window_y, window_xy, err_window_x, err_window_y, err_window_xy, data_folder_path, file_name)          # calculates the time windows again with the new, accepted events


        # ------------------------
        # Z-class-analysis
        # ------------------------
                if z_class_analysis[0] == 'yes':
                    temp_len = len(spot_positions_rec)
                    spot_positions_z = z_class_reduction(spot_positions_rec, z_class_analysis)
                    print 'Z-class reduction ({}-{})  ::  {}/{}'.format(str(int(z_class_analysis[1])), str(int(z_class_analysis[2])), len(spot_positions_z), temp_len)
                else:
                    spot_positions_z = spot_positions_rec

                counts_info_dict['{}'.format(file_name)]['#z-spots'] = len(spot_positions_z)      # number of reconstructed spots
                counts_info_dict['{}'.format(file_name)]['#z-spots-distribution'] = count_distribution(spot_positions_z)


        # ------------------------
        # ToF fit/cut
        # ------------------------


                if tof_cut_on == 0:
                    tof_range_fit = root_1D_unbinned_gauss_fit([x[3] for x in spot_positions_rec], file_name)
                    print 'ToF window fit           ::  (%3.2f +/- %3.2f) us' % (tof_range_fit[0], 2.5 * tof_range_fit[2])
                    tof_min = tof_range_fit[0] - 2.5 * tof_range_fit[2]
                    tof_max = tof_range_fit[0] + 2.5 * tof_range_fit[2]
                    tof_cut_on = 1

                spot_positions = []
                for spot in spot_positions_z:
                    if spot[3] < tof_max and spot[3] > tof_min:
                        spot_positions.append(spot)

                tof_plot(spot_positions, 50, file_name, dont_show)


                counts_info_dict['{}'.format(file_name)]['#tof-spots'] = len(spot_positions)      # number of reconstructed spots
                counts_info_dict['{}'.format(file_name)]['#tof-spots-distribution'] = count_distribution(spot_positions)
                counts_info_dict['{}'.format(file_name)]['percent-contamination'] = (float(len(spot_positions_z))-float(len(spot_positions)))/float(len(spot_positions_z))      # number of reconstructed spots


        # ------------------------
        # convert coordinates also to polar coordinates and save spot positions
        # ------------------------
                x_cartesian = []
                y_cartesian = []
                # for i in range(len(spot_positions)):
                #     x_cartesian.append(spot_positions[i][1] - origin_cartesian[0])
                #     y_cartesian.append(spot_positions[i][2] - origin_cartesian[1])
                # r_polar, phi_polar = convert_cartesian_to_polar(x_cartesian, y_cartesian)
                # spot_positions_polar = [['Event #', 'Radius / mm', 'Angle / radian', 'ToF / ns']]
                # for j in range(len(r_polar)):
                #     spot_positions_polar.append([spot_positions[j][0], float(r_polar[j]) * axis_scaling_factor, phi_polar[j], spot_positions[j][3]])
                # write_csv(data_folder_path, spot_positions_polar, '%s_spot_positions_polar' % (file_name))

        # # ------------------------
        # # Z-class-analysis
        # # ------------------------
        #         if z_class_analysis[0] == 'yes':
        #             temp_len = len(spot_positions)
        #             spot_positions = z_class_reduction(spot_positions, z_class_analysis)
        #             print 'Z-class reduction ({}-{})  ::  {}/{}'.format(str(int(z_class_analysis[1])), str(int(z_class_analysis[2])), len(spot_positions), temp_len)

        # ------------------------
        # Ellipse position selection
        # ------------------------
                if ellipse_selection == True:
                    if one_or_two_peaks_to_fit == 'one':
                        es = Ellipse_selector()
                        analysis_dict = es.main(spot_positions, file_name, [])
                        spot_positions = analysis_dict['cut_data']
                        write_csv(data_folder_path, spot_positions, '%s_spot_positions' % (file_name))
                        counts_info_dict['{}'.format(file_name)]['#manual-spots'] = len(spot_positions)      # number of reconstructed spots
                        counts_info_dict['{}'.format(file_name)]['#manual-spots-distribution'] = count_distribution(spot_positions)

                    elif one_or_two_peaks_to_fit == 'two':
                        es_1 = Ellipse_selector()
                        analysis_dict = es_1.main(spot_positions, file_name, [])
                        spot_positions_dominant = analysis_dict['cut_data']      # dominant spot
                        write_csv(data_folder_path, spot_positions_dominant, '%s_spot_positions_dominant' % (file_name))
                        counts_info_dict['{}'.format(file_name)]['#manual-spots-ground'] = len(spot_positions_dominant)      # number of reconstructed spots
                        counts_info_dict['{}'.format(file_name)]['#manual-spots-ground-distribution'] = count_distribution(spot_positions_dominant)

                        time.sleep(0.1)
                        es_2 = Ellipse_selector()
                        analysis_dict = es_2.main(spot_positions, file_name, [])
                        spot_positions_recessive = analysis_dict['cut_data']      # recessive spot
                        write_csv(data_folder_path, spot_positions_recessive, '%s_spot_positions_recessive' % (file_name))
                        counts_info_dict['{}'.format(file_name)]['#manual-spots-isomer'] = len(spot_positions_recessive)      # number of reconstructed spots
                        counts_info_dict['{}'.format(file_name)]['#manual-spots-isomer-distribution'] = count_distribution(spot_positions_recessive)
                        counts_info_dict['{}'.format(file_name)]['percent-ground'] = float(len(spot_positions_dominant))/(float(len(spot_positions_recessive))+float(len(spot_positions_dominant)))      # number of reconstructed spots
                        counts_info_dict['{}'.format(file_name)]['percent-isomer'] = float(len(spot_positions_recessive))/(float(len(spot_positions_recessive))+float(len(spot_positions_dominant)))      # number of reconstructed spots


                else:
                    write_csv(data_folder_path, spot_positions, '%s_spot_positions' % (file_name))


        # ------------------------
        # measurement window calculation from data
        # ------------------------
                # print 'Calculation:    X-Window = %s (%s) ; Y-Window = %s (%s) ; X-Y-Window = %s (%s)\n' % (window_x, round(err_window_x), window_y, round(err_window_y), window_xy, round(err_window_xy))

                window_x = COPY.copy(window_x)
                window_y = COPY.copy(window_y)
                window_xy = COPY.copy(window_xy)
                err_window_x = round(COPY.copy(err_window_x))
                err_window_y = round(COPY.copy(err_window_y))
                err_window_xy = round(COPY.copy(err_window_xy))
        # ------------------------
        # calculation of informative properties of the data
        # ------------------------
                min_event_number = 0
                max_event_number = 0
                missing_event_counter = 0
                zero_event_counter = 0

                for i in range(0, len(results), 5):                     # check each ejection
                    min_event_number += min(results[i][0], results[i+1][0], results[i+2][0], results[i+3][0])
                for i in range(0, len(results), 5):
                    if results[i][0] != results[i+1][0] or results[i+1][0] != results[i+2][0] or results[i+2][0] != results[i+3][0]:
                        missing_event_counter += 1
                for i in range(0, len(results), 5):
                    max_event_number += max(results[i][0], results[i+1][0], results[i+2][0], results[i+3][0])
                for i in range(0, len(results), 5):
                    if results[i][0] == 0 or results[i+1][0] == 0 or results[i+2][0] == 0 or results[i+3][0] == 0:
                        zero_event_counter += 1
        # ------------------------
        # Try to reconstruct missing events
        # ------------------------
                # spot_positions_reconstructed, check_event_taken_rec, double_taken_event_counter_rec1 = rec_1(results, trigger_on_off, tof_cut_on, tof_min, tof_max, window_x, window_y, err_window_x, err_window_y, window_xy, err_window_xy, sigma, check_event_taken)
                # spot_positions_reconstructed_2, check_event_taken_rec_2, double_taken_event_counter_rec2 = rec_2(results, trigger_on_off, tof_cut_on, tof_min, tof_max, window_x, window_y, err_window_x, err_window_y, window_xy, err_window_xy, sigma, check_event_taken_rec)
                # print '\nClear events: %s ; Reconstruction I: %s ; Reconstruction II: %s ; Maximum possible events: %s\n' % (len(spot_positions), len(spot_positions_reconstructed), len(spot_positions_reconstructed_2), max_event_number)
        # ------------------------
        # Histogram plot normal + zoomed
        # ------------------------
                nice_plots_but_problems_in_fitting = 1
                # H, xedges, yedges, xs, ys = histo2D(spot_positions, file_name, pattern, 1, hist2d_min_max, nice_plots_but_problems_in_fitting, bins, hist2d_spines_off, axis_scaling_factor, dont_show, color_map)
                nice_plots_but_problems_in_fitting = 0  # This procedure allows plotting both a 2Dhistogram of the detector picture as well as a zoomed image for better fitting.
                # H, xedges, yedges, xs, ys = histo2D(spot_positions, file_name, pattern, 1, hist2d_min_max, nice_plots_but_problems_in_fitting, bins, hist2d_spines_off, axis_scaling_factor, dont_show, color_map)

                if trigger_on_off == 1 and one_or_two_peaks_to_fit == 'one':
                   # tof_plot(spot_positions, 50, file_name, dont_show)
                   ion_cloud_3d(spot_positions, axis_scaling_factor, dont_show, file_name)

        # ------------------------
        # Simple plot if wanted
        # ------------------------
                # if simple_plot_on == 1:
                #     simple_plot(xs, ys, file_name, pattern, 1, dont_show)
        # ------------------------
        # Projections + Fit
        # ------------------------
                # hist_x, hist_y, hist_x_weights, hist_y_weights = hist1D(H, bins, 1, x_fit_min, x_fit_max, y_fit_min, y_fit_max, dont_show, xedges, yedges)
                # all_peak_positions, histo_number_of_peaks_x, histo_number_of_peaks_y, histo_peak_height_x, histo_peak_height_y, histo_peak_pos_x, histo_peak_pos_y, x_fit_no_success, y_fit_no_success, x_counts, y_counts = projection_fitter(hist_x, hist_y, hist_x_weights, hist_y_weights, x_fit_min, x_fit_max, y_fit_min, y_fit_max, xedges, yedges, bins, peak_finder_plot, file_name, cut_counter, dont_show, all_peak_positions, axis_scaling_factor, x_fit_no_success, y_fit_no_success)
        # ------------------------
        # 3D
        # ------------------------
                if surface3D_plot_on == 1 and one_or_two_peaks_to_fit == 'one':
                    surface3D_fit(xedges, yedges, H, histo_number_of_peaks_x, histo_number_of_peaks_y, histo_peak_height_x, histo_peak_height_y, histo_peak_pos_x, histo_peak_pos_y, file_name, dont_show, bins)

        # ------------------------
        # Mean of the wanted spots
        # ------------------------

                if one_or_two_peaks_to_fit == 'one':

                    hilf_xx = []
                    hilf_yy = []
                    if fit_range == [-750, 750, -750, 750]:     # calculate the mean of the wanted spot
                        for i in range(len(spot_positions)):
                            if np.sqrt(spot_positions[i][1] ** 2 + spot_positions[i][2] ** 2) < 650:
                                hilf_xx.append(spot_positions[i][1])
                                hilf_yy.append(spot_positions[i][2])
                        mean_x = np.mean(hilf_xx)
                        mean_y = np.mean(hilf_yy)
                    else:
                        for i in range(len(spot_positions)):
                            if spot_positions[i][1] > fit_range[0] and spot_positions[i][1] < fit_range[1]:
                                hilf_xx.append(spot_positions[i][1])
                            if spot_positions[i][2] > fit_range[2] and spot_positions[i][2] < fit_range[3]:
                                hilf_yy.append(spot_positions[i][2])
                        mean_x = np.mean(hilf_xx)
                        mean_y = np.mean(hilf_yy)


                elif one_or_two_peaks_to_fit == 'two':
                    for spot_pos in [spot_positions_dominant, spot_positions_recessive]:
                        hilf_xx = []
                        hilf_yy = []
                        for i in range(len(spot_pos)):
                            if np.sqrt(spot_pos[i][1] ** 2 + spot_pos[i][2] ** 2) < 650:
                                hilf_xx.append(spot_pos[i][1])
                                hilf_yy.append(spot_pos[i][2])
                        if spot_pos == spot_positions_dominant:
                            mean_x_dom = np.mean(hilf_xx)
                            mean_y_dom = np.mean(hilf_yy)
                        else:
                            mean_x_rec = np.mean(hilf_xx)
                            mean_y_rec = np.mean(hilf_yy)

        # ------------------------
        # Interative fitting
        # ------------------------

                if one_or_two_peaks_to_fit == 'one':
                    if fit_range == [-750, 750, -750, 750]:     # cuts away the stattered ions
                        fit_range_new = [mean_x-200, mean_x+200, mean_y-200, mean_y+200]

                        first_fit_info = spot_fit_root(file_name, spot_positions, all_peak_positions, mean_x, mean_y, fit_range_new, 'no', 'no', 1, FWHM_parameter_fix_dict)
                        second_fit_info = spot_fit_root(file_name, spot_positions, all_peak_positions, first_fit_info[1], first_fit_info[3], [first_fit_info[1] - 1.8 * first_fit_info[5], first_fit_info[1] + 1.8 * first_fit_info[5], first_fit_info[3] - 1.8 * first_fit_info[7], first_fit_info[3] + 1.8 * first_fit_info[7]], 'no', 'no', 2, FWHM_parameter_fix_dict)
                        all_peak_positions.append(spot_fit_root(file_name, spot_positions, all_peak_positions, second_fit_info[1], second_fit_info[3], [second_fit_info[1] - 1.8 * second_fit_info[5], second_fit_info[1] + 1.8 * second_fit_info[5], second_fit_info[3] - 1.8 * second_fit_info[7], second_fit_info[3] + 1.8 * second_fit_info[7]], 'no', 'yes', 3, FWHM_parameter_fix_dict))
                    else:
                        all_peak_positions.append(spot_fit_root(file_name, spot_positions, all_peak_positions, mean_x, mean_y, fit_range, nll_plot, 'yes', 0, FWHM_parameter_fix_dict))


                elif one_or_two_peaks_to_fit == 'two':

                    first_fit_info = two_spot_fit_root(file_name, all_peak_positions, mean_x_dom, mean_y_dom, mean_x_rec, mean_y_rec, 'no', 1, FWHM_parameter_fix_dict)

                    second_fit_info = two_spot_fit_root(file_name, all_peak_positions, first_fit_info[1], first_fit_info[3], first_fit_info[5], first_fit_info[7], 'no', 2, FWHM_parameter_fix_dict)

                    all_peak_positions.append(two_spot_fit_root(file_name, all_peak_positions, second_fit_info[1], second_fit_info[3], second_fit_info[5], second_fit_info[7], 'yes', 3, FWHM_parameter_fix_dict))

                    all_peak_positions[0] = ['File name', 'X-Pos-dom', 'X-Pos-dom-error', 'Y-Pos-dom', 'Y-Pos--domerror', 'X-Pos-rec', 'X-Pos-rec-error', 'Y-Pos-rec', 'Y-Pos-rec-error', 'X-FWHM', 'X-FWHM-error', 'Y-FWHM', 'Y-FWHM-error', 'Counts']
                    # print all_peak_positions

        # ------------------------
        # Conversion to polar coordinates
        # ------------------------

                if one_or_two_peaks_to_fit == 'one':
                    for i in range(len(spot_positions)):        # convert to polar coordinates
                        x_cartesian.append(spot_positions[i][1] - all_peak_positions[-1][1])
                        y_cartesian.append(spot_positions[i][2] - all_peak_positions[-1][3])
                    r_polar, phi_polar = convert_cartesian_to_polar(x_cartesian, y_cartesian)
                    spot_positions_polar = [['Event #', 'Radius / mm', 'Angle / radian', 'ToF / ns'], [all_peak_positions[-1][1], all_peak_positions[-1][3]]]
                    for j in range(len(r_polar)):
                        spot_positions_polar.append([spot_positions[j][0], float(r_polar[j]) * axis_scaling_factor, phi_polar[j], spot_positions[j][3]])
                    write_csv(data_folder_path, spot_positions_polar, '%s_spot_positions_polar' % (file_name))




















# ------------------------
# Option: cutting one file into smaller pieces
# ------------------------



            elif trigger_splitter == 1:
                if trigger_numbers > int(len(results)/5):
                    hilf_file_name = file_name
                    pass
                else:
                    cut_counter = 0
                    if pattern == 1 or pattern == 2:
                        trigger_numbers = trigger_numbers/2
                    for i in range(0, (len(results)-trigger_numbers*5), trigger_numbers*5):
                        # all_peak_positions = []
                        cut_counter += 1
                        hilf_file_name = file_name + str('_part%s') % (cut_counter)
                        cut_results = []
                        for j in range(0, ((trigger_numbers)*5), 1):
                            cut_results.append(results[i+j])
                        # print len(cut_results)

                        input_list['trigger_on_off'] = trigger_on_off
                        input_list['tof_cut_on'] = tof_cut_on
                        input_list['tof_min'] = tof_min
                        input_list['tof_max'] = tof_max
                        input_list['window_x'] = window_x
                        input_list['window_y'] = window_y
                        input_list['err_window_x'] = err_window_x
                        input_list['err_window_y'] = err_window_y
                        input_list['window_xy'] = window_xy
                        input_list['err_window_xy'] = err_window_xy
                        input_list['sigma'] = sigma

                        input_list['results'] = cut_results


                        spot_positions_rec = multiprocess_rec_0(input_list)

                        if ellipse_selection == True:
                            es = Ellipse_selector()
                            analysis_dict = es.main(spot_positions_rec, file_name, [])
                            spot_positions_rec = analysis_dict['cut_data']

                        # spot_positions, check_event_taken, double_taken_event_counter, temp_results = rec_0(cut_results, trigger_on_off, tof_cut_on, tof_min, tof_max, window_x, window_y, err_window_x, err_window_y, window_xy, err_window_xy, sigma)
                        # window_x, window_y, window_xy, err_window_x, err_window_y, err_window_xy, sum_x, sum_y = window(1, temp_results, window_x, window_y, window_xy, err_window_x, err_window_y, err_window_xy, data_folder_path, file_name)          # calculates the time windows again with the new, accepted events

                # ------------------------
                # ToF fit/cut
                # ------------------------
                        if tof_cut_on == 0:
                            tof_range_fit = root_1D_unbinned_gauss_fit([x[3] for x in spot_positions_rec], file_name)
                            print 'Auto-ToF-window: (%3.2f +/- %3.2f) us' % (tof_range_fit[0], 2.2 * tof_range_fit[2])
                            tof_min = tof_range_fit[0] - 2.2 * tof_range_fit[2]
                            tof_max = tof_range_fit[0] + 2.2 * tof_range_fit[2]
                            tof_cut_on = 1

                        spot_positions = []
                        for spot in spot_positions_rec:
                            if spot[3] < tof_max and spot[3] > tof_min:
                                spot_positions.append(spot)

                        tof_plot(spot_positions, 50, file_name, dont_show)

                # ------------------------
                # Z-class-analysis
                # ------------------------
                        if z_class_analysis[0] == 'yes':
                            spot_positions = z_class_reduction(spot_positions, z_class_analysis)

                # ------------------------
                # measurement window calculation from data
                # ------------------------
                        write_csv(data_folder_path, spot_positions, '%s_spot_positions' % (hilf_file_name))
                        # with open('%s_spot_positions.csv' % (hilf_file_name), 'wb') as f:   # save all accepted positions into csv
                        #     writer = csv.writer(f)
                        #     writer.writerows(spot_positions)
                        # print 'Calculation:    X-Window = %s (%s) ; Y-Window = %s (%s) ; X-Y-Window = %s (%s)\n' % (window_x, round(err_window_x), window_y, round(err_window_y), window_xy, round(err_window_xy))

                        window_x = COPY.copy(window_x)
                        window_y = COPY.copy(window_y)
                        window_xy = COPY.copy(window_xy)
                        err_window_x = round(COPY.copy(err_window_x))
                        err_window_y = round(COPY.copy(err_window_y))
                        err_window_xy = round(COPY.copy(err_window_xy))
                # ------------------------
                # calculation of informative properties of the data
                # ------------------------
                        min_event_number = 0
                        max_event_number = 0
                        missing_event_counter = 0
                        zero_event_counter = 0

                        for i in range(0, len(results), 5):                     # check each ejection
                            min_event_number += min(results[i][0], results[i+1][0], results[i+2][0], results[i+3][0])
                        for i in range(0, len(results), 5):
                            if results[i][0] != results[i+1][0] or results[i+1][0] != results[i+2][0] or results[i+2][0] != results[i+3][0]:
                                missing_event_counter += 1
                        for i in range(0, len(results), 5):
                            max_event_number += max(results[i][0], results[i+1][0], results[i+2][0], results[i+3][0])
                        for i in range(0, len(results), 5):
                            if results[i][0] == 0 or results[i+1][0] == 0 or results[i+2][0] == 0 or results[i+3][0] == 0:
                                zero_event_counter += 1
                # ------------------------
                # Try to reconstruct missing events
                # ------------------------
                        # spot_positions_reconstructed, check_event_taken_rec, double_taken_event_counter_rec1 = rec_1(cut_results, trigger_on_off, tof_cut_on, tof_min, tof_max, window_x, window_y, err_window_x, err_window_y, window_xy, err_window_xy, sigma, check_event_taken)
                        # spot_positions_reconstructed_2, check_event_taken_rec_2, double_taken_event_counter_rec2 = rec_2(cut_results, trigger_on_off, tof_cut_on, tof_min, tof_max, window_x, window_y, err_window_x, err_window_y, window_xy, err_window_xy, sigma, check_event_taken_rec)
                        # print '\nClear events: %s ; Reconstruction I: %s ; Reconstruction II: %s ; Maximum possible events: %s\n' % (len(spot_positions), len(spot_positions_reconstructed), len(spot_positions_reconstructed_2), max_event_number)
                # ------------------------
                # Histogram plot normal + zoomed
                # ------------------------
                        nice_plots_but_problems_in_fitting = 1
                        # H, xedges, yedges, xs, ys = histo2D(spot_positions, hilf_file_name, pattern, 1, hist2d_min_max, nice_plots_but_problems_in_fitting, bins, hist2d_spines_off, axis_scaling_factor, dont_show, color_map)
                        nice_plots_but_problems_in_fitting = 0  # This procedure allows plotting both a 2Dhistogram of the detector picture as well as a zoomed image for better fitting.
                        # H, xedges, yedges, xs, ys = histo2D(spot_positions, hilf_file_name, pattern, 1, hist2d_min_max, nice_plots_but_problems_in_fitting, bins, hist2d_spines_off, axis_scaling_factor, dont_show, color_map)

                        if trigger_on_off == 1:
                            # tof_plot(spot_positions, bins, hilf_file_name, dont_show)
                            ion_cloud_3d(spot_positions, axis_scaling_factor, dont_show, hilf_file_name)
                # ------------------------
                # Simple plot if wanted
                # ------------------------
                        # if simple_plot_on == 1:
                        #     simple_plot(xs, ys, hilf_file_name, pattern, 1, dont_show)
                # ------------------------
                # Projections + Fit
                # ------------------------
                        # hist_x, hist_y, hist_x_weights, hist_y_weights = hist1D(H, bins, 1, x_fit_min, x_fit_max, y_fit_min, y_fit_max, dont_show, xedges, yedges)
                        # all_peak_positions, histo_number_of_peaks_x, histo_number_of_peaks_y, histo_peak_height_x, histo_peak_height_y, histo_peak_pos_x, histo_peak_pos_y, x_fit_no_success, y_fit_no_success, x_counts, y_counts = projection_fitter(hist_x, hist_y, hist_x_weights, hist_y_weights, x_fit_min, x_fit_max, y_fit_min, y_fit_max, xedges, yedges, bins, peak_finder_plot, hilf_file_name, cut_counter, dont_show, all_peak_positions, axis_scaling_factor, x_fit_no_success, y_fit_no_success)
                # ------------------------
                # 3D + Fit
                # ------------------------
                        # if surface3D_plot_on == 1:
                        #     surface3D_fit(xedges, yedges, H, histo_number_of_peaks_x, histo_number_of_peaks_y, histo_peak_height_x, histo_peak_height_y, histo_peak_pos_x, histo_peak_pos_y, hilf_file_name, dont_show, bins)
                        hilf_xx = []
                        hilf_yy = []
                        x_cartesian = []
                        y_cartesian = []


                        if one_or_two_peaks_to_fit == 'one':

                            if fit_range == [-750, 750, -750, 750]:     # calculate the mean of the wanted spot
                                for i in range(len(spot_positions)):
                                    if np.sqrt(spot_positions[i][1] ** 2 + spot_positions[i][2] ** 2) < 650:
                                        hilf_xx.append(spot_positions[i][1])
                                        hilf_yy.append(spot_positions[i][2])
                                mean_x = np.mean(hilf_xx)
                                mean_y = np.mean(hilf_yy)
                            else:
                                for i in range(len(spot_positions)):
                                    if spot_positions[i][1] > fit_range[0] and spot_positions[i][1] < fit_range[1]:
                                        hilf_xx.append(spot_positions[i][1])
                                    if spot_positions[i][2] > fit_range[2] and spot_positions[i][2] < fit_range[3]:
                                        hilf_yy.append(spot_positions[i][2])
                                mean_x = np.mean(hilf_xx)
                                mean_y = np.mean(hilf_yy)


                            if fit_range == [-750, 750, -750, 750]:     # cuts away the stattered ions
                                fit_range_new = [mean_x-200, mean_x+200, mean_y-200, mean_y+200]

                                first_fit_info = spot_fit_root(file_name, spot_positions, all_peak_positions, mean_x, mean_y, fit_range_new, 'no', 'no', 1, FWHM_parameter_fix_dict)
                                second_fit_info = spot_fit_root(file_name, spot_positions, all_peak_positions, first_fit_info[1], first_fit_info[3], [first_fit_info[1] - 1.8 * first_fit_info[5], first_fit_info[1] + 1.8 * first_fit_info[5], first_fit_info[3] - 1.8 * first_fit_info[7], first_fit_info[3] + 1.8 * first_fit_info[7]], 'no', 'no', 2, FWHM_parameter_fix_dict)
                                all_peak_positions.append(spot_fit_root(file_name, spot_positions, all_peak_positions, second_fit_info[1], second_fit_info[3], [second_fit_info[1] - 1.8 * second_fit_info[5], second_fit_info[1] + 1.8 * second_fit_info[5], second_fit_info[3] - 1.8 * second_fit_info[7], second_fit_info[3] + 1.8 * second_fit_info[7]], 'no', 'yes', 3, FWHM_parameter_fix_dict))
                            else:
                                all_peak_positions.append(spot_fit_root(file_name, spot_positions, all_peak_positions, mean_x, mean_y, fit_range, nll_plot, 'yes', 0, FWHM_parameter_fix_dict))

                        elif one_or_two_peaks_to_fit == 'two':

                            sys.exit('This is not yet fully implemented.')

                            first_fit_info = two_spot_fit_root(file_name, spot_positions, all_peak_positions, mean_x_1, mean_y_1, mean_x_2, mean_y_2, 'no', 1, FWHM_parameter_fix_dict)

                            second_fit_info = two_spot_fit_root(file_name, spot_positions, all_peak_positions, first_fit_info[1], first_fit_info[3], first_fit_info[5], first_fit_info[7], 'no', 2, FWHM_parameter_fix_dict)

                            all_peak_positions.append(two_spot_fit_root(file_name, spot_positions, all_peak_positions, second_fit_info[1], second_fit_info[3], second_fit_info[5], second_fit_info[7], 'yes', 3, FWHM_parameter_fix_dict))






                        for i in range(len(spot_positions)):
                            x_cartesian.append(spot_positions[i][1] - all_peak_positions[-1][1])
                            y_cartesian.append(spot_positions[i][2] - all_peak_positions[-1][3])
                        r_polar, phi_polar = convert_cartesian_to_polar(x_cartesian, y_cartesian)
                        spot_positions_polar = [['Event #', 'Radius / mm', 'Angle / radian', 'ToF / ns'], [all_peak_positions[-1][1], all_peak_positions[-1][3]]]
                        for j in range(len(r_polar)):
                            spot_positions_polar.append([spot_positions[j][0], float(r_polar[j]) * axis_scaling_factor, phi_polar[j], spot_positions[j][3]])
                        write_csv(data_folder_path, spot_positions_polar, '%s_spot_positions_polar' % (file_name))


    # ------------------------
    # Error output
    # ------------------------
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            exception_counter += 1
            print '\nError appeared in %s.txt, file will be skipped and batch process continues.' % file_name
            print 'Python error: %s \n"%s" (line %s)' % (exc_type, e, exc_tb.tb_lineno)
            pass
    # ------------------------
    # Save batch results
    # ------------------------
    if file_or_folder == 'folder':
        save_fit_info(fit_range_save, pattern, folder_name, all_peak_positions, cut_counter, data_folder_path)
    elif file_or_folder == 'file':
        save_fit_info(fit_range_save, pattern, file_name, all_peak_positions, cut_counter, data_folder_path)
    if analysis_old_new == 'new':
        write_csv(data_folder_path, time_info, '_time_info_%s' % (folder_name))

    with open('counts_info_dict.json', 'w') as f:
        json.dump(counts_info_dict,f)

    # print 'Number of files which have been skipped due to an error: %s' % exception_counter
    # print 'Number of unsuccessfull fits: X = %s ; Y = %s' %(x_fit_no_success, y_fit_no_success)
    # ------------------------
    # change the folder back to the folder of the script
    # ------------------------
    os.chdir(dname)

    return(all_peak_positions)


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
            window_x = COPY.copy(window_x)
            window_y = COPY.copy(window_y)
            window_xy = COPY.copy(window_xy)
            err_window_x = COPY.copy(err_window_x)
            err_window_y = COPY.copy(err_window_y)
            err_window_xy = COPY.copy(err_window_xy)
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
    check_event_taken_rec = COPY.copy(check_event_taken)
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
    check_event_taken_rec_2 = COPY.copy(check_event_taken_rec)
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
    if platform.system() == 'Darwin':
        tree = root.TTree( 'tree', 'tree' )

        x = ARRAY.array('d', [ 0. ])
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
        return([x_pos, x_pos_err, x_sigma, x_sigma_err])
    else:
        return([mean_x,0.001,std_x,0.001])


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
    ### ALWAYS CHANGE. READ "READ-ME-first" FIRST!
    upper_folder_path = '/Volumes/dfs/DATA/2017/2017-06-Cd-run_prep/PI-ICR_Scans/Ret1_vs_Ret2/133Cs/p1p2' # 'E:\\Jonas\\2017-05-13_85Rb_88Sr_88Rb_600ms_002\\85Rb\\p1p2' # '/Users/jonaskarthein/cernbox/Python/2017-05-05-Harmonic_and_long_run' # 'E:\\Welker\\Harmonic_and_long_run' # Windoof: 'F:\\2017-05-09_85Rb_Ret1_UT3' # Macintosh: '/Volumes/Macintosh HD/Users/jonaskarthein/cernbox/Python/8587Rb_over_night_26.04.17_failed/85Rb/p1p2'
    latex_name = '^{133}$Cs'

    ### ONLY CHANGE IF YOU KNOW WHAT YOU DO!
    axis_scaling_factor = 0.031746032
    dont_show = 1           # if =1 the plots are not displayed, for showing =0! file_name = '133Cs_c_091'  # if only one file should be looked at, uncomment! if peak_finder_plot = 1 this number will automatically set to 0 since they don't work together.
    pattern = 3  # chose if you want pattern one (='1'; p1+p2 in file) or pattern two (='2'; p1+p2 in file) should be calculated. If only one pattern is saved in the file (e.g. the center) just put a '3' there
    sigma = 15               # changes the measurement window; eqal the number of sigma of the gaussian distribution calculated for the totaltime
    bins = 60  # 85         # changes the bin number per axis of the histogram plot
    fit_range = [-750, 750, -750, 750]	# vector containing [0] = cut_x_min, [1] = cut_x_max, [2] = cut_y_min, [3] = cut_y_max, default: fit_range = [-750, 750, -750, 750]
    cut_counter = 1         # adds a nametag to certain files when you enable cutting
    tof_cut_on = 0          # if you want a TOF cut with TOF_min/max, set it to '1', if you want to see the whole data and the whole TOF range set it to '0'
    tof_min = 47            # only valid if tof_cut_on == 1! unit is microseconds
    tof_max = 57            # only valid if tof_cut_on == 1! unit is microseconds
    trigger_on_off = 1      # when triggered (-> MCP time != 0) then trigger = 1, else set it to zero (e.g. for darkcount tests)
    clear_one = 0           # when there are engough events with only 1 event at each delayline set it to 1, elso to zero (if 1 time window gets calculated)
    peak_finder_plot = 0    # when neither 1 nor 3 peaks for both could be found an error will appear (e.g. in line 485 no parameter found). To better check and correct for that (adjust the peakfinder function in width or in height) set peak_finder_plot to one and you will get a plot with the found peaks. Normally you don't need that, so set it then to 0
    all_peak_positions = [] # list to save all fit positions. Ordering: file_name, X-Fit-Pos, X-Fit-Pos-err, Y-Fit-Pos, Y-Fit-Pos-Err
    hist2d_spines_off = 1   # if you want the axis lines for the 2d Histogram to be switched off set this variable to '1', if not set it to '0'
    hist2d_min_max = 750    # this variable sets the +- axis range for the 2d Histogram. Set 750 for the whole MCP
    simple_plot_on = 0      # if you want an additional plot that dispays the whole detector area set it to '1', otherwise to '0'
    surface3D_plot_on = 0   # if you want a 3D plot set it to '1', otherwise to '0'
    trigger_splitter = 0    # if = '0' it will go through all files in one calculation. if = '1' it will cut the raw data in portions of #trigger_number and calculates for each
    trigger_numbers = 25   # if trigger_splitter = '1' it will cut the raw data in portions of #trigger_numbers trigger.
    color_map = 2           # 1 for blue, 2 for multicolor
    plot_type = 'scatter'   # 'scatter' or 'histogram'
    file_or_folder = 'folder'               # in case you just want to analyze only one single file e.g. for cutting an isomer out
    file_name_input = '85Rb_002.txt'    # if file_or_folder == 'file' input here the fil name WITH the ending .txt
    mode = 'c'  #['c', 'p1', 'p2']
    analysis = 'new'
    z_class_analysis = ['no', 0, 5]     # list with three entries: string z-class-analysis 'yes'/'no', min number of ions and max number of ions
    # data_folder_path = '%s/%s' % (upper_folder_path, mode)  # Macintosh /Volumes/dfs/DATA/2017/2017_Feb-PI-ICR/2017-02-cross-checks/P1-P2/  # Windows C://Users/jkarthei/cernbox/Python/test2      #  G://Experiments/ISOLTRAP/DATA/2017/2017_Feb-PI-ICR/2017-02-long-run
    data_folder_path = upper_folder_path
    path, folder_name = os.path.split(data_folder_path)  # saves the folder name
    os.chdir(data_folder_path)


    all_peak_positions_c = []
    all_peak_positions_c = piicr_daq(data_folder_path, axis_scaling_factor, dont_show, pattern, sigma, bins, fit_range, cut_counter, tof_cut_on, tof_min, tof_max, trigger_on_off, clear_one, peak_finder_plot, all_peak_positions, hist2d_spines_off, hist2d_min_max, simple_plot_on, surface3D_plot_on, trigger_splitter, trigger_numbers, file_or_folder, file_name_input, color_map, analysis, z_class_analysis)
    # number_of_patterns = 1
    # all_positions_andree = [['File name', 'X-Pos', 'X-Pos-error', 'Y-Pos', 'Y-Pos-error', 'X-FWHM', 'X-FWHM-error', 'Y-FWHM', 'Y-FWHM-error', 'Counts']]
    # for i in range(1, number_of_patterns+1):
    #     all_peak_positions_c = []
    #     all_peak_positions_c = piicr_daq(data_folder_path, axis_scaling_factor, dont_show, [i, number_of_patterns], sigma, bins, fit_range, cut_counter, tof_cut_on, tof_min, tof_max, trigger_on_off, clear_one, peak_finder_plot, all_peak_positions, hist2d_spines_off, hist2d_min_max, simple_plot_on, surface3D_plot_on, trigger_splitter, trigger_numbers, file_or_folder, file_name_input, color_map, analysis)
    #     all_positions_andree.append(all_peak_positions_c[-1])

    # write_csv(upper_folder_path, all_positions_andree, 'all_positions_andree_slow')

