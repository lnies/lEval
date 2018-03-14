
# ---------------------------------------------------------------------------
# Written by Jonas Karthein in 2016/2017. Questions to jonas.karthein@cern.ch
# ---------------------------------------------------------------------------

import copy as CpY
import csv
import glob
import sys, os
from piicr_daq_functions import *
from read_write_functions import *
import numpy
from ellipse_selector import *
import time
import json


def piicr_daq(data_folder_path, axis_scaling_factor, dont_show, pattern, sigma, bins, fit_range, cut_counter, tof_cut_on, tof_min, tof_max, trigger_on_off, clear_one, peak_finder_plot, all_peak_positions, hist2d_spines_off, hist2d_min_max, simple_plot_on, surface3D_plot_on, trigger_splitter, trigger_numbers, file_or_folder, file_name_input, color_map, analysis_old_new, z_class_analysis, FWHM_parameter_fix, one_or_two_peaks_to_fit):
    axis_scaling_factor = 0.031746032       # convert a.u. to mm
    origin_cartesian = [0, -31.5]       # in case the origin for a ring shape analysis is not at 0/0 change it here
    single_analysis = True
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
            print('\n### Loading\n\nFile name                ::  %s') % (file_name)
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
                if single_analysis == True:
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

                window_x = CpY.copy(window_x)
                window_y = CpY.copy(window_y)
                window_xy = CpY.copy(window_xy)
                err_window_x = round(CpY.copy(err_window_x))
                err_window_y = round(CpY.copy(err_window_y))
                err_window_xy = round(CpY.copy(err_window_xy))
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
                        mean_x = numpy.mean(hilf_xx)
                        mean_y = numpy.mean(hilf_yy)
                    else:
                        for i in range(len(spot_positions)):
                            if spot_positions[i][1] > fit_range[0] and spot_positions[i][1] < fit_range[1]:
                                hilf_xx.append(spot_positions[i][1])
                            if spot_positions[i][2] > fit_range[2] and spot_positions[i][2] < fit_range[3]:
                                hilf_yy.append(spot_positions[i][2])
                        mean_x = numpy.mean(hilf_xx)
                        mean_y = numpy.mean(hilf_yy)


                elif one_or_two_peaks_to_fit == 'two':
                    for spot_pos in [spot_positions_dominant, spot_positions_recessive]:
                        hilf_xx = []
                        hilf_yy = []
                        for i in range(len(spot_pos)):
                            if np.sqrt(spot_pos[i][1] ** 2 + spot_pos[i][2] ** 2) < 650:
                                hilf_xx.append(spot_pos[i][1])
                                hilf_yy.append(spot_pos[i][2])
                        if spot_pos == spot_positions_dominant:
                            mean_x_dom = numpy.mean(hilf_xx)
                            mean_y_dom = numpy.mean(hilf_yy)
                        else:
                            mean_x_rec = numpy.mean(hilf_xx)
                            mean_y_rec = numpy.mean(hilf_yy)

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

                        if single_analysis == True:
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

                        window_x = CpY.copy(window_x)
                        window_y = CpY.copy(window_y)
                        window_xy = CpY.copy(window_xy)
                        err_window_x = round(CpY.copy(err_window_x))
                        err_window_y = round(CpY.copy(err_window_y))
                        err_window_xy = round(CpY.copy(err_window_xy))
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
                                mean_x = numpy.mean(hilf_xx)
                                mean_y = numpy.mean(hilf_yy)
                            else:
                                for i in range(len(spot_positions)):
                                    if spot_positions[i][1] > fit_range[0] and spot_positions[i][1] < fit_range[1]:
                                        hilf_xx.append(spot_positions[i][1])
                                    if spot_positions[i][2] > fit_range[2] and spot_positions[i][2] < fit_range[3]:
                                        hilf_yy.append(spot_positions[i][2])
                                mean_x = numpy.mean(hilf_xx)
                                mean_y = numpy.mean(hilf_yy)


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

