# ---------------------------------------------------------------------------
# Written by Jonas Karthein in 2016/2017. Questions to jonas.karthein@cern.ch
# ---------------------------------------------------------------------------

import os
import sys
import math
import numpy as np
from scipy.stats import chisquare
import pandas as pd
import datetime
import xlsxwriter
import xlrd
from get_ame_data import Ame
from read_write_functions import *
from python_plotter_functions import *
import platform


def GetNextHigh(temp, templist):
    """Outputs the next higher value in the list."""
    templist = (int(t) for t in templist if t != '')
    templist = [t for t in templist if t > int(temp)]
    if templist: return min(templist)
    else: return None

def GetNextLow(temp, templist):
    """Outputs the next lower value in the list."""
    templist = (int(t) for t in templist if t != '')
    templist = [t for t in templist if t < int(temp)]
    if templist: return max(templist)
    else: return None


def folder_creation_cross_checks(folder_name, ref, meas):

    if platform.system() == 'Windows':
        if os.path.isdir(folder_name+'\\cross-checks') == False:
            os.mkdir(folder_name+'\\cross-checks')
        if os.path.isdir(folder_name+'\\cross-checks\\'+ref+'-'+meas+'-'+ref) == False:
            os.mkdir(folder_name+'\\cross-checks\\'+ref+'-'+meas+'-'+ref)

    elif platform.system() == 'Darwin': # MacOS
        if os.path.isdir(folder_name+'/cross-checks/') == False:
            os.mkdir(folder_name+'/cross-checks/')
        if os.path.isdir(folder_name+'/cross-checks/'+ref+'-'+meas+'-'+ref+'/') == False:
            os.mkdir(folder_name+'/cross-checks/'+ref+'-'+meas+'-'+ref+'/')

def cross_checks(upper_folder_path, reference_isotope, measurement_isotope, isomer_state__empty_str_if_not_isotope_analysis, mass_excess_meas_import, mass_excess_meas_err_import):
    '''The function performs cross checks.'''

    # ---------------------------------------------------------------------------
    # ADVANCED SETTINGS -- Only change if you know what you do!!!
    # ---------------------------------------------------------------------------
    if platform.system() == 'Windows':
        upper_folder_path_reference_ion = upper_folder_path+'\\'+reference_isotope
        upper_folder_path_measurement = upper_folder_path+'\\'+measurement_isotope
        upper_folder_path_cross_checks = upper_folder_path+'\\cross-checks\\'+reference_isotope+'-'+measurement_isotope+'-'+reference_isotope

    elif platform.system() == 'Darwin': # MacOS
        upper_folder_path_reference_ion = upper_folder_path+'/'+reference_isotope
        upper_folder_path_measurement = upper_folder_path+'/'+measurement_isotope
        upper_folder_path_cross_checks = upper_folder_path+'/cross-checks/'+reference_isotope+'-'+measurement_isotope+'-'+reference_isotope

    analysis = 'new'
    modes = ['c', 'p1', 'p2']
    atomic_mass_unit_in_keV = 931494.0023
    atomic_mass_unit_in_keV_err = 0.0007
    electron_mass_in_unit = 5.4857990946E-04
    electron_mass_in_unit_err = 2.2E-13
    magnetic_field_disp = 6.35E-11      # unit = [min^{-1}]
    mass_dependent_shift_correction = 0    #1.6E-10
    systematic_error = 0 #8E-9
    isomer = isomer_state__empty_str_if_not_isotope_analysis      # '_ground' or '_isomer' or ''


    if platform.system() == 'Windows':
        if os.path.isfile(upper_folder_path_measurement+'\\freq_config/piicr_excitation.csv'):
            piicr_excitation_import = read_csv(upper_folder_path_measurement+'\\freq_config', 'piicr_excitation')

    elif platform.system() == 'Darwin': # MacOS
        if os.path.isfile(upper_folder_path_measurement+'/freq_config/piicr_excitation.csv'):
            piicr_excitation_import = read_csv(upper_folder_path_measurement+'/freq_config/', 'piicr_excitation')
    phase_accu_time = float(piicr_excitation_import[1][3])*0.001    #t_acc in ms

    # ---------------------------------------------------------------------------
    # get mass excess and atomic mass from AME
    # ---------------------------------------------------------------------------
    # data_folder_path = '%s/%s' % (upper_folder_path, mode)  # Macintosh /Volumes/dfs/DATA/2017/2017_Feb-PI-ICR/2017-02-cross-checks/P1-P2/  # Windows C://Users/jkarthei/cernbox/Python/test2      #  G://Experiments/ISOLTRAP/DATA/2017/2017_Feb-PI-ICR/2017-02-long-run
    path, folder_name_ref = os.path.split(upper_folder_path_reference_ion)  # saves the folder name
    path, folder_name_meas = os.path.split(upper_folder_path_measurement)  # saves the folder name
    path, folder_name = os.path.split(upper_folder_path_cross_checks)  # saves the folder name
    latex_name_ref = '$^{%s}$%s' % (filter(str.isdigit, folder_name_ref), filter(str.isalpha, folder_name_ref))
    latex_name_meas = '$^{%s}$%s' % (filter(str.isdigit, folder_name_meas), filter(str.isalpha, folder_name_meas))
    # latex_name_meas = '$^{%s}$%s' % ('129', 'Cd')


    ame_table = Ame()
    symbol_ref, atom_mass_ref, atom_mass_unc_ref = ame_table.get_ame_mass('%s' % (folder_name_ref))
    if isomer == '_ground' or isomer == '_isomer':
        symbol_meas, atom_mass_meas, atom_mass_unc_meas = ame_table.get_ame_mass('%s' % (folder_name_meas.replace('m', '').replace('g', '')))
    elif isomer == '':
        symbol_meas, atom_mass_meas, atom_mass_unc_meas = ame_table.get_ame_mass('%s' % (folder_name_meas))
    # symbol_meas, atom_mass_meas, atom_mass_unc_meas = ame_table.get_ame_mass('%s' % ('129Cd'))
    a_ref = int(filter(str.isdigit, symbol_ref))
    a_meas = int(filter(str.isdigit, symbol_meas))
    # a_meas = 129
    mass_excess_ref = (atom_mass_ref * 0.000001 - a_ref) * atomic_mass_unit_in_keV
    mass_excess_ref_err = np.sqrt((atom_mass_unc_ref * 0.000001 * atomic_mass_unit_in_keV) ** 2 + ((atom_mass_ref * 0.000001 - a_ref) * atomic_mass_unit_in_keV_err) ** 2)
    mass_excess_meas = (atom_mass_meas * 0.000001 - a_meas) * atomic_mass_unit_in_keV
    # if isomer == '_isomer' and symbol_meas == '123Cd':
    #     mass_excess_meas = -77271.2
    #     mass_excess_meas_err = 4.8
    mass_excess_meas_err = np.sqrt((atom_mass_unc_meas * 0.000001 * atomic_mass_unit_in_keV) ** 2 + ((atom_mass_meas * 0.000001 - a_meas) * atomic_mass_unit_in_keV_err) ** 2)

    if mass_excess_meas_import != 0:
        mass_excess_meas = mass_excess_meas_import
        mass_excess_meas_err = mass_excess_meas_err_import

    # mass_excess_meas = -68743.4     # titan 127gCd
    # mass_excess_meas_err = 56       # titan
    # mass_excess_meas = -68460.1     # titan 127mCd
    # mass_excess_meas_err = 47       # titan

    print '\n### Cross checks\n'
    print 'Reference isotope        :: {}: m = ({:.3f} +/- {:.3f}) uu --> ME = ({:.3f} +/- {:.3f}) keV'.format(symbol_ref, atom_mass_ref, atom_mass_unc_ref, mass_excess_ref, mass_excess_ref_err)
    print 'Measurement isotope      :: {}: m = ({:.3f} +/- {:.3f}) uu --> ME = ({:.3f} +/- {:.3f}) keV'.format(symbol_meas, atom_mass_meas, atom_mass_unc_meas, mass_excess_meas, mass_excess_meas_err)
    # ---------------------------------------------------------------------------
    # import and sort frequency data
    # ---------------------------------------------------------------------------

    ref_freq_import = read_csv(upper_folder_path_reference_ion, 'piicr_cyc_freq_%s' % (folder_name_ref))
    meas_freq_import = read_csv(upper_folder_path_measurement, 'piicr_cyc_freq_%s%s' % (folder_name_meas, isomer))

    header_ref = ref_freq_import.pop(0)
    ref_freq_import.sort(key=lambda x: x[0])
    ref_freq_import.insert(0, header_ref)

    header_meas = meas_freq_import.pop(0)
    meas_freq_import.sort(key=lambda x: x[0])
    meas_freq_import.insert(0, header_meas)

    # ---------------------------------------------------------------------------
    # get calculated frequencies for ref and meas ion
    # ---------------------------------------------------------------------------

    meas_freq_import.pop(0)     # delete header
    ref_freq_import.pop(0)

    hilf_meas_nr = []
    for i in range(0, len(meas_freq_import), 1):       # get the file number
        hilf_meas_nr.append(int(meas_freq_import[i][0][-3:]))
    hilf_ref_nr = []
    for i in range(0, len(ref_freq_import), 1):
        hilf_ref_nr.append(int(ref_freq_import[i][0][-3:]))


    hilf_meas = []
    hilf_ref = []
    for meas_nr in hilf_meas_nr:        # sort and combine always REF-MEAS-REF and consider that always a REF must be at the beginning and end!
        if GetNextLow(meas_nr, hilf_ref_nr) != None and GetNextHigh(meas_nr, hilf_ref_nr) != None:

            # print 'Next lower to %s: '%(meas_nr), GetNextLow(meas_nr, hilf_ref_nr)
            # print 'Next higher to %s: ' %(meas_nr), GetNextHigh(meas_nr, hilf_ref_nr)

            for i in range(len(meas_freq_import)):       # create list for all needed MEAS
                if meas_nr == int(meas_freq_import[i][0][-3:]):
                    hilf_meas.append([meas_freq_import[i][0], float(meas_freq_import[i][1]), float(meas_freq_import[i][2]), int(meas_freq_import[i][0][-3:])])
            for i in range(len(ref_freq_import)):         # create list for all needed REF
                if GetNextLow(meas_nr, hilf_ref_nr) == int(ref_freq_import[i][0][-3:]):
                    hilf_ref.append([ref_freq_import[i][0], float(ref_freq_import[i][1]), float(ref_freq_import[i][2]), GetNextLow(meas_nr, hilf_ref_nr)])
                if GetNextHigh(meas_nr, hilf_ref_nr) == int(ref_freq_import[i][0][-3:]):
                    hilf_ref.append([ref_freq_import[i][0], float(ref_freq_import[i][1]), float(ref_freq_import[i][2]), GetNextHigh(meas_nr, hilf_ref_nr)])

    # ---------------------------------------------------------------------------
    # get time information for ref and meas ion
    # ---------------------------------------------------------------------------

    if analysis == 'new':
        time_info_p2_ref_import = read_csv('%s/p1p2' % (upper_folder_path_reference_ion), '_time_info_p1p2')
        time_info_p2_meas_import = read_csv('%s/p1p2' % (upper_folder_path_measurement), '_time_info_p1p2')
    else:
        time_info_p2_ref_import = read_csv('%s/p2' % (upper_folder_path_reference_ion), '_time_info_p2')
        time_info_p2_meas_import = read_csv('%s/p2' % (upper_folder_path_measurement), '_time_info_p2')
    time_info_p2_ref_begin = []
    time_info_p2_ref_end = []
    time_info_p2_meas_begin = []
    time_info_p2_meas_end = []
    for number in range(len(hilf_ref)):
        for row in range(len(time_info_p2_ref_import)):
            if int(time_info_p2_ref_import[row][0][-3:]) == int(hilf_ref[number][3]):
                time_info_p2_ref_begin_hilf = []
                time_info_p2_ref_end_hilf = []
                time_info_p2_ref_begin_hilf = [time_info_p2_ref_import[row][0], datetime.datetime.strptime(time_info_p2_ref_import[row][1]+' '+time_info_p2_ref_import[row][2], '%m/%d/%Y %H:%M:%S')]
                time_info_p2_ref_end_hilf = [time_info_p2_ref_import[row][0], datetime.datetime.strptime(time_info_p2_ref_import[row][3]+' '+time_info_p2_ref_import[row][4], '%m/%d/%Y %H:%M:%S')]
                time_info_p2_ref_begin.append(time_info_p2_ref_begin_hilf)
                time_info_p2_ref_end.append(time_info_p2_ref_end_hilf)
    for number1 in range(len(hilf_meas)):
        for row1 in range(len(time_info_p2_meas_import)):
            if int(time_info_p2_meas_import[row1][0][-3:]) == int(hilf_meas[number1][3]):
                time_info_p2_meas_begin_hilf = []
                time_info_p2_meas_end_hilf = []
                time_info_p2_meas_begin_hilf = [time_info_p2_meas_import[row1][0], datetime.datetime.strptime(time_info_p2_meas_import[row1][1]+' '+time_info_p2_meas_import[row1][2], '%m/%d/%Y %H:%M:%S')]
                time_info_p2_meas_end_hilf = [time_info_p2_meas_import[row1][0], datetime.datetime.strptime(time_info_p2_meas_import[row1][3]+' '+time_info_p2_meas_import[row1][4], '%m/%d/%Y %H:%M:%S')]
                time_info_p2_meas_begin.append(time_info_p2_meas_begin_hilf)
                time_info_p2_meas_end.append(time_info_p2_meas_end_hilf)

    time_info_ref = []
    for i in range(len(time_info_p2_ref_begin)):
        time_info_ref.append([symbol_ref+'-'+time_info_p2_ref_begin[i][0][-3:], (time_info_p2_ref_begin[i][1] + (time_info_p2_ref_end[i][1]-time_info_p2_ref_begin[i][1])/2)])
    time_info_meas = []
    for i in range(len(time_info_p2_meas_begin)):
        time_info_meas.append([symbol_ref+'-'+time_info_p2_meas_begin[i][0][-3:], (time_info_p2_meas_begin[i][1] + (time_info_p2_meas_end[i][1]-time_info_p2_meas_begin[i][1])/2)])

    # ---------------------------------------------------------------------------
    # initialize final vector
    # ---------------------------------------------------------------------------

    cross_checks = [['File name', '#', 'Ref/Meas', 'Time / y-m-d h:m:s', 'PI-ICR cyc. freq. / Hz',
                     'Cyc. freq. err. / Hz', 'Cyc. freq. ref. interpolated / Hz',
                     'Cyc. freq. ref. interpolated err. / Hz', 'Ratio r', 'r err.',
                     'Corr. r err.', '1/(r err.)^2', 'r/(r err.)^2', '[(mean-r)/(r err.)]^2',
                     'Calc. mass / u', 'Calc. mass err. / u', 'Mass excess / keV',
                     'Mass excess err. / keV', 'Diff. to AME in mass excess / keV',
                     'Diff. within 1 sigma', 'Diff. within 2 sigma',
                     'Corr. diff. to AME in mass excess / keV',
                     'Corr. diff. to AME within 1 sigma',
                     'Corr. diff. to AME within 2 sigma']]

    for i in range(len(hilf_meas)+len(hilf_ref)):       # creating an array of zeros
        cross_checks.append([0, i+1, 0, '0', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    ref_counter = 0
    meas_counter = 0
    for i in range(1, len(cross_checks), 3):
        # REF before MEAS
        cross_checks[i][0] = symbol_ref+'-'+hilf_ref[ref_counter][0][-3:]                     # file name     R1
        cross_checks[i][2] = 'REF'                                                            # ref/meas
        cross_checks[i][3] = time_info_ref[ref_counter][1].strftime("%Y-%m-%d %H:%M:%S")      # time          R1
        cross_checks[i][4] = hilf_ref[ref_counter][1]                                         # cyc freq      R1
        cross_checks[i][5] = hilf_ref[ref_counter][2]                                         # err cyc freq  R1
        ref_counter += 1

        # MEAS
        cross_checks[i+1][0] = symbol_meas+'-'+hilf_meas[meas_counter][0][-3:]                # file name     M
        cross_checks[i+1][2] = 'MEAS'                                                         # ref/meas
        cross_checks[i+1][3] = time_info_meas[meas_counter][1].strftime("%Y-%m-%d %H:%M:%S")  # time          M
        cross_checks[i+1][4] = hilf_meas[meas_counter][1]                                     # cyc freq      M
        cross_checks[i+1][5] = hilf_meas[meas_counter][2]                                     # err cyc freq  M
        meas_counter += 1

        # REF after MEAS
        cross_checks[i+2][0] = symbol_ref+'-'+hilf_ref[ref_counter][0][-3:]                   # file name     R2
        cross_checks[i+2][2] = 'REF'                                                          # ref/meas
        cross_checks[i+2][3] = time_info_ref[ref_counter][1].strftime("%Y-%m-%d %H:%M:%S")    # time          R2
        cross_checks[i+2][4] = hilf_ref[ref_counter][1]                                       # cyc freq      R2
        cross_checks[i+2][5] = hilf_ref[ref_counter][2]                                       # err cyc freq  R2
        ref_counter += 1

    # ---------------------------------------------------------------------------
    # calculation of interpolated freq
    # ---------------------------------------------------------------------------

    for i in range(1, len(cross_checks), 3):
        # REF before MEAS
        cross_checks[i][6] = ''
        cross_checks[i][7] = ''

        # MEAS
        cross_checks[i+1][6] = ((cross_checks[i+2][4] - cross_checks[i][4])
                                * ((datetime.datetime.strptime(cross_checks[i+1][3], "%Y-%m-%d %H:%M:%S")
                                    - datetime.datetime.strptime(cross_checks[i][3], "%Y-%m-%d %H:%M:%S")).total_seconds())
                                / ((datetime.datetime.strptime(cross_checks[i+2][3], "%Y-%m-%d %H:%M:%S")
                                    - datetime.datetime.strptime(cross_checks[i][3], "%Y-%m-%d %H:%M:%S")).total_seconds())
                                + cross_checks[i][4])
        cross_checks[i+1][7] = (np.sqrt(((cross_checks[i][5]) ** 2)
                                        * (((datetime.datetime.strptime(cross_checks[i+2][3], "%Y-%m-%d %H:%M:%S")
                                             - datetime.datetime.strptime(cross_checks[i+1][3], "%Y-%m-%d %H:%M:%S")).total_seconds()) ** 2)
                                        + (((cross_checks[i+2][5]) ** 2)
                                        * (((datetime.datetime.strptime(cross_checks[i+1][3], "%Y-%m-%d %H:%M:%S")
                                             - datetime.datetime.strptime(cross_checks[i][3], "%Y-%m-%d %H:%M:%S")).total_seconds()) ** 2)))
                                / ((datetime.datetime.strptime(cross_checks[i+2][3], "%Y-%m-%d %H:%M:%S")
                                    - datetime.datetime.strptime(cross_checks[i][3], "%Y-%m-%d %H:%M:%S")).total_seconds()))

        # REF after MEAS
        cross_checks[i+2][6] = ''
        cross_checks[i+2][7] = ''

    # ---------------------------------------------------------------------------
    # calculation of ratio r
    # ---------------------------------------------------------------------------

    for i in range(1, len(cross_checks), 3):
        # REF before MEAS
        cross_checks[i][8] = ''
        cross_checks[i][9] = ''
        cross_checks[i][10] = ''
        cross_checks[i][11] = ''
        cross_checks[i][12] = ''

        # MEAS
        cross_checks[i+1][8] = (cross_checks[i+1][6] / cross_checks[i+1][4])
        cross_checks[i+1][9] = (cross_checks[i+1][8]
                                * np.sqrt(((cross_checks[i+1][7] / cross_checks[i+1][6]) ** 2)
                                          + ((cross_checks[i+1][5] / cross_checks[i+1][4]) ** 2)))
        cross_checks[i+1][10] = np.sqrt((cross_checks[i+1][9] ** 2)
                                        + (((datetime.datetime.strptime(cross_checks[i+2][3], "%Y-%m-%d %H:%M:%S")
                                             - datetime.datetime.strptime(cross_checks[i][3], "%Y-%m-%d %H:%M:%S")).total_seconds())
                                        / 60 * magnetic_field_disp
                                        * cross_checks[i+1][8]) ** 2)
        cross_checks[i+1][11] = 1 / (cross_checks[i+1][10] ** 2)
        cross_checks[i+1][12] = cross_checks[i+1][8] * cross_checks[i+1][11]

        # REF after MEAS
        cross_checks[i+2][8] = ''
        cross_checks[i+2][9] = ''
        cross_checks[i+2][10] = ''
        cross_checks[i+2][11] = ''
        cross_checks[i+2][12] = ''

    mean_1_sigma = 0
    mean_r_sigma = 0
    for i in range(1, len(cross_checks), 3):
        mean_1_sigma += cross_checks[i+1][11]
        mean_r_sigma += cross_checks[i+1][12]
    mean_ratio_r = mean_r_sigma / mean_1_sigma      # weighted
    print 'Mean ratio R             :: {:.3f} +/- {:.3f}'.format(mean_ratio_r, mean_r_sigma)
    for i in range(1, len(cross_checks), 3):
        # REF before MEAS
        cross_checks[i][13] = ''

        # MEAS
        cross_checks[i+1][13] = ((mean_ratio_r - cross_checks[i+1][8]) / cross_checks[i+1][10]) ** 2

        # REF after MEAS
        cross_checks[i+2][13] = ''

    # ---------------------------------------------------------------------------
    # calculation of the mass
    # ---------------------------------------------------------------------------
    for i in range(1, len(cross_checks), 3):
        # REF before MEAS
        cross_checks[i][14] = ''
        cross_checks[i][15] = ''

        # MEAS
        cross_checks[i+1][14] = cross_checks[i+1][8] * (atom_mass_ref * 0.000001 - electron_mass_in_unit) + electron_mass_in_unit
        cross_checks[i+1][15] = np.sqrt(((atom_mass_ref * 0.000001 - electron_mass_in_unit) ** 2) * (cross_checks[i+1][9] ** 2)     # shouldn't it be here cross_checks[i+1][10] instead of cross_checks[i+1][9] ???!?!?!
                                        + ((cross_checks[i+1][8] ** 2) * ((atom_mass_unc_ref * 0.000001) ** 2)))  # uncertainty of electron mass not included

        # REF after MEAS
        cross_checks[i+2][14] = ''
        cross_checks[i+2][15] = ''

    # ---------------------------------------------------------------------------
    # calculation of the mass excess
    # ---------------------------------------------------------------------------

    for i in range(1, len(cross_checks), 3):
        # REF before MEAS
        cross_checks[i][16] = ''
        cross_checks[i][17] = ''
        cross_checks[i][18] = ''
        cross_checks[i][19] = ''
        cross_checks[i][20] = ''

        # MEAS
        cross_checks[i+1][16] = ((cross_checks[i+1][14] - a_meas) * atomic_mass_unit_in_keV)
        cross_checks[i+1][17] = (cross_checks[i+1][15] * atomic_mass_unit_in_keV)
        cross_checks[i+1][18] = mass_excess_meas - cross_checks[i+1][16] # difference to AME
        if abs(cross_checks[i+1][18]) > (cross_checks[i+1][17] + mass_excess_meas_err):
            cross_checks[i+1][19] = True
        else:
            cross_checks[i+1][19] = False
        if abs(cross_checks[i+1][18]) > (2*(cross_checks[i+1][17] + mass_excess_meas_err)):
            cross_checks[i+1][20] = True
        else:
            cross_checks[i+1][20] = False

        # REF after MEAS
        cross_checks[i+2][16] = ''
        cross_checks[i+2][17] = ''
        cross_checks[i+2][18] = ''
        cross_checks[i+2][19] = ''
        cross_checks[i+2][20] = ''

    # ---------------------------------------------------------------------------
    # correction for mass dependent shift
    # ---------------------------------------------------------------------------


    mean_ratio_r_corrected = mean_ratio_r + mass_dependent_shift_correction * mean_ratio_r * (atom_mass_meas * 0.000001 - atom_mass_ref * 0.000001)
    mass_isoltrap_final = mean_ratio_r_corrected * (atom_mass_ref * 0.000001 - electron_mass_in_unit) + electron_mass_in_unit
    mass_isoltrap_uncorrected = mean_ratio_r * (atom_mass_ref * 0.000001 - electron_mass_in_unit) + electron_mass_in_unit
    mass_excess_meas_isoltrap_final = (mass_isoltrap_final - a_meas) * atomic_mass_unit_in_keV
    mass_excess_meas_isoltrap_uncorrected = (mass_isoltrap_uncorrected - a_meas) * atomic_mass_unit_in_keV

    print 'Final ISOLTRAP mass      :: {:.3f}'.format(mass_isoltrap_final*1000000)
    print 'Final ISOLTRAP ME        :: {:.3f} keV'.format(mass_excess_meas_isoltrap_final)



    for i in range(1, len(cross_checks), 3):
        # REF before MEAS
        cross_checks[i][21] = ''
        cross_checks[i][22] = ''
        cross_checks[i][23] = ''

        # MEAS
        cross_checks[i+1][21] = cross_checks[i+1][18] - (mass_excess_meas_isoltrap_final - mass_excess_meas_isoltrap_uncorrected) # difference to AME corrected
        if abs(cross_checks[i+1][21]) > (cross_checks[i+1][17] + mass_excess_meas_err):
            cross_checks[i+1][22] = False
        else:
            cross_checks[i+1][22] = True
        if abs(cross_checks[i+1][21]) > (2*(cross_checks[i+1][17] + mass_excess_meas_err)):
            cross_checks[i+1][23] = False
        else:
            cross_checks[i+1][23] = True

        # REF after MEAS
        cross_checks[i+2][21] = ''
        cross_checks[i+2][22] = ''
        cross_checks[i+2][23] = ''

    # ---------------------------------------------------------------------------
    # uncertainties
    # ---------------------------------------------------------------------------
    avg_red_chi_sq = 0
    counter_chi = 0
    for i in range(1, len(cross_checks), 3):
        counter_chi += 1
        avg_red_chi_sq += cross_checks[i+1][13]
    avg_red_chi_sq = avg_red_chi_sq / (counter_chi - 1)

    abs_statistical_unc = 0
    for i in range(1, len(cross_checks), 3):
        abs_statistical_unc += cross_checks[i+1][11]
    if avg_red_chi_sq < 1:
        abs_statistical_unc = np.sqrt(1 / abs_statistical_unc)
    else:
        abs_statistical_unc = np.sqrt(1 / abs_statistical_unc) * np.sqrt(avg_red_chi_sq)

    abs_mass_dep_shift = abs(mass_dependent_shift_correction * (atom_mass_ref * 0.000001 - atom_mass_meas * 0.000001)*mean_ratio_r)
    abs_resid_sys_shift = mean_ratio_r * systematic_error

    sum__abs_stat__mass_dep__syst_uncertainty = np.sqrt(abs_statistical_unc ** 2 + abs_mass_dep_shift ** 2 + abs_resid_sys_shift ** 2)
    rel_stat__mass_dep__syst_uncertainty = sum__abs_stat__mass_dep__syst_uncertainty / mean_ratio_r

    mass_isoltrap_final_err = np.sqrt(((sum__abs_stat__mass_dep__syst_uncertainty * (atom_mass_ref * 0.000001-electron_mass_in_unit)) ** 2) + ((atom_mass_unc_ref * 0.000001 * mean_ratio_r_corrected) ** 2)) # uncertainty of electron mass not included
    mass_excess_meas_isoltrap_final_err = mass_isoltrap_final_err * atomic_mass_unit_in_keV

    diff_to_AME = (mass_excess_meas - mass_excess_meas_isoltrap_final)
    diff_to_AME_err = (mass_excess_meas_isoltrap_final_err + mass_excess_ref_err)

    print 'Avg. red. Chi. sq.       :: {:.3f}'.format(avg_red_chi_sq)
    print 'Difference to AME        :: ({:.3f} +/- {:.3f}) keV'.format(diff_to_AME, diff_to_AME_err)

    # ---------------------------------------------------------------------------
    # plot
    # ---------------------------------------------------------------------------

    os.chdir(upper_folder_path_cross_checks)

    difference_to_AME_final = [['#', 'Difference to AME / keV', 'Err / keV']]
    difference_to_AME_final_time = [['time', 'Difference to AME / keV', 'Err / keV']]
    counter = 0
    for i in range(1, len(cross_checks), 3):
        counter += 1
        difference_to_AME_final.append([counter, cross_checks[i+1][21], cross_checks[i+1][17]])
        difference_to_AME_final_time.append([datetime.datetime.strptime(cross_checks[i+1][3], "%Y-%m-%d %H:%M:%S"), cross_checks[i+1][21], cross_checks[i+1][17]])
    python_plot(difference_to_AME_final, 'Cross checks %s: %s-%s-%s' % (isomer, latex_name_ref, latex_name_meas, latex_name_ref), 'cross_checks_%s%s' % (folder_name, isomer), 'no', '', 'no', 'linear', 'Utopia', 'red', 'green', 'full', 1, 5, 1, 'scatter', 20, 5, 'on', 'on', [avg_red_chi_sq, phase_accu_time, mass_dependent_shift_correction, systematic_error], diff_to_AME, diff_to_AME_err)
    python_plot(difference_to_AME_final_time, 'Cross checks %s: %s-%s-%s' % (isomer, latex_name_ref, latex_name_meas, latex_name_ref), 'cross_checks_%s%s_time' % (folder_name, isomer), 'no', '', 'no', 'linear', 'Utopia', 'red', 'green', 'full', 1, 5, 1, 'scatter', 20, 5, 'on', 'on', [avg_red_chi_sq, phase_accu_time, mass_dependent_shift_correction, systematic_error], diff_to_AME, diff_to_AME_err)
    # ---------------------------------------------------------------------------
    # save
    # ---------------------------------------------------------------------------
    cross_checks.insert(0, ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''])
    cross_checks.insert(0, ['Avg. red. chi sq.', '%9.4f' % avg_red_chi_sq, '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''])
    cross_checks.insert(0, ['AME-ISOLTRAP', '%9.4f (%2.4f) keV' % (diff_to_AME, diff_to_AME_err), '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''])
    cross_checks.insert(0, ['ISOLTRAP final mass excess', '%9.4f (%2.4f) keV' % (mass_excess_meas_isoltrap_final, mass_excess_meas_isoltrap_final_err), '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''])
    cross_checks.insert(0, ['ISOLTRAP final mass', '%9.7f (%2.7f) unit' % (mass_isoltrap_final, mass_isoltrap_final_err), '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''])

    write_csv(upper_folder_path_cross_checks, cross_checks, 'cross_checks_%s%s' % (folder_name, isomer))
    write_excel(upper_folder_path_cross_checks, cross_checks, 'cross_checks_%s%s' % (folder_name, isomer))

    # more info: http://xlsxwriter.readthedocs.io/format.html#format
    workbook = xlsxwriter.Workbook('cross_checks_%s%s_with_pic.xlsx' % (folder_name, isomer))
    worksheet = workbook.add_worksheet()
    workbook.set_properties({
        'title':    'Cross-checks-%s' % (folder_name),
        'subject':  'Penning trap mass spectrometry (PI-ICR)',
        'author':   'Jonas Karthein',
        'manager':  'Prof. Dr. Klaus Blaum',
        'company':  'MPIK / CERN',
        'category': 'cross-checks',
        'keywords': 'cross-checks, %s, PI-ICR' % (folder_name),
        'comments': 'Created with Python and XlsxWriter',
        'status':   'Quo',
    })
    worksheet.set_column('A:X', 15)
    worksheet.set_row(5, 20)
    red = workbook.add_format({'font_color': 'red'})
    bold2 = workbook.add_format({'bold': True})
    yellow = workbook.add_format({'bold': True, 'bg_color': 'yellow'})
    middle = workbook.add_format({'bold': True, 'font_size': 15, 'top': 2, 'bottom': 2})
    left = workbook.add_format({'bold': True, 'font_size': 15, 'top': 2, 'bottom': 2, 'left': 2})
    right = workbook.add_format({'bold': True, 'font_size': 15, 'top': 2, 'bottom': 2, 'right': 2})
    row = 0
    col = 0
    for a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x in (cross_checks):
        if row == 2:
            worksheet.write(row, col,     a, bold2)
            worksheet.write(row, col + 1, b, yellow)
            worksheet.write(row, col + 2, c)
            worksheet.write(row, col + 3, d)
            worksheet.write(row, col + 4, e)
            worksheet.write(row, col + 5, f)
            worksheet.write(row, col + 6, g)
            worksheet.write(row, col + 7, h)
            worksheet.write(row, col + 8, i)
            worksheet.write(row, col + 9, j)
            worksheet.write(row, col + 10, k)
            worksheet.write(row, col + 11, l)
            worksheet.write(row, col + 12, m)
            worksheet.write(row, col + 13, n)
            worksheet.write(row, col + 14, o)
            worksheet.write(row, col + 15, p)
            worksheet.write(row, col + 16, q)
            worksheet.write(row, col + 17, r)
            worksheet.write(row, col + 18, s)
            worksheet.write(row, col + 19, t)
            worksheet.write(row, col + 20, u)
            worksheet.write(row, col + 21, v)
            worksheet.write(row, col + 22, w)
            worksheet.write(row, col + 23, x)
        elif row == 5:
            worksheet.write(row, col,     a, left)
            worksheet.write(row, col + 1, b, middle)
            worksheet.write(row, col + 2, c, middle)
            worksheet.write(row, col + 3, d, middle)
            worksheet.write(row, col + 4, e, middle)
            worksheet.write(row, col + 5, f, middle)
            worksheet.write(row, col + 6, g, middle)
            worksheet.write(row, col + 7, h, middle)
            worksheet.write(row, col + 8, i, middle)
            worksheet.write(row, col + 9, j, middle)
            worksheet.write(row, col + 10, k, middle)
            worksheet.write(row, col + 11, l, middle)
            worksheet.write(row, col + 12, m, middle)
            worksheet.write(row, col + 13, n, middle)
            worksheet.write(row, col + 14, o, middle)
            worksheet.write(row, col + 15, p, middle)
            worksheet.write(row, col + 16, q, middle)
            worksheet.write(row, col + 17, r, middle)
            worksheet.write(row, col + 18, s, middle)
            worksheet.write(row, col + 19, t, middle)
            worksheet.write(row, col + 20, u, middle)
            worksheet.write(row, col + 21, v, middle)
            worksheet.write(row, col + 22, w, middle)
            worksheet.write(row, col + 23, x, right)
        elif c == 'MEAS':
            worksheet.write(row, col,     a, red)
            worksheet.write(row, col + 1, b, red)
            worksheet.write(row, col + 2, c, red)
            worksheet.write(row, col + 3, d, red)
            worksheet.write(row, col + 4, e, red)
            worksheet.write(row, col + 5, f, red)
            worksheet.write(row, col + 6, g, red)
            worksheet.write(row, col + 7, h, red)
            worksheet.write(row, col + 8, i, red)
            worksheet.write(row, col + 9, j, red)
            worksheet.write(row, col + 10, k, red)
            worksheet.write(row, col + 11, l, red)
            worksheet.write(row, col + 12, m, red)
            worksheet.write(row, col + 13, n, red)
            worksheet.write(row, col + 14, o, red)
            worksheet.write(row, col + 15, p, red)
            worksheet.write(row, col + 16, q, red)
            worksheet.write(row, col + 17, r, red)
            worksheet.write(row, col + 18, s, red)
            worksheet.write(row, col + 19, t, red)
            worksheet.write(row, col + 20, u, red)
            worksheet.write(row, col + 21, v, red)
            worksheet.write(row, col + 22, w, red)
            worksheet.write(row, col + 23, x, red)
        else:
            worksheet.write(row, col,     a)
            worksheet.write(row, col + 1, b)
            worksheet.write(row, col + 2, c)
            worksheet.write(row, col + 3, d)
            worksheet.write(row, col + 4, e)
            worksheet.write(row, col + 5, f)
            worksheet.write(row, col + 6, g)
            worksheet.write(row, col + 7, h)
            worksheet.write(row, col + 8, i)
            worksheet.write(row, col + 9, j)
            worksheet.write(row, col + 10, k)
            worksheet.write(row, col + 11, l)
            worksheet.write(row, col + 12, m)
            worksheet.write(row, col + 13, n)
            worksheet.write(row, col + 14, o)
            worksheet.write(row, col + 15, p)
            worksheet.write(row, col + 16, q)
            worksheet.write(row, col + 17, r)
            worksheet.write(row, col + 18, s)
            worksheet.write(row, col + 19, t)
            worksheet.write(row, col + 20, u)
            worksheet.write(row, col + 21, v)
            worksheet.write(row, col + 22, w)
            worksheet.write(row, col + 23, x)
        row += 1
    worksheet.insert_image('A%s' % (len(cross_checks)+2), 'cross_checks_%s%s.PNG' % (folder_name, isomer))
    workbook.close()
    # ---------------------------------------------------------------------------
    # convert cross_checks into pandas DataFrame (nice to read in terminal)
    # ---------------------------------------------------------------------------

    cross_checks.pop(0)
    cross_checks.pop(0)
    cross_checks.pop(0)
    pd.options.display.width = 180
    pd.set_option('display.max_column', 25)
    cross_checks_df = pd.DataFrame(cross_checks, columns=cross_checks.pop(0))
    # print cross_checks_df


def plot_all_cross_checks(upper_folder_paths):

    combined_cross_checks = [['run #', 'difference to AME / keV', 'Err / keV']]
    run_names = ''
    avg_red_chi_sq_all = []
    for upper_folder_path in upper_folder_paths:
        os.chdir(upper_folder_path)

        folder_name = os.path.dirname(os.path.abspath(__file__))
        current_run = folder_name.split('_run')[-1]

        if run_names == '':
            run_names += ('run'+current_run)
        else:
            run_names += '_' + 'run'+current_run

        if platform.system() == 'Windows':
            upper_folder_path_cross_checks = upper_folder_path+'\\cross-checks\\'
        elif platform.system() == 'Darwin': # MacOS
            upper_folder_path_cross_checks = upper_folder_path+'/cross-checks/'

        isotope_directory_list = []
        isotope_directory_list = [x[0] for x in os.walk(upper_folder_path_cross_checks)]

        for i in isotope_directory_list:
            os.chdir(i)
            csv_files = []
            for file in glob.glob("*.csv"):
                csv_files.append(file)
            for j in csv_files:
                input_list = []
                input_list = read_csv(i, os.path.splitext(j)[0])
                combined_temp = []
                avg_red_chi_sq_all.extend([float(input_list[3][1]) for x in range(((len(input_list)-6)/3))])
                for k in input_list:
                    if k[18] != '' and k[18] != 'Diff. to AME in mass excess / keV':
                        combined_temp.append([k[0]+'_'+current_run, float(k[18]), float(k[17])])
                combined_cross_checks.extend(combined_temp)

    # print combined_cross_checks
    weighted_avg = float(np.sum([x[1]/(x[2] ** 2) for x in combined_cross_checks[1:]])) / float(np.sum([1/(x[2] ** 2) for x in combined_cross_checks[1:]]))
    inner_unc_diff = float(np.sqrt(float(1/np.sum([1/(x[2] ** 2) for x in combined_cross_checks[1:]]))))
    outer_unc_diff = float(np.sqrt(np.sum([ (x[1] - weighted_avg) ** 2 / (x[2] ** 2) for x in combined_cross_checks[1:]]) / (np.sum([1/(x[2] ** 2) * (len(combined_cross_checks)-2) for x in combined_cross_checks[1:]]))))
    total_unc = max([inner_unc_diff, outer_unc_diff])
    birge_ratio = outer_unc_diff / inner_unc_diff

    # print inner_unc_diff, outer_unc_diff, birge_ratio
    os.chdir(upper_folder_paths[0])
    # print upper_folder_paths[0], run_names

    write_csv(upper_folder_paths[0], combined_cross_checks, 'combined_cross_checks_%s' % (run_names))
    python_plot(combined_cross_checks, 'Combined cross checks of %s' % (run_names), 'combined_cross_checks_%s' % (run_names), 'yes', '', 'no', 'linear', 'Utopia', 'red', 'green', 'full', 1, 5, 1, 'scatter', 20, 5, 'on', 'on', [birge_ratio, run_names, -1, -1], weighted_avg, total_unc)


def avg_freq(upper_folder_paths):
    avg_freq = [['Run # & isotope', 'Weighted avg. $\\nu_c$ / Hz', 'Err / Hz']]
    run_names = ''
    for upper_folder_path in upper_folder_paths:
        os.chdir(upper_folder_path)

        folder_name = os.path.dirname(os.path.abspath(__file__))
        current_run = folder_name.split('_run')[-1]
        if run_names == '':
            run_names += ('run'+current_run)
        else:
            run_names += '_' + ('run'+current_run)

        if platform.system() == 'Windows':
            upper_folder_path_cross_checks = upper_folder_path+'\\cross-checks\\'
        elif platform.system() == 'Darwin': # MacOS
            upper_folder_path_cross_checks = upper_folder_path+'/cross-checks/'

        isotope_directory_list = []
        isotope_directory_list = [x[0] for x in os.walk(upper_folder_path)]

        for i in isotope_directory_list:
            os.chdir(i)
            csv_files = []
            for file in glob.glob("piicr_cyc*.csv"):
                csv_files.append(file)
            for j in csv_files:
                input_list = []
                input_list = read_csv(i, os.path.splitext(j)[0])
                weighted_avg = float(np.sum([float(x[1])/(float(x[2]) ** 2) for x in input_list[1:]])) / float(np.sum([1/(float(x[2]) ** 2) for x in input_list[1:]]))
                inner_unc_diff = float(np.sqrt(float(1/np.sum([1/(float(x[2]) ** 2) for x in input_list[1:]]))))
                outer_unc_diff = float(np.sqrt(np.sum([ (float(x[1]) - weighted_avg) ** 2 / (float(x[2]) ** 2) for x in input_list[1:]]) / (np.sum([1/(float(x[2]) ** 2) * (len(input_list)-2) for x in input_list[1:]]))))
                total_unc = max([inner_unc_diff, outer_unc_diff])
                birge_ratio = outer_unc_diff / inner_unc_diff
                folder_name_i = os.path.dirname(os.path.abspath(__file__))
                current_isotope = folder_name_i.split('/')[-1]

                avg_freq.append([str(current_run+'_'+current_isotope), weighted_avg, total_unc])

    write_csv(upper_folder_paths[0], avg_freq, 'avg_freq_%s' % (run_names))
    python_plot(avg_freq, 'Average weighted $\\nu_c$ for %s' % (run_names), 'avg_freq_%s' % (run_names), 'yes', '', 'no', 'linear', 'Utopia', 'red', 'green', 'full', 1, 5, 1, 'scatter', 20, 5, 'on', 'on', '', weighted_avg, total_unc)


def weighted_mean(input_list):
    '''
    Takes a list with [[header, header_unc], [x, x_unc.],...] and calculates the weighted mean, inner and outer unc, total unc and the birge ratio.

    The ''internal'' uncertainty includes only measurement uncertainty, that is, uncertainty on N, in computing the uncertainty on the exposure age t. The internal uncertainty could be sometimes called the ''measurement uncertainty,'' although that is less clear and could mean several different things. The ''external'' uncertainty includes all input uncertainties - uncertainties on both N and P in the example above - and is also often called the ''total uncertainty.''

    taken from: https://cosmognosis.wordpress.com/2013/10/25/internal-vs-external-uncertainties-again/
    '''
    weighted_avg = float(np.sum([float(x[0])/(float(x[1]) ** 2) for x in input_list[1:]])) / float(np.sum([1/(float(x[1]) ** 2) for x in input_list[1:]]))
    internal_unc = float(np.sqrt(float(1/np.sum([1/(float(x[1]) ** 2) for x in input_list[1:]]))))
    external_unc = float(np.sqrt(np.sum([ (float(x[0]) - weighted_avg) ** 2 / (float(x[1]) ** 2) for x in input_list[1:]]) / (np.sum([1/(float(x[1]) ** 2) * (len(input_list)-2) for x in input_list[1:]]))))

    total_unc = max([internal_unc, external_unc])
    birge_ratio = external_unc / internal_unc

    return(weighted_avg, internal_unc, external_unc, total_unc, birge_ratio)


if __name__ == '__main__':
    # ---------------------------------------------------------------------------
    # Please change before use to the correct folder -- Instructions in REAM-ME-first
    # ---------------------------------------------------------------------------

    upper_folder_path = '/Volumes/dfs/DATA/2017/2017-06-Cd-run_prep/Cross-checks/test'
    reference_isotope = '85Rb'
    measurement_isotope = '87Rb'

    # cross_checks(upper_folder_path, reference_isotope, measurement_isotope)
    avg_freq(['/Users/jonaskarthein/cernbox/Python/2017-07-16-85Rb-nuamp-133Cs_run20'])
