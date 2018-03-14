import pandas as pd
import numpy as np
import os
import ConfigParser
import csv
import itertools
from fortranformat import FortranRecordReader
from python_plotter_functions2 import *
import glob


def read_excel(file_path, file_name):
    """
    The function reads an excel sheet and returns it as np.array.

    It reads only the first sheet inside the excel sheet. It expects a header
    in the first line and returns a list of lists with the header in the first
    list saved as str and the following as floats.
    """
    os.chdir(file_path)
    xl = pd.ExcelFile('%s.xlsx' % (file_name))
    xl.sheet_names
    df = xl.parse(xl.sheet_names[0])
    df.head()
    output_header = np.array(df.columns)
    output_unicode = np.array(df)
    output_hilf1 = output_header.tolist()
    output_hilf2 = [i.encode('ascii', 'ignore') for i in output_hilf1]       # convert header from unicode to str
    output = output_unicode.tolist()
    output.insert(0, output_hilf2)      # add header
    return(output)


def read_csv(file_path, file_name):
    """
    The function reads a csv file and returns a list of lists.

    The file name is just the name without the .csv ending. All the output
    is str().
    """
    os.chdir(file_path)

    with open('%s.csv' % (file_name), 'rU') as f:   # opens PW file
        reader = csv.reader(f)
        output = list(list(rec) for rec in csv.reader(f, delimiter=','))  # reads csv into a list of lists
    return(output)


def write_excel(file_path, output_list, output_file_name):
    """The function writes a list in an excel sheet."""
    os.chdir(file_path)
    df = pd.DataFrame(output_list)
    df.to_excel('%s.xlsx' % (output_file_name), header=False, index=False)


def write_csv(file_path, output_list, output_file_name):
    """The function writes a list in a csv file."""
    os.chdir(file_path)
    with open('%s.csv' % (output_file_name), 'wb') as f:   # save the calculated cyc_acc_time into csv
        writer = csv.writer(f)
        writer.writerows(output_list)


def read_freq_config(file_path, file_name):
    '''The function reads the Config files with the freq information for PI-ICR'''
    config = ConfigParser.RawConfigParser(allow_no_value=True)   # also allow empty entries
    os.chdir(file_path)
    out_files = []
    for file in glob.glob("*.out"):
        if str(file[:len(file_name)]) == str(file_name):
            out_files.append(file)
    print out_files
    config.read(out_files[0])
    # get frequencies and calculate the cyc acc time = (n_plus_rounds / nu_plus) + (n_minus_rounds / nu_minus)
    cyc_freq = config.getfloat('UT_P2', 'setfreq2')
    red_cyc_freq = config.getfloat('UT_P2', 'setfreq1')
    mag_freq = config.getfloat('UT_Mag', 'setfreq1')
    cyc_acc_time = config.getfloat('Delays', 'accu_p') * 1000000    # in micro seconds
    n_acc = config.getfloat('UT_P2', 'n_acc')

    # save info in a csv
    piicr_excitation = [['Cyclotron frequency / Hz', 'Reduced cyclotron frequency / Hz',
                         'Magnetron frequency / Hz', 'Cyclotron accumulation time / microseconds', 'Applied accumulation rounds']]
    piicr_excitation.append([cyc_freq, red_cyc_freq, mag_freq, cyc_acc_time, n_acc])
    write_csv(file_path, piicr_excitation, 'piicr_excitation')

    return(piicr_excitation)

# if __name__ == '__main__':
#     file_path = '/Volumes/Macintosh HD/Users/jonaskarthein/cernbox/Python/8587Rb_over_night_26.04.17_failed/87Rb/freq_config'
#     file_name = '87Rb'
#     read_freq_config(file_path, file_name)


def read_mpant(file_path, file_name):
    '''
    Read MPANT files and returns for all sections list of (name, value) in a dict.

    The function uses ConfigParser. Find very helpful information at:
    https://docs.python.org/2/library/configparser.html
    '''
    config = ConfigParser.RawConfigParser(allow_no_value=True)     # also allow empty entries
    os.chdir(file_path)
    config.read(file_name)          # read MPANT file to process it
    sections = config.sections()    # get list of sections
    mpant_dict = {}                 # dictionary with sections and lists of their (name, value)
    for i in sections:
        mpant_dict['%s' % (i)] = config.items(i)
    return mpant_dict


def mpant_cyc_acc_time(file_or_folder, folder_path, mpant_file_name):
    '''
    The function imports the data from MPANT and returns the cyc acc time / ns.

    It is possible to either read just one MPANT file (-> file_or_folder ==
    'file') or to process all MPANT files in the according 'data_folder_path'.

    The function uses ConfigParser. Find very helpful information at:
    https://docs.python.org/2/library/configparser.html
    '''
    config = ConfigParser.RawConfigParser(allow_no_value=True)     # also allow empty entries
    mpant_files = []

    if file_or_folder == 'file':                # search for mpant files and save filename(s)
        os.chdir(folder_path)
        mpant_files.append(mpant_file_name)
    elif file_or_folder == 'folder':            # in case of only one specific file, give name as 'mpant_file_name'
        os.chdir(folder_path)
        for file in glob.glob("*.mpa"):
            mpant_files.append(file)
    print mpant_files

    cyc_acc_time = [['Filename', 'Cyclotron accumulation time', 'Err']]
    time_offset = 0     # in ns

    for i in mpant_files:
        begin_line = []
        lookup = 'TDAT'
        with open(i) as myFile:
            for num, line in enumerate(myFile, 1):      # save histogram + header in list
                if lookup in line:
                    begin_line.append(num)
        print begin_line
        load_file = []
        config.readfp(open('%s' % (i)))          # read MPANT file to process it
        bin_range = config.getint('MPA4A', 'range')     # get bin range
        which_begin_line = 0    # sometimes 2 channels are active. here I decide which histogram to fit
        time_offset = 0
        terminate = ''
        for j in ['CHN1', 'CHN2', 'CHN3', 'CHN4', 'CHN5', 'CHN6']:
            if config.getint(j, 'active') == 1 and config.getint(j, 'TOTALSUM') > 10:
                time_offset = float(config.getfloat(j, 'caloff'))      # in nanoseconds
            elif config.getint(j, 'active') == 1 and time_offset == 0:
                which_begin_line += 1
                if which_begin_line == len(begin_line):
                    terminate = 'now'
        if terminate == '':
            histogram_data = [['Ch. / 100 ps', '#']]
            histogram_data_reduced = [['Ch. / 100 ps', '#']]
            with open(i, 'rb') as infile:           # delete header from imported list
                load_file = [[str(h) for h in line.strip().split()] for line in infile]
            maxi = 0

            for k in range(begin_line[which_begin_line]+1, begin_line[which_begin_line]+bin_range, 1):  # get histogram data from file (not possible with configparser with the non-standard mpant files -.-)
                help_load_file = [float(k-begin_line[which_begin_line])]
                help_load_file.extend([float(l) for l in load_file[k]])
                if help_load_file[1] > maxi:        # get the maximum for reducing fit range
                    maxi = help_load_file[1]
                    maxIndex = k-begin_line[which_begin_line]
                histogram_data.append(help_load_file)

            for o in range(maxIndex-150, maxIndex+150, 1):
                help_hist = histogram_data[o]
                histogram_data_reduced.append(help_hist)

            python_plot(histogram_data, 'Cyc.acc.time histogram', os.path.splitext(i)[0], 'no', '', 'yes', 'gauss', 'Utopia', 'red', 'green', 'partly', int(maxIndex-150), int(maxIndex+150), maxIndex, 'scatter', 20, 2, 'off', 'off', '', '', '')
            xss, yss, yerr, parameter_mlf, parameter_mlf_err = python_plot(histogram_data_reduced, 'Cyc.acc.time histogram', '%s_zoom' % os.path.splitext(i)[0], 'no', '', 'yes', 'gauss', 'Utopia', 'red', 'green', 'full', int(maxIndex-150), int(maxIndex+150), maxIndex, 'scatter', 20, 3, 'off', 'off', '', '', '')
            hilf_cyc_acc_time = []
            hilf_cyc_acc_time = [str(os.path.splitext(i)[0]), float(parameter_mlf[1]), float(parameter_mlf_err[2]), time_offset]
            cyc_acc_time.append(hilf_cyc_acc_time)
    if terminate == '':
        for b in range(1, len(cyc_acc_time), 1):
            cyc_acc_time[b][1] = cyc_acc_time[b][1] * 0.1 + cyc_acc_time[b][3]      # 1 channel = 100ps and offset in ns
            cyc_acc_time[b][2] = cyc_acc_time[b][2] * 0.1                    # 1 channel = 100ps and offset in ns

        write_csv(folder_path, cyc_acc_time, 'cyc_acc_time')
        write_excel(folder_path, cyc_acc_time, 'cyc_acc_time')

        return(cyc_acc_time) # if you just need the histogram data, simply add histogram_data to return


def tmp_red_cyc_freq(file_or_folder, folder_path, tmp_file_name):
    '''
    Imports the data from tmpxxx.dat and returns the red cyc freq / Hz.

    It is possible to either read just one tmp.dat file (-> file_or_folder ==
    'file') or to process all tmp.dat files in the according 'data_folder_path'.

    The function uses ConfigParser. Find very helpful information at:
    https://docs.python.org/2/library/configparser.html
    '''
    config = ConfigParser.RawConfigParser(allow_no_value=True)     # also allow empty entries

    tmp_files = []

    if file_or_folder == 'file':                # save filename(s)
        os.chdir(folder_path)
        tmp_files.append(tmp_file_name)
    elif file_or_folder == 'folder':
        os.chdir(folder_path)
        for file in glob.glob("*.dat"):
            if 'corrected_' in file:
                pass
            else:
                tmp_files.append(file)

    excitation_frequencies = [['File name', 'Cyclotron frequency / Hz', 'Reduced cyclotron frequency / Hz', 'Magnetron frequency / Hz']]

    for i in tmp_files:
        cheap_line_counter = 0
        replacements = {',': '\n', ' ': ''}          # correct tmp for proper config standard
        with open('%s' % i) as infile, open('corrected_%s' % i, 'w') as outfile:
            for line in infile:
                cheap_line_counter += 1
                if cheap_line_counter < 133:
                    for src, target in replacements.iteritems():
                        line = line.replace(src, target)
                    outfile.write(line)

        load_file = []
        fh = open('corrected_%s' % i, 'rb')
        fh.readline()                                           # discard first line
        config.readfp(fh)
        cyc_freq_help = config.getfloat('Excit', 'Freq')          # get cyc freq
        mag_freq_help = config.getfloat('Cooling1', 'Freq')       # get cyc freq
        red_cyc_freq_help = cyc_freq_help - mag_freq_help
        hilf_23 = [i, cyc_freq_help, red_cyc_freq_help, mag_freq_help]
        excitation_frequencies.append(hilf_23)

    write_csv(folder_path, excitation_frequencies, 'excitation_frequencies')
    write_excel(folder_path, excitation_frequencies, 'excitation_frequencies')

    return(excitation_frequencies)


def read_ame_table():
    '''
    The function reads the AME table and returns a list of lists and a DataFrame.

    Input name:     'AME2012.mas12.ff'
    File info:      Before starting please remove all '*' by ' '
                    and all '#' by '.'
    Needed package: https://pypi.python.org/pypi/fortranformat
    Created by:     Jonas Karthein (jonas.karthein@cern.ch)
    Changelog:      V1.0 - 11.04.2017 - Jonas Karthein
    '''
    # fortran list structure, you can find structure code in header of AME-table
    fline=FortranRecordReader('(a1, i3, i5, i5, i5, 1x, a3, a4, 1x, f13.5, f11.5,\
                                f11.3, f9.3, 1x, a2, f11.3, f9.3, 1x, i3, 1x, f12.5,\
                                f11.5)')
    # AME-table header, you can find it in AME-table file
    header = ['cc', 'NZ', 'N', 'Z', 'A', 'el', 'o', 'mass excess / keV',\
              'mass excess unc / keV', 'binding energy / keV',\
              'binding energy unc / keV', 'B', 'beta decay energy / keV',\
              'beta decay energy unc / keV', 'atomic mass (int) / Dalton',\
              'atomic mass / micro Dalton', 'atomic mass unc']
    ame_table = []

    with open('AME16.txt') as f:
        for line in itertools.islice(f, 39, None):  # skips header (length see AME table)
            nucl = fline.read(line)                 # using fortran structure to read
            ame_table.append(nucl)

    for i in range(len(ame_table)):  # entry 16 was only the float part in micro Dalton
        ame_table[i][15] = ame_table[i][14] * 1E6 + ame_table[i][15]  # now full atomic mass

    ame_table_df = pd.DataFrame(ame_table, columns=header)  # create DataFrame
    ame_table.insert(0, header)     # add header to list of lists

    return(ame_table, ame_table_df)

    # ### needs to be tested first!!!!

    # ame_import = np.genfromtxt('mass16.txt', skip_header=39, dtype=['a1', 'int', 'int', 'int', 'int', 'a4', 'a4', 'float', 'float', 'float', 'float', 'a3', 'float', 'float', 'int', 'float', 'float'], delimiter=[1,3,5,5,5,4,4,16,12,12,6,3,11,9,4,13,10])
    # ame_table = [list(elem) for elem in ame_import.tolist()]    # convert list of tuples to list of lists
    # for i in range(len(ame_table)):
    #     ame_table[i][15] = ame_table[i][14] * 1E6 + ame_table[i][15]    # calculate full atomic mass (int*1000000+float)
    #     for j in range(len(ame_table[i])):
    #         if type(ame_table[i][j]) == str:
    #             ame_table[i][j] = ame_table[i][j].replace(" ", "")      # delete all spaces --> makes searching easier
    #         ame_table[i][j] = str(ame_table[i][j])                      # convert all entries to int (needed for Dinkos class)
    # ame_table_df = pd.DataFrame(ame_table, columns=['cc', 'NZ', 'N', 'Z', 'A', 'el', 'o', 'mass excess / keV', 'mass excess unc / keV', 'binding energy / keV', 'binding energy unc / keV', 'B', 'beta decay energy / keV', 'beta decay energy unc / keV', 'atomic mass (int) / Dalton', 'atomic mass / micro Dalton', 'atomic mass unc'])
    # print ame_table_df


if __name__ == '__main__':
    print read_csv('E:\\Jonas\\Harmonisierung', 'Efield_scan_1V45406_p3_part1_spot_positions')

