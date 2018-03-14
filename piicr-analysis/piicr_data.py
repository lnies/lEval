import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import csv
import os
import ConfigParser
import glob
from array import array
import numpy as np
import time
import multiprocessing
import ROOT as root
from collections import Counter
import json
import pandas as pd


class PIICR_data:
    def __init__(self):
        # ---------------------------------------------------------------------------
        # Needed variables / input
        # ---------------------------------------------------------------------------
        self.axis_scaling_factor = 0.031746032  # converts channels to mm
        self.file_name = ''
        self.txt_files = []
        self.folder_path = os.getcwd()
        self.script_path = os.getcwd()
        self.general = {}
        self.out_files = []
        self.cyc_freq = 0.0
        self.red_cyc_freq = 0.0
        self.cyc_acc_time = 0.0
        self.piicr_excitation = [['Cyclotron frequency / Hz',
                                  'Reduced cyclotron frequency / Hz',
                                  'Magnetron frequency / Hz',
                                  'Cyclotron accumulation time / microseconds']]
        self.number_cpu_cores = 2
        self.win_x_low = 1500
        self.win_x_high = 1700
        self.win_y_low = 1500
        self.win_y_high = 1700
        self.win_xy_low = 3100
        self.win_xy_high = 3300
        self.spot_positions = {} # 1st column = number of ejection, 2nd column x position, 3rd column = y position, 4th column ToF
        self.check_event_taken = []
        self.load_file = []
        self.raw_data_p1 = []  # now create a list of lists depending on the pattern.
        self.raw_data_p2 = []  # now create a list of lists depending on the pattern.
        self.data_patterns = {}
        self.tree_dict = {}
        self.trees = {}
        self.z_classes = [[1, 1], [2, 2], [3, 3], [4, 4], [5, 5],
                          [1, 2], [1, 3], [1, 4], [1, 5]]
        self.z_class_analysis = [False, 0, 100]
        self.tof_range = []
        mpl.rc('font', family='serif', serif='Utopia')      # Utopia LaTeX font!!
        mpl.rc('text', usetex=False)




    def batch(self, folder_path_input, z_class, tof_range):
        self.set_folder_path(folder_path_input)
        self.set_z_class(z_class)
        self.set_tof_range(tof_range)
        self.get_txt()
        self.get_out()
        print '\n### Batch running'
        for txt_file in self.txt_files:
            self.set_file_name(txt_file)
            print 'Reconstruction started   :: ', self.file_name
            self.load_raw_data()
            self.z_class_reduction()
            self.tof_cut()
# ellipse cut
            self.fill_tree()
            self.save_spot_positions()
# fit
        return(self.spot_positions)


    def single(self, file_path_input):
        return(False)



    def get_txt(self):
        for file in glob.glob('*.txt'):
            self.txt_files.append(file)
        print 'Data files in the batch  :: ', self.txt_files

    def set_file_name(self, name):
        self.file_name = name.split('.')[0]

    def set_z_class(self, z_class):
        if z_class != []:
            self.z_class_analysis[0] = True
            self.z_class_analysis[1] = z_class[0]
            self.z_class_analysis[2] = z_class[1]

    def set_tof_range(self, tof_range):
        if tof_range != []:
            self.tof_range = tof_range

    def set_folder_path(self, path):
        self.folder_path = path
        os.chdir(path)

    def get_out(self):
        self.read_freq_config()
        return self.general

    def read_freq_config(self):
        '''Reads the config files with the freq information for PI-ICR.'''
        config = ConfigParser.RawConfigParser(allow_no_value=True)  # also allow empty entries
        os.chdir(self.folder_path)
        for file in glob.glob('*.out'):
            self.out_files.append(file)
        print 'Config files found       :: ', self.out_files

        for out in self.out_files:
            config.read(out)
            self.general.update({'cyc_freq' : config.getfloat('UT_P2', 'setfreq2')})
            self.general.update({'red_cyc_freq' : config.getfloat('UT_P2', 'setfreq1')})
            self.general.update({'mag_freq': config.getfloat('UT_Mag', 'setfreq1')})
            self.general.update({'cyc_acc_time': 1E6 * config.getfloat('Delays', 'accu_p')})
            #mcp_ps_win  = [SumX, SumY, SumX-Sigma, SumY-Sigma, SumXY, SumXY-Sigma, Sigma-acceptance-factor]
            self.general.update({'sumx':1600})
            self.general.update({'sumy': 1600})
            self.general.update({'sumx-sigma': 30})
            self.general.update({'sumy-sigma': 30})
            self.general.update({'sumxy': 3200})
            self.general.update({'sumxy-sigma': 50})
            self.general.update({'factor': 7})


            # get frequencies and calculate the cyc acc time = (n_plus_rounds / nu_plus) + (n_minus_rounds / nu_minus)
            self.cyc_freq = config.getfloat('UT_P2', 'setfreq2')
            self.red_cyc_freq = config.getfloat('UT_P2', 'setfreq1')
            self.mag_freq = config.getfloat('UT_Mag', 'setfreq1')
            self.cyc_acc_time = config.getfloat('Delays', 'accu_p') * 1E6  # in micro seconds
            # save info in a csv
            self.piicr_excitation.append([config.getfloat('UT_P2', 'setfreq2'), config.getfloat('UT_P2', 'setfreq1'),
                                          config.getfloat('UT_Mag', 'setfreq1'), 1E6 * config.getfloat('Delays', 'accu_p')])
        self.write_csv('piicr_excitation', self.piicr_excitation)
        self.window_edges()


    def write_csv(self, name, list_of_lists):
        with open('{}.csv'.format(name), 'w') as f:  # save the calculated cyc_acc_time into csv
            writer = csv.writer(f)
            writer.writerows(list_of_lists)


    def load_raw_data(self):
        """
        Data loading function.

        The function loads the given raw data file and saves the data depending on the
        way the raw data was saved in the .txt (e.g. p1 and p2 alternating). In the
        example case it would only save p1 OR p2.
        """

        with open('{}.txt'.format(self.file_name), 'r') as infile:
            self.load_file = [[str(i) for i in line.strip().split()] for line in infile]

        time_info_hilf = self.load_file[-1]  # The time info is written in the last line.
        time_info = [i for i in time_info_hilf]

        self.load_file.pop(-1)  # delete time stamp
        self.load_file.pop(-1)  # when the time stamps are added, 5 empty lines are also added (we don't know why)
        self.load_file.pop(-1)
        self.load_file.pop(-1)
        self.load_file.pop(-1)
        self.load_file.pop(-1)
        # check if there was a bad triggering at the beginning.
        # This would lead to a swap in P1/P2, therefore the first ejection must be deleted.
        if not len(self.load_file) % 10 == 0:
            self.load_file.pop(0)
            self.load_file.pop(0)
            self.load_file.pop(0)
            self.load_file.pop(0)
            self.load_file.pop(0)
            print('File info                ::  Bad trigger found, first ejection will be deleted!')

        self.load_file = [map(int, x) for x in self.load_file]  # convert list of lists of str's to int
        self.raw_data_p1 = []
        self.raw_data_p2 = []
        for i in range(0, len(self.load_file), 5):
            for j in range(5):
                if (i%2) == 0:
                    self.raw_data_p1.append(self.load_file[i + j])
                else:
                    self.raw_data_p2.append(self.load_file[i + j])

        self.data_patterns['raw_p1'] = self.raw_data_p1
        self.data_patterns['raw_p2'] = self.raw_data_p2
        self.multiprocess_reconstruction()

    def save_spot_positions(self):
        with open('{}_spot_positions.json'.format(self.file_name), 'w') as f:
            json.dump(self.spot_positions,f)

    def fill_tree(self):
        f = root.TFile.Open('{}.root'.format(self.file_name), 'recreate')
        for option in self.spot_positions:
            tree_name = self.file_name.split(os.path.sep)[-1]
            tree = root.TTree(option, tree_name)
            # tmp_data = np.zeros((3, 1), dtype=float)
            tmp_data = np.empty((4,1), dtype=float)
            x_list = [x[1] for x in self.spot_positions[option]['raw_data']]
            y_list = [x[2] for x in self.spot_positions[option]['raw_data']]
            t_list = [x[3] for x in self.spot_positions[option]['raw_data']]
            ej_list = [x[0] for x in self.spot_positions[option]['raw_data']]
            tree.Branch('raw_data', tmp_data, 'x/D:y/D:tof/D:ej/D')
            for i in range(len(x_list)):
                tmp_data[0][0] = x_list[i]
                # print 'tmp_data[0][0] = ', tmp_data[0][0], 'x_list[i] = ', x_list[i]
                tmp_data[1][0] = y_list[i]
                tmp_data[2][0] = t_list[i]
                tmp_data[3][0] = ej_list[i]
                tree.Fill()
            if self.z_class_analysis[0]:
                tmp_data_2 = np.empty((4,1), dtype=float)
                x_list = [x[1] for x in self.spot_positions[option]['Z{}-{}'.format(self.z_class_analysis[1], self.z_class_analysis[2])]]
                y_list = [x[2] for x in self.spot_positions[option]['Z{}-{}'.format(self.z_class_analysis[1], self.z_class_analysis[2])]]
                t_list = [x[3] for x in self.spot_positions[option]['Z{}-{}'.format(self.z_class_analysis[1], self.z_class_analysis[2])]]
                ej_list = [x[0] for x in self.spot_positions[option]['Z{}-{}'.format(self.z_class_analysis[1], self.z_class_analysis[2])]]
                tree.Branch('cut_data', tmp_data_2, 'x/D:y/D:tof/D:ej/D')
                for i in range(len(x_list)):
                    tmp_data_2[0][0] = x_list[i]
                    # print 'tmp_data[0][0] = ', tmp_data[0][0], 'x_list[i] = ', x_list[i]
                    tmp_data_2[1][0] = y_list[i]
                    tmp_data_2[2][0] = t_list[i]
                    tmp_data_2[3][0] = ej_list[i]
                    tree.Fill()
                # tree.Print()
                tree.Write()
            else:
                tmp_data_2 = np.empty((4,1), dtype=float)
                x_list = [x[1] for x in self.spot_positions[option]['tof_cut']]
                y_list = [x[2] for x in self.spot_positions[option]['tof_cut']]
                t_list = [x[3] for x in self.spot_positions[option]['tof_cut']]
                ej_list = [x[0] for x in self.spot_positions[option]['tof_cut']]
                tree.Branch('cut_data', tmp_data_2, 'x/D:y/D:tof/D:ej/D')
                for i in range(len(x_list)):
                    tmp_data_2[0][0] = x_list[i]
                    # print 'tmp_data[0][0] = ', tmp_data[0][0], 'x_list[i] = ', x_list[i]
                    tmp_data_2[1][0] = y_list[i]
                    tmp_data_2[2][0] = t_list[i]
                    tmp_data_2[3][0] = ej_list[i]
                    tree.Fill()
                # tree.Print()
                tree.Write()
        f.Close()

    def divide(self, list_in, parts, pattern):
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
            dict_parts["raw_p{}_part{}".format(pattern, i+1)]=hilf
        hilf = []
        for j in range((parts-1)*part_length, length):
            hilf.append(list_in[j])
        dict_parts["raw_p{}_part{}".format(pattern, parts)]=hilf
        return dict_parts

    def multiprocess_reconstruction(self):
        '''Function to prepare multiprocess reconstruction'''
        self.number_cpu_cores = multiprocessing.cpu_count()
        print 'Parallelization          :: ', self.number_cpu_cores, 'CPU cores in use'

        for i in [1, 2]:
            data_p_part = self.divide(self.data_patterns['raw_p{}'.format(i)], self.number_cpu_cores, i)

            m = multiprocessing.Manager()
            q = m.Queue()        # shared queue (= shared memory)

            t1 = time.time()
            if self.number_cpu_cores == 2:
                proc1 = multiprocessing.Process(target=self.reconstruction, args=([data_p_part['raw_p{}_part1'.format(i)], 1], q))
                proc2 = multiprocessing.Process(target=self.reconstruction, args=([data_p_part['raw_p{}_part2'.format(i)], 2], q))

                proc1.start()
                proc2.start()

                proc1.join()
                proc2.join()
            elif self.number_cpu_cores == 4:
                proc1 = multiprocessing.Process(target=self.reconstruction, args=([data_p_part['raw_p{}_part1'.format(i)], 1], q))
                proc2 = multiprocessing.Process(target=self.reconstruction, args=([data_p_part['raw_p{}_part2'.format(i)], 2], q))
                proc3 = multiprocessing.Process(target=self.reconstruction, args=([data_p_part['raw_p{}_part3'.format(i)], 3], q))
                proc4 = multiprocessing.Process(target=self.reconstruction, args=([data_p_part['raw_p{}_part4'.format(i)], 4], q))

                proc1.start()
                proc2.start()
                proc3.start()
                proc4.start()

                proc1.join()
                proc2.join()
                proc3.join()
                proc4.join()
            elif self.number_cpu_cores == 8:
                proc1 = multiprocessing.Process(target=self.reconstruction, args=([data_p_part['raw_p{}_part1'.format(i)], 1], q))
                proc2 = multiprocessing.Process(target=self.reconstruction, args=([data_p_part['raw_p{}_part2'.format(i)], 2], q))
                proc3 = multiprocessing.Process(target=self.reconstruction, args=([data_p_part['raw_p{}_part3'.format(i)], 3], q))
                proc4 = multiprocessing.Process(target=self.reconstruction, args=([data_p_part['raw_p{}_part4'.format(i)], 4], q))
                proc5 = multiprocessing.Process(target=self.reconstruction, args=([data_p_part['raw_p{}_part5'.format(i)], 5], q))
                proc6 = multiprocessing.Process(target=self.reconstruction, args=([data_p_part['raw_p{}_part6'.format(i)], 6], q))
                proc7 = multiprocessing.Process(target=self.reconstruction, args=([data_p_part['raw_p{}_part7'.format(i)], 7], q))
                proc8 = multiprocessing.Process(target=self.reconstruction, args=([data_p_part['raw_p{}_part8'.format(i)], 8], q))

                proc1.start()
                proc2.start()
                proc3.start()
                proc4.start()
                proc5.start()
                proc6.start()
                proc7.start()
                proc8.start()

                proc1.join()
                proc2.join()
                proc3.join()
                proc4.join()
                proc5.join()
                proc6.join()
                proc7.join()
                proc8.join()
            temp = []
            while q.empty() is False:       # save all reconstructed spots which are put into a queue
                temp.append(q.get())

            self.spot_positions['p{}'.format(i)] = {'raw_data': temp}


    def zerolistmaker(self, n):
        '''Creates a list of certain length filled with zeros.'''
        return [0] * (n + 1)

    def tof_check(self, var):
        if var > 250000:
            return True
        else:
            return False

    def window_edges(self):
        self.win_x_low = self.general['sumx'] - self.general['sumx-sigma'] * self.general['factor']
        self.win_x_high = self.general['sumx'] + self.general['sumx-sigma'] * self.general['factor']
        self.win_y_low = self.general['sumy'] - self.general['sumy-sigma'] * self.general['factor']
        self.win_y_high = self.general['sumy'] + self.general['sumy-sigma'] * self.general['factor']
        self.win_xy_low = self.general['sumxy'] - self.general['sumxy-sigma'] * self.general['factor']
        self.win_xy_high = self.general['sumxy'] + self.general['sumxy-sigma'] * self.general['factor']

    def reconstruction(self, input_list, que):
        """
        Event reconstruction function.

        The function reconstructs all possible events in one ejection with 1 time
        stamp in all 5 lines --> event only possible if times are within the
        measurement window(= on the detector) and 1 time per X1, X2, Y1, Y2 and MCP
        channel is found. If data is lost in one or more channels the event is not
        reconstructed.
        Input_list consists of the raw_data in the first entry and a dictionary with further values.
        """
        # ---------------------------------------------------------------------------
        # Reconstruction loop
        # ---------------------------------------------------------------------------

        for i in range(0, len(input_list[0]), 5):  # check each ejection
            # Matrix with size of results to check if the event was used for a match
            self.check_event_taken.append(self.zerolistmaker(input_list[0][i][0]))
            self.check_event_taken.append(self.zerolistmaker(input_list[0][i + 1][0]))
            self.check_event_taken.append(self.zerolistmaker(input_list[0][i + 2][0]))
            self.check_event_taken.append(self.zerolistmaker(input_list[0][i + 3][0]))
            self.check_event_taken.append(self.zerolistmaker(input_list[0][i + 4][0]))
            if input_list[0][i][0] > 0 and input_list[0][i + 1][0] > 0 and input_list[0][i + 2][0] > 0 and input_list[0][i + 3][0] > 0:  # make sure that there is at least one hit at the detector, then build every possible combination (even when a event is missing) and store all reasonable values
                self.x1_loop(input_list[0], que, input_list[1], i)

    def x1_loop(self, elements, que, proc_num, i):
        for j in range(1, elements[i][0] + 1):  # x1
            if self.tof_check(elements[i][j]):  # excludes trigger-errors: saves time
                self.x2_loop(elements, que, proc_num, i, j)

    def x2_loop(self, elements, que, proc_num, i, j):
        for k in range(1, elements[i + 1][0] + 1):  # x2
            if self.tof_check(elements[i + 1][k]):  # excludes trigger-errors: saves time
                self.y1_loop(elements, que, proc_num, i, j, k)

    def y1_loop(self, elements, que, proc_num, i, j, k):
        for l in range(1, elements[i + 2][0] + 1):  # y1
            if self.tof_check(elements[i + 2][l]):  # excludes trigger-errors: saves time
                self.y2_loop(elements, que, proc_num, i, j, k, l)

    def y2_loop(self, elements, que, proc_num, i, j, k, l):
        for m in range(1, elements[i + 3][0] + 1):  # y2
            if self.tof_check(elements[i + 3][m]):  # excludes trigger-errors: saves time
                self.mcp_loop(elements, que, proc_num, i, j, k, l, m)

    def mcp_loop(self, elements, que, proc_num, i, j, k, l, m):
        # input_list[0]
        # short_save = [0., 0., 0.]   # nur ein Zwischenspeicher (buffer variable)
        # short_save_2 = [1., 0.]     # nur ein Zwischenspeicher (buffer variable)
        # double_taken_event_counter = 0
        # temp_results = []
        for n in range(1,  elements[i + 4][0] + 1):  # MCP
            if self.tof_check(elements[i + 4][n]):  # excludes trigger-errors: saves time
                if (elements[i][j] + elements[i + 1][k] - 2 * elements[i + 4][n]) < self.win_x_high and \
                                (elements[i][j] + elements[i + 1][k] - 2 * elements[i + 4][n]) > self.win_x_low and \
                                (elements[i + 2][l] + elements[i + 3][m] - 2 * elements[i + 4][n]) < self.win_y_high and \
                                (elements[i + 2][l] + elements[i + 3][m] - 2 * elements[i + 4][n]) > self.win_y_low and \
                                (elements[i][j] + elements[i + 1][k] + elements[i + 2][l] + elements[i + 3][m] - 4 * elements[i + 4][n]) < self.win_xy_high and \
                                (elements[i][j] + elements[i + 1][k] + elements[i + 2][l] + elements[i + 3][m] - 4 * elements[i + 4][n]) > self.win_xy_low:  # check for x pairs and y pairs in their x/y-measurement window
                    short_save = [0., 0., 0., 0.]
                    short_save[0] = (float(i) + len(elements) * ( proc_num - 1)) / 5
                    short_save[1] = (float(elements[i][j]) - float(elements[i + 1][k])) / 2
                    short_save[2] = (float(elements[i + 2][l]) - float(elements[i + 3][m])) / 2
                    short_save[3] = float(elements[i + 4][n]) * 0.000025  # ch to micro seconds
                    if short_save[1] < self.general['sumx'] and short_save[1] > -self.general['sumx'] and \
                                    short_save[2] < self.general['sumy'] and short_save[2] > -self.general['sumy']:  # check if the matched signal is also on the mcp
                        que.put(short_save)
                        self.check_event_taken[i][j] += 1
                        self.check_event_taken[i + 1][k] += 1
                        self.check_event_taken[i + 2][l] += 1
                        self.check_event_taken[i + 3][m] += 1
                        self.check_event_taken[i + 4][n] += 1

    def get_spot_positions(self):
        return self.spot_positions

    def z_class_reduction(self):
        if self.z_class_analysis[0]:
            if [self.z_class_analysis[1], self.z_class_analysis[2]] not in self.z_classes:
                self.z_classes.append([self.z_class_analysis[1], self.z_class_analysis[2]])
            for option in ['p1', 'p2']:
                for z_class in self.z_classes:
                    self.spot_positions[option]['Z{}-{}'.format(z_class[0], z_class[1])] = self.z_class_reducer(self.spot_positions[option]['raw_data'], z_class)

            print 'Z-class reduction ({}-{})  ::  {}/{}'.format(str(int(self.z_class_analysis[1])), str(int(self.z_class_analysis[2])), len(self.spot_positions[option]['Z{}-{}'.format(z_class[0], z_class[1])]), len(self.spot_positions[option]['raw_data']))
        else:
            'No Z-class reduction performed..'


    def z_class_reducer(self, spot_positions, z_class):
        '''
        :param: spot_positions      = list of lists with entries of ejection number, x-pos and y-pos of reconstructed spot as well as its ToF
        :param: z_class_analysis    = list with three entries: string z-class-analysis 'yes'/'no', min number of ions and max number of ions

        Function which reduces the spot positions list to a list which contains only spots with given number of ions in the trap.
        '''
        ejections = []
        for entry in spot_positions:
            ejections.append(entry[0])

        ejections_dict = dict(Counter(ejections))   # count how many ions per ejection were found

        reduced_spot_positions = []
        for key in ejections_dict:
            if ejections_dict[key] >= z_class[0] and ejections_dict[key] <= z_class[1]:
                for entry in spot_positions:
                    if key == entry[0]:
                        reduced_spot_positions.append(entry)

        return(reduced_spot_positions)


    def tof_cut(self):
        if self.tof_range == []:
            for option in ['p1']:
                tmp_range = self.get_tof_range([x[3] for x in self.spot_positions[option]['raw_data']])
                self.tof_range = [tmp_range[0] - 2.5 * tmp_range[2],
                                  tmp_range[0] + 2.5 * tmp_range[2]]
        print 'ToF window (fit)         ::  [{:.2f} .. {:.2f}] us'.format(self.tof_range[0], self.tof_range[1])
        for option in ['p1', 'p2']:
            tmp_key_list = []
            for key in self.spot_positions[option]:
                tmp_key_list.append(key)
            for key in tmp_key_list:
                if key != 'raw_data':
                    self.spot_positions[option][key+'-tof'] = []
                    for spot in self.spot_positions[option][key]:
                        if spot[3] < self.tof_range[1] and spot[3] > self.tof_range[0]:
                            self.spot_positions[option][key+'-tof'].append(spot)
        self.general['tof_range'] = self.tof_range
        self.tof_range = []     # reset to calculate new for each file



    def get_tof_range(self, liste):
        '''Function takes (1D) list with unbinned data and returns a list with gaussian fit parameters: x, x_unc., sigma, sigma_unc.'''
        root.gErrorIgnoreLevel = root.kInfo
        root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
        import_list = [float(x) for x in liste]

        # get the x where most of the ions arrive
        hist, hilf_edges = np.histogram(import_list, bins=2000, range=(0,200))

        fig = plt.figure(1110002, figsize=(9,6))

        plt.plot(hilf_edges[1:], hist)

        plt.xlabel('time of flight / $\mu$s', fontsize=22)
        plt.ylabel('counts', fontsize=22)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)

        x_edges = []
        for i in range(len(hilf_edges)-1):
            x_edges.append(np.mean([hilf_edges[i], hilf_edges[i+1]]))
        x_max = x_edges[self.index_max_value(hist)]

        plt.plot([x_max], [max(hist)], 'r *')
        plt.tight_layout()
        plt.savefig('{}_tof_total.pdf'.format(self.file_name))
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


    def index_max_value(self, liste):
    ########## find index of maximum in list
        index_max = 0
        max_value = -99999999999.
        for i in range(len(liste)):
            if liste[i] > max_value:
                max_value = liste[i]
                index_max = i
        return(index_max)

if __name__ == '__main__':
    PIICR_data = PIICR_data()
    spots = PIICR_data.batch('/Volumes/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/PI-ICR/129Cd/test-run7', [2, 2], [])
    spots_name_dict = {0: 'ej', 1:'x', 2:'y', 3:'tof'}
    for i in ['p1', 'p2']:
        for j in spots[i].keys():
            for k in range(max([len(spots['p1']['raw_data']), len(spots['p2']['raw_data'])])-len(spots[i][j])):
                spots[i][j].append([0., 0., 0., 0.])
    print [x[0] for x in spots['p1']['Z2-2']]

    for i in spots.keys():
        for j in spots[i].keys():
            for k in range(max([len(spots['p1']['raw_data']), len(spots['p2']['raw_data'])])):
                for l in [0,1,2,3]:
                    print spots[i][j][k][l]

    df = pd.DataFrame({(i,j,spots_name_dict[l]): spots[i][j][k][l]
                           for i in spots.keys()
                           for j in spots[i].keys()
                           for k in range(max([len(spots['p1']['raw_data']), len(spots['p2']['raw_data'])]))
                           for l in [0,1,2,3]}, index=range(max([len(spots['p1']['raw_data']), len(spots['p2']['raw_data'])])))

    # das scheiß Ding funktioniert nicht so wie es soll, vgl letzte for loop mit dem ergebnis des df
    print df
    print df[('p1', 'Z2-2', 'ej')].tolist()
