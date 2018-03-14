import csv
import os
import ConfigParser
import glob
from array import array
import numpy as np
import time
import multiprocessing
import ROOT


class PiIcr:
    def __init__(self):
        # ---------------------------------------------------------------------------
        # Needed variables / input
        # ---------------------------------------------------------------------------
        self.axis_scaling_factor = 0.031746032  # converts channels to mm
        self.file_name = ''
        self.file_path = os.getcwd()
        self.general = {}
        self.out_files = []
        self.cyc_freq = 0.0
        self.red_cyc_freq = 0.0
        self.cyc_acc_time = 0.0
        self.piicr_excitation = [['Cyclotron frequency / Hz', 'Reduced cyclotron frequency / Hz',
                                  'Magnetron frequency / Hz', 'Cyclotron accumulation time / microseconds']]
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

    def set_file_name(self, name):
        self.file_name = name.split('.')[0]

    def set_file_path(self, path):
        self.file_path = path

    def load_out(self):
        self.read_freq_config()
        return self.general

    def read_freq_config(self):
        '''Reads the config files with the freq information for PI-ICR.'''
        config = ConfigParser.RawConfigParser(allow_no_value=True)  # also allow empty entries
        os.chdir(self.file_path)
        for file in glob.glob('*.out'):
            print file, str(file[:len(self.file_name)]), str(self.file_name)
            # if str(file[:len(self.file_name)]) == str(self.file_name.split('/')[-1]):
            self.out_files.append(file)
        print 'PiIcr.read_freq_config:: out_files = ', self.out_files

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
            self.write_csv(out)
            self.window_edges()

    def write_csv(self, out):
        """The function writes a list in a csv file."""
        os.chdir(self.file_path)
        with open('%s.csv' % (out), 'w') as f:  # save the calculated cyc_acc_time into csv
            writer = csv.writer(f)
            writer.writerows(self.piicr_excitation)

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
        print time_info

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
            print('BAD TRIGGERING, first 5 lines deleted!!!')

        self.load_file = [map(int, x) for x in self.load_file]  # convert list of lists of str's to int
        print self.load_file
        for i in range(0, len(self.load_file), 5):
            for j in range(5):
                if (i%2) == 0:
                    self.raw_data_p1.append(self.load_file[i + j])
                else:
                    self.raw_data_p2.append(self.load_file[i + j])
        self.data_patterns.update({'raw_p1': self.raw_data_p1})
        self.data_patterns.update({'raw_p2': self.raw_data_p2})
        self.multiprocess_reconstruction()
        self.fill_tree()

    def fill_tree(self):
        f = ROOT.TFile.Open('{}.root'.format(self.file_name), 'recreate')
        for option in self.spot_positions:
            tree_name = self.file_name.split(os.path.sep)[-1]
            tree = ROOT.TTree(option, tree_name)
            print 'option = ', option
            # tmp_data = np.zeros((3, 1), dtype=float)
            tmp_data = np.empty((3,1), dtype=float)
            x_list = [x[1] for x in self.spot_positions[option]]
            y_list = [x[2] for x in self.spot_positions[option]]
            t_list = [x[3] for x in self.spot_positions[option]]
            tree.Branch('data', tmp_data, 'x/D:y/D:tof/D')
            for i in range(len(x_list)):
                tmp_data[0][0] = x_list[i]
                # print 'tmp_data[0][0] = ', tmp_data[0][0], 'x_list[i] = ', x_list[i]
                tmp_data[1][0] = y_list[i]
                tmp_data[2][0] = t_list[i]
                tree.Fill()
            tree.Print()
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
        print '\n', self.number_cpu_cores, 'CPU cores found and used for reconstruction.\n'

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
            self.spot_positions.update({'p{}'.format(i):temp})

    def zerolistmaker(self, n):
        '''Creates a list of certain length filled with zeros.'''
        return [0] * (n + 1)

    def tof_check(self, var):
        if var > 1E6:
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
                        # self.spot_positions.append(short_save)
                        que.put(short_save)
                        self.check_event_taken[i][j] += 1
                        self.check_event_taken[i + 1][k] += 1
                        self.check_event_taken[i + 2][l] += 1
                        self.check_event_taken[i + 3][m] += 1
                        self.check_event_taken[i + 4][n] += 1
                        # short_save_2 = [1., 0.]  # next 15 lines save an accepted event into temp_results so that the function window can calculate the time windows again
                        # short_save_2[1] = elements[i][j]
                        # temp_results.append(short_save_2)
                        # short_save_2 = [1., 0.]
                        # short_save_2[1] = elements[i + 1][k]
                        # temp_results.append(short_save_2)
                        # short_save_2 = [1., 0.]
                        # short_save_2[1] = elements[i + 2][l]
                        # temp_results.append(short_save_2)
                        # short_save_2 = [1., 0.]
                        # short_save_2[1] = elements[i + 3][m]
                        # temp_results.append(short_save_2)
                        # short_save_2 = [1., 0.]
                        # short_save_2[1] = elements[i + 4][n]
                        # temp_results.append(short_save_2)
                        # if self.check_event_taken[i][j] > 1 or \
                        #                 self.check_event_taken[i + 1][k] > 1 or \
                        #                 self.check_event_taken[i + 2][l] > 1 or \
                        #                 self.check_event_taken[i + 3][m] > 1 or \
                        #                 self.check_event_taken[i + 4][n] > 1:
                            # print 'Double counting in ejection %s' % ((i+1)/5)
                            # double_taken_event_counter += 1

    def get_spot_positions(self):
        return self.spot_positions