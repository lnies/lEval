import os
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import platform
import time
try:
    import ROOT as root
    root_available = True
except:
    root_available = False

class Jonas():
    def __init__(self):
        if platform.system() == 'Darwin': # MacOS
            self.file_path = 'C:'
        elif platform.system() == 'Windows':
            self.file_path = '/'
        self.header_lines = 1
        self.import_data = []
        self.zoomed_data = []
        self.header = True
        self.columns = ['x', 'y']
        self.usecols = range(len(self.columns))
        self.delimiter = ';'
        self.plot_type = 'scatter'
        self.plot_type_options = {'line' : 'line plot (default)',
                                  'bar' : 'vertical bar plot',
                                  'barh' : 'horizontal bar plot',
                                  'hist' : 'histogram',
                                  'box' : 'boxplot',
                                  'kde' : 'Kernel Density Estimation plot',
                                  'density' : 'same as "kde"',
                                  'area' : 'area plot',
                                  'pie' : 'pie plot',
                                  'scatter' : 'scatter plot',
                                  'hexbin' : 'hexbin plot'}
        self.fit_function_dict = {'x^2': '<A> * x^2 + <B> * x + <C>',
                                'linear': '<M> * x + <N>',
                                'exp': '<A> * exp( <B> * (x - <C>)) + <D>',
                                'sin': '<A> * sin( <B> * x ) + <C>',
                                'sqrt': '<A> * sqrt( <B> * x ) + <C>'}
        self.fit_parameter_num = {'x^2': 3,
                                'linear': 2,
                                'exp': 4,
                                'sin': 3,
                                'sqrt': 3}

        self.figsize = (9,6)
        self.tick_fontsize = 16
        self.png_on_off = True
        self.x_lim = [0, 999999999]
        self.x_lim_zoom = [0, 999999999]
        self.y_lim = [0, 999999999]
        self.y_lim_zoom = [0, 999999999]
        self.fit_option = False
        self.fit_flag = True
        self.fit_function = ''
        self.root_available = root_available
        self.lim_index_list = []
        self.model = 0
        self.popt = []
        self.perr = []
        self.isotope = ''
        self.fig_num = 0
        self.fit_bounds = []
        mpl.rc('font', family='serif', serif='Utopia')
        mpl.rc('text', usetex=False)
        mpl.rc('axes', labelsize=22)    # fontsize of the x and y labels


    def load(self, new_file_path_or_df_or_list, delimiter=self.delimiter):
        '''Loads columns from text (.csv, .txt, ...) or .xlsx files or a dataframe or a list and returns a dataframe'''
        if type(new_file_path_or_df_or_list)==str:
            self.file_path = new_file_path_or_df_or_list
            if os.path.splitext(self.file_path)[-1] == '.xlsx':
                self.load_xlsx()
            else:
                self.load_data_df()

            print('\nSuccessfully loaded :: {}'.format(os.path.split(self.file_path)[-1]))
            print('Location            :: {}\n'.format(os.path.split(self.file_path)[0]))


        elif type(new_file_path_or_df_or_list) == pd.core.frame.DataFrame:
            self.import_data = new_file_path_or_df_or_list
            print('\nSuccessfully loaded :: list\n')

        elif type(new_file_path_or_df_or_list) == list:
            self.import_data = pd.DataFrame(new_file_path_or_df_or_list[1:], columns=new_file_path_or_df_or_list[0])
            print('\nSuccessfully loaded :: DataFrame\n')

        return(self.import_data)


    def load_data_df(self):
        '''Imports data from raw text file (.csv, .txt, ...) to dataframe'''
        self.import_data = pd.read_table(filepath_or_buffer=self.file_path,
                                         delimiter=self.delimiter,
                                         header=range(self.header_lines)[-1],
                                         # usecols=range(len(self.columns)),
                                         # usecols=self.usecols,#range(len(self.columns)),
                                         )
        if self.header:
            self.header = [str(x) for x in list(self.import_data.columns)]


    def load_xlsx(self):
        '''Imports data from an .xlsx file as dataframe'''
        self.import_data = pd.ExcelFile(self.file_path).parse(pd.ExcelFile(self.file_path).sheet_names[0])
        if self.header:
            self.header = [str(x) for x in list(self.import_data.columns)]


    def df_to_dict(self):
        for i in self.header:
            self.df_dict[i] = self.import_data[i].tolist()


    def plot_piicr(self, new_file_path_or_df_or_list, columns, isotope='', delimiter=';', header=''):
        '''Plot import_data'''
        # plt.ion()
        self.delimiter = delimiter
        if header != '':
            self.header = header
        self.columns=columns
        if usecolumns == '':
            self.usecols = range(len(self.columns))
        else:
            self.usecols = usecolumns
        if isotope != '':
            self.isotope=isotope
        self.set_font()
        self.load(new_file_path_or_df_or_list)


        self.fig_num += 1
        fig = plt.figure(self.fig_num, figsize=self.figsize)
        ax = fig.add_subplot(111)

        self.import_data.plot(x=self.import_data.columns[1],
                              y=self.import_data.columns[2],
                              kind='scatter',
                              figsize=(9,9),
                              fontsize=self.tick_fontsize,
                              ax=ax,
                              color='k',
                              ms=3
                              )

        plt.axis([-700,700,-700,700])


if __name__ == '__main__':
    jonas = Jonas()
    df = jonas.load('/Volumes/Python-DAQ/129mCd_002_p2_spot_positions.csv')
    print df
    # import root_pandas as RP
    # RP.to_root(df, 'df.root', key='dftree')
    jonas.plot_piicr('/Users/jonaskarthein/cernbox/Python/Code/129mCd_002_p2_spot_positions.csv', ['x', 'y'], isotope='98Mo', delimiter=';', header=['bla / $\mu$s', 'bla'])





