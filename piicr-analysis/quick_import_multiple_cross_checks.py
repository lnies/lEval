import os
import pandas as pd
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import datetime



class Import_multi_cross_checks():
    def __init__(self):
        self.isotopes = {'REF': ['133Cs', '129Cs'],
                         'MEAS': ['129gCd', '129mCd']}
        self.modes = ['ground', 'isomer']
        self.naming = {'129gCdground': '129Cd-ground-ov',
                       '129gCdisomer': '129Cd-isomer-non-ov',
                       '129mCdground': '129Cd-ground-non-ov',
                       '129mCdisomer': '129Cd-isomer-ov',}
        self.wd = ''
        self.run_paths = []


    def set_font(self):
        '''Set Utopia font and label font size'''
        mpl.rc('font', family='serif', serif='Utopia')
        mpl.rc('text', usetex=False)
        mpl.rc('axes', labelsize=22)    # fontsize of the x and y labels


    def load_cross(self, upper_folder_path):
        self.wd = upper_folder_path
        self.get_runs()
        df = self.load_all_csv()
        df.to_csv(self.wd+'/all_fit_data.csv')
        return(df)


    def get_runs(self):
        for folder_name in os.listdir(self.wd):
            if folder_name[0] == '.':       # skip folders like '.DS_Store', ...
                pass
            else:
                self.run_paths.append(self.wd+'/'+folder_name)


    def load_cross_csv(self, run, mode, meas):
        for folder_name in os.listdir(run+'/cross-checks'):
            if meas in folder_name:
                temp_dir = run+'/cross-checks/'+folder_name+'/cross_checks_{}_{}.csv'.format(folder_name, mode)
                temp_df = pd.read_csv(temp_dir, delimiter=',', skiprows=5)
                temp_df.insert(0, 'Run', run[-3:]+'-'+folder_name[:5])
                try:
                    df_zs = df_zs.append(temp_df, ignore_index=True)
                except NameError:
                    df_zs = temp_df
        return(df_zs)


    def load_all_cross_csv(self):
        for meas in self.isotopes['MEAS']:
            for mode in self.modes:
                for run in self.run_paths:
                    temp_df_in = self.load_csv(run, mode, meas)
                    temp_df_in.insert(0, 'isotope', self.naming[meas+mode])
                    try:
                        df = df.append(temp_df_in, ignore_index=True)
                    except NameError:
                        df = temp_df_in
        return(df)


    def load_all_spot_positions(self, folder_path):
        for path, subdirs, files in os.walk(folder_path):
            for filename in files:
                if 'positions.csv' in filename and 'Cd' in filename and not '_c_' in filename:
                    temp_df = self.load_spot_positions(path+'/'+filename)
                    try:
                        df = pd.concat([df, temp_df], axis=1)
                    except NameError:
                        df = temp_df.copy()
        df.to_csv(folder_path+'/all_129Cd_spot_positions_x.csv')
        import matplotlib.pyplot as plt
        print 'Total sum of candmium counts per ejection', df.sum(axis=1)
        df.sum(axis=1).plot(kind='bar')
        plt.xlabel('cadmium counts per ejection')
        plt.ylabel('total sum of counts')
        plt.savefig(folder_path+'/all_129Cd_spot_positions_x.pdf')
        plt.show()


    def load_spot_positions(self, file_path):
        '''Loads spot_positions by file_path and returns pandas series with number of counts per ej'''
        temp_df = pd.read_csv(file_path, delimiter=',', names=['ej', 'x', 'y', 'tof'])
        series = pd.cut(temp_df['ej'], bins=range(0,int(temp_df.max(0)['ej']+1))).value_counts()
        counts_bin = pd.cut(series, bins=range(-1,10+1), labels=[str(i) for i in range(10+1)]).value_counts().sort_index()
        counts_bin_df = counts_bin.to_frame('{}-{}'.format(file_path.split('/')[8], file_path.split('/')[11].split('.')[0]))
        return(counts_bin_df)
        # import matplotlib.pyplot as plt
        # counts_bin.plot(kind='bar')
        # plt.show()


    def load_plot_ISOLTRAP_temp(self, file_path):
        df = pd.read_csv(file_path, delimiter=';', parse_dates={'Date_Time': [0,1]}, index_col='Date_Time')
        self.set_font()
        df.ix[datetime.datetime(year=2017, month=7, day=2, hour=8):datetime.datetime(year=2017, month=7, day=2, hour=11)].plot(use_index=True, figsize=(4,6), y='UT_lower', fontsize=16)


        plt.tight_layout()
        plt.savefig('{}/ut_lower_temperature.pdf'.format(os.path.dirname(file_path)))
        plt.show()


    def read_counts_info_dict(self, folder_path):
        counts_info_dict = {}
        count_rate = {}

        for i in os.walk(folder_path).next()[1]:
            file_path = folder_path+'/'+i+'/129gCd/p1p2/counts_info_dict.json'
            if os.path.isfile(file_path):
                df_in = pd.read_json(file_path)
                for name in df_in.ix['#rec-spots-distribution'].axes[0]:
                    counts_info_dict[i+'/'+str(name)] = df_in.ix['#rec-spots-distribution'][name]
                    count_rate[i+'/'+str(name)] = sum(df_in.ix['#rec-spots-distribution'][name])

        df = pd.DataFrame(counts_info_dict)
        df_2 = pd.DataFrame(count_rate, index=['counts'])
        self.set_font()
        fig = plt.figure(1, figsize=(9,6))
        ax = fig.add_subplot(111)
        df.sum(axis=1).plot(kind='bar', axes=ax, fontsize=16, title='IS574 $^{129g}$Cd')
        print df_2.sum(axis=1)['counts']
        ax.set_xlabel('Z-class')
        ax.set_ylabel('$\Sigma$ counts')
        ax.title.set_size(24)
        plt.text(8.4, 3200, '{} entries'.format(df_2.sum(axis=1)['counts']), backgroundcolor='white', fontsize=16)
        plt.savefig('{}/z-class-distribution-all-129gCd.pdf'.format(folder_path))
        plt.show()



if __name__ == '__main__':
    im = Import_multi_cross_checks()
    # df = im.load('/Volumes/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/PI-ICR/129Cd/Z1-2')
    # print(df[df.isotope=='129Cd-ground-ov'])
    # im.load_all_spot_positions('/Volumes/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/PI-ICR/129Cd/Z1-1')
    # im.load_plot_ISOLTRAP_temp('/Volumes/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/Info/temperature/temperature-26-06--03-07.csv')
    im.read_counts_info_dict('/Volumes/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/PI-ICR/129Cd/Z2-2-before-tof-cut')
