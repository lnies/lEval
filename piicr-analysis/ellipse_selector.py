import matplotlib as mpl
# mpl.use('GTKAgg')
# mpl.use('Qt5Agg')
mpl.use('Qt5Agg')
from matplotlib.widgets import  EllipseSelector
from pylab import *
import matplotlib.pylab as plt
from matplotlib.patches import Ellipse
import json
import re
import time



class Ellipse_selector():
    def __init__(self):
        self.cut_list = []
        self.import_list = []
        mpl.rc('font', family='serif', serif='Utopia')      # Utopia LaTeX font!!
        mpl.rc('text', usetex=False)
        self.fig = figure((int(time.time()*100)), figsize=(8,8))
        self.ax = subplot(111)
        self.toggle_selector_ES = EllipseSelector(self.ax, self.onselect, drawtype='box')
        connect('key_press_event', self.toggle_selector)
        self.axis_scaling_factor = 0.031746032
        self.event_counter = 0
        self.analysis_dict = {}
        self.file_name = ''
        self.pattern = ''
        self.isotope = ''


    def main(self, import_list, file_name, freq_info):
        '''main function to execute'''

        self.split_file_name(file_name)
        self.load_list(import_list)
        self.plot_list()
        self.create_save_dict(import_list, freq_info)
        return(self.analysis_dict)

    def create_save_dict(self, import_list, freq_info):
        '''function to save the original list, the cut list and some info in dict and json file'''

        cut_list_output = []
        for i in self.cut_list:
            cut_list_output.append([i[0],
                             i[1]/self.axis_scaling_factor,
                             i[2]/self.axis_scaling_factor,
                             i[3]])

        self.analysis_dict = {'file name' : self.file_name,
                              'isotope': self.isotope,
                              'pattern' : self.pattern,
                              'frequency info' : freq_info,
                              'data' : import_list,
                              'cut_data' : cut_list_output}

        json_file = json.dumps(self.analysis_dict)
        f = open('{}_{}.json'.format(self.file_name, 'analysis_dict'),'w')
        f.write(json_file)
        f.close()


    def split_file_name(self, file_name):
        '''extract information from file name'''
        self.file_name = file_name
        if '_p1' in self.file_name:
            self.pattern = 'p1'
        elif '_p2' in self.file_name:
            self.pattern = 'p2'
        else:
            self.pattern = 'c'
        self.isotope = self.file_name.split('_')[0]


    def plot_list(self):
        '''function to plot the original list'''

        self.ax.plot([i[1] for i in self.import_list],
                     [i[2] for i in self.import_list],
                     'o', color='#1f77b4', markersize=5)
        self.set_plot_properties()
        # plt.ion()
        # plt.show()
        plt.pause(6)
        plt.savefig('{}_mcp_cut_t{}.pdf'.format(self.file_name, time.strftime('%M%S', time.localtime())))
        plt.close()

    def set_plot_properties(self):
        '''function to set some nice plotting parameters'''

        plt.xlabel('x / mm', fontsize=26)
        plt.ylabel('y / mm', fontsize=26)
        plt.xticks([-18,-12,-6,0,6,12,18], fontsize=22)
        plt.yticks([-18,-12,-6,0,6,12,18], fontsize=22)

        plt.axis([-23.809524, 23.809524, -23.809524, 23.809524])
        ax2 = plt.gca()
        mcp = Ellipse(xy=(0, 0), width=45.71428608, height=45.71428608, edgecolor='k', fc='None', lw=2)
        ax2.add_patch(mcp)

        plt.grid()
        plt.axhline(y=0, linewidth=1, color = 'k')
        plt.axvline(x=0, linewidth=1, color='k')
        plt.tight_layout()


    def load_list(self, import_list):
        '''function to load spot_positions with conversion ch. --> mm'''

        for i in import_list:
            self.import_list.append([i[0],
                                     i[1]*self.axis_scaling_factor,
                                     i[2]*self.axis_scaling_factor,
                                     i[3]])


    def onselect(self, eclick, erelease):
        '''function to handle click events with ellipse spans and cutting down the list depending on the selection'''
        # eclick and erelease are matplotlib events at press and release

        center, r_x, r_y = self.get_ellipse([eclick.xdata, eclick.ydata], [erelease.xdata, erelease.ydata])

        for i in self.import_list:
            if self.check_point_within_ellipse([i[1], i[2]], center, r_x, r_y):
                already_taken_flag = False
                for j in self.cut_list:
                    if i == j:
                        already_taken_flag = True
                if already_taken_flag == False and eclick.button == 1:
                    self.cut_list.append(i)
                elif already_taken_flag == True and eclick.button == 3:
                    self.cut_list.remove(i)

        if self.event_counter == 0:
            del self.ax.lines[3]        # avoids over-plotting of data points
        else:
            del self.ax.lines[-1]
            del self.ax.lines[-1]
        self.event_counter += 1

        self.ax.plot([i[1] for i in self.import_list],
                     [i[2] for i in self.import_list],
                     'o', color='#1f77b4', label='ions not considered for fit', markersize=5)
        self.ax.plot([i[1] for i in self.cut_list],
                     [i[2] for i in self.cut_list],
                     'x r', label='ions considered for fit', markersize=5)
        legend = plt.legend(title=str(self.file_name), fontsize=14)
        plt.setp(legend.get_title(),fontsize=18)
        plt.text(-23, -23, '{}/{} entries'.format(len(self.cut_list), len(self.import_list)), backgroundcolor='white')
        print 'Ions considered for fit  ::  {}/{}'.format(len(self.cut_list), len(self.import_list))


    def toggle_selector(self, event):
        '''function to activate/deactivate the ellipse functionality'''
        print(' Key pressed.')
        if event.key in ['Q', 'q'] and self.toggle_selector_ES.active:
            print(' EllipseSelector deactivated.')
            self.toggle_selector_ES.set_active(False)
        if event.key in ['A', 'a'] and not self.toggle_selector_ES.active:
            print(' EllipseSelector activated.')
            self.toggle_selector_ES.set_active(True)


    def get_ellipse(self, start, end):
        '''function calculates center, semi-major and semi-minor axis of an ellipse spanned by a sqare box'''

        x1 = start[0]
        y1 = start[1]
        x2 = end[0]
        y2 = end[1]
        center = [(x1+x2)/2, (y1+y2)/2]   # center of ellipse
        r_x = (x2-x1)/2                   # semi-major/minor axis in x direction
        r_y = (y2-y1)/2                   # semi-major/minor axis in y direction
        return(center, r_x, r_y)


    def check_point_within_ellipse(self, point, center_ellipse, r_x, r_y):
        '''function check whether a point is within a ellipse'''

        x = point[0]
        y = point[1]
        h = center_ellipse[0]
        k = center_ellipse[1]
        if ( (((x-h)**2) / (r_x**2)) + (((y-k)**2) / (r_y**2)) ) <= 1:
            return True
        else:
            return False



if __name__ == '__main__':
    list_x = [[1.0,-130.0,-36.5,2402457.0],[3.0,151.0,186.0,2482932.0],[5.0,239.0,220.5,2486431.0],[8.0,186.5,16.5,2485881.0],[9.0,17.5,209.5,2474744.0],[25.0,179.5,222.0,2485512.0],[30.0,196.0,156.5,2469798.0],[31.0,258.5,137.0,2474085.0],[35.0,299.0,193.0,2479816.0],[36.0,267.0,200.0,2478193.0],[40.0,-113.0,219.0,2475943.0],[53.0,0.0,154.5,2496631.0],[58.0,186.0,-38.0,2490658.0],[67.0,188.5,135.5,2504323.0],[69.0,-117.5,141.0,2489333.0],[73.0,216.5,219.0,2473443.0],[84.0,183.5,138.0,2496904.0],[97.0,-94.0,159.5,2466104.0],[99.0,-80.5,39.0,2489934.0],[109.0,203.0,182.0,2490617.0],[110.0,209.0,136.0,2494844.0],[111.0,289.5,208.5,2486257.0],[114.0,170.0,186.0,2489714.0],[118.0,325.0,103.0,2500314.0],[120.0,-95.0,139.5,2496229.0],[121.0,-85.0,166.0,2489965.0],[123.0,-131.5,213.5,2479931.0],[130.0,292.0,232.5,2491278.0],[131.0,211.5,218.5,2475633.0],[132.0,-153.0,175.0,2475701.0],[136.0,164.0,200.5,2490194.0],[139.0,144.5,125.5,2488658.0],[149.0,-153.0,189.0,2472362.0],[150.0,-151.5,179.5,2496810.0],[151.0,264.0,140.0,2471638.0],[156.0,189.5,192.5,2490653.0],[157.0,267.5,189.0,2481933.0],[226.0,218.0,200.5,2483474.0]]
    es = Ellipse_selector()
    a, temp_fig = es.main(list_x, '129Cd', 'p2', [])
    print a['cut_data']
