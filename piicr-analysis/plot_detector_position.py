import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy
from mpl_toolkits.mplot3d import axes3d
import mle
import pylab as P
import math
import datetime
from matplotlib.patches import Ellipse
from matplotlib.offsetbox import AnchoredText
from matplotlib import rcParams, gridspec
import time
import sys
from read_write_functions import *
from cross_checks_functions import *


def plot_detector_position(file_path, xlsx_file_name):
    axis_scaling_factor = 0.031746032       # convert a.u. to mm

    x_xerr_y_yerr = read_excel(file_path, xlsx_file_name)
    xss = []
    xerr = []
    yss = []
    yerr = []

    for j in range(1, len(x_xerr_y_yerr), 1):
        xss.append(float(x_xerr_y_yerr[j][1]))
        xerr.append(float(x_xerr_y_yerr[j][2]))
        yss.append(float(x_xerr_y_yerr[j][3]))
        yerr.append(float(x_xerr_y_yerr[j][4]))

################## Plot

    mpl.rc('font', family='serif', serif='Utopia')      # Utopia LaTeX font!!
    mpl.rc('text', usetex=False)

    plt.figure(1, figsize=(9, 9), facecolor='white')
    plt.errorbar(xss, yss, xerr=xerr, yerr=yerr, fmt='g o', markersize=5)

    plt.xlabel('x / mm', fontsize=22)
    plt.ylabel('y / mm', fontsize=22)

    plt.axis([-6, 6, -6, 6])
    ax = plt.gca()
    # mcp = Ellipse(xy=(0, 0), width=1440, height=1440, edgecolor='k', fc='None', lw=2)
    # ax.add_patch(mcp)


    # center_x = xss[-9]
    # center_y = yss[-9]
    # center_x = numpy.mean(xss[:-10])
    # center_y = numpy.mean(yss[:-10])

    info = []
    for i in range(0, len(xss), 10):
        all_r = [['Radius / mm', 'Unc. / mm']]
        center_x = numpy.mean(xss[i:i+10])
        center_y = numpy.mean(yss[i:i+10])
        for j in range(10):
            all_r.append([np.sqrt((xss[i+j]-center_x)**2+(yss[i+j]-center_y)**2), np.sqrt(((xss[i+j]-center_x)*xerr[i+j])**2 + ((yss[i+j]-center_y)*yerr[i+j])**2) / np.sqrt((xss[i+j]-center_x)**2+(yss[i+j]-center_y)**2)])
        weighted_avg, inner_unc_diff, outer_unc_diff, total_unc, birge_ratio = weighted_mean(all_r)
        info.append([weighted_avg, total_unc, center_x, center_y])
        avg_radius = Ellipse(xy=(center_x, center_y), width=weighted_avg*2, height=weighted_avg*2, edgecolor='r', fc='None', lw=1)
        ax.add_patch(avg_radius)

    print len(info)
    anchored_text = AnchoredText('W.avg. radius @ mag.amp=1.35V: %3.2f $\pm$ %3.2f mm (c_x = %3.2f mm, c_y = %3.2f mm)\nW.avg. radius @ mag.amp=1.85V: %3.2f $\pm$ %3.2f mm (c_x = %3.2f mm, c_y = %3.2f mm)\nW.avg. radius @ mag.amp=2.35V: %3.2f $\pm$ %3.2f mm (c_x = %3.2f mm, c_y = %3.2f mm)\nW.avg. radius @ mag.amp=2.85V: %3.2f $\pm$ %3.2f mm (c_x = %3.2f mm, c_y = %3.2f mm)\nW.avg. radius @ mag.amp=0.85V: %3.2f $\pm$ %3.2f mm (c_x = %3.2f mm, c_y = %3.2f mm)' % (info[0][0], info[0][1],info[0][2], info[0][3],info[1][0], info[1][1],info[1][2], info[1][3], info[2][0], info[2][1], info[2][2], info[2][3], info[3][0], info[3][1], info[3][2], info[3][3], info[4][0], info[4][1], info[4][2], info[4][3]), loc=1)
    # anchored_text = AnchoredText('W.avg. radius @ mag.amp=1.35V: %3.2f $\pm$ %3.2f mm\nW.avg. radius @ mag.amp=1.85V: %3.2f $\pm$ %3.2f mm\nW.avg. radius @ mag.amp=2.35V: %3.2f $\pm$ %3.2f mm\nW.avg. radius @ mag.amp=2.85V: %3.2f $\pm$ %3.2f mm' % (info[0][0], info[0][1],info[1][0], info[1][1], info[2][0], info[2][1], info[3][0], info[3][1]), loc=1)
    anchored_text.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    anchored_text.patch.set_alpha(0.8)
    ax.add_artist(anchored_text)


    plt.grid()
    plt.title('%s' % str('Ejection phase scan for $^{85}$Rb'), fontsize=26)
    plt.tight_layout()
    plt.savefig('%s.pdf' % xlsx_file_name)
    plt.savefig('%s.PNG' % xlsx_file_name)
    plt.show()


if __name__ == '__main__':
    file_path = '/Volumes/dfs/DATA/2017/Screenshots/2017-07-Xe-Cs-prep/'
    xlsx_file_name = 'piicr_detector_position_scan'
    plot_detector_position(file_path, xlsx_file_name)
