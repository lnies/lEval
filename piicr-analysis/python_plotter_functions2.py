# ---------------------------------------------------------------------------
# Written by Jonas Karthein in 2016/2017. Questions to jonas.karthein@cern.ch
# ---------------------------------------------------------------------------
# --Andree Welker added x^6, x^8 fit 23.01.2018

import matplotlib as mpl
mpl.use('Qt5Agg')

import matplotlib.pyplot as plt
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
import os


def python_plot(x_y_yerr, title, file_name, show_plot, legend_text, fit, fit_function, font, color_fit, color_points, fit_range, x_fit_min, x_fit_max, fit_b_parameter_start, plot_type, bins, markersize, grid_on_off, png_on_off, chi_sq_in_plot, diff_to_AME, diff_to_AME_err):
    """
    Function to plot (+ fit) a dataset.

    The sheet should consist of a header ('name / unit') and at least data
    for x and y. Possible formats are:
    - x,y
    - x,y,yerr

    Possible fit functions are:
    'linear', 'gauss', 'x^2', 'x^4', 'x^6', 'exp', 'sin'
    """
    timestamp_plot = (int(time.time()*100))
    xss = []
    xss_str = []
    yss = []
    yerr = []
    yss2 = []
    color = {'green': 'g', 'blue': 'b', 'yellow': 'y', 'black': 'k', 'red': 'r', 'cyan': 'c'}
    for j in range(1, len(x_y_yerr), 1):
        if type(x_y_yerr[1][0]) == str:
            if x_y_yerr[1][0].isdigit() == False:
                xss = [x for x in range(len(x_y_yerr)-1)]
                xss_str.append(str(x_y_yerr[j][0]))
        else:
            xss.append(x_y_yerr[j][0])

        # xss.append(datetime.datetime.strptime(str(x_y_yerr[j][0]), '%m/%d/%Y %H:%M:%S'))
        if plot_type == 'scatter' or plot_type == '2Dhistogram' or plot_type == '2dhistogram-mcp' or plot_type == 'polar' or plot_type == '2lines-no-error' or plot_type == 'mcp':
            yss.append(x_y_yerr[j][1])
        if len(x_y_yerr[0]) == 3 and not plot_type == '2lines-no-error':
            yerr.append(x_y_yerr[j][2])
        elif len(x_y_yerr[0]) == 3 and plot_type == '2lines-no-error':
            yss2.append(x_y_yerr[j][2])

    if plot_type == '2Dhistogram' or plot_type == 'mcp' or plot_type == '2dhistogram-mcp':
        plt.figure(timestamp_plot, figsize=(8, 8), facecolor='white')
    else:
        if type(x_y_yerr[1][0]) == str:
            if x_y_yerr[1][0].isdigit() == False:
                plt.figure(timestamp_plot, figsize=(18, 8), facecolor='white')
        else:
            plt.figure(timestamp_plot, figsize=(9, 6), facecolor='white')

    if font == 'Utopia':
        mpl.rc('font', family='serif', serif='Utopia')      # Utopia LaTeX font!!
        mpl.rc('text', usetex=False)
    elif font == 'LaTeX':
        mpl.rc('font', family='serif')      # Standard LaTeX font!!
        mpl.rc('text', usetex=True)

    if plot_type == 'scatter':
        range_x = 2 - 0 #max(xss) - min(xss) #usually use these settings
        range_y = 10 - 0 #max(yss) - min(yss)
        if yerr == []:
            plt.plot(xss, yss, '%s o' % color[color_points], label='%s' % legend_text, markersize=markersize)
        else:
            plt.errorbar(xss, yss, yerr=yerr, fmt='%s o' % color[color_points], label='%s' % legend_text, markersize=markersize)  # , range=[43, 58]
        plt.ylabel(r'$\frac{\nu_- - \nu_{-,ref}}{\nu_{-,ref}}$', fontsize=22)        #'%s' % x_y_yerr[0][1], fontsize=22    was delete
        plt.xlabel(r'$\rho_- [mm]$', fontsize=20)
        if type(xss[0]) == datetime.datetime:
            range_time = (max(xss)-min(xss)).total_seconds()
            if yerr == []:
                plt.axis([((min(xss) - datetime.timedelta(seconds=0.05*range_time))), ((max(xss) + datetime.timedelta(seconds=0.05*range_time))), min(yss) - 0.1 * range_y, max(yss) + 0.1 * range_y])
            else:
                plt.axis([((min(xss) - datetime.timedelta(seconds=0.05*range_time))), ((max(xss) + datetime.timedelta(seconds=0.05*range_time))), min(yss) - max(yerr) - 0.1 * range_y, max(yss) + max(yerr) + 0.1 * range_y])
        else:
            if yerr == []:
                plt.axis([min(xss) - 0.01 * range_x, max(xss) + 0.01 * range_x, min(yss) - 0.0001 * range_y, max(yss) + 0.0001 * range_y])
            else:
                plt.axis([min(xss) - 0.01 * range_x, max(xss) + 0.01 * range_x, min(yss) - max(yerr) +0.0000035, max(yss) + max(yerr) -0.000001])
                plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        if type(x_y_yerr[1][0]) == str:
            if x_y_yerr[1][0].isdigit() == False:
                plt.xticks(range(len(xss_str)), [str(i) for i in xss_str], rotation=80, size='small')
                plt.axis([min(xss) - 1, max(xss) + 1, min(yss) - max(yerr) - 0.1 * range_y, max(yss) + max(yerr) + 0.1 * range_y])
    elif plot_type == 'histogram':
        hist_x, xedges, patches = plt.hist(xss, bins=bins, histtype='step', stacked=True, fill=False)
        plt.ylabel('#', fontsize=22)
        plt.xlabel('%s' % x_y_yerr[0][0], fontsize=22)
        if type(xss[1]) == 'float' or type(xss[1]) == 'int':
            range_x = max(xedges) - min(xedges)
            range_y = max(hist_x) - min(hist_x)
            plt.axis([min(xedges) - 0.05 * range_x, max(xedges) + 0.05 * range_x, 0, max(hist_x) + 0.05 * range_y])
    elif plot_type == '2Dhistogram' or plot_type == '2dhistogram-mcp':
        colors = [(1,1,1), (0.10196, 0.42745, 0), (0.1294, 0.50588,0), (0.22745, 0.8039, 0.1843), (0.87058, 1, 0), (0.9882, 0.996, 0.1961), (0.9686, 0.8, 0.19215), (0.9529, 0.5922, 0.0118), (0.9451, 0.3882, 0.0157), (0.9333, 0.0314, 0.0157), (0.6078, 0.0118, 0), (0.1882, 0, 0)]  # easy to see background
        # colors = [(1,1,1), (0.76078, 0.9254, 0.7254), (0.10196, 0.42745, 0), (0.1294, 0.50588,0), (0.22745, 0.8039, 0.1843), (0.87058, 1, 0), (0.9882, 0.996, 0.1961), (0.9686, 0.8, 0.19215), (0.9529, 0.5922, 0.0118), (0.9451, 0.3882, 0.0157), (0.9333, 0.0314, 0.0157), (0.6078, 0.0118, 0), (0.1882, 0, 0)]   # low background
        n_bins = 1000
        cmap_name = 'my_name'

        cm = mpl.colors.LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins, gamma=1.5)
# no normalization:
        # counts, xedges, yedges, Image = plt.hist2d(xss, yss, bins=[bins/2.5, bins], cmap=cm)
# no normalization but maximum number of counts per bin:
        # counts, xedges, yedges, Image = plt.hist2d(xss, yss, bins=[bins/2.5, bins], cmap=cm, cmax=15)
# logarithmic normalization:
        counts, xedges, yedges, Image = plt.hist2d(xss, yss, bins=[bins/2.5, bins], cmap=cm, norm=mpl.colors.LogNorm(vmin=1, vmax=100))
        plt.ylabel('%s' % x_y_yerr[0][1], fontsize=22)
        plt.xlabel('%s' % x_y_yerr[0][0], fontsize=22)
        range_x = max(xedges) - min(xedges)
        range_y = max(yedges) - min(yedges)
        plt.axis([min(xedges) - 0.05 * range_x, max(xedges) + 0.05 * range_x, min(yedges) - 0.05 * range_x, max(yedges) + 0.05 * range_x])
    elif plot_type == 'polar':
        plt.polar(xss, yss, 'g o', alpha=.15)
        plt.xlabel('%s' % x_y_yerr[0][0], fontsize=22)
        plt.gca().set_rlim(0, 20)
    elif plot_type == '2lines-no-error':
        range_x = max(xss) - min(xss)
        range_y = max(yss) - min(yss2)
        plt.gca().set_xticklabels([])
        if yerr == []:
            plt.figure(timestamp_plot)
            plt.plot(xss, yss, '%s' % color[color_points[0]], label='%s' % legend_text[0], markersize=markersize)
            plt.figure(timestamp_plot)
            plt.plot(xss, yss2, '%s' % color[color_points[1]], label='%s' % legend_text[1], markersize=markersize)
        plt.ylabel('%s' % x_y_yerr[0][1], fontsize=22)
        plt.xlabel('%s' % x_y_yerr[0][0], fontsize=22)
        if type(xss[0]) == datetime.datetime:
            range_time = (sum(int(x) * 60 ** i for i,x in enumerate(reversed(str(max(xss)-min(xss)).split(":")))))
            plt.axis([((min(xss) - datetime.timedelta(seconds=0.05*range_time))), ((max(xss) + datetime.timedelta(seconds=0.05*range_time))), min([min(yss), min(yss2)]) - 0.1 * range_y, max([max(yss), max(yss2)]) + 0.1 * range_y])
        else:
            plt.axis([min(xss) - 0.05 * range_x, max(xss) + 0.05 * range_x, min([min(yss), min(yss2)]) - 0.1 * range_y, max([max(yss), max(yss2)]) + 0.1 * range_y])
        plt.legend(fontsize=16)
    if plot_type == 'mcp':
        if markersize > 2:
            plt.plot(xss, yss, '%s o' % color[color_points], label='%s' % legend_text, markersize=markersize, alpha=0.08)
        else:
            plt.plot(xss, yss, '%s o' % color[color_points], label='%s' % legend_text, markersize=markersize)
    if (plot_type == 'mcp' or plot_type == '2dhistogram-mcp') and diff_to_AME_err != -99999:
        plt.ylabel('y / mm', fontsize=26)
        plt.xlabel('x / mm', fontsize=26)
        plt.xticks([-566.999995464,-377.999996976,-188.999998488,0,188.999998488,377.999996976,566.999995464], [u'\u221218', u'\u221212', u'\u22126',0,6,12,18], fontsize=22)
        plt.yticks([-566.999995464,-377.999996976,-188.999998488,0,188.999998488,377.999996976,566.999995464], [u'\u221218', u'\u221212', u'\u22126',0,6,12,18], fontsize=22)
        plt.axhline(y=0, linewidth=1, color = 'k')
        plt.axvline(x=0, linewidth=1, color='k')
        plt.text(-724.499994204, -724.499994204, '{} entries'.format(len(xss)))

        plt.axis([-750, 750, -750, 750])
        ax = plt.gca()
        mcp = Ellipse(xy=(0, 0), width=1440, height=1440, edgecolor='k', fc='None', lw=2)
        ax.add_patch(mcp)
        if diff_to_AME_err != -99999:
            if len(diff_to_AME_err) == 4:
                fwhm = Ellipse(xy=(diff_to_AME_err[0], diff_to_AME_err[2]), width=diff_to_AME[0]*2, height=diff_to_AME[2]*2, edgecolor='r', fc='None', lw=2)
                ax.add_patch(fwhm)
                anchored_text = AnchoredText('Fit X pos. = %3.2f (%3.2f)\nFit Y pos. = %3.2f (%3.2f)\nFit X FWHM = %3.2f (%3.2f)\nFit Y FWHM = %3.2f (%3.2f)' %(float(diff_to_AME_err[0]), float(diff_to_AME_err[1]), float(diff_to_AME_err[2]), float(diff_to_AME_err[3]), float(diff_to_AME[0]), float(diff_to_AME[1]), float(diff_to_AME[2]), float(diff_to_AME[3])), loc=1)
                ax.add_artist(anchored_text)
            elif len(diff_to_AME_err) > 4:
                fwhm_1 = Ellipse(xy=(diff_to_AME_err[0], diff_to_AME_err[2]), width=diff_to_AME[0]*2, height=diff_to_AME[2]*2, edgecolor='r', fc='None', lw=2)
                ax.add_patch(fwhm_1)
                fwhm_2 = Ellipse(xy=(diff_to_AME_err[4], diff_to_AME_err[6]), width=diff_to_AME[0]*2, height=diff_to_AME[2]*2, edgecolor='r', fc='None', lw=2)
                ax.add_patch(fwhm_2)
                anchored_text = AnchoredText('Fit X pos. dom. = %3.2f (%3.2f)\nFit Y pos. dom. = %3.2f (%3.2f)\nFit X pos. rec. = %3.2f (%3.2f)\nFit Y rec. dom. = %3.2f (%3.2f)\nFit X FWHM = %3.2f (%3.2f)\nFit Y FWHM = %3.2f (%3.2f)\nStatus: %s' %(float(diff_to_AME_err[0]), float(diff_to_AME_err[1]), float(diff_to_AME_err[2]), float(diff_to_AME_err[3]), float(diff_to_AME_err[4]), float(diff_to_AME_err[5]), float(diff_to_AME_err[6]), float(diff_to_AME_err[7]),float(diff_to_AME[0]), float(diff_to_AME[1]), float(diff_to_AME[2]), float(diff_to_AME[3]), diff_to_AME_err[-1]), loc=1)
                ax.add_artist(anchored_text)

    if plot_type == '2dhistogram-mcp' and diff_to_AME_err == -99999:
        plt.ylabel('%s' % x_y_yerr[0][1], fontsize=22)
        plt.xlabel('%s' % x_y_yerr[0][0], fontsize=22)
        plt.xticks([-18,-12,-6,0,6,12,18], fontsize=22)
        plt.yticks([-18,-12,-6,0,6,12,18], fontsize=22)

        plt.axis([-23.81, 23.81, -23.81, 23.81])
        ax = plt.gca()
        mcp = Ellipse(xy=(0, 0), width=45.714, height=45.714, edgecolor='k', fc='None', lw=2)
        ax.add_patch(mcp)


    if title == '':
        pass
    elif plot_type == 'polar':
        plt.title('%s' % str(title), fontsize=26, y=1.03)
    else:
        plt.title('%s' % str(title), fontsize=26)
        if ('dominant' in title or 'recessive' in title) and 'vector' not in title:
            plt.text(0, min(yss)-max(yerr), '{}'.format(os.getcwd().split(os.sep)[-2]))    # get folder name of run



    if fit == 'yes':
        x = mle.var('x', observed=True, vector=True)
        y = mle.var('y', observed=True, vector=True)

        a = mle.var('a')
        b = mle.var('b')
        c = mle.var('c')
        d = mle.var('d')
        e = mle.var('e')
        sigma = mle.var('sigma')


        if fit_function == 'gauss':
            model = mle.Normal(y, a * numpy.exp(-(x - b)**2.0 / (2 * c**2)), sigma)
            def func(x, a, b, c):
                return a * numpy.exp(-(x - b)**2.0 / (2 * c**2))
        elif fit_function == 'x^2':
            model = mle.Normal(y, a * (x**2.0) + b * x + c, sigma)
            def func(x, a, b, c):
                return a * (x**2.0) + b * x + c
        elif fit_function == 'x^4':
            model = mle.Normal(y, a * (x**2.0) + b * (x**4.0) + c, sigma)
            def func(x, a, b, c):
                return a * (x**2.0) + b * (x**4.0) + c
        elif fit_function == 'x^6':
            model = mle.Normal(y, a * (x**2.0) + b * (x**4.0) + c * (x**6.0) + d, sigma)
            def func(x, a, b, c, d):
                return a * (x**2.0) + b * (x**4.0) + c * (x**6.0) + d
        elif fit_function == 'x^8':
            model = mle.Normal(y, a * (x**2.0) + b * (x**4.0) + c * (x**6.0) + d * (x**8.0) + e, sigma)
            def func(x, a, b, c, d, e):
                return a * (x**2.0) + b * (x**4.0) + c * (x**6.0) + d * (x**8.0) + e
        elif fit_function == 'linear':
            model = mle.Normal(y, a * x + b + c, sigma)
            def func(x, a, b, c):
                return a * x + b + c
        elif fit_function == 'sin':
            model = mle.Normal(y, a * numpy.sin(b*x) + c, sigma)
            def func(x, a, b, c):
                return a * numpy.sin(b*x) + c
        elif fit_function == 'exp':
            model = mle.Normal(y, a * numpy.exp(b * x) + c, sigma)
            def func(x, a, b, c):
                return a * numpy.exp(b * x) + c
        elif fit_function == 'log':
            model = mle.Normal(y, a * numpy.log(b * x) + c, sigma)
            def func(x, a, b, c):
                return a * numpy.log(b * x) + c

        if fit_range == 'full':
            if plot_type == 'scatter':
                xs = numpy.array(xss)
                ys = numpy.array(yss)
            elif plot_type == 'histogram':
                xs_hilf = numpy.array(xedges)
                xs = xs_hilf[:-1]
                ys = numpy.array(hist_x)
        elif fit_range == 'partly':
            if plot_type == 'scatter':
                xs_list = []
                ys_list = []
                xs_hilf = numpy.array(xss)
                ys_hilf = numpy.array(yss)
                for kko in range(len(xss)):
                    if xs_hilf[kko] >= x_fit_min and xs_hilf[kko] <= x_fit_max:
                        xs_list.append(xs_hilf[kko])
                        ys_list.append(ys_hilf[kko])
                xs = numpy.array(xs_list)
                ys = numpy.array(ys_list)
            elif plot_type == 'histogram':
                xs = numpy.array(xedges[int(x_fit_min*bins):int(x_fit_max*bins)])
                xs = xs_hilf[:-1]
                ys = numpy.array(hist_x)
        result_mlf = model.fit({'x': xs, 'y': ys}, {'a': 1, 'b': 1, 'c': 1, 'd': 1, 'e': 1, 'sigma': 1})         # FIT Startbedingungen
        parameter_mlf = []
        parameter_mlf.append(result_mlf.x['a'])
        parameter_mlf.append(result_mlf.x['b'])
        parameter_mlf.append(result_mlf.x['c'])
        #parameter_mlf.append(result_mlf.x['d'])
        #parameter_mlf.append(result_mlf.x['e'])    #For the plotter it is somehow an issue if a letter is not returned.
        #parameter_mlf.append(result_mlf.x['f'])
        parameter_mlf_err = []
        parameter_mlf_err = numpy.sqrt(numpy.diag(result_mlf['hess_inv']))      # one standard deviation error of parameters
        print parameter_mlf, parameter_mlf_err
        # print parameter_mlf, parameter_mlf_err

        if plot_type == 'scatter':
            hilfx = numpy.linspace(min(xss), max(xss), 10000)
        elif plot_type == 'histogram':
            hilfx = numpy.linspace(min(xedges[:bins]), max(xedges[:bins]), 1000)

        plt.figure(timestamp_plot, figsize=(8, 8), facecolor='white')
        plt.plot(hilfx, func(hilfx, *parameter_mlf), '%s-' % color[color_fit], linewidth=2, label='Max. likelihood fit')
        plt.tick_params(labelsize=16)

        if legend_text == '' or plot_type == '2lines-no-error':
            pass
        else:
            plt.legend(fontsize=16, numpoints=1)

        if fit_function == 'gauss':
            plt.text(min(xss), max(yss) - 0.2 * range_y, 'Max. likelihood estimation: \na = %3.2f (%3.2f) \nb = %3.2f (%3.2f) \nc = %3.2f (%3.2f)' % (parameter_mlf[0], parameter_mlf_err[1], parameter_mlf[1], parameter_mlf_err[2], parameter_mlf[2], parameter_mlf_err[3]), fontsize=10)
            plt.text(min(xss), max(yss) - 0.05 * range_y, 'Fit-Funktion: \n$a*exp\\left(\\frac{-(x - b)^2}{2 * c^2}\\right)$', fontsize=10)
        elif fit_function == 'x^2':
            plt.text(min(xss), max(yss) - 0.3 * range_y, 'Max. likelihood estimation: \na = %3.2f (%3.2f) \nb = %3.2f (%3.2f) \nc = %3.2f (%3.2f) \nextremum = %3.2f' % (parameter_mlf[0], parameter_mlf_err[1], parameter_mlf[1], parameter_mlf_err[2], parameter_mlf[2], parameter_mlf_err[3], (-parameter_mlf[1]/2/parameter_mlf[0])), fontsize=10)
            plt.text(min(xss), max(yss) - 0.05 * range_y, 'Fit-Funktion: \n$a*x^2 + b*x + c$', fontsize=10)
        elif fit_function == 'x^4':
            plt.text(min(xss), max(yss) - 0.3 * range_y, 'Max. likelihood estimation: \na = %3.8f (%3.8f) \nb = %3.8f (%3.8f) \nc = %3.8f (%3.8f)' % (parameter_mlf[0], parameter_mlf_err[1], parameter_mlf[1], parameter_mlf_err[2], parameter_mlf[2], parameter_mlf_err[3]), fontsize=10)
            plt.text(min(xss), max(yss) - 0.05 * range_y, 'Fit-Funktion: \n$a*x^2 + b*x^4 + c$', fontsize=10)
        elif fit_function == 'x^6':
            plt.text(1, 0 - 0.3 * range_y, 'Max. likelihood estimation: \na = %3.8f (%3.8f) \nb = %3.8f (%3.8f) \nc = %3.8f (%3.8f) \nd = %3.8f (%3.8f) ' % (parameter_mlf[0], parameter_mlf_err[1], parameter_mlf[1], parameter_mlf_err[2], parameter_mlf[2], parameter_mlf_err[3], parameter_mlf[3], parameter_mlf_err[4]), fontsize=10)
            plt.text(min(xss), max(yss) - 0.05 * range_y, 'Fit-Funktion: \n$a*x^2 + b*x^4 + c*x^6 + d$', fontsize=10)
        elif fit_function == 'x^8':
            plt.text(1, 0 - 0.3 * range_y, 'Max. likelihood estimation: \na = %3.8f (%3.8f) \nb = %3.8f (%3.8f) \nc = %3.8f (%3.8f) \nd = %3.8f (%3.8f) \ne = %3.8f (%3.8f) ' % (parameter_mlf[0], parameter_mlf_err[1], parameter_mlf[1], parameter_mlf_err[2], parameter_mlf[2], parameter_mlf_err[3], parameter_mlf[3], parameter_mlf_err[4], parameter_mlf[4], parameter_mlf_err[5]), fontsize=10)
            plt.text(min(xss), max(yss) - 0.05 * range_y, 'Fit-Funktion: \n$a*x^2 + b*x^4 + c*x^6 + d*x^8 + e$', fontsize=10)
        elif fit_function == 'linear':
            plt.text(min(xss), max(yss) - 0.2 * range_y, 'Max. likelihood estimation: \na = %3.2f (%3.2f) \nb = %3.2f (%3.2f)' % (parameter_mlf[0], parameter_mlf_err[1], parameter_mlf[1] + parameter_mlf[2], parameter_mlf_err[2] + parameter_mlf_err[3]), fontsize=10)
            plt.text(min(xss), max(yss) - 0.05 * range_y, 'Fit-Funktion: \n$a*x+b$', fontsize=10)
        elif fit_function == 'sin':
            plt.text(min(xss), max(yss) - 0.2 * range_y, 'Max. likelihood estimation: \na = %3.2f (%3.2f) \nb = %3.2f (%3.2f) \nc = %3.2f (%3.2f)' % (parameter_mlf[0], parameter_mlf_err[1], parameter_mlf[1], parameter_mlf_err[2], parameter_mlf[2], parameter_mlf_err[3]), fontsize=10)
            plt.text(min(xss), max(yss) - 0.05 * range_y, 'Fit-Funktion: \n$a*sin\\left(b*x\\right) + c$', fontsize=10)
        elif fit_function == 'exp':
            plt.text(min(xss), max(yss) - 0.2 * range_y, 'Max. likelihood estimation: \na = %3.2f (%3.2f) \nb = %3.2f (%3.2f) \nc = %3.2f (%3.2f)' % (parameter_mlf[0], parameter_mlf_err[1], parameter_mlf[1], parameter_mlf_err[2], parameter_mlf[2], parameter_mlf_err[3]), fontsize=10)
            plt.text(min(xss), max(yss) - 0.05 * range_y, 'Fit-Funktion: \n$a*exp\\left(b*x\\right) + c$', fontsize=10)
        elif fit_function == 'log':
            plt.text(min(xss), max(yss) - 0.2 * range_y, 'Max. likelihood estimation: \na = %3.2f (%3.2f) \nb = %3.2f (%3.2f) \nc = %3.2f (%3.2f)' % (parameter_mlf[0], parameter_mlf_err[1], parameter_mlf[1], parameter_mlf_err[2], parameter_mlf[2], parameter_mlf_err[3]), fontsize=10)
            plt.text(min(xss), max(yss) - 0.05 * range_y, 'Fit-Funktion: \n$a*log\\left(b*x\\right) + c$', fontsize=10)

    if grid_on_off == 'on':
        plt.grid(True)
        # gridlines = ax.get_xgridlines() + ax.get_ygridlines()
        # for line in gridlines:
        #     line.set_linestyle('-.')
    # print chi_sq_in_plot
    if not chi_sq_in_plot == '' and not plot_type == 'mcp' and not plot_type == '2dhistogram-mcp':
        # print 'fucking here'
        # print chi_sq_in_plot[2]
        plt.tick_params(labelsize=16)
        plt.axhline(y=0, color='k', linestyle='-', linewidth=2)

        plt.axhline(y=diff_to_AME, color='r', linestyle='-')
        plt.axhline(y=diff_to_AME+diff_to_AME_err, color='r', linestyle='--', alpha=0.75)
        plt.axhline(y=diff_to_AME-diff_to_AME_err, color='r', linestyle='--', alpha=0.75)
        if chi_sq_in_plot[2] >= 0:
            plt.axhline(y=diff_to_AME-diff_to_AME_err*math.sqrt(chi_sq_in_plot[0]), color='r', linestyle='-.', alpha=0.65)
            plt.axhline(y=diff_to_AME+diff_to_AME_err*math.sqrt(chi_sq_in_plot[0]), color='r', linestyle='-.', alpha=0.65)
            ax = P.gca()
            anchored_text = AnchoredText('$\Delta_{AME}$ = %3.2f ($\pm$%3.2f) keV\nAvg. $\chi^2_{red}$ = %3.2f\nAvg. unc. = %3.2f keV\nt$_{acc}$ = %3.2f ms\nMass dep. shift = %s\nSys. unc. = %s' %(diff_to_AME, diff_to_AME_err, chi_sq_in_plot[0], numpy.mean(yerr), chi_sq_in_plot[1], chi_sq_in_plot[2], chi_sq_in_plot[3]), loc=1)
        else:
            ax = P.gca()
            # plt.axhline(y=diff_to_AME-diff_to_AME_err*math.sqrt(chi_sq_in_plot[0]), color='r', linestyle='-.', alpha=0.65)
            # plt.axhline(y=diff_to_AME+diff_to_AME_err*math.sqrt(chi_sq_in_plot[0]), color='r', linestyle='-.', alpha=0.65)
            anchored_text = AnchoredText('$\Delta_{AME}$ = %3.2f ($\pm$%3.2f) keV\nAvg. $\chi^2_{red}$ = %3.2f \n%s' %(diff_to_AME, diff_to_AME_err, chi_sq_in_plot[0], chi_sq_in_plot[1]), loc=1)
            # print 'or here?!?'
        anchored_text.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        anchored_text.patch.set_alpha(0.8)
        ax.add_artist(anchored_text)

        if type(xss[0]) == datetime.datetime:
            range_time = (max(xss)-min(xss)).total_seconds()
            interval_time = [((min(xss) - datetime.timedelta(seconds=0.05*range_time))), ((max(xss) + datetime.timedelta(seconds=0.05*range_time)))]
            plt.fill_between(interval_time, y1=diff_to_AME+diff_to_AME_err, y2=diff_to_AME-diff_to_AME_err, color='r', alpha=0.1)
            plt.fill_between(interval_time, y1=diff_to_AME+diff_to_AME_err*math.sqrt(chi_sq_in_plot[0]), y2=diff_to_AME-diff_to_AME_err*math.sqrt(chi_sq_in_plot[0]), color='r', alpha=0.1)
        else:
            plt.fill_between([min(xss) - 0.05 * range_x, max(xss) + 0.05 * range_x], y1=diff_to_AME+diff_to_AME_err, y2=diff_to_AME-diff_to_AME_err, color='r', alpha=0.1)
            plt.fill_between([min(xss) - 0.05 * range_x, max(xss) + 0.05 * range_x], y1=diff_to_AME+diff_to_AME_err*math.sqrt(chi_sq_in_plot[0]), y2=diff_to_AME-diff_to_AME_err*math.sqrt(chi_sq_in_plot[0]), color='r', alpha=0.1)

    plt.tight_layout()

    plt.savefig('%s.pdf' % file_name)
    if png_on_off == 'on':
        plt.savefig('%s.PNG' % file_name)

    if show_plot == 'yes':
        plt.show()
    elif show_plot == 'no':
        plt.close()

    if fit == 'yes':
        return(xss, yss, yerr, parameter_mlf, parameter_mlf_err)
    elif plot_type == 'histogram':
        return(xss, yss, yerr, hist_x, xedges)
    else:
        return(xss, yss, yerr)
