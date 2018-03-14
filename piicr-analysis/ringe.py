from read_write_functions import *
import os
from piicr_daq_functions import *
import glob
import math
import numpy as np
import matplotlib.pyplot as plt


def ringe(folder_path):
    '''bla'''
    csv_files = []
    os.chdir(folder_path)
    for file in glob.glob("*polar.csv"):        # searches for all ...polar.csv files
        csv_files.append(file)
    for i in csv_files:
    	import_list = read_csv(folder_path, os.path.splitext(i)[0])

        fit_x_y = import_list.pop(1) # beinhaltet fit_x und fit_y, damit laesst sich ein winkel berechnen um den anfangspunkt um winkel und um 180 grad zu verschieben. probiere modulo
        angle_fit = math.atan2(float(fit_x_y[1]), float(fit_x_y[0]))
        spot_positions = []
        for j in range(1, len(import_list), 1):       # convert to float
            spot_positions.append([float(import_list[j][0]), float(import_list[j][1]), float(import_list[j][2]) + math.pi, float(import_list[j][3])])   # make angle from 0 to 2pi

        plot_r = []
        plot_phi = []
        for xx in range(len(spot_positions)):
            plot_r.append(spot_positions[xx][1])
            plot_phi.append(spot_positions[xx][2])
        ax = plt.subplot(111, projection='polar')
        ax.set_ylim(0, 20)
        c = ax.scatter(plot_phi, plot_r, s=0.5, alpha=0.15)
        plt.savefig('0polar_plot_'+i+'.pdf')
        plt.close()



        spot_positions_temp_1 = []
        spot_positions_temp_2 = []
        spot_positions_temp_3 = []
        for l in range(len(spot_positions)):
            if spot_positions[l][1] > 0 and spot_positions[l][1] < 0.33*2.5 and spot_positions[l][3] > 47/0.000025 and spot_positions[l][3] < 55/0.000025:
                spot_positions_temp_1.append(spot_positions[l]) 
            elif spot_positions[l][1] > 0.33*2.5 and spot_positions[l][1] < 0.66*2.5 and spot_positions[l][3] > 47/0.000025 and spot_positions[l][3] < 55/0.000025:
                spot_positions_temp_2.append(spot_positions[l]) 
            elif spot_positions[l][1] > 0.66*2.5 and spot_positions[l][1] < 1*2.5 and spot_positions[l][3] > 47/0.000025 and spot_positions[l][3] < 55/0.000025:
                spot_positions_temp_3.append(spot_positions[l]) 
        tof_plot(spot_positions_temp_1, 60, '1'+i, 1)
        tof_plot(spot_positions_temp_2, 60, '2'+i, 1)
        tof_plot(spot_positions_temp_3, 60, '3'+i, 1)

        circle_parts = 8
        intervals = [(x*2*math.pi/circle_parts + angle_fit)%(2*math.pi) for x in range(circle_parts+1)]     # divides the circle into param:circle_parts pieces but starts at 
        hist_2d_1 = []
        hist_2d_2 = []
        hist_2d_3 = []
        for m in range(len(intervals)-1):
            for n in range(len(spot_positions_temp_1)):
                if intervals[m] < intervals[m+1]:
                    if spot_positions_temp_1[n][2] > intervals[m] and spot_positions_temp_1[n][2] < intervals[m+1]:
                        hist_2d_1.append([((intervals[m]-intervals[m+1])/2+intervals[m]), spot_positions_temp_1[n][3]])
                else:
                    if (spot_positions_temp_1[n][2] > intervals[m] and spot_positions_temp_1[n][2] < 2*math.pi) or (spot_positions_temp_1[n][2] < intervals[m+1] and spot_positions_temp_1[n][2] > 0):
                        hist_2d_1.append([((intervals[m-1]-intervals[m])/2+intervals[m])%(2*math.pi), spot_positions_temp_1[n][3]])

            for n in range(len(spot_positions_temp_2)):
                if intervals[m] < intervals[m+1]:
                    if spot_positions_temp_2[n][2] > intervals[m] and spot_positions_temp_2[n][2] < intervals[m+1]:
                        hist_2d_2.append([((intervals[m]-intervals[m+1])/2+intervals[m]), spot_positions_temp_2[n][3]])
                else:
                    if (spot_positions_temp_2[n][2] > intervals[m] and spot_positions_temp_2[n][2] < 2*math.pi) or (spot_positions_temp_2[n][2] < intervals[m+1] and spot_positions_temp_2[n][2] > 0):
                        hist_2d_2.append([((intervals[m-1]-intervals[m])/2+intervals[m])%(2*math.pi), spot_positions_temp_2[n][3]])
            for n in range(len(spot_positions_temp_3)):
                if intervals[m] < intervals[m+1]:
                    if spot_positions_temp_3[n][2] > intervals[m] and spot_positions_temp_3[n][2] < intervals[m+1]:
                        hist_2d_3.append([((intervals[m]-intervals[m+1])/2+intervals[m]), spot_positions_temp_3[n][3]])
                else:
                    if (spot_positions_temp_3[n][2] > intervals[m] and spot_positions_temp_3[n][2] < 2*math.pi) or (spot_positions_temp_3[n][2] < intervals[m+1] and spot_positions_temp_3[n][2] > 0):
                        hist_2d_3.append([((intervals[m-1]-intervals[m])/2+intervals[m])%(2*math.pi), spot_positions_temp_3[n][3]])

        xs1 = []
        ys1 = []
        for o in range(len(hist_2d_1)):
            xs1.append(hist_2d_1[o][0])
            ys1.append(hist_2d_1[o][1])
        xs2 = []
        ys2 = []
        for p in range(len(hist_2d_2)):
            xs2.append(hist_2d_2[p][0])
            ys2.append(hist_2d_2[p][1])
        xs3 = []
        ys3 = []
        for q in range(len(hist_2d_3)):
            xs3.append(hist_2d_3[q][0])
            ys3.append(hist_2d_3[q][1])

        # H1, xedges1, yedges1 = np.histogram2d(xs1, ys1, bins=(60, 60))#, range=np.array([(-hist2d_min_max, hist2d_min_max), (-hist2d_min_max, hist2d_min_max)]))
        # H2, xedges2, yedges2 = np.histogram2d(xs2, ys2, bins=(60, 60))#, range=np.array([(-hist2d_min_max, hist2d_min_max), (-hist2d_min_max, hist2d_min_max)]))
        # H3, xedges3, yedges3 = np.histogram2d(xs3, ys3, bins=(60, 60))#, range=np.array([(-hist2d_min_max, hist2d_min_max), (-hist2d_min_max, hist2d_min_max)]))

        colors = [(1,1,1),(0.10196, 0.42745, 0), (0.1294, 0.50588,0), (0.22745, 0.8039, 0.1843), (0.87058, 1, 0), (0.9882, 0.996, 0.1961), (0.9686, 0.8, 0.19215), (0.9529, 0.5922, 0.0118), (0.9451, 0.3882, 0.0157), (0.9333, 0.0314, 0.0157), (0.6078, 0.0118, 0), (0.1882, 0, 0)]    # 4&5 (0.6078, 0.8078, 0.19215), (0.8, 0.99215, 0.3842), 
        n_bins = 1000
        cmap_name = 'my_name'
        cm = mpl.colors.LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins, gamma=1.5)

        hist1 = plt.figure(1, figsize=(9,9))
        plt.hist2d(xs1, ys1, bins=50, cmap=cm)
        plt.xlabel('Angle / radian', fontsize=22)
        plt.ylabel('Time of Flight / ns', fontsize=22)
        plt.axis([-0.5, 2*math.pi+0.5, 1880000, 2200000])
        hist1.savefig('1stacked_hist2d_'+i+'.pdf')
        plt.close()
        hist2 = plt.figure(2, figsize=(9,9))
        plt.hist2d(xs2, ys2, bins=50, cmap=cm)
        plt.xlabel('Angle / radian', fontsize=22)
        plt.ylabel('Time of Flight / ns', fontsize=22)
        plt.axis([-0.5, 2*math.pi+0.5, 1880000, 2200000])
        hist2.savefig('2stacked_hist2d_'+i+'.pdf')
        plt.close()
        hist3 = plt.figure(3, figsize=(9,9))
        plt.hist2d(xs2, ys2, bins=50, cmap=cm)
        plt.xlabel('Angle / radian', fontsize=22)
        plt.ylabel('Time of Flight / ns', fontsize=22)
        plt.axis([-0.5, 2*math.pi+0.5, 1880000, 2200000])
        hist3.savefig('3stacked_hist2d_'+i+'.pdf')
        plt.close()


        plot_r_1 = []
        plot_phi_1 = []
        for xx1 in range(len(spot_positions_temp_1)):
            plot_r_1.append(spot_positions_temp_1[xx1][1])
            plot_phi_1.append(spot_positions_temp_1[xx1][2])
        ax = plt.subplot(111, projection='polar')
        ax.set_ylim(0, 20)
        c = ax.scatter(plot_phi_1, plot_r_1, s=0.5, alpha=0.15)
        plt.savefig('1polar_plot_'+i+'.pdf')
        plt.close()

        plot_r_2 = []
        plot_phi_2 = []
        for xx2 in range(len(spot_positions_temp_2)):
            plot_r_2.append(spot_positions_temp_2[xx2][1])
            plot_phi_2.append(spot_positions_temp_2[xx2][2])
        ax = plt.subplot(111, projection='polar')
        ax.set_ylim(0, 20)
        c = ax.scatter(plot_phi_2, plot_r_2, s=0.5, alpha=0.15)
        plt.savefig('2polar_plot_'+i+'.pdf')
        plt.close()

        plot_r_3 = []
        plot_phi_3 = []
        for xx3 in range(len(spot_positions_temp_3)):
            plot_r_3.append(spot_positions_temp_3[xx3][1])
            plot_phi_3.append(spot_positions_temp_3[xx3][2])
        ax = plt.subplot(111, projection='polar')
        ax.set_ylim(0, 20)
        c = ax.scatter(plot_phi_3, plot_r_3, s=0.5, alpha=0.15)
        plt.savefig('3polar_plot_'+i+'.pdf')
        plt.close()


if __name__ == '__main__':
	folder_path = 'E:\\Jonas\\harm-new3'
	ringe(folder_path)















# def histo2D(spot_pos, file_name, pattern, histo_number, hist2d_min_max, nice_plots_but_problems_in_fitting, bins, hist2d_spines_off, axis_scaling_factor, dont_show, color_map):
#     """
#     2D histogram plotting function.

#     The function plots a 2D histogram for a given data set vector. It is
#     executed 2 times. The first time the scale range is set so that the
#     histogram shows the whole delayline detector including a circle which
#     shows the round MCP. The second time is auto-scales the axis to provide
#     a zoomed image and to easify the projection fits (they take the data from
#     this histogram function).
#     """

#     xs = [x[1] for x in spot_pos]
#     ys = [x[2] for x in spot_pos]

# # start with a rectangular Figure
#     mainFig = plt.figure(((histo_number-1)*3+1), figsize=(8, 6.3), facecolor='white')

# # define some gridding.
#     if color_map == 1:
#         axHist2d = plt.subplot2grid((9, 9), (1, 0), colspan=8, rowspan=8)
#         axHistx = plt.subplot2grid((9, 9), (0, 0), colspan=8)
#         axHisty = plt.subplot2grid((9, 9), (1, 8), rowspan=8)
#     elif color_map == 2:
#         axHisty = plt.subplot2grid((9, 9), (0, 0), rowspan=1)
#         axHist2d = plt.subplot2grid((9, 9), (0, 1), colspan=9, rowspan=8)


# # # Hist2D colormap def
#     if color_map == 1:      # blue
#         colors = [(1, 1, 1), (0, 0, 1), (0, 0, 0)]
#         n_bins = 1000
#         cmap_name = 'my_name'
#         cm = mpl.colors.LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins, gamma=1.5)
#     elif color_map == 2:    # multi color
#         colors = [(1,1,1),(0.10196, 0.42745, 0), (0.1294, 0.50588,0), (0.22745, 0.8039, 0.1843), (0.87058, 1, 0), (0.9882, 0.996, 0.1961), (0.9686, 0.8, 0.19215), (0.9529, 0.5922, 0.0118), (0.9451, 0.3882, 0.0157), (0.9333, 0.0314, 0.0157), (0.6078, 0.0118, 0), (0.1882, 0, 0)]    # 4&5 (0.6078, 0.8078, 0.19215), (0.8, 0.99215, 0.3842), 
#         n_bins = 1000
#         cmap_name = 'my_name'
#         cm = mpl.colors.LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins, gamma=1.5)
#     elif color_map == 3:    # blau pink
#         colors = [(1, 1, 1), (0.32, 0.9726, 0.9961), (0.2539, 0.7734, 0.8711), (0.6836, 0.4219, 0.4453), (0.9609, 0.2969, 0.4063), (0.9570, 0.1484, 0.3203)]    # , (0.3125, 0.5625, 0.6094)
#         n_bins = 1000
#         cmap_name = 'my_name'
#         cm = mpl.colors.LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins, gamma=1.5)

# # the 2D Histogram:
#     if nice_plots_but_problems_in_fitting == 1:
#         H, xedges, yedges = np.histogram2d(xs, ys, bins=(bins, bins), range=np.array([(-hist2d_min_max, hist2d_min_max), (-hist2d_min_max, hist2d_min_max)]))
#     else:
#         H, xedges, yedges = np.histogram2d(xs, ys, bins=(bins, bins))
#     cax = axHist2d.imshow(H.T, interpolation='nearest', aspect='auto', cmap=cm)

# # make histograms for x and y seperately.
#     if color_map == 1:
#         axHistx.hist(xs, bins=xedges, facecolor='blue', alpha=0.5, edgecolor='None')
#         axHisty.hist(ys, bins=yedges, facecolor='blue', alpha=0.5, orientation='horizontal', edgecolor='None')
# # print some correlation coefficients at the top of the image.
# # mainFig.text(0.05,.90,'r='+str(round(np.corrcoef( xs, ys )[1][0],2))+'; rho='+str(round(spearmanr( xs, ys )[0],2)), style='italic', fontsize=10 )
#     mainFig.text(0.20, .95, 'Event position on MCP-PS', style='normal', fontsize=26)  # '$^{85}$Rb event position on MCP-PS'

# # set axes
#     if color_map == 1:
#         axHistx.set_xlim([xedges.min(), xedges.max()])
#         axHisty.set_ylim([yedges.min(), yedges.max()])
#     axHist2d.set_ylim([axHist2d.get_ylim()[1], axHist2d.get_ylim()[0]])

# # remove some labels
#     nullfmt = mpl.ticker.NullFormatter()
#     if color_map == 1:
#         axHistx.xaxis.set_major_formatter(nullfmt)
#         axHistx.yaxis.set_major_formatter(nullfmt)
#     axHisty.xaxis.set_major_formatter(nullfmt)
#     axHisty.yaxis.set_major_formatter(nullfmt)

# # remove some axes lines
#     if hist2d_spines_off == 1:
#         axHist2d.spines['top'].set_visible(False)
#         axHist2d.spines['right'].set_visible(False)
#         axHist2d.spines['left'].set_visible(False)
#         axHist2d.spines['bottom'].set_visible(False)
#     if color_map == 1:
#         axHistx.spines['top'].set_visible(False)
#         axHistx.spines['right'].set_visible(False)
#         axHistx.spines['left'].set_visible(False)
#     axHisty.spines['top'].set_visible(False)
#     axHisty.spines['bottom'].set_visible(False)
#     axHisty.spines['right'].set_visible(False)
#     if color_map == 2:
#         axHisty.spines['left'].set_visible(False)


#     # remove some ticks
#     if color_map == 1:
#         axHistx.set_xticks([])
#         axHistx.set_yticks([])
#     axHisty.set_xticks([])
#     axHisty.set_yticks([])

# # label 2d hist axes
#     myTicks = np.arange(0, bins, bins/6)
#     axHist2d.set_xticks(myTicks)
#     axHist2d.set_yticks(myTicks)
#     axHist2d.set_xticklabels(np.round([x*axis_scaling_factor for x in xedges[myTicks]]), fontsize=18)
#     axHist2d.set_yticklabels(np.round([x*axis_scaling_factor for x in yedges[myTicks]]), fontsize=18)
#     axHist2d.grid()
#     if hist2d_min_max == 750 and nice_plots_but_problems_in_fitting == 1 and bins == 60:
#         circle = plt.Circle((30.8, 30.5), 28, color='k', linewidth=2, fill=False)
#         axHist2d.add_patch(circle)
#     if hist2d_min_max == 750 and nice_plots_but_problems_in_fitting == 1 and bins == 90:
#         circle = plt.Circle((45, 45), 40, color='k', linewidth=2, fill=False)
#         axHist2d.add_patch(circle)
#     if hist2d_min_max == 750 and nice_plots_but_problems_in_fitting == 1 and bins == 120:
#         circle = plt.Circle((60, 60), 50, color='k', linewidth=2, fill=False)
#         axHist2d.add_patch(circle)
#     for a in [axHist2d]:
#         plt.setp(a.get_xticklabels()[0], visible=False)
#         plt.setp(a.get_yticklabels()[0], visible=False)

#     if color_map == 2:
#         cbar = mainFig.colorbar(cax)  # , ticks=[0, 3, 6, 9, 12])
#         cbar.set_label('# of ions', fontsize=22, labelpad=10, rotation=90)
#         cbar.ax.tick_params(labelsize=18)
#         # cbar.ax.set_yticklabels(['0', '3', '6', '9', '12'])

# # set titles
#     axHist2d.set_xlabel('x / mm', fontsize=22)
#     axHist2d.set_ylabel('y / mm', fontsize=22)
# # axHistx.set_title('x / a.u.', fontsize=10)
# # axHisty.yaxis.set_label_position("right")
# # axHisty.set_ylabel('y / a.u.', fontsize=10, rotation=-90, verticalalignment='top', horizontalalignment='center' )

# # save figure
#     if nice_plots_but_problems_in_fitting == 1:
#         if color_map == 1:
#             plt.savefig('%s_hist_%s_blue.pdf' % (file_name, histo_number))
#         elif color_map == 2:
#             plt.savefig('%s_hist_%s_multicolor.pdf' % (file_name, histo_number))
#     else:
#         if color_map == 1:
#             plt.savefig('%s_hist_%s_zoom_blue.pdf' % (file_name, histo_number))
#         elif color_map == 2:
#             plt.savefig('%s_hist_%s_zoom_multicolor.pdf' % (file_name, histo_number))


#     if dont_show == 1:
#         plt.close()
#     return(H, xedges, yedges, xs, ys)