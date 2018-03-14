from python_plotter_functions import *
from read_write_functions import *
import math

import matplotlib.pyplot as plt


file_path = '/Users/jonaskarthein/cernbox/google_sync/CERN/Paper/129Cd/physical-review-letters/figures/'
file_name = 'Z2-run4-129gCd-run7-p2'
title = '' # '$^{85}$Rb$^+$ mag amp in PI-ICR'      # if you don't need it, leave it empty like ''
legend_text = ''    # # if you don't need it, leave it empty like ''; Otherwise give the points a name!
show_plot = 'yes'   # only 'yes' or 'no' are possible
fit = 'no'  # only 'yes' or 'no' are possible
fit_function = 'x^2'   # chose from: 'linear', 'gauss', 'x^2', 'x^4', 'exp', 'sin'
font = 'Utopia'        # chose from: 'LaTeX' = Computer Modern Serif or 'Utopia' = \\usepackage{fourier}
color_points = 'green'      # chose from: 'green', 'blue', 'yellow', 'black', 'red'
color_fit = 'red'        # chose from: 'green', 'blue', 'yellow', 'black', 'red'
fit_range = 'full'      # chose from: 'full', 'partly'
x_fit_min = 400
x_fit_max = 500
fit_b_parameter_start = 1
plot_type = '2dhistogram-mcp'
bins = 15
markersize = 5
grid_on_off = 'on'
png_on_off = 'off'
diff_to_AME_err = -99999
axis_scaling_factor = 0.031746032       # convert a.u. to mm


if font == 'Utopia':
    mpl.rc('font', family='serif', serif='Utopia')      # Utopia LaTeX font!!
    mpl.rc('text', usetex=False)
elif font == 'LaTeX':
    mpl.rc('font', family='serif')      # Standard LaTeX font!!
    mpl.rc('text', usetex=True)


### import normal
import_list = read_csv(file_path, file_name)
print import_list

### import additional points
import_c = read_csv(file_path, 'Z2-run4-129gCd-run7-c')
# import_list_2 = read_csv(file_path, '85Rb_c_216_spot_positions')
# import_list = read_csv(file_path, '85Rb_c_016_spot_positions')
# import_list_2 = read_csv(file_path, '85Rb_c_026_spot_positions')

### manipulate position in ch.
# for i in import_list[1:]:
#     i[2] = float(i[2]) + 100
# for i in import_list_2:
#     i[1] = float(i[1]) - 89.61
#     i[2] = float(i[2]) + 229.06

### manipulate in polar coordinates
# radius_manipulator = -3      # increase or decrease in mm
# angle_manipulator = 1.5       # rotate in rad

# import_list_polar = [[float(entry[1])*axis_scaling_factor/2, float(entry[2])*axis_scaling_factor/2] for entry in import_list[len(import_list)/3:4*len(import_list)/5]]
# list_polar = [[math.sqrt(entry[0]**2+entry[1]**2) + radius_manipulator, math.atan2(entry[1], entry[0]) + angle_manipulator] for entry in import_list_polar]

# x_y_yerr = [[float(entry[0])*math.cos(entry[1]), float(entry[0])*math.sin(entry[1])*0.9] for entry in list_polar]



### convert data to mm
x_y_yerr = [[float(entry[1])*axis_scaling_factor, float(entry[2])*axis_scaling_factor] for entry in import_list[1:]]
            # if float(entry[1])*axis_scaling_factor<3 and float(entry[1])*axis_scaling_factor>-3
            # and float(entry[2])*axis_scaling_factor<3 and float(entry[2])*axis_scaling_factor>-3]

### manipulate position in mm
x_y_yerr_c = [[float(j)*axis_scaling_factor for j in i[1:3]] for i in import_c[1:]]
for i in x_y_yerr_c:
    i[0] = float(i[0]) + 1.9
#     i[1] = float(i[1]) + 2.6
# x_y_yerr = [i for i in x_y_yerr_1 if (i[0] > 6.2 and i[0] < 11 and i[1] < 1.3 and i[1] > -3.9) or (i[0] > -7.7 and i[0] < -3 and i[1] < -7.4 and i[1] > -10.7)]

### additional isomer stuff
# count_isomer = [i for i in x_y_yerr_1 if (i[0] > 6.2 and i[0] < 11 and i[1] < 1.3 and i[1] > -3.9)]
# count_ground = [i for i in x_y_yerr_1 if (i[0] > -7.7 and i[0] < -3 and i[1] < -7.4 and i[1] > -10.7)]
# print 'Isomeric state :: ',len(count_isomer)
# print 'Ground state   :: ',len(count_ground)
# print 'Ratio          :: ',float(len(count_isomer))/float(len(count_ground))
x_y_yerr.extend(x_y_yerr_c)

### insert axis labels
x_y_yerr.insert(0, ['x (mm)', 'y (mm)'])

plt.figure(figsize=(8,8))
plt.plot([x[0] for x in x_y_yerr[1:]], [y[1] for y in x_y_yerr[1:]], 'o', alpha=0.18, markersize=7, markeredgewidth=0.0)
plt.xlabel('x (mm)', fontsize=22)
plt.ylabel('y (mm)', fontsize=22)
plt.xticks([-18,-12,-6,0,6,12,18], fontsize=22)
plt.yticks([-18,-12,-6,0,6,12,18], fontsize=22)

plt.axis([-23.81, 23.81, -23.81, 23.81])
ax = plt.gca()
mcp = Ellipse(xy=(0, 0), width=45.714, height=45.714, edgecolor='k', fc='None', lw=2)
ax.add_patch(mcp)

plt.tight_layout()

plt.savefig('/Users/jonaskarthein/cernbox/google_sync/CERN/Paper/129Cd/physical-review-letters/figures/IS574_129Cd.pdf')
plt.show()

### plot
# xss, yss, yerr = python_plot(x_y_yerr, title, file_name, show_plot, legend_text, fit, fit_function, font, color_fit, color_points, fit_range, x_fit_min, x_fit_max, fit_b_parameter_start, plot_type, bins, markersize, grid_on_off, png_on_off, '', '', diff_to_AME_err)
