from python_plotter_functions import *
from read_write_functions import *


file_path = 'G:\\Experiments\\ISOLTRAP\\USER\\Karthein_Jonas\\Master\\2017-05-05-Harmonic_and_long_run\\Efield' #E:\\Welker\\test
file_name = 'Efield_scan_1V35406_spot_positions'
title = '' # '$^{85}$Rb$^+$ mag amp in PI-ICR'      # if you don't need it, leave it empty like ''
legend_text = ''    # # if you don't need it, leave it empty like ''; Otherwise give the points a name!
show_plot = 'yes'   # only 'yes' or 'no' are possible
fit = 'no'  # only 'yes' or 'no' are possible
fit_function = 'x^2'   # chose from: 'linear', 'gauss', 'x^2', 'x^4', 'exp', 'sin'
font = 'Utopia'        # chose from: 'LaTeX' = Computer Modern Serif or 'Utopia' = \\usepackage{fourier}
color_points = 'blue'      # chose from: 'green', 'blue', 'yellow', 'black', 'red'
color_fit = 'red'        # chose from: 'green', 'blue', 'yellow', 'black', 'red'
fit_range = 'full'      # chose from: 'full', 'partly'
x_fit_min = 400
x_fit_max = 500
fit_b_parameter_start = 1
plot_type = '2dhistogram-mcp'
bins = 250 # 1200 is good for small signatures, 300 is a good one for publication
markersize = 5 #usual it is 5
grid_on_off = 'off'
png_on_off = 'on'
diff_to_AME_err = -99999
axis_scaling_factor = 0.031746032*0.5       # convert a.u. to mm


import_list = read_csv(file_path, file_name)
# import_c = read_csv(file_path, '127Cd_c_2017')
# import_list_2 = read_csv(file_path, '85Rb_c_216_spot_positions')
# import_list = read_csv(file_path, '85Rb_c_016_spot_positions')
# import_list_2 = read_csv(file_path, '85Rb_c_026_spot_positions')
# for i in import_list[1:]:
#     i[2] = float(i[2]) + 100
# for i in import_list_2:
#     i[1] = float(i[1]) - 89.61
#     i[2] = float(i[2]) + 229.06
x_y_yerr = [[float(j)*axis_scaling_factor*2 for j in i[1:3]] for i in import_list[1:]]
# x_y_yerr_c = [[float(j)*axis_scaling_factor for j in i[1:3]] for i in import_c[1:]]
for i in x_y_yerr:
    # i[0] = float(i[0]) + 1.4
    i[1] = float(i[1]) + 2.6
# x_y_yerr = [i for i in x_y_yerr_1 if (i[0] > 6.2 and i[0] < 11 and i[1] < 1.3 and i[1] > -3.9) or (i[0] > -7.7 and i[0] < -3 and i[1] < -7.4 and i[1] > -10.7)]
# count_isomer = [i for i in x_y_yerr_1 if (i[0] > 6.2 and i[0] < 11 and i[1] < 1.3 and i[1] > -3.9)]
# count_ground = [i for i in x_y_yerr_1 if (i[0] > -7.7 and i[0] < -3 and i[1] < -7.4 and i[1] > -10.7)]
# print 'Isomeric state :: ',len(count_isomer)
# print 'Ground state   :: ',len(count_ground)
# print 'Ratio          :: ',float(len(count_isomer))/float(len(count_ground))
# x_y_yerr.extend(x_y_yerr_c[:100])
x_y_yerr.insert(0, ['x / mm', 'y / mm'])

xss, yss, yerr = python_plot(x_y_yerr, title, file_name, show_plot, legend_text, fit, fit_function, font, color_fit, color_points, fit_range, x_fit_min, x_fit_max, fit_b_parameter_start, plot_type, bins, markersize, grid_on_off, png_on_off, '', '', diff_to_AME_err)
