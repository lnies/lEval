from python_plotter_functions import *
from read_write_functions import *


# ------------------------
# Procedure to start:
# ------------------------
# Shell:
# 	cd /
# 	cd \\\\cerndfs13.cern.ch\\dfs\\Experiments\\ISOLTRAP\\Software\\PI-ICR\\Python-DAQ
# Script:
# 	create excel file with x data in first column, y data in second column and optional yerr data in third. The first line
#	change file_name
# ------------------------
# Variables to change:
# ------------------------
file_path = 'G:\\Experiments\\ISOLTRAP\\DATA\\2017\\Screenshots\\2017-08-ArKr-run-prep' # 'G:\\Experiments\\ISOLTRAP\\DATA\\2017\\2017-07-Xe-Cs-run-prep\\analyzed\\2017-07-16-85Rb-nuamp-133Cs_run20' # 'G:\\Experiments\\ISOLTRAP\\DATA\\2017\\Screenshots\\2017-06' # '/Volumes/dfs/USER/Karthein_Jonas/harm-new' # 'G:\\Experiments\\ISOLTRAP\\DATA\\2017\\2017-06-tests' #'/Volumes/dfs/DATA/2017/2017-04-tuning/' # 'G:\\Experiments\\ISOLTRAP\\DATA\\2017\\2017-05_ArKr_run'
data_folder_path = file_path # 'G:\\Experiments\\ISOLTRAP\\DATA\\2017\\2017-05_ArKr_run'		# MAC '/Volumes/Macintosh HD/Users/jonaskarthein/cernbox/Python/test10'
file_name = 'UTmagPhase_85Rb'  # name of the excel sheet without the ending .xlsx
title = '$^{85}$Rb UTmagPhase scan' # '$^{85}$Rb$^+$ mag amp in PI-ICR'      # if you don't need it, leave it empty like ''
legend_text = ''    # # if you don't need it, leave it empty like ''; Otherwise give the points a name!
show_plot = 'yes'   # only 'yes' or 'no' are possible
fit = 'yes'  # only 'yes' or 'no' are possible
fit_function = 'x^2'   # chose from: 'linear', 'gauss', 'x^2', 'x^4', 'exp', 'sin'
font = 'Utopia'        # chose from: 'LaTeX' = Computer Modern Serif or 'Utopia' = \\usepackage{fourier}
color_points = 'green'      # chose from: 'green', 'blue', 'yellow', 'black', 'red'
color_fit = 'red'        # chose from: 'green', 'blue', 'yellow', 'black', 'red'
fit_range = 'partly'      # chose from: 'full', 'partly'
x_fit_min = 370
x_fit_max = 820
fit_b_parameter_start = 1
plot_type = 'scatter'
bins = 20
markersize = 5
grid_on_off = 'on'
png_on_off = 'on'

x_y_yerr = read_excel(file_path, file_name)
if fit == 'yes':
	xss, yss, yerr, parameter_mlf, parameter_mlf_err = python_plot(x_y_yerr, title, file_name, show_plot, legend_text, fit, fit_function, font, color_fit, color_points, fit_range, x_fit_min, x_fit_max, fit_b_parameter_start, plot_type, bins, markersize, grid_on_off, png_on_off, '', '', '')
else:
	xss, yss, yerr = python_plot(x_y_yerr, title, file_name, show_plot, legend_text, fit, fit_function, font, color_fit, color_points, fit_range, x_fit_min, x_fit_max, fit_b_parameter_start, plot_type, bins, markersize, grid_on_off, png_on_off, '', '', '')
