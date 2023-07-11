# -*- coding: utf-8 -*-
"""
Created on Mon 11 January 2022
Modified by Lukas.Nies@cern.ch on 04/02/2021
@author: Lukas Nies
@contact: Lukas.Nies@cern.ch
@license: MIT
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy as sc
from scipy.signal import find_peaks
from scipy.optimize import least_squares
from scipy import stats
import re
import math 
import time
import mmap

from utilities import FitToDict

def simple_error_plt(y, y_err, x='', x_labels='', \
					 label = ["ISOLTRAP"], x_label='', y_label=[''], title='', \
					 ref_value=None, ref_err=None, ref_legend_label='AME20 Error', ref_axis=0,
					 with_lines = False, h_lines = [],
					 x_share = False, figsize = (4.5*1.5, 4.5),
					 ):
	'''
	Simple scatter plot with y error bars.
	Parameters:
	- y: y-values
	- y_err: y-errors
	- x: x-data (exclusive with x_labels)
	- x_labels: array of strings to be used as x-labels
	- x_label: x-axis labeling
	- y_label: y-axis labeling
	- with_lines: if True, connects scatter data points with lines
	- h_lines: array-like. Draws in hlines
	- title: plot title
	'''
	colors = custom_colors()
	colors_keys = list(colors['Jo-s_favs'].keys())
	mpl.rc('text', usetex=False)
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize)
	
	if len(x_labels) != 0 and len(x) == 0:
		x = np.arange(0,len(x_labels),1)
		plt.xticks(x, x_labels, size = 10, rotation = 0)
	elif len(x_labels) == '' and len(x) != '':
		x = x
	else:
		print("Wrong x-data passed. Can only pass either x or x_labels")
		return 0

	# Loop through list of arrays passed
	twins = []
	for i in np.arange(0,len(y),1):
		# If x_labels are passed
		if len(x_labels) != 0:
			x_plot = np.arange(0,len(x_labels),1)
			ax.set_xticklabels(labels=x_labels)

		else:
			x_plot = x[i]
		#
		if i == 0:
			ax.errorbar(x_plot, y[i], y_err[i],
				fmt='o', color=colors['Jo-s_favs'][colors_keys[i]], zorder = 2, 
				label=label[i], markeredgewidth=2, mec = 'black',
				fillstyle='full', mfc="black", linewidth=2, ms =10
			)
			if with_lines:
				ax.plot(x_plot, y[i], "-", color=colors['Jo-s_favs'][colors_keys[i]], zorder = 2, linewidth=2)
			#
			ax.set_ylabel(y_label[i], size=18) #fontweight='bold')
			ax.tick_params(direction='out')
			#
		else:
			if x_share:
				twins.append(ax.twinx())
				twins[i-1].spines['right'].set_position(("axes", 0.8+i*0.2))
				color = colors['Jo-s_favs'][colors_keys[i]]
				twins[i-1].set_ylabel(y_label[i], color=color, fontsize=18)  # we already handled the x-label with ax1
				twins[i-1].tick_params(axis = 'y', direction='in', color=color)
				twins[i-1].yaxis.label.set_color(color)
				ax.spines["right"].set_color(color)
				twins[i-1].spines["right"].set_color(color)
			else:
				twins.append(ax)
				color = colors['Jo-s_favs'][colors_keys[i]]

			twins[i-1].errorbar(x_plot, y[i], y_err[i],
				   fmt='o', color=color, zorder = 2, 
				   label=label[i], markeredgewidth=2, mec = 'black',
				   fillstyle='full', mfc=color, linewidth=2, ms =10
			)
			if with_lines:
				twins[i-1].plot(x_plot, y[i], color=color, zorder = 2, linewidth=2)


			# ax.set_zorder(twins[i-1].get_zorder()+1)
			# ax.patch.set_visible(False)


	#

	# ax_t = ax.secondary_xaxis('top')
	# ax_t.set_xticks(x)

	ax.set_xlabel(x_label, size=18) #fontweight='bold')
	# ax_r = ax.secondary_yaxis('right')
	# # ax.set_xlabel("Mercury Isotopes", size=14, fontweight='bold')
	# ax_t.tick_params(axis='x', direction='in', labeltop=False)
	# ax_r.tick_params(axis='y', direction='in', labelright=False)


	if ref_value is not None and ref_err is not None:
		# Error band
		if len(x_labels) != 0:
			# 
			x_plot = np.arange(0,len(x_labels),1)
			if ref_axis == 0:
				ax.fill_between(x_plot, ref_value-ref_err, ref_value+ref_err, facecolor='0.5', alpha=0.5,
						label = ref_legend_label, zorder = 1)
				ax.plot(x_plot, [ref_value for value in x_plot], color='grey', zorder = 1)
			else:
				twins[ref_axis-1].fill_between(x_plot, ref_value-ref_err, ref_value+ref_err, facecolor='0.5', alpha=0.5,
						label = ref_legend_label, zorder = 1)
				twins[ref_axis-1].plot(x_plot, [ref_value for value in x_plot], color='grey', zorder = 1)


	# h_lines
	# if len(hlines) != 0:
	# 	for line in h_lines:



	# plt.axhspan(ref_value-ref_err, ref_value+ref_err, facecolor='0.5', alpha=0.5)

	# handles, labels = ax.get_legend_handles_labels()
	# handles = [h[0] for h in handles] # Remove errorbars from legend
	
	plt.legend()#handles, labels, fontsize=12, frameon=False, ncol=1, loc="upper left")
	
	plt.title(title, fontsize=20)

	plt.tight_layout()
	# plt.savefig("./Mercury_Masses_Comparison.pdf", dpi=300)
	plt.show()

def simple_fit_plot(df, fit_files, bins = 1, log=False, file_out=False, legend=True, style='errorbar',
					labels=[]
					):
		"""  
		Wrapper for plotting fits and data
			- bins: number of bins to rebin. Defaults to 1, e.g. no rebinning
			- log: plot y-scale in log. Defaults to false
			- file_out: name to save plot as .pdf. Defaults to False, e.g. no saving
			- legend: Plot legend if
		"""
		#
		colors = custom_colors()
		colors_keys = list(colors['Jo-s_favs'].keys())
		#
		if len(labels) == 0:
			labels = ['Fit' for i in range(len(fit_files))]
		# Load and plot the data
		xdata = df.tof
		tof_mean = np.mean(xdata)
		binning = round((xdata.max()-xdata.min())/0.8/bins)
		n, xe = np.histogram(xdata, bins=binning)
		cx = 0.5 * (xe[1:] + xe[:-1]) 
		dx = np.diff(xe)
		# Plot data
		if style == 'errorbar':
			plt.errorbar(cx - tof_mean,
					 n, n ** 0.5, fmt="ok", zorder=1, label=f"Data (bins={bins})")
		elif style == 'hist':
			plt.hist((xdata - tof_mean), bins=binning, color='black', label=f"Data (bins={bins})")

		# Loop through all passed fit files for xdata
		i = 0
		for file in fit_files:
			#
			fitfromfile = FitToDict(file)
			Var_dict = {}
			#
			if not 'FIT-VALUES' in fitfromfile.fit:
				print(f"Fit file {file} has no fit values that were exported.")
				return 0
			#
			xm = fitfromfile.fit['FIT-VALUES'].tof
			y_val = fitfromfile.fit['FIT-VALUES'].fit
			#
			for idx,row in fitfromfile.fit['RESULTS-TABLE'].iterrows():
				Var_dict[row['var']] = row['value']
			Var_dict['xmin'] = float(fitfromfile.fit['META-DATA']['xmin'])
			Var_dict['xmax'] = float(fitfromfile.fit['META-DATA']['xmax'])

			# Normalize values
			integral_cut = sum(y_val) * np.diff(xm)[0]
			left_n_cut = len(xe[xe<Var_dict['xmin']])
			right_n_cut = len(xe[xe<Var_dict['xmax']])
			n_cut = n[left_n_cut:right_n_cut]        
			y_val = y_val / integral_cut * sum(n_cut) * dx[0]

			# Plot fit	
			plt.plot(xm - tof_mean, 
					 y_val, label=f"{labels[i]}", c=colors['Jo-s_favs'][colors_keys[i+1]], zorder=3, linewidth=3)

			#
			i += 1
		
		# Get y axis limits
		ylims = plt.ylim()
		if log:
			plt.yscale("log")
			plt.ylim(0.5,2*ylims[1])

		# Add axis labels
		plt.xlabel(f'Time-of-Flight [ns] - {tof_mean:.1f}ns', fontsize=20)
		plt.ylabel(f'Counts per bin', fontsize=20)

		# Format Legend
		if legend:
			plt.legend(fontsize=20)

		# Save plot
		if file_out != False:
			print(f"Plot fit save as {file_out}")
			plt.savefig(file_out, dpi=300)

def plot_laser_on_off(df, bins = 10, n_per_laser = 100):
	"""  
	Wrapper for plotting laser-on-off data (shot-to-shot basis)
		- bins: number of bins to rebin. Defaults to 10
		- n_per_laser: number of laser on / laser on shots per slice (MCS slice)
	"""
	
	# Get binning of raw data
	# Get min and max tof from data frame
	minn = df.tof.min()
	maxx = df.tof.max()
	# Avoid having empty binning when maxx is equal to minn
	if minn == maxx:
		maxx += 1
	#
	binning = round((maxx-minn)/0.8/bins)	

	fs_labels = 20
	fs_ticks = 20
	
	df['sweeps_floored'] = np.floor(df.sweep/n_per_laser)

	xdata = df.tof
	n, xe = np.histogram(xdata, bins=binning, range=(df.tof.min(), df.tof.max()))
	cx = 0.5 * (xe[1:] + xe[:-1])
	dx = np.diff(xe)

	xdata1 = df.tof[df['sweeps_floored']%2==0] #df.tof[df['sweeps_floored']%2!=0]
	n1, xe1 = np.histogram(xdata1, bins=binning, range=(df.tof.min(), df.tof.max()))
	cx1 = 0.5 * (xe1[1:] + xe1[:-1])
	dx1 = np.diff(xe1)

	xdata2 = df.tof[df['sweeps_floored']%2!=0] #df.tof[df['sweeps_floored']%2!=0]
	n2, xe2 = np.histogram(xdata2, bins=binning, range=(df.tof.min(), df.tof.max()))
	cx2 = 0.5 * (xe2[1:] + xe2[:-1])
	dx2 = np.diff(xe2)

	#use sqrt(n) as error, if n==1 use smaller error to avoid having inifite long error bars in log-scale
	n_sub = n1-n2
	plt.step(cx1, n_sub, #[val ** 0.5 if val != 1 else 0.75 for val in n] ,
		color='r', label=f"Laser On - Laser off (bins={bins})",
		 linewidth=4, zorder = 2,
		#fmt="ok", zorder=2, label=f"Laser On - Laser off (bins={bins})")
	   )
	plt.step(cx, n, #[val ** 0.5 if val != 1 else 0.75 for val in n] ,
		color='black', label=f"Laser On",
		 linewidth=4, zorder = 1,
		#fmt="ok", zorder=2, label=f"Laser On - Laser off (bins={bins})")
	   )


	# plt.errorbar(cx, n, [val ** 0.5 if val != 1 else 0.75 for val in n] ,
	#         ecolor='b', elinewidth=1, color='b', 
	#         fmt="ok", zorder=1, label=f"Laser on (bins={bins})")
	# plt.errorbar(cx1, n_sub, [val ** 0.5 if val != 1 else 0.75 for val in n] ,
	#         ecolor='r', elinewidth=1, color='r', mfc='r',
	#         fmt="ok", zorder=2, label=f"Laser On - Laser off (bins={bins})")

	# plt.errorbar(cx, n, [val ** 0.5 if val != 1 else 0.75 for val in n] ,
	#         ecolor='b', elinewidth=1, color='b', 
	#         fmt="ok", zorder=1, label=f"Laser on (bins={bins})")
	# plt.plot(xdata, np.zeros_like(xdata)-5, "|", alpha=0.1, label = "ToF Data", zorder = 3)

	# # plt.hist((xdata), bins=get_binning(df = df, bins=bins), 
	# #          color='grey', edgecolor='black', linewidth=0.1, label=f"Data (bins={bins})")


	# # plt.errorbar(cx, n, n ** 0.5, fmt="ok", zorder=1)
	# #
	# xm = np.linspace(xe[0], xe[-1], num=1000)
	plt.legend();
	# # plt.xlim(peaks.pos[0]-300, peaks.pos[0]+300)

	# # Add axis labels
	plt.xlabel(f'Time-of-Flight [ns]', fontsize=fs_labels)
	plt.ylabel(f'Counts per bin', fontsize=fs_labels)

	# # Set ticks size 
	# plt.tick_params(axis='both', which='major', labelsize=fs_ticks)

	plt.tight_layout()