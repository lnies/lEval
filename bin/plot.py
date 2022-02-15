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

from utilities import FitToDict, Peaks, softCool

def custom_colors():
	# Color definition from coolors.co
	colors = {
		# https://coolors.co/d1495b-edae49-abe188-1a6fdf-0c0a3e
		'Fall_rgb': {'red':'#D1495B', 'orange':'#EDAE49', 'green':'#ABE188','blue':'#1A6FDF','purple':'#0C0A3E'},
		'Jo-s_favs': 
		{				
						'black': "#000000",
						'red': "#FF2D55", 
						'blue': "#00A2FF",
						'orange': "#FFCC00", 
						'green': "#61D935", 
						'grey': "#C0C0C0", 
						'purple': "#C177DA", 
						'lightblue': "#6FF1E9",
		}
	}
	return colors

def simple_error_plt(y, y_err, x='', x_labels='', \
					 label = ["ISOLTRAP"], x_label='', y_label=[''], title='', \
					 ref_value=None, ref_err=None, ref_legend_label='AME20 Error', ref_axis=0,
					 x_share = False,
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
	- title: plot title
	'''
	colors = custom_colors()
	colors_keys = list(colors['Jo-s_favs'].keys())
	mpl.rc('text', usetex=False)
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.5*1.5, 4.5))
	
	if len(x_labels) != 0 and len(x) == 0:
		x = np.arange(0,len(x_labels),1)
		plt.xticks(x, x_labels, size = 10, rotation = 0)
	elif len(x_labels) == '' and len(x) != '':
		x = x
	else:
		print("Wrong x-data passed. Can only pass either x or x_labels")

	# Loop through list of arrays passed
	twins = []
	for i in np.arange(0,len(y),1):
		# If x_labels are passed
		if len(x_labels) != 0:
			x_plot = np.arange(0,len(x_labels),1)
		else:
			x_plot = x[i]
		#
		if i == 0:
			ax.errorbar(x_plot, y[i], y_err[i],
					   fmt='o', color=colors['Jo-s_favs'][colors_keys[i]], zorder = 2, 
					   label=label[i], 
					   fillstyle='full', mfc="black", linewidth=2, ms =10
			)
			ax.set_ylabel(y_label[i], size=18) #fontweight='bold')
			ax.tick_params(direction='out')
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
				   label=label[i], 
				   fillstyle='full', mfc=color, linewidth=2, ms =10
			)

			ax.set_zorder(twins[i-1].get_zorder()+1)
			ax.patch.set_visible(False)


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
						label = ref_legend_label)
				ax.plot(x_plot, [ref_value for value in x_plot], color='grey')
			else:
				twins[ref_axis-1].fill_between(x_plot, ref_value-ref_err, ref_value+ref_err, facecolor='0.5', alpha=0.5,
						label = ref_legend_label)
				twins[ref_axis-1].plot(x_plot, [ref_value for value in x_plot], color='grey')


	# plt.axhspan(ref_value-ref_err, ref_value+ref_err, facecolor='0.5', alpha=0.5)

	# handles, labels = ax.get_legend_handles_labels()
	# handles = [h[0] for h in handles] # Remove errorbars from legend
	
	plt.legend()#handles, labels, fontsize=12, frameon=False, ncol=1, loc="upper left")
	
	plt.tight_layout()
	plt.title(title, fontsize=20)
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