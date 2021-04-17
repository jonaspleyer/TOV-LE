#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 06:07:29 2020

@author: jonas
"""

import pymongo
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import tikzplotlib
from fractions import Fraction

class Plotter():
	def __init__(self):
		# Initialise database access
		self.client = pymongo.MongoClient()
		self.db = self.client["TOV-LE"]
		self.coll = self.db["Exponents-polytropic-EOS"]
		
		# Create plotstyles
# 		self.curve_styles = ['-', '--', '-.', ':',(0, (1, 10)), (0, (1, 1))]
		self.curve_styles = [
			(0, (1, 2)),
			(0, (7, 1, 4, 1)),
			(0, (1, 4)),
			(0, (5, 2, 1, 2)),
			(0, (1, 6)),
			(0, (3, 2, 1, 2, 1, 2))
			]
		self.curve_colours = ['k','r', 'g', 'b', 'orangered', 'lime','cornflowerblue']
	
	def plotInitialVals(self, var1_name, var1_vals, var1_fixed_val, var2_name, var2_vals, var2_fixed_val, terms, ylim_low=None, ylim_up=None, xlim_low=None, xlim_up=None, noExport=True):
		# Think about var1_name and var2_name as something like 'p_init'
		# Think about var1_vals as something like [0.001,0.01,0.1,1,10]
		# Think about var2_fixed_val as something like 10
		# var1 is the evolving Variable for which different results will be plotted
		# var2 is fixed to allow comparison between the results
		# terms is a list [0,1] that indicates which terms will be plotted
		
		if xlim_up==None:
 			xlim_up=99
		if xlim_low==None:
 			xlim_low=-99
		
		# Get a dict with all LE results for different var1_vals
		results_LE_1 = {}
		for var1_val in var1_vals:
			results_LE_1[var1_val] = list(self.coll.find({var1_name:var1_val, var2_name:var2_fixed_val,'succ':True,'equ_type':'LE', 'terms':0,'n':{"$lt":xlim_up, "$gt":xlim_low}}))
		
		results_LE_2 = {}
		for var2_val in var2_vals:
			results_LE_2[var2_val] = list(self.coll.find({var2_name:var2_val, var1_name:var1_fixed_val,'succ':True,'equ_type':'LE', 'terms':0,'n':{"$lt":xlim_up, "$gt":xlim_low}}))
		
		# Create a list of dicts with all TOV results for different var1_vals
		results_TOV_1 = [{}]*(len(terms)+1)
		for term in terms:
			for var1_val in var1_vals:
				results_TOV_1[term][var1_val] = list(self.coll.find({var1_name:var1_val, var2_name:var2_fixed_val,'succ':True,'equ_type':'TOV', 'terms':term, 'n':{"$lt":xlim_up, "$gt":xlim_low}}))
		
		results_TOV_2 = [{}]*(len(terms)+1)
		for term in terms:
			for var2_val in var2_vals:
				results_TOV_2[term][var2_val] = list(self.coll.find({var2_name:var2_val, var1_name:var1_fixed_val,'succ':True,'equ_type':'TOV', 'terms':term, 'n':{"$lt":xlim_up, "$gt":xlim_low}}))
		
		# Finally get the plot results for LE
		plot_results_LE_1 = {}
		for var1_val in var1_vals:
			arr = np.array([[res_LE['n'], res_LE['r_max']] for res_LE in results_LE_1[var1_val]])
			plot_results_LE_1[var1_val] = arr[arr[:,0].argsort()]
		
		plot_results_LE_2 = {}
		for var2_val in var2_vals:
			arr = np.array([[res_LE['n'], res_LE['r_max']] for res_LE in results_LE_2[var2_val]])
			plot_results_LE_2[var2_val] = arr[arr[:,0].argsort()]
		
		# and for TOV
		plot_results_TOV_1 = [{}]*(len(terms)+1)
		for term in terms:
			for var1_val in var1_vals:
				arr = np.array([[res_TOV['n'], res_TOV['r_max']] for res_TOV in results_TOV_1[term][var1_val]])
				plot_results_TOV_1[term][var1_val] = arr[arr[:,0].argsort()]
		
		plot_results_TOV_2 = [{}]*(len(terms)+1)
		for term in terms:
			for var2_val in var2_vals:
				arr = np.array([[res_TOV['n'], res_TOV['r_max']] for res_TOV in results_TOV_2[term][var2_val]])
				plot_results_TOV_2[term][var2_val] = arr[arr[:,0].argsort()]
		
		if var1_name=="p_init":
			var1_name=r'$p_0$'
		else:
			var1_name=r'A'
			
		if var2_name=="p_init":
			var2_name=r'$p_0$'
		else:
			var2_name=r'A'
		
		# Initialise plot with right size
		cm = 1/2.54
		plt.figure(figsize=[16*cm,18*cm])
		
		gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[1, 1]) 
		ax1 = plt.subplot(gs[0])
		ax1.xaxis.set_label_position('top') 
		ax2 = plt.subplot(gs[1])
		
		for i, var1_val in enumerate(var1_vals):
			ax1.plot(plot_results_LE_1[var1_val][:,0], plot_results_LE_1[var1_val][:,1], linestyle=self.curve_styles[2*i], c=self.curve_colours[0],  label="LE for " + var1_name + "=" + str(var1_val))
			for term in terms:
				if len(terms)==1:
					label = "TOV for " + var1_name + "=" + str(var1_val)
				else:
					label = "TOV with term=" + str(term) + " for " + var1_name + "=" + str(var1_val)
				ax1.plot(plot_results_TOV_1[term][var1_val][:,0], plot_results_TOV_1[term][var1_val][:,1], linestyle=self.curve_styles[2*i+1], c=self.curve_colours[0], label=label)
		
		ax1.set_yscale('log')
		ax1.set_ylabel(r'Radial Zero Value $r_0$')
		ax1.set_ylim([ylim_low,ylim_up])
		ax1.set_xlabel(r'$\gamma=1+n^{-1}$')
		ax1.legend(title=r'Fix ' + var2_name + "=" + str(var2_fixed_val))
		
		for i, var2_val in enumerate(var2_vals):
			ax2.plot(plot_results_LE_2[var2_val][:,0], plot_results_LE_2[var2_val][:,1], linestyle=self.curve_styles[2*i], c=self.curve_colours[0],  label="LE for " + var2_name + "=" + str(var2_val))
			for term in terms:
				if len(terms)==1:
					label = "TOV for " + var2_name + "=" + str(var2_val)
				else:
					label = "TOV with term=" + str(term) + " for " + var2_name + "=" + str(var2_val)
				ax2.plot(plot_results_TOV_2[term][var2_val][:,0], plot_results_TOV_2[term][var2_val][:,1], linestyle=self.curve_styles[2*i+1], c=self.curve_colours[0], label=label)
		
		ax2.set_yscale('log')
		ax2.set_ylim([ylim_low,ylim_up])
		ax2.set_xlabel(r'Polytropic Index $n=\frac{1}{\gamma-1}$')
		ax2.legend(title=r'Fix ' + var1_name + "=" + str(var1_fixed_val))
		
		locs_1 = ax2.get_xticks()[1:-1]
		f = lambda x: Fraction((x+1)/x).limit_denominator() if x != 0 else r'$\infty$'
		labels_2 = [f(loc1) for loc1 in locs_1]
		ax1.xaxis.tick_top()
		ax1.set_xticks(locs_1)
		ax1.set_xticklabels(labels_2)
		
		ax2.set_ylabel(r'Radial Zero Value $r_0$')
		plt.tight_layout()
		
		filename = 'pictures/TOV-Exponents-LESubs-InitialVals-Database-PlotResults-Combo'
		
		plt.savefig(filename + '.svg')
		plt.subplots_adjust(hspace=.0)
		plt.show()
		
		matplotlib.use("pgf")
		matplotlib.rcParams.update({
# 		    "pgf.texsystem": "pdflatex",
		    'font.family': 'serif',
		    'text.usetex': True,
		    'pgf.rcfonts': False,
		})
		if noExport==False:
 			plt.savefig(filename + ".pgf", dpi=1000, bbox_inches='tight')
 			print(filename + ".pgf")
 			tikzplotlib.save(filename + ".tex")
 			print(filename + ".tex")

Plodda = Plotter()

Plodda.plotInitialVals(var1_name='p_init', var1_vals=[0.1,1], var1_fixed_val=1, var2_name='A_init', var2_vals=[0.1,1,8], var2_fixed_val=1, terms=[0], xlim_low=0, ylim_up=None, noExport=False)