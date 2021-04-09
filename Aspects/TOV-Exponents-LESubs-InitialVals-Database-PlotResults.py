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
	
	def plotInitialVals(self, var1_name, var1_vals, var2_name, var2_fixed_val, terms, ylim_low=None, ylim_up=None, xlim_low=None, xlim_up=None, noExport=True):
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
		results_LE = {}
		for var1_val in var1_vals:
			results_LE[var1_val] = list(self.coll.find({var1_name:var1_val, var2_name:var2_fixed_val,'succ':True,'equ_type':'LE', 'terms':0,'n':{"$lt":xlim_up},'n':{"$gt":xlim_low}}))
		# Create a list of dicts with all TOV results for different var1_vals
		results_TOV = [{}]*(len(terms)+1)
		for term in terms:
			for var1_val in var1_vals:
				results_TOV[term][var1_val] = list(self.coll.find({var1_name:var1_val, var2_name:var2_fixed_val,'succ':True,'equ_type':'TOV', 'terms':term, 'n':{"$lt":xlim_up},'n':{"$gt":xlim_low}}))
		
		# Finally get the plot results for LE
		plot_results_LE = {}
		for var1_val in var1_vals:
			arr = np.array([[res_LE['n'], res_LE['r_max']] for res_LE in results_LE[var1_val]])
			plot_results_LE[var1_val] = arr[arr[:,0].argsort()]
		# and for TOV
		plot_results_TOV = [{}]*(len(terms)+1)
		for term in terms:
			for var1_val in var1_vals:
				arr = np.array([[res_TOV['n'], res_TOV['r_max']] for res_TOV in results_TOV[term][var1_val]])
				plot_results_TOV[term][var1_val] = arr[arr[:,0].argsort()]
		
		if var1_name=="p_init":
			var1_name=r'$p_0$'
		else:
			var1_name=r'A'
		
		# Initialise plot with right size
		cm = 1/2.54
		plt.figure(figsize=[16*cm,9*cm])
		
		# Create the plot axis
		fig, axs = plt.subplots()
		for i, var1_val in enumerate(var1_vals):
			axs.plot(plot_results_LE[var1_val][:,0], plot_results_LE[var1_val][:,1], linestyle=self.curve_styles[2*i], c=self.curve_colours[0],  label="LE for " + var1_name + "=" + str(var1_val))
			for term in terms:
				if len(terms)==1:
					label = "TOV for " + var1_name + "=" + str(var1_val)
				else:
					label = "TOV with term=" + str(term) + " for " + var1_name + "=" + str(var1_val)
				axs.plot(plot_results_TOV[term][var1_val][:,0], plot_results_TOV[term][var1_val][:,1], linestyle=self.curve_styles[2*i+1], c=self.curve_colours[0], label=label)
		axs.set_yscale('log')
		
		axs.set_ylabel(r'$r_ 0$')
		axs.set_ylim([ylim_low,ylim_up])
# 		axs.set_xlim([xlim_low,xlim_up])
		axs.set_xlabel(r'$n=\frac{1}{\gamma-1}$')
		
		# Resize the whole plot to fit the legend next to it.
# 		box = axs.get_position()
# 		axs.set_position([box.x0, box.y0, box.width * 0.8, box.height])
		# Create the legend right of plot
# 		axs.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		axs.legend()
		plt.tight_layout()
		
		if var2_name=="p_init":
			var2_name=r'$p_0$'
		else:
			var2_name=r'A'
		plt.title("$r_0$ for fixed " + var2_name + r'$=' + str(var2_fixed_val) + r'$')
# 		plt.title("$xi_0$ for different initial parameters")
		plt.savefig('pictures/TOV-Exponents-LESubs-InitialVals-Database-PlotResults.svg')
		plt.show()
		
		matplotlib.use("pgf")
		matplotlib.rcParams.update({
# 		    "pgf.texsystem": "pdflatex",
		    'font.family': 'serif',
		    'text.usetex': True,
		    'pgf.rcfonts': False,
		})
		if noExport==False:
 			plt.savefig("pictures/TOV-Exponents-LESubs-InitialVals-Database-PlotResults.pgf", dpi=1000, bbox_inches='tight')
 			print("Saved Plot under pictures/TOV-Exponents-LESubs-InitialVals-Database-PlotResults.pgf \n")
 			tikzplotlib.save("pictures/TOV-Exponents-LESubs-InitialVals-Database-PlotResults.tex")
 			print("Saved Plot under pictures/TOV-Exponents-LESubs-InitialVals-Database-PlotResults.tex \n")

Plodda = Plotter()

Plodda.plotInitialVals(var1_name='p_init', var1_vals=[0.1,1], var2_name='A_init', var2_fixed_val=1, terms=[0], xlim_low=0, ylim_up=None, noExport=False)