#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import tikzplotlib

# Only use this to import the package DiffEqSolver 
# from ../Solvers/ correctly
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from Solvers.Solver import DiffEqSolver

class Plotter(DiffEqSolver):
	def solveAndPlotResults(self, r0, u0, p0, R, rend, dr, noExport=True):
		# Print information
		print("===== Parameter Values =====")
		print("gamma     = " + str(self.gamma))
		print("A         = " + str(self.factor))
		print("p0        = " + str(p0))
		
		# Solve the equations seperately and store the result
		results_TOV, results_TOV_small, succ, r_max_TOV = self.solveTOV(r0, u0, p0, R, rend, dr)
		results_LE, succ, r_max_LE = self.convertSolveLE(r0, u0, p0, R, rend, dr)

		# Plot the results
		# Create a new figure with the right dimesions
		cm = 1/2.54
		plt.figure(figsize=[16*cm,12*cm])
		
		# Create a subplot for the different pressures and for the mass and density
		# Subplot for pressure
		plt.subplot(2,2,1)
# 		plt.title(r'Pressure $p(r)$')
		plt.plot(results_TOV[:, 0], results_TOV[:, 2], label=r'$p_{TOV}$', linestyle='-', c='black')
		plt.plot(results_LE[:, 0], results_LE[:, 2], label=r'$p_{LE}$', linestyle=':', c='black')
		plt.ylabel(r'Pressure $p$')
		plt.xlabel(r'Radius $r$')
		plt.legend(loc='upper right')
		
		# Next subplot for Density
		plt.subplot(2,2,2)
# 		plt.title(r'Density $\rho(r)$')
		y_vals_TOV = [results_TOV[:,1][i]/results_TOV[:,0][i]**3/4/np.pi*3 if i >= 1 else np.nan for i in range(0,len(results_TOV[:,0]))]
		y_vals_LE  = [results_LE[:,1][i]/results_LE[:,0][i]**3/4/np.pi*3 if i >= 1 else np.nan  for i in range(0,len(results_LE[:,0]))]
		plt.plot(results_TOV[:,0], results_TOV[:,3], label=r'$\rho_{TOV}$', linestyle='-', c='black')
		plt.plot(results_TOV[:,0], y_vals_TOV, label=r'$\bar{\rho}_{TOV}$', linestyle='--', c='black')
		plt.plot(results_LE[:,0], results_LE[:,3], label=r'$\rho_{LE}$', linestyle=':', c='black')
		plt.plot(results_LE[:,0], y_vals_LE, label=r'$\bar{\rho}_{LE}$', linestyle=(0,(1,5)), c='black')
		plt.ylabel(r'Density $\rho$')
		plt.xlabel(r'Radius $r$')
		plt.legend(loc='upper right')
		
		# Next subplot for Mass
		plt.subplot(2,2,3)
# 		plt.title(r'Mass $m(r)$')
		plt.plot(results_TOV[:, 0], results_TOV[:, 1], label=r'$m_{TOV}$', linestyle='-', c='black')
		plt.plot(results_LE[:,0],results_LE[:,1], label='$m_{LE}$', linestyle=':', c='black')
		plt.ylabel(r'Mass $m$')
		plt.xlabel(r'Radius $r$')
		plt.legend(loc='upper left')
		
		# Next subplot for m(r)/r**3
		plt.subplot(2,2,4)
# 		plt.title(r'Mass radius ratio $m/r$')
		# Calculate new vals and normalise with self.eos(p0,0)
		y_vals_TOV = [results_TOV[:,1][i]/results_TOV[:,0][i] if i >= 1 else np.nan for i in range(0,len(results_TOV[:,0]))]
		y_vals_LE  = [results_LE[:,1][i]/results_LE[:,0][i] if i >= 1 else np.nan  for i in range(0,len(results_LE[:,0]))]
		plt.plot(results_TOV[:,0], y_vals_TOV, label=r'$m_{TOV}/r$', linestyle='-', c='black')
		plt.plot(results_LE[:,0], y_vals_LE, label=r'$m_{LE}/r$', linestyle=':', c='black')
		plt.ylabel(r'Mass radius ratio $m/r$')
		plt.xlabel(r'Radius $r$')
		plt.legend(loc='upper right')

		# Save the total picture
		plt.tight_layout()
		plt.subplots_adjust(hspace=0.3)
		plt.savefig('pictures/TOV-LE-Combi.svg')
		print("Saved Plot under pictures/TOV-LE-Combi.svg \n")
		plt.show()
		
		matplotlib.use("pgf")
		matplotlib.rcParams.update({
# 		    "pgf.texsystem": "pdflatex",
		    'font.family': 'serif',
		    'text.usetex': True,
		    'pgf.rcfonts': False,
		})
		if noExport==False:
			plt.savefig("pictures/TOV-LE-Combi.pgf", dpi=1000, bbox_inches='tight')
			print("Saved Plot under pictures/TOV-LE-Combi.pgf \n")
			tikzplotlib.save("pictures/TOV-LE-Combi.tex")
			print("Saved Plot under pictures/TOV-LE-Combi.tex \n")
		
# Values of interest are mainly between 1.5 and 3 since those correspond to gamma = 1+1/n = 4/3, 5/3
n = 3
gamma = 1+1/n
A = 2

# Create instance of solver with exponent
Solver = Plotter(gamma,A)

# Set initial values
r0 = 0
u0 = 0
p0 = 0.5
R = 2.5
rend = R
dr = 0.01

# Solve and plot results
Solver.solveAndPlotResults(r0,u0,p0,R,rend,dr,noExport=False)