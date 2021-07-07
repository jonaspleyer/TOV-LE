#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Only use this to import the package DiffEqSolver 
# from ../Solvers/ correctly
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from Solvers.Solver import DiffEqSolver
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import matplotlib

# We create a subclass if the DiffEqSolver Class to still have all functions needed
# Furthermore we add a function to solve and plot the TOV equation for different given terms
# and the Lane-Emden equation as well.
class Plotter(DiffEqSolver):
	def solveAndPlotResults(self, r0, u0, p0, R, rend, dr):
		# Define the linestyles to be used for plotting
		linestyles = ['-', '--', '-.', ':']
		
		# Define result arrays for different amount of terms present in equation
		results = [[]]*4
		results_small = [[]]*4
		succ = [[]]*4
		r_max = [[]]*4
		
		# First we need to solve the equations
		for i in range(4):
			results[i], results_small[i], succ[i], r_max[i] = self.solveTOV(r0, u0, p0, R, rend, dr, terms=i)
		# Check if the solving was successful
		
		# Initialise plot with right size
		cm = 1/2.54
		plt.figure(figsize=[14.755060*cm,11.066295*cm])
		
		for i in range(4):
			if succ[i] == True:
				print("Solving TOV for i=" + str(i) + " terms")
				
				# Plot pressure
				plt.subplot(2,2,1)
				if i == 1:
					plt.title(r'Pressure $p(r)$')
				plt.plot(results[i][:, 0], results[i][:, 2], label=r'$i=$' + str(3-i), linestyle=linestyles[i], c='black')
				
				# Plot mass
				plt.subplot(2,2,3)
				if i == 1:
					plt.title(r'Mass $m(r)$')
				plt.plot(results[i][:, 0], results[i][:, 1], label=r'$i=$' + str(3-i), linestyle=linestyles[i], c='black')
				
				# Plot density
				plt.subplot(2,2,2)
				if i == 1:
					plt.title(r'Density $\rho(r)$')
					
				# Calculate the density via the eos
				plt.plot(results[i][:,0], [self.eos(x[2], x[0]) for x in results[i]], label=r'$i=$' + str(3-i), linestyle=linestyles[i], c='black')
			else:
				print("Solving was not possible for i=" + str(i))
		for i in range(3):
			plt.subplot(2,2,i+1)
			plt.legend()
		
		# Comppare With Lane Emden
		plt.subplot(2,2,4)
		
		# This contains all the results for the LE results
		# in form [[r,m,p,rho],...]
		results_LE, succ, r_max = self.convertSolveLE(r0, u0, p0, R, rend, dr)
		
		print("Solving LE")
		# results_LE, succ_LE, xi_max = self.solveLE(0, 1, 0 ,xi_end, dxi)
		
		# Plot pressure LE
		# plt.plot(p_results_LE_transf[:, 0], p_results_LE_transf[:, 1], label=r'$p_{LE}$', linestyle=linestyles[3], c='red')
		
		# Get reduced range of TOV results (LE will not be solvable for as long as TOV Result)
		results_reduced_TOV = np.array([res for res in results[3] if res[0]<= r_max*1.2])
		
		# Interpolate for the pressure of TOV and LE
		# to be able to calculate the difference between p_TOV and P_LE
		p_TOV_interpolate = interp1d(results_reduced_TOV[:,0],results_reduced_TOV[:,1], kind='cubic')
		p_LE_interpolate = interp1d(results_LE[:, 0], results_LE[:, 1], kind='cubic')
		
		# Create a linspace for plotting
		r_range = np.linspace(r0,min(results_reduced_TOV[:,0][-1],results_LE[:,0][-1]))
		# Calculate the difference
		difference_p_TOV_LE = [p_LE_interpolate(r)-p_TOV_interpolate(r) for r in r_range]
		# Plot the difference
		plt.plot(r_range, difference_p_TOV_LE, label=r'$p_{LE}-p_{TOV,0}$', linestyle=linestyles[1], c='black')
		
		plt.legend()
		plt.title("TOV$_0$ and LE results")
		
		# General configuration of the total plot
		plt.tight_layout()
		plt.subplots_adjust(hspace=0.3)
		plt.savefig('pictures/TOV-Terms.svg')
		plt.show()
		matplotlib.use("pgf")
		matplotlib.rcParams.update({
# 		    "pgf.texsystem": "pdflatex",
		    'font.family': 'serif',
		    'text.usetex': True,
		    'pgf.rcfonts': False,
		})
		plt.savefig("pictures/TOV-Terms.pgf", dpi=1000, bbox_inches='tight')

# Define initial values
r0 = 0
u0 = 0
p0 = 0.5
R  = 2.5
rend = R
dr = 0.001

# Create an instance of the Solver with a particular EOS
# n=0 Rocky Planets
# n=0.5,1 Neutron stars
# n=1.5 fully convective star cores (red giants, brown dwarfs, giant gas planets)
# n=3 wight dwarfs
# n=5 stellar system
# For more see https://en.wikipedia.org/wiki/Polytrope
n = 3
gamma = 1 + 1/n
A = 2
Solver = Plotter(gamma, A)

# Solve and plot the results
Solver.solveAndPlotResults(r0,u0,p0,R,rend,dr)

# Also calulate the mass
m, r_max, succ = Solver.getMass(r0,u0,p0,R,rend,dr)

# Print results
if succ == True:
	print("The Total Mass is " + str(round(m,4)))
	print("The ratio M/R is " + str(round(m/r_max,4)))
	print("The theoretical limit is " + str(round(4/9,4)))
