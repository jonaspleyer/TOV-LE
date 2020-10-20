#!/bin/spyder

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.misc import derivative
import scipy.integrate as integrate

# Only use this to import the package DiffEqSolver 
# from ../Solvers/ correctly
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from Solvers.Solver import DiffEqSolver


class Plotter(DiffEqSolver):
	def solveAndPlotResults(self, r0, u0, p0, R, rend, dr):
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
		plt.figure(figsize=[6.4,8])
		
		# Create a subplot for the different pressures and for the mass and density
		plt.subplot(2,1,1)
		
		# Plot pressure TOV and LE
		plt.plot(results_TOV[:, 0], results_TOV[:, 2], label=r'Pressure $p_{TOV}$', linestyle='-', c='black')
		plt.plot(results_LE[:, 0], results_LE[:, 2], label=r'Pressure $p_{LE}$', linestyle=':', c='black')
		
		# Create a legend
		plt.legend()
		
		# Create the next subplot
		plt.subplot(2,1,2)
		
		# Plot the density for the LE and TOV solution
		plt.plot(results_LE[:,0], results_LE[:,3], label=r'$\rho_{TOV}$', linestyle='--', c='black')
		plt.plot(results_TOV[:,0], results_TOV[:,3], label=r'$\rho_{LE}$', linestyle='-.', c='black')
		
		# Plot m(r) TOV and LE
		plt.plot(results_TOV[:, 0], results_TOV[:, 1], label=r'$m_{TOV}$', linestyle='-', c='black')
		plt.plot(results_LE[:,0],results_LE[:,1], label='$m_{LE}$', linestyle=':', c='black')

		# Plot the legend
		plt.legend()

		# Save the total picture
		plt.savefig('pictures/TOV-LE-Combi.svg')
		print("Saved Plot under pictures/TOV-LE-Combi.svg \n")
		plt.show()
		
# Values of interest are mainly between 1.5 and 3 since those correspond to gamma = 1+1/n = 4/3, 5/3
n = 3
gamma = 1+1/n
A = 2

# Create instance of solver with exponent
Solver = Plotter(gamma,A)

# Set initial values
r0 = 0
u0 = 0
p0 = 1
R = 1.2
rend = R
dr = 0.01

# Solve and plot results
Solver.solveAndPlotResults(r0,u0,p0,R,rend,dr)