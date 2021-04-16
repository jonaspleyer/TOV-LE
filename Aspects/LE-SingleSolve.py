#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Only use this to import the package DiffEqSolver 
# from ../Solvers/ correctly
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from Solvers.Solver import DiffEqSolver
import matplotlib.pyplot as plt
import matplotlib

class Plotter(DiffEqSolver):
	def solveAndPlotResults(self, xi0, T0, dT0, xi_max, dxi, exponents, suppressWarning=False):
		# Define linestyles for later plotting
		linestyles = [
			(0, (2, 1)),
			(0, (2, 4)),
			(0, (2, 8))]
		
		# Define arrays that store the information for different runs
		results = [[]]*len(exponents)
		succ = [[]]*len(exponents)
		xi_end = [[]]*len(exponents)
		
		for i, exponent in enumerate(exponents):
			results[i], succ[i], xi_end[i] = self.solveLE(xi0, T0, dT0, xi_max, dxi, exponent=exponent, suppressWarning=suppressWarning)
			# Check if the solving was successful
			if succ[i] == True:
				# Plot Theta
				plt.plot(results[i][:, 0], results[i][:, 1], label=r'$n=$' + str(exponent), linestyle=linestyles[i], c='black')
			else:
				print("Solving was not possible for $n=$"+str(exponent))
		plt.legend()
		plt.ylabel(r'LE Solution $\theta$')
		plt.xlabel(r'Radial Coordinate $\xi$')
		plt.savefig('pictures/LE-SingleSolve.svg')
		plt.show()
		matplotlib.use("pgf")
		matplotlib.rcParams.update({
# 		    "pgf.texsystem": "pdflatex",
		    'font.family': 'serif',
		    'text.usetex': True,
		    'pgf.rcfonts': False,
		})
		plt.savefig("pictures/LE-SingleSolve.pgf", dpi=1000, bbox_inches='tight')


# Define inital values
xi0 = 0
T0 = 1
dT0 = 0
xi_max = 10
dxi = 0.01

# Create an instance of the Solver with polytropic EOS
# Both of these values do not matter
n = 2
gamma = 1+1/n
A = 1

# Define a range of exponents to plot functions for (maximum 4)
exponents = [0,1,5]

# Initialise Solver
Solver = Plotter(gamma, A)

Solver.solveAndPlotResults(xi0, T0, dT0, xi_max, dxi, exponents, suppressWarning=False)