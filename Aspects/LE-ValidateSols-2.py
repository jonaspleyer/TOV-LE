#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Only use this to import the package DiffEqSolver 
# from ../Solvers/ correctly
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from Solvers.Solver import DiffEqSolver
from Standards import PlottingStandards as standards
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

class Plotter(DiffEqSolver):
	def solveAndPlotResults(self, xi0, T0, dT0, xi_max, dxi, suppressWarning=False):
		# Known exact solutions are for n=0,1,5
		exponents = [0,1,5]
		solutions = [lambda x:1-x**2/6, lambda x:np.sin(x)/x if x>0 else 1, lambda x: 1/(1+1/3*x**2)**(1/2)]
		
# 		steps = [0.5,0.1,0.05,0.01,0.005]
		steps = [1,0.75,0.5,0.25,0.1,0.075,0.05,0.025,0.01,0.0075,0.005,0.0025,0.001,0.00075]
		
		# Define arrays that store the information for different runs
		diffs = [[[]]*len(steps)]*len(exponents)
		
		# Initialise plot with right size
		cm = 1/2.54
		plt.figure(figsize=[16*cm,9*cm])
		
		for i, step in enumerate(steps):
			for k, exponent in enumerate(exponents):
				results, succ, xi_end = self.solveLE(xi0, T0, dT0, xi_max, step, exponent=exponent, suppressWarning=suppressWarning)
				# Check if the solving was successful
				if succ == True:
# 					diffs[k].append(
 					diffs[k][i] = max([abs(results[:,1][j]-solutions[k](results[:,0][j])) for j in range(len(results[:,0]))])
				else:
					print("Solving was not possible for $n=$"+str(exponent))
		plt.loglog()
		for i in range(len(exponents)):
			plt.plot(steps, diffs[i], c='k', linestyle=standards.linestyles[i], label="$\Delta$ for $n=$"+str(i))
		plt.legend()
		plt.savefig('pictures/LE-ValidateSols-2.svg')
		plt.show()
		matplotlib.use("pgf")
		matplotlib.rcParams.update({
# 		    "pgf.texsystem": "pdflatex",
		    'font.family': 'serif',
		    'text.usetex': True,
		    'pgf.rcfonts': False,
		})
		plt.savefig("pictures/LE-ValidateSols-2.pgf", dpi=1000, bbox_inches='tight')


# Define inital values
xi0 = 0
T0 = 1
dT0 = 0
xi_max = 100
dxi = 0.03

# Create an instance of the Solver with polytropic EOS
# Both of these values do not matter
n = 2
gamma = 1+1/n
A = 1

# Initialise Solver
Solver = Plotter(gamma, A)

Solver.solveAndPlotResults(xi0, T0, dT0, xi_max, dxi, suppressWarning=False)