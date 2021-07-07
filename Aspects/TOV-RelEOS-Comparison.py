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
from scipy.special import kv
from scipy.interpolate import interp1d
import numpy as np

class Plotter(DiffEqSolver):
	def solveAndPlotResults(self, r0, u0, p0, R, rend, dr):
		# First we need to define the function that we want to invert
# 		f_alpha = lambda a: self.factor/(kv(2,a)*a**2)*np.exp(-a*(kv(1,a)+kv(3,a))/(2*kv(2,a)))
		p_alpha = lambda a: p0*self.factor/(kv(2,a)*a**2)*np.exp(-a*(kv(1,a)+kv(3,a))/(2*kv(2,a))) if a > 0 else self.factor/(np.exp(1)**2*2)
		# Now we need to try and estimate what the starting value will be.
		# First check if the current initial value is even obtainable. Since f_alpha is monotously decreasing
		# We can check the value at a very small value. If this value is not larger than 
		
		a_max = 1
		while p_alpha(a_max) > 0.001*p0:
			a_max += 1
		a_max += -1
		
		# Define the x- and y-vals for interpolation
		a_vals = np.arange(0,a_max,0.1)
		p_vals = [p_alpha(x) for x in a_vals]
		# interpolate the inverse by switching y and x 
		alpha = interp1d(p_vals, a_vals)
		# Define the EOS
		eos_new = lambda p, r:p*(1+alpha(p)*(kv(1,alpha(p))+kv(3,alpha(p)))/(2*kv(2,alpha(p))))
		
		# Choose the multiplication factor such that polytropic and rel eos results
		# can be compared later on
		
		A = eos_new(p0,0)/p0**(1/self.gamma)
		self.eos = lambda p,r: A*p**(1/self.gamma) if p > 0 else 0
		
		# Solve the equations for the polytropic eos
		results_1, results_small_1, succ_1, r_max_1 = self.solveTOV(r0, u0, p0, R, rend, dr)
		
		# Solve the equations for the rel EOS
		self.eos = eos_new
		results_2, results_small_2, succ_2, r_max_2 = self.solveTOV(r0, u0, p0, R, rend, dr)
		
		# Plot the results
		# Create a new figure with the right dimesions
		cm = 1/2.54
		plt.figure(figsize=[14.755060*cm,11.066295*cm])
		
		# Check if the solving was successful
		if succ_1 == True and succ_2 == True:
			# Plot pressure
			plt.subplot(2,2,1)
			plt.plot(results_1[:, 0], results_1[:, 2], label=r'$p(r)$ pol', linestyle='-', c='black')
			plt.plot(results_2[:, 0], results_2[:, 2], label=r'$p(r)$ rel', linestyle='--', c='black')
			plt.legend()
# 			plt.title("Pressure")
			plt.ylabel(r'Pressure $p$')
			plt.xlabel(r'Radius $r$')
			
			# Plot u(r)
			plt.subplot(2,2,2)
			plt.plot(results_1[:, 0], results_1[:, 1], label=r'$m(r)$ pol', linestyle='-', c='black')
			plt.plot(results_2[:, 0], results_2[:, 1], label=r'$m(r)$ rel', linestyle='--', c='black')
			plt.legend()
# 			plt.title("Mass")
			plt.ylabel(r'Mass $m$')
			plt.xlabel(r'Radius $r$')
			
			# Plot density
			plt.subplot(2,2,3)
			plt.plot(results_1[:,0], results_1[:,3], label=r'$\rho (r)$ pol', linestyle='-', c='black')
			plt.plot(results_2[:,0], results_2[:,3], label=r'$\rho (r)$ rel', linestyle='--', c='black')
			plt.legend()
# 			plt.title("Density")
			plt.ylabel(r'Density $\rho$')
			plt.xlabel(r'Radius $r$')
			
			# Plot EOS
			plt.subplot(2,2,4)
			f = lambda p,r: A*p**(1/self.gamma) if p > 0 else 0
			p_vals = [p_val for p_val in p_vals[1:] if p_val <= p0]
			plt.plot(p_vals,[f(p_val,0) for p_val in p_vals], label=r'$\rho(p)$ pol',linestyle='-', c='black')
			plt.plot(p_vals,self.eos(p_vals,0), label=r'$\rho(p)$ rel',linestyle='--', c='black')
			
			locs, labels = plt.xticks()
			m = max(locs[1:-1])
			labels = [str(round(loc/m,2)) for loc in locs[1:-1]]
			plt.xticks(locs[1:-1], labels)
			plt.xlabel(r'Pressure $p/p_0$')
			
			locs, labels = plt.yticks()
			m = max(locs[1:-1])
			labels = [str(round(loc/m,2)) for loc in locs[1:-1]]
			plt.yticks(locs[1:-1], labels)
			plt.ylabel(r'Density $\rho/\rho_0$')
			
			plt.legend()
# 			plt.title("Equation of State")
			
# 			plt.subplots_adjust(wspace=0.35, hspace=0.35)
			plt.tight_layout()
			plt.savefig("pictures/TOV-RelEOS-Comparison.svg", dpi=1000, bbox_inches='tight')
			plt.show()
			matplotlib.use("pgf")
			matplotlib.rcParams.update({
			    "pgf.texsystem": "pdflatex",
			    'font.family': 'serif',
			    'text.usetex': True,
			    'pgf.rcfonts': False,
			})
			plt.savefig("pictures/TOV-RelEOS-Comparison.pgf", dpi=1000, bbox_inches='tight')
		else:
			print("Solving was not possible.")

# Define initial values
r0 = 0
u0 = 0
p0 = 0.4
R  = 1.5
rend = R
dr = 0.001

# Create an instance of the Solver with polytropic EOS
n = 1
gamma = 1+1/n
A = 20

# Initialise Solver
Solver = Plotter(gamma, A)

Solver.solveAndPlotResults(r0, u0, p0, R, rend, dr)
