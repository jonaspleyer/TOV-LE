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
from scipy.special import kv
from scipy.interpolate import interp1d
import numpy as np

class Plotter(DiffEqSolver):
	def solveAndPlotResults(self, r0, u0, p0, R, rend, dr):
		# First we need to define the function that we want to invert
		f_alpha = lambda a: self.factor/(kv(2,a)*a**2)*np.exp(-a*(kv(1,a)+kv(3,a))/(2*kv(2,a)))
		# Now we need to try and estimate what the starting value will be.
		# First check if the current initial value is even obtainable. Since f_alpha is monotously decreasing
		# We can check the value at a very small value. If this value is not larger than 
		a0 = 0.000000001
# 		if f_alpha(a0) < p0:
#  			return False
		
		# Define the x- and y-vals for interpolation
		y_vals = np.linspace(a0,100)
		x_vals = f_alpha(y_vals)
		# interpolate the inverse by switching y and x 
		alpha = interp1d(y_vals, x_vals)
		# Define the EOS
		eos_new = lambda p,r:p*(1+alpha(p)*(kv(1,alpha(p))+kv(3,alpha(p)))/(2*kv(2,alpha(p))))
		
		# Choose the multiplication factor such that polytropic and rel eos results
		# can be compared later on
		
		self.factor = eos_new(p0,0)/p0**(1/self.gamma)
		self.eos = lambda p,r: self.factor*p**(1/gamma) if p > 0 else 0
		
		# Solve the equations for the polytropic eos
		results_1, results_small_1, succ_1, r_max_1 = self.solveTOV(r0, u0, p0, R, rend, dr)
		
		# Solve the equations for the rel EOS
		self.eos = eos_new
		results_2, results_small_2, succ_2, r_max_2 = self.solveTOV(r0, u0, p0, R, rend, dr)
		
		# Plot the results
		# Create a new figure with the right dimesions
		cm = 1/2.54
		plt.figure(figsize=[16*cm,12*cm])
		
		# Check if the solving was successful
		if succ_1 == True and succ_2 == True:
			# Plot pressure
			plt.subplot(2,2,1)
			plt.plot(results_1[:, 0], results_1[:, 2], label=r'$p(r)$ polEOS', linestyle='-', c='black')
			plt.plot(results_2[:, 0], results_2[:, 2], label=r'$p(r)$ rel', linestyle='--', c='black')
			plt.legend()
			plt.title("Pressure")
			
			# Plot u(r)
			plt.subplot(2,2,2)
			plt.plot(results_1[:, 0], results_1[:, 1], label=r'$m(r)$ pol', linestyle='-', c='black')
			plt.plot(results_2[:, 0], results_2[:, 1], label=r'$m(r)$ rel', linestyle='--', c='black')
			plt.legend()
			plt.title("Mass")
			
			# Plot density
			plt.subplot(2,2,3)
			plt.plot(results_1[:,0], results_1[:,3], label=r'$\rho (r)$ pol', linestyle='-', c='black')
			plt.plot(results_2[:,0], results_2[:,3], label=r'$\rho (r)$ rel', linestyle='--', c='black')
			plt.legend()
			plt.title("Density")
			
			# Plot EOS
			plt.subplot(2,2,4)
			p_vals = np.linspace(0.001*p0,p0)
			plt.plot(p_vals,self.eos(p_vals,0), label=r'$\rho(p)$',c='black')
			plt.legend()
			plt.title("Equation of State")
			
			plt.subplots_adjust(wspace=0.5, hspace=0.5)
			plt.show()
		else:
			print("Solving was not possible.")

# Define initial values
r0 = 0
u0 = 0
p0 = 0.4
R  = 1.5
rend = R
dr = 0.05

# Create an instance of the Solver with polytropic EOS
n = 2
gamma = 1+1/n
A = 0.01

# Initialise Solver
Solver = Plotter(gamma, A)

Solver.solveAndPlotResults(r0, u0, p0, R, rend, dr)