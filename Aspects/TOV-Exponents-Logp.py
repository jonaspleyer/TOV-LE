#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Only use this to import the package DiffEqSolver 
# from ../Solvers/ correctly
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from Solvers.SolverLogp import DiffEqSolverLogp
import matplotlib.pyplot as plt
import time
import numpy as np
import multiprocessing

class Plotter(DiffEqSolverLogp):
	def solveMultiprocExponents(self, n_0, n_max, n_step, r0, u0, p0, R, rend, dr, suppressOutput=False, N_threads = 13):
		#
		self.exponent_vals = np.arange(n_0,n_max+n_step,n_step)
		
		# max_r_results = multiprocessing.Array('l',[[0,False,0,0] for i in range(0,len(self.exponent_vals))])
		self.r_maxes = multiprocessing.Array('f',[0 for i in range(0,len(self.exponent_vals))])
		self.xi_maxes = multiprocessing.Array('f',[0 for i in range(0,len(self.exponent_vals))])
		self.r_Bool = multiprocessing.Array('l',[0 for i in range(0,len(self.exponent_vals))])
		self.xi_Bool = multiprocessing.Array('l',[0 for i in range(0,len(self.exponent_vals))])
		self.times = multiprocessing.Array('f',[0 for i in range(0,len(self.exponent_vals))])
		
		# Initialise threads
		processes = []
		
		# Split process in several subthreads
		# and distribute load equally over all threads
		for i in range(0,N_threads):
			# Helper array to store the threads
			array = []
			j = 0
			# This distribution is chosen such that all the values with high exponent n
			# (for which the solving takes longer) are distributed equally about the different threads
			while N_threads*j+i < len(self.exponent_vals):
				array.append(N_threads*j + i)
				j += 1
			# Generate the process
			x = multiprocessing.Process(target=self.solveForexponent_vals, args=(array, r0, u0, p0, R, rend, dr,))
			# Store it in the x array
			processes.append(x)
			x.start()
		
		# Join the processes to have them end continue the method at the time when all processes have finished
		for process in processes:
			process.join()
		
		# Create lists to store values for plotting
		# We differentiatie between values for which the solution 
		# was successfull and ones where it was not
		r_plot_zero_values = []
		xi_plot_zero_values = []
		r_plot_exponent_values = []
		xi_plot_exponent_values = []
		
		for i in range(0,len(self.exponent_vals)):
			if suppressOutput == False:
				print(str(self.exponent_vals[i]) + "  " + str(self.r_maxes[i]) + "  " + str(self.r_Bool[i]) + "  " + str(self.xi_maxes[i]) + " " + str(self.xi_Bool[i]))
			if self.r_Bool[i] == 1:
				r_plot_zero_values.append(self.r_maxes[i])
				r_plot_exponent_values.append(self.exponent_vals[i])
			if self.xi_Bool[i] == 1:
				xi_plot_zero_values.append(self.xi_maxes[i])
				xi_plot_exponent_values.append(self.exponent_vals[i])
		
		# Stop the timer
		end = time.time()
		
		# Print Duration
		print("===== Finished Process =====")
		print("Duration was " + str(round(end-start,3)) + " Seconds for " + str(len(self.exponent_vals)) + " Test Samples with " + str(N_threads) + " threads")
		
		# Plot the results
		# Create figure with right dimensions
		plt.figure(figsize=[6.4,4])
		plt.plot(r_plot_exponent_values, r_plot_zero_values, label=r'$r_0$', linestyle="-", c='black')
		plt.plot(xi_plot_exponent_values, xi_plot_zero_values, label=r'$\xi_0[r]$', linestyle="-.", c='black')
		plt.legend()
		plt.title(r'$r_0$ where $p(r_0)=0$')
		plt.xlabel(r"Exponent $n=\frac{1}{\gamma-1}$")
		# plt.ylabel(r"$\xi_0")
		plt.yscale('log')
		# Plot a vertical line at xi=5 with the correct height
		plt.vlines(5,0,max(r_plot_zero_values), colors='k', linestyle="--")
		# Save the plot to a file
		plt.savefig("pictures/TOV-Exponents-Logp.svg")
		plt.show()
		
	def solveForexponent_vals(self, array, r0, u0, p0, R, rend, dr):
		for i in array:
			self.times[i] = time.time()
			results_TOV, results_TOV_small, succ_TOV, r_max = Solver.solveTOV(r0, u0, p0, R, rend, dr, terms=0, exponent=Solver.exponent_vals[i], suppressWarning=True)
			results_LE, succ_LE, xi_max = Solver.convertSolveLE(r0, u0, p0, R, rend, dr, exponent=Solver.exponent_vals[i], suppressWarning=True, suppressOutput=True)
			self.r_maxes[i] = r_max
			self.xi_maxes[i] = xi_max
			self.r_Bool[i] = 1 if r_max < rend else 0
			self.xi_Bool[i] = 1 if xi_max < rend else 0
			self.times[i] = time.time()-self.times[i]
		
# Keep track of total time spent
start = time.time()

# Create an instance of the Solver with polytropic EOS
n = 2
gamma = 1+1/n
A = 5

# Initialise Solver with arbitrary values for A and gamma (will not be used)
Solver = Plotter(gamma, A)

# Define initial values
r0 = 0
u0 = 0.0
p0 = 1
R  = 50
rend = R
dr = 0.01

# Define the range of exponents to solve for
n_0 = 0.01
n_max = 5.01
n_step = 0.01

# HINT for choosing N_threads
# CPU: AMD Ryzen 3700X (8 Cores, 16Threads)
# r0=0, u0=0, p0=1, R=100, rend=5000, dr=0.01, n0 = 0.01, n_max = 4.01, n_step = 0.01
# Threads Time R=100  Time R=500   Time R=5000
# 1       45.88s      231.37s      
# 2       23.57s      123.12s
# 4       12.75s       73.09s
# 8        6.74s       64.60s      
# 16       5.95s      108.62s
# REMARK: The following variables have the most effect:
# R  > 100
# dr < 0.05
# REMARK2: For good performance/effort choose
# R  = 100
# dr = 0.01
# N_threads = maximum-1
# Every integer value >=1 is possible
print("===== Starting Process =====")
Solver.solveMultiprocExponents(n_0, n_max, n_step, r0, u0, p0, R, rend, dr, suppressOutput=False, N_threads=8)