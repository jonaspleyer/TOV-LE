#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Only use this to import the package DiffEqSolver 
# from ../Solvers/ correctly
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from Solvers.SolverLESubs import DiffEqSolverLESubs
import matplotlib.pyplot as plt
import time
import numpy as np
import multiprocessing
import ctypes

class Plotter(DiffEqSolverLESubs):
	def solveMultiprocExponents(self, n_0, n_max, n_step, r0, u0, p0, R, rend, dr, suppressOutput=False, N_threads = 2, N_terms=1):
		# Define linestyles for plotting later
		linestyles = ['-', '--', '-.', ':']
		
		# Define the array which contains all exponents n = 1/(gamma-1)
		self.exponent_vals = np.arange(n_0,n_max+n_step,n_step)
		N_exponents = len(self.exponent_vals)
		
		# Define multiprocessing arrays that can be accessed even if processes run in parallel
		self.r_maxes_all = multiprocessing.Array(ctypes.c_float, N_exponents*(N_terms+1))
		self.r_Bool_all = multiprocessing.Array(ctypes.c_float, N_exponents*(N_terms+1))
		
		# for j in range(0,N_terms):
		# 	self.r_maxes.append(multiprocessing.Array('f',[0 for i in range(0,len(self.exponent_vals))]))
		#	self.r_Bool.append(multiprocessing.Array('l',[0 for i in range(0,len(self.exponent_vals))]))
		self.xi_maxes = multiprocessing.Array('f',[0 for i in range(0,N_exponents)])
		self.xi_Bool = multiprocessing.Array('l',[0 for i in range(0,N_exponents)])
		self.times = multiprocessing.Array('f',[0 for i in range(0,N_exponents)])
		
		# Initialise threads; Create list that contains all processes
		processes = []
		
		# Split process in several subthreads
		# and distribute load equally over all threads
		for i in range(0,N_threads):
			# Helper array to store the threads
			array = []
			j = 0
			# This distribution is chosen such that all the values with high exponent n
			# (for which the solving takes longer) are distributed equally about the different threads
			while N_threads*j+i < N_exponents:
				if j%2 == 0:
					array.append(N_threads*j + i)
				if j%2 == 1 and N_threads*(j+1)-(i+1) < N_exponents:
					array.append(N_threads*(j+1)-(i+1))
				j += 1
			# Generate the process
			x = multiprocessing.Process(target=self.solveForexponent_vals, args=(array, r0, u0, p0, R, rend, dr, N_terms,))
			# Store it in the x array
			processes.append(x)
			x.start()
		
		# Join the processes to have them end continue the method at the time when all processes have finished
		for process in processes:
			process.join()
		
		# Create lists to store values for plotting
		# We differentiatie between values for which the solution 
		# was successfull and ones where it was not
		r_plot_zero_values = [[]]*(N_terms+1)
		xi_plot_zero_values = []
		r_plot_exponent_values = [[]]*(N_terms+1)
		xi_plot_exponent_values = []
		
		for i in range(0,N_exponents):
			if suppressOutput == False:
					# print(str(round(self.exponent_vals[i],3)) + "  " + str(round(self.r_maxes[i],3)) + "  " + str(self.r_Bool[i]) + "  " + str(round(self.xi_maxes[i],3)) + " " + str(round(self.xi_Bool[i],3)))
					pass
		
		for j in range(0,N_terms+1):
			r_plot_zero_values[j] = [self.r_maxes_all[i] for i in range(j*N_exponents,(j+1)*N_exponents) if self.r_Bool_all[i] == 1]
			r_plot_exponent_values[j] = [self.exponent_vals[i] for i in range(0, N_exponents) if self.r_Bool_all[j*N_exponents + i] == 1]
																			   
		xi_plot_zero_values = [self.xi_maxes[i] for i in range(0,N_exponents) if self.xi_Bool[i] == 1]
		xi_plot_exponent_values = [self.exponent_vals[i] for i in range(0,N_exponents) if self.xi_Bool[i] == 1]
		
		# Stop the timer
		end = time.time()
		
		# Print Duration
		print("===== Finished Process =====")
		print("Duration was " + str(round(end-start,3)) + " Seconds for " + str(N_exponents*(N_terms+1)) + " Test Samples with " + str(N_threads) + " threads")
		
		# Plot the results
		# Create figure with right dimensions
		plt.figure(figsize=[6.4,4])
		for j in range(0,N_terms+1):
			plt.plot(r_plot_exponent_values[j], r_plot_zero_values[j], label=r'$r_0$ $j=$' + str(j), linestyle=linestyles[j+1], c='black')
		plt.plot(xi_plot_exponent_values, xi_plot_zero_values, label=r'$\xi_0[r]$', linestyle="-", c='black')
		# Calculate exact results for LE
		A = self.factor
		alpha = lambda n: np.sqrt((n+1)*A**(-2)*p0**((n-1)/(n+1))/(4*np.pi))
		LE_exact_x_vals =  [0,1]
		LE_exact_y_vals = [alpha(0)*np.sqrt(6), alpha(1)*np.pi]
		# Plot them
		plt.scatter(LE_exact_x_vals, LE_exact_y_vals, label=r'$\xi_0[r]$ exact', c='r', marker="P")
		
		plt.legend()
		plt.title(r'$r_0$ where $p(r_0)=0$')
		plt.xlabel(r"Exponent $n=\frac{1}{\gamma-1}$")
		# plt.ylabel(r"$\xi_0")
		plt.yscale('log')
		# Plot a vertical line at xi=5 with the correct height
		plt.vlines(5,0,max(xi_plot_zero_values), colors='k', linestyle="-")
		# Save the plot to a file
		plt.savefig("pictures/TOV-Exponents-LESubs.svg")
		plt.show()
		
	def solveForexponent_vals(self, array, r0, u0, p0, R, rend, dr, N_terms):
		N_exponents = len(self.exponent_vals)
		
		# r_maxes_proc = np.frombuffer(self.r_maxes_all.get_obj())
		# r_maxes = r_maxes_proc.reshape((N_terms, N_exponents))
		# r_Bool_proc = np.frombuffer(self.r_Bool_all.get_obj())
		# r_Bool = r_Bool_proc.reshape((N_terms, N_exponents))
		
		for i in array:
			self.times[i] = time.time()
			results_TOV = [[]]*(N_terms+1)
			results_TOV_small = [[]]*(N_terms+1)
			succ_TOV = [[]]*(N_terms+1)
			r_max = [[]]*(N_terms+1)
			for j in range(0,N_terms+1):
				results_TOV[j], results_TOV_small[j], succ_TOV[j], r_max[j] = Solver.solveTOV(r0, u0, p0, R, rend, dr, terms=j, exponent=Solver.exponent_vals[i], suppressWarning=True)
				self.r_maxes_all[j*N_exponents + i] = r_max[j]
				self.r_Bool_all[j*N_exponents + i] = 1 if r_max[j] < rend else 0
			results_LE, succ_LE, xi_max = Solver.convertSolveLE(r0, u0, p0, R, rend, dr, exponent=Solver.exponent_vals[i], suppressWarning=True, suppressOutput=True, nointerpolate=True)
			self.xi_maxes[i] = xi_max				
			self.xi_Bool[i] = 1 if xi_max < rend else 0
			self.times[i] = time.time()-self.times[i]
		
# Keep track of total time spent
start = time.time()

# Create an instance of the Solver with polytropic EOS
n = 2
gamma = 1+1/n
A = 1

# Initialise Solver with arbitrary values for A and gamma (will not be used)
Solver = Plotter(gamma, A)

# Define initial values
r0 = 0
u0 = 0.0
p0 = 1
R  = 100
rend = R
dr = 0.01

# Define the range of exponents to solve for
n_0 = 0.01
n_max = 5.01
n_step = 0.1

# HINT for choosing N_threads
# CPU: AMD Ryzen 3700X (8 Cores, 16Threads)
# r0=0, u0=0, p0=1, rend=R, dr=0.01, n0 = 0.01, n_max = 4.01, n_step = 0.01
# Threads Time R=100  Time R=500   Time R=1000
# 1       45.88s      231.37s      
# 2       23.57s      123.12s
# 4       12.75s       73.09s
# 8        6.74s       64.60s      455.7s
# 16       5.95s      108.62s
# REMARK: The following variables have the most effect:
# R  > 100
# dr < 0.05
# REMARK2: For good performance/effort choose
# R  = 50
# dr = 0.01
# n_0 = 0.01
# n_max = 5.01
# n_step = 0.02
# N_threads = maximum-1

print("===== Starting Process =====")
Solver.solveMultiprocExponents(n_0, n_max, n_step, r0, u0, p0, R, rend, dr, suppressOutput=False, N_threads=14, N_terms=2)