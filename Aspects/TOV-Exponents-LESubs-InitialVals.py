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
	def solveMultiprocExponents(self, n_0, n_max, n_step, A_initials, p0_initials, r0, u0, p0, R, rend, dr, suppressOutput=False, N_threads = 2, N_terms=1):
		# Define linestyles for plotting later
		linestyles = ['-', '--', '-.', ':']
		
		# Define the array which contains all exponents n = 1/(gamma-1)
		self.exponent_vals = np.arange(n_0,n_max+n_step,n_step)
		self.A_initial_vals = np.array(A_initials)
		self.p0_initial_vals = np.array(p0_initials)
		N_exponents = len(self.exponent_vals)
		N_A_initial = len(self.A_initial_vals)
		N_p0_initial = len(self.p0_initial_vals)
		
		# Define multiprocessing arrays that can be accessed even if processes run in parallel
		self.r_maxes_all = multiprocessing.Array(ctypes.c_float, N_exponents*(N_terms+1)*N_A_initial*N_p0_initial)
		self.r_Bool_all = multiprocessing.Array(ctypes.c_float, N_exponents*(N_terms+1)*N_A_initial*N_p0_initial)
		self.exponents_all = multiprocessing.Array(ctypes.c_float, N_exponents*(N_terms+1)*N_A_initial*N_p0_initial)
		self.A_initials_all = multiprocessing.Array(ctypes.c_float, N_exponents*(N_terms+1)*N_A_initial*N_p0_initial)
		self.p0_initials_all = multiprocessing.Array(ctypes.c_float, N_exponents*(N_terms+1)*N_A_initial*N_p0_initial)
		self.xi_maxes_all = multiprocessing.Array(ctypes.c_float, N_exponents*N_A_initial*N_p0_initial)
		self.xi_Bool_all = multiprocessing.Array(ctypes.c_float, N_exponents*N_A_initial*N_p0_initial)
		
		self.xi_maxes = multiprocessing.Array('f',[0 for i in range(0,N_exponents)])
		self.xi_Bool = multiprocessing.Array('l',[0 for i in range(0,N_exponents)])
		self.times = multiprocessing.Array('f',[0 for i in range(0,N_exponents)])
		
		# Initialise threads; Create list that contains all processes
		processes = []
		
		# Split process in several subthreads
		# and distribute load equally over all threads
		for i in range(0,N_threads):
			# Helper array to store the threads
			
			exponents = []
			j = 0
			# This distribution is chosen such that all the values with high exponent n
			# (for which the solving takes longer) are distributed equally about the different threads
			while N_threads*j+i < N_exponents:
				exponents.append(N_threads*j + i)
				j += 1
			
			# Generate the process
			x = multiprocessing.Process(target=self.solveForexponent_vals, args=(exponents, r0, u0, p0, R, rend, dr, N_terms,))
			# Store x in the processes array
			processes.append(x)
			# Actually start solving the equations in different processes
			x.start()
			# Join the processes to have them end continue the method at the time when all processes have finished
			x.join()
		
		# Create lists to store values for plotting
		# We differentiatie between values for which the solution 
		# was successfull and ones where it was not
		r_plot_zero_values_all = [[]]*(N_terms+1)
		r_plot_A_initial_values_all = [[]]*(N_terms+1)
		r_plot_p0_initial_values_all = [[]]*(N_terms+1)
		r_plot_exponent_values_all = [[]]*(N_terms+1)
		
		xi_plot_zero_values = []
		xi_plot_exponent_values = []
		r_plot_exponent_values = [[]]*(N_terms+1)
		r_plot_zero_values = [[]]*(N_terms+1)
		
		for j in range(0,N_terms+1):
			r_plot_zero_values_all[j] = np.array([[[self.r_maxes_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] for i in range(0,N_exponents)
							  if self.r_Bool_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] == 1]
							 for k in range(0,N_A_initial)] for l in range(0,N_p0_initial)])
			r_plot_exponent_values_all[j] = np.array([[[self.exponents_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] for i in range(0,N_exponents)
							  if self.r_Bool_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] == 1]
							 for k in range(0,N_A_initial)] for l in range(0,N_p0_initial)])
			r_plot_A_initial_values_all[j] = np.array([[[self.A_initials_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] for i in range(0,N_exponents)
							  if self.r_Bool_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] == 1]
							 for k in range(0,N_A_initial)] for l in range(0,N_p0_initial)])
			r_plot_p0_initial_values_all[j] = np.array([[[self.p0_initials_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] for i in range(0,N_exponents)
							  if self.r_Bool_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] == 1]
							 for k in range(0,N_A_initial)] for l in range(0,N_p0_initial)])
			
			r_plot_zero_values[j] = r_plot_zero_values_all[j][0][0]
			r_plot_exponent_values[j] = r_plot_exponent_values_all[j][0][0]
			
# 			r_plot_zero_values[j] = [[self.r_maxes_all[i] for i in range(j*N_exponents,(j+1)*N_exponents) if self.r_Bool_all[i] == 1]]
# 			r_plot_exponent_values[j] = [self.exponent_vals[i] for i in range(0, N_exponents) if self.r_Bool_all[j*N_exponents + i] == 1]
		
		for i, val in enumerate(self.r_maxes_all):
			if suppressOutput == False:
				print("r_end=" + str(round(val,3)) + "   at n=" + str(round(self.exponents_all[i],3)) + "   with A=" + str(round(self.A_initials_all[i],3)) + "   with p0=" + str(round(self.p0_initials_all[i],3)))
		
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

		# Define colors for different graphs (maximum 6)
		colors = ['r', 'g', 'b', 'orangered', 'lime','cornflowerblue']
		# Now set j=0 (corresponds to plain TOV solution)
		# and l=0 corresponds to the first p0 initial value
		j = 0
		l = 0
		# For all different initial values of A we plot 2 graphs: one for TOV and one for LE
		for k in range(0,N_A_initial):
			# Define X and Y values for TOV and LE respectively.
			X_vals = np.array([self.exponents_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] for i in range(0,N_exponents)])
			Y_vals_TOV = np.array([self.r_maxes_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] for i in range(0,N_exponents)])
			Y_vals_LE = np.array([self.xi_maxes_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] for i in range(0,N_exponents)])
			# Also create a mask to hide values which exceeded the integration limit
			Y_mask_TOV = [self.r_Bool_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i]==1 for i in range(0,N_exponents)]
			Y_mask_LE = [self.xi_Bool_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i]==1 for i in range(0,N_exponents)]
			# Plot the graphs
			plt.plot(X_vals[Y_mask_TOV], Y_vals_TOV[Y_mask_TOV], label=r"$r_0$ TOV for $A=$" + str(self.A_initial_vals[k]), color=colors[k], linestyle='-.')
			plt.plot(X_vals[Y_mask_LE], Y_vals_LE[Y_mask_LE], label=r"$r_0$ LE  for $A=$" + str(self.A_initial_vals[k]), color=colors[k], linestyle='-')

		# Formatting of the newly generated plot		
		plt.legend()
		plt.title(r'$r_0$ where $p(r_0)=0$')
		plt.xlabel(r"Exponent $n=\frac{1}{\gamma-1}$")
		plt.ylabel(r"$r_0$")
		plt.yscale('log')
		# Plot a vertical line at xi=5 with the correct height
# 		plt.vlines(5,0,max(xi_plot_zero_values), colors='k', linestyle="-")
		# Save the plot to a file
		plt.savefig("pictures/TOV-Exponents-LESubs-InitialVals.svg")
		plt.show()
		
	def solveForexponent_vals(self, exponents_index, r0, u0, p0, R, rend, dr, N_terms):
		N_exponents = len(self.exponent_vals)
		N_A_initial = len(self.A_initial_vals)
		N_p0_initial = len(self.p0_initial_vals)
		
		# r_maxes_proc = np.frombuffer(self.r_maxes_all.get_obj())
		# r_maxes = r_maxes_proc.reshape((N_terms, N_exponents))
		# r_Bool_proc = np.frombuffer(self.r_Bool_all.get_obj())
		# r_Bool = r_Bool_proc.reshape((N_terms, N_exponents))
		
		results_TOV = [[]]*(N_terms+1)
		results_TOV_small = [[]]*(N_terms+1)
		succ_TOV = [[]]*(N_terms+1)
		r_max = [[]]*(N_terms+1)
		
		for i in exponents_index:
			self.times[i] = time.time()
			for k, A_init in enumerate(self.A_initial_vals):
				for l, p_init in enumerate(self.p0_initial_vals):
					self.factor = A_init
					for j in range(0,N_terms+1):
						results_TOV[j], results_TOV_small[j], succ_TOV[j], r_max[j] = Solver.solveTOV(r0, u0, p0, R, rend, dr, terms=j, exponent=Solver.exponent_vals[i], suppressWarning=True)
						self.r_maxes_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = r_max[j]
						self.r_Bool_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = 1 if r_max[j] < rend else 0
						self.exponents_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = self.exponent_vals[i]
						self.A_initials_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = A_init
						self.p0_initials_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = p_init
					results_LE, succ_LE, xi_max = Solver.convertSolveLE(r0, u0, p0, R, rend, dr, exponent=Solver.exponent_vals[i], suppressWarning=True, suppressOutput=True, nointerpolate=True)
					self.xi_maxes_all[k*N_exponents*N_p0_initial + l*N_exponents + i] = xi_max
					self.xi_Bool_all[k*N_exponents*N_p0_initial + l*N_exponents + i] = 1 if xi_max < rend else 0
			self.times[i] = time.time()-self.times[i]
		
# Keep track of total time spent
start = time.time()

# Create an instance of the Solver with polytropic EOS
n = 2
gamma = 1+1/n
A = 3

# Initialise Solver with arbitrary values for A and gamma (will not be used)
Solver = Plotter(gamma, A)

# Define initial values
r0 = 0
u0 = 0.0
p0 = 1
R  = 1000
rend = R
dr = 0.01

# Define the range of exponents to solve for
n_0 = 0.01
n_max = 5.01
n_step = 0.1

# Define the range for initial values
A_initials = [0.01,0.1,1,10]

p0_initials = [1]

print("===== Starting Process =====")
Solver.solveMultiprocExponents(n_0, n_max, n_step, A_initials, p0_initials, r0, u0, p0, R, rend, dr, suppressOutput=True, N_threads=12, N_terms=2)