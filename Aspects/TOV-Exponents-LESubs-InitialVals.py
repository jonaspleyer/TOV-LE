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
		
		self.counter = multiprocessing.Value(ctypes.c_int, 0)
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
			x = multiprocessing.Process(target=self.solveForexponent_vals, args=(exponents, r0, u0, p0, R, rend, dr, N_terms, ))
			# Store x in the processes array
			processes.append(x)
			# Actually start solving the equations in different processes
			x.start()
		
		# Join the processes to have them end continue the method at the time when all processes have finished
		for x in processes:
			x.join()
		
		sys.stdout.write("\n")
		
		for i, val in enumerate(self.r_maxes_all):
			if suppressOutput == False:
				print("r_end=" + str(round(val,3)) + "   at n=" + str(round(self.exponents_all[i],3)) + "   with A=" + str(round(self.A_initials_all[i],3)) + "   with p0=" + str(round(self.p0_initials_all[i],3)))
		
		# Stop the timer
		end = time.time()
		
		# Print Duration
		print("===== Finished Process =====")
		print("Duration was " + str(round(end-start,3)) + " Seconds for " + str(N_exponents*(N_terms+1)*N_A_initial*N_p0_initial) + " Test Samples with " + str(N_threads) + " threads")
		
		# Plot the results
		# Create figure with right dimensions
		fig = plt.figure(figsize=[7,4])
		ax = plt.subplot(1,1,1)

		# Define colors for different graphs (maximum 6)
		colors = ['k','r', 'g', 'b', 'orangered', 'lime','cornflowerblue']
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
			ax.plot(X_vals[Y_mask_TOV], Y_vals_TOV[Y_mask_TOV], label=r"$r_0$ TOV for $A=$" + str(self.A_initial_vals[k]), color=colors[k], linestyle='-.')
			ax.plot(X_vals[Y_mask_LE], Y_vals_LE[Y_mask_LE], label=r"$r_0$ LE  for $A=$" + str(self.A_initial_vals[k]), color=colors[k], linestyle='-')

		# Formatting of the newly generated plot		
		ax.set_title(r'$r_0$ where $p(r_0)=0$')
		ax.set_xlabel(r"Exponent $n=\frac{1}{\gamma-1}$")
		ax.set_ylabel(r"$r_0$")
		ax.set_yscale('log')
		
		# Resize the whole plot to fit the legend next to it.
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
		# Create the legend right of plot
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		
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
			A = time.time()
			for k, A_init in enumerate(self.A_initial_vals):
				for l, p_init in enumerate(self.p0_initial_vals):
					self.factor = A_init
					for j in range(0,N_terms+1):
						if i>=1 and 0 in [self.r_Bool_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + h] for h in range(0,i-1)] and 0 in [int(self.r_maxes_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + h]>0) for h in range(0,i-1)]:
							self.r_maxes_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = rend
							self.r_Bool_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = 0
						else:
							results_TOV[j], results_TOV_small[j], succ_TOV[j], r_max[j] = Solver.solveTOV(r0, u0, p0, R, rend, dr, terms=j, exponent=Solver.exponent_vals[i], suppressWarning=True)
							self.r_maxes_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = r_max[j]
							self.r_Bool_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = 1 if r_max[j] < rend else 0
						# Do this nevertheless if solving was successfull or not
						self.exponents_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = self.exponent_vals[i]
						self.A_initials_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = A_init
						self.p0_initials_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = p_init
					results_LE, succ_LE, xi_max = Solver.convertSolveLE(r0, u0, p0, R, rend, dr, exponent=Solver.exponent_vals[i], suppressWarning=True, suppressOutput=True, nointerpolate=True)
					self.xi_maxes_all[k*N_exponents*N_p0_initial + l*N_exponents + i] = xi_max
					self.xi_Bool_all[k*N_exponents*N_p0_initial + l*N_exponents + i] = 1 if xi_max < rend else 0
			self.counter.value = self.counter.value + 1
			self.times[i] = time.time()-A
			print("\rFinished " + str(self.counter.value) + "/" + str(N_exponents) + " exponents with n=" + str(round(self.exponent_vals[i],3)) + " and time t=" + str(round(self.times[i],3)), end='')
		
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
R  = 100000
rend = R
dr = 0.5

# Define the range of exponents to solve for
n_0 = 0.01
n_max = 4.91
n_step = 0.1

# Define the range for initial values
A_initials = [0.1]

p0_initials = [1]

print("===== Starting Process =====")
Solver.solveMultiprocExponents(n_0, n_max, n_step, A_initials, p0_initials, r0, u0, p0, R, rend, dr, suppressOutput=True, N_threads=14, N_terms=2)