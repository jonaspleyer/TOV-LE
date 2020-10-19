import numpy as np
# Only use this to import the package DiffEqSolver 
# from ../Solvers/ correctly
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from Solvers.Solver import DiffEqSolver
import multiprocessing
import time
import matplotlib.pyplot as plt

class Plotter(DiffEqSolver):
	def solveMultiprocExponents(self, n_0, n_max, n_step, suppressOutput=False, N_threads = 13):
		#
		self.exponent_vals = np.arange(n_0,n_max+n_step,n_step)
		
		# max_r_results = multiprocessing.Array('l',[[0,False,0,0] for i in range(0,len(self.exponent_vals))])
		self.xi_maxes = multiprocessing.Array('f',[0 for i in range(0,len(self.exponent_vals))])
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
			x = multiprocessing.Process(target=self.solveForexponent_vals, args=(array, xi0, T0, dT0, xi_end, dxi,))
			# Store it in the x array
			processes.append(x)
			x.start()
		
		# Join the processes to have them end continue the method at the time when all processes have finished
		for process in processes:
			process.join()
		
		# Create lists to store values for plotting
		# We differentiatie between values for which the solution 
		# was successfull and ones where it was not
		plot_zero_values = []
		plot_zero_values_Not = []
		plot_exponent_values = []
		plot_exponent_values_Not = []
		
		for i in range(0,len(self.exponent_vals)):
			if suppressOutput == False:
				print(str(self.xi_maxes[i]) + "  " + str(self.xi_Bool[i]) + "  " + str(self.exponent_vals[i]) + "  " + str(self.times[i]))
			if self.xi_Bool[i] == 1:
				plot_zero_values.append(self.xi_maxes[i])
				plot_exponent_values.append(self.exponent_vals[i])
			else:
				plot_zero_values_Not.append(self.xi_maxes[i])
				plot_exponent_values_Not.append(self.exponent_vals[i])
		
		# Stop the timer
		end = time.time()
		
		# Print Duration
		print("===== Finished Process =====")
		print("Duration was " + str(round(end-start,3)) + " Seconds for " + str(len(self.exponent_vals)) + " Test Samples with " + str(N_threads) + " threads")
		
		# Plot the results
		# Create figure with right dimensions
		plt.figure(figsize=[6.4,4])
		plt.plot(plot_exponent_values, plot_zero_values, label=r'$\xi_0$', linestyle="-", c='black')
		plt.legend()
		plt.title(r'$\xi_0$ where $\theta(\xi_0)=0$')
		plt.xlabel(r"Exponent $n$")
		# plt.ylabel(r"$\xi_0")
		plt.yscale('log')
		# Plot a vertical line at xi=5 with the correct height
		plt.vlines(5,0,max(plot_zero_values), colors='k', linestyle="--")
		# Save the plot to a file
		plt.savefig("pictures/LE-Exponents.svg")
		plt.show()
		
	def solveForexponent_vals(self, array, xi0, T0, dT0, xi_end, dxi):
		for i in array:
			self.times[i] = time.time()
			result, succ, xi_max = Solver.solveLE(xi0, T0, dT0, xi_end, dxi, self.exponent_vals[i], suppressWarning=True)
			self.xi_maxes[i] = xi_max
			self.xi_Bool[i] = 1 if xi_max < xi_end else 0
			self.times[i] = time.time()-self.times[i]
		
# Keep track of total time spent
start = time.time()

# Initialise Solver with arbitrary values for A and gamma (will not be used)
Solver = Plotter(1, 1)

# Define initial values
# These values should not be altered
xi0 = 0
T0 = 1
dT0 = 0

# Define stepsize and endpoint for solving
# These may be changed to fit your convenience
xi_end = 5000
dxi = 0.1

# Define the range of exponents to solve for
n_0 = 0.01
n_max = 5.01
n_step = 0.01

# HINT for choosing N_threads
# CPU: AMD Ryzen 3700X (8 Cores, 16Threads)
# xi0=0, T0=1, dT0=0, xi_end=5000, dxi=0.1, 
# Threads Time (Average over 5 samples)
# 1       7.27s
# 2       3.96s
# 4       2.79s
# 8       2.43s
# 16      2.36
# Every integer value >=1 is possible
print("Start Solving")
Solver.solveMultiprocExponents(n_0, n_max, n_step, suppressOutput=True, N_threads=4)
