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
import pymongo

class TOVLECalculator(DiffEqSolverLESubs):
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
# 		self.A_initials_all = multiprocessing.Array(ctypes.c_float, N_exponents*(N_terms+1)*N_A_initial*N_p0_initial)
# 		self.p0_initials_all = multiprocessing.Array(ctypes.c_float, N_exponents*(N_terms+1)*N_A_initial*N_p0_initial)
		self.xi_maxes_all = multiprocessing.Array(ctypes.c_float, N_exponents*N_A_initial*N_p0_initial)
		self.xi_Bool_all = multiprocessing.Array(ctypes.c_float, N_exponents*N_A_initial*N_p0_initial)
		
		self.xi_maxes = multiprocessing.Array('f',[0 for i in range(0,N_exponents)])
		self.xi_Bool = multiprocessing.Array('l',[0 for i in range(0,N_exponents)])
		self.times = multiprocessing.Array('f',[0 for i in range(0,N_exponents)])
		
		self.counter = multiprocessing.Value(ctypes.c_int, 0)
		self.update_counter = multiprocessing.Value(ctypes.c_int, 0)
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
		
		print("\n")
		
		for i, val in enumerate(self.r_maxes_all):
			if suppressOutput == False:
				print("r_end=" + str(round(val,3)) + "   at n=" + str(round(self.exponents_all[i],3)) + "   with A=" + str(round(self.A_initials_all[i],3)) + "   with p0=" + str(round(self.p0_initials_all[i],3)))
		
		# Stop the timer
		end = time.time()
		
		# Print Duration
		print("===== Finished Process =====")
		print("Duration was " + str(round(end-start,3)) + " Seconds for " + str(N_exponents*(N_terms+2)*N_A_initial*N_p0_initial) + " Test Samples with " + str(N_threads) + " threads")
		
	def solveForexponent_vals(self, exponents_index, r0, u0, p0, R, rend, dr, N_terms):
		# Create a client to talk to mongodb
		client = pymongo.MongoClient()
		db = client["TOV-LE"]
		exponents_collection = db["Exponents-polytropic-EOS"]
		
		# Define useful parameters
		N_exponents = len(self.exponent_vals)
		N_A_initial = len(self.A_initial_vals)
		N_p0_initial = len(self.p0_initial_vals)
		
		# Define the array that holds the temporary results
		results_TOV = [[]]*(N_terms+1)
		results_TOV_small = [[]]*(N_terms+1)
		succ_TOV = [[]]*(N_terms+1)
		r_max = [[]]*(N_terms+1)
		
		# Iterate over all exponent indices given to the thread
		for i in exponents_index:
			# Start the timer
			A = time.time()
			# Iterate over initial Values for A
			for k, A_init in enumerate(self.A_initial_vals):
				# Iterate over initial Values for p0
				for l, p_init in enumerate(self.p0_initial_vals):
					# Set the value of A_init in the Solver class
					self.factor = A_init
					# Iterate over the given terms
					for j in range(0,N_terms+1):
						# Check if in this particular calculation previously an integration was running over the limit of r_end
						# If this is the case we no longer need to solve since we know that the zero value will only increase.
						if self.checkCalculateResult(exponents_collection, Solver.exponent_vals[i], j, A_init, p_init, 'TOV', rend)==True:
							# Checks if the solving routine before encountered any results that could not be solved.
							if i>=2 and 2 == sum([0==self.r_Bool_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + h] for h in range(0,i-1)]) and 3 == sum([0==int(self.r_maxes_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + h]>0) for h in range(0,i-1)]):
								self.r_maxes_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = rend
								self.r_Bool_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = 0
								# Write here that solving was not successfull with the current parameters
								self.updateDBResults(exponents_collection, Solver.exponent_vals[i], j, A_init, p_init, rend, rend, dr, equ_type='TOV', succ=False)
							# If the previous case was not fulfilled we need to calculate the result.
							else:
								results_TOV[j], results_TOV_small[j], succ_TOV[j], r_max[j] = Solver.solveTOV(r0, u0, p_init, R, rend, dr, terms=j, exponent=Solver.exponent_vals[i], suppressWarning=True)
								self.r_maxes_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = r_max[j]
								self.r_Bool_all[j*N_exponents*N_A_initial*N_p0_initial + k*N_exponents*N_p0_initial + l*N_exponents + i] = 1 if r_max[j] < rend else 0
							
								# Update the results already present
								self.updateDBResults(exponents_collection, Solver.exponent_vals[i], j, A_init, p_init, r_max[j], rend, dr, equ_type='TOV', succ=bool(r_max[j]<rend))
								
					if self.checkCalculateResult(exponents_collection, Solver.exponent_vals[i], 0, A_init, p_init, 'LE', rend) == True:
						# Get the LE results. We do not need to iterate over different terms.
						results_LE, succ_LE, xi_max = Solver.convertSolveLE(r0, u0, p_init, R, rend, dr, exponent=Solver.exponent_vals[i], suppressWarning=True, suppressOutput=True, nointerpolate=True)
						self.xi_maxes_all[k*N_exponents*N_p0_initial + l*N_exponents + i] = xi_max
						self.xi_Bool_all[k*N_exponents*N_p0_initial + l*N_exponents + i] = 1 if xi_max < rend else 0
						self.updateDBResults(exponents_collection, Solver.exponent_vals[i], 0, A_init, p_init, xi_max, rend, dr, equ_type='LE', succ=xi_max<rend)
						
					# Set the counter value for display plus one and append the time it took to calculate.
					self.counter.value = self.counter.value + 1
					self.times[i] = time.time()-A
					print("\rFinished " + str(self.counter.value) + "/" + str(N_exponents*N_A_initial*N_p0_initial) + " and updated " + str(round(self.update_counter.value/(N_exponents*(N_terms+2)*N_A_initial*N_p0_initial)*100,1)) + "% exponents. Last solved exponent was n=" + str(round(self.exponent_vals[i],3)) + " and time t=" + str(round(self.times[i],3)), end='')
			 
	def updateDBResults(self, collection, n, term, A_init, p_init, r_max, r_end, dr, equ_type, succ):
		# Create the details dict with additional information about the solving process
		details_dict = {'r_end':float(r_end), 'dr':float(dr)}
		# Create the succ variable that defines if the solving was successfull in determining the zero value of the function
		succ = bool(r_max<r_end)
		# Create the save dict. The structure with all necessary information that will be saved
		save_dict = {"n":float(n),"terms":int(term),"A_init":float(A_init),'p_init':float(p_init),'r_max':float(r_max), 'equ_type':equ_type, 'succ':succ,'details':details_dict}
		# Search the database for already existing entries
		presents = list(collection.find({"n":float(n),"terms":int(term),"A_init":float(A_init),'p_init':float(p_init),'r_max':float(r_max), 'equ_type':equ_type}))
		# If there are no entries, just save the newly calculated result
		if len(presents)<1:
			collection.insert_one(save_dict)
			self.update_counter.value = self.update_counter.value + 1
		# If there is one or more entries present with the same saving values, we will replace them. 
		# It was previously decided if we need to recalculate these entries or not.
		else:
			for pres in presents:
				collection.delete_one(pres)
			collection.insert_one(save_dict)
			self.update_counter.value = self.update_counter.value + 1
	
	def checkCalculateResult(self, collection, n, term, A_init, p_init, equ_type, rend):
		# Check if a result should be newly calculated
		# Returns TRUE if result should be calculated. Return False otherwise
		results = list(collection.find({"n":float(n),"terms":int(term),"A_init":float(A_init),'p_init':float(p_init), 'equ_type':equ_type}))
		if len(results)==0:
			return True
		# There should not be more than one result. If there are more than one, delete them all.
		elif len(results) >1:
			for res in results:
				collection.delete_one(res)
				# Also return False statement to calculate new result.
				return True
		# In this case we can assume len(results)=1
		else:
			# Check if current solving method would go further in solving
			if results[0]['succ'] == False and results[0]['r_max']<rend:
				return True
			# See if the current solving method would be more accurate (dr is smaller)
			elif dr < results[0]['details']['dr'] and results[0]['r_max'] < rend:
				return True
			else:
				# Return False if the method previously used solved over a bigger 
				# range and was more precise
				return False
		
# Keep track of total time spent
start = time.time()

# Create an instance of the Solver with polytropic EOS
n = 2
gamma = 1+1/n
A = 3

# Initialise Solver with arbitrary values for A and gamma (will not be used)
Solver = TOVLECalculator(gamma, A)

# Define initial values
r0 = 0
u0 = 0.0
p0 = 1
R  = 100
rend = R
dr = 0.001

# Define the range of exponents to solve for
n_0 = 0.01
n_max = 5.00
n_step = 0.01

# Define the range for initial values
A_initials = [0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2]

p0_initials = [0.05,1.5,500]

print(time.strftime('%X %x %Z'))
print("===== Starting Process =====")
Solver.solveMultiprocExponents(n_0, n_max, n_step, A_initials, p0_initials, r0, u0, p0, R, rend, dr, suppressOutput=True, N_threads=6, N_terms=0)