#!/bin/spyder

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.misc import derivative
import scipy.integrate as integrate

# Ignore MatplotlibDeprecationWarning
import warnings
import matplotlib.cbook
from Solvers.Solver import DiffEqSolver
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)


# NOTICE: We inherit every method from the DiffEqSolver Class previously defined
# but overwrite the solving method for the TOV equation
class DiffEqSolverLogp(DiffEqSolver):
	# Solves the Tov equation (also for different additional terms present in equation)
	def solveTOV(self, r0, u0, p0, R, rend, dr, terms=0):
		# We have to adjust the EOS to accomodate for the change of variables.
		self.eos_new = lambda q,r: self.factor*np.exp(q*1/self.gamma)
		
		# Define the function that returns the increment
		# The parameter "terms" defines how many additional terms are present when comparing
		# the LE and TOV equation
		if terms==0:
			du = self.RK4(lambda r, u, q: [4 * np.pi * self.eos_new(q,r) * r ** 2,
										   -np.exp(-q)*(np.exp(q) + self.eos_new(q,r)) * (4 * np.pi * r ** 3 * np.exp(q) + u) / (r**2 - 2 * u * r) if r>0 else 0])
		elif terms==1:
			du = self.RK4(lambda r, u, q: [4 * np.pi * self.eos_new(q,r) * r ** 2,
										   - np.exp(-q)/r**2 *(self.eos_new(q,r) + np.exp(q)) * (4 * np.pi * r ** 3 * np.exp(q) + u) if r>0 else 0])
		elif terms==2:
			du = self.RK4(lambda r, u, q: [4 * np.pi * self.eos_new(q,r) * r ** 2,
										   - u*np.exp(-q)/r**2 *(np.exp(q) + self.eos_new(p,r)) if r>0 else 0])
		elif terms==3:
			du = self.RK4(lambda r, u, q: [4 * np.pi * self.eos_new(q,r) * r ** 2,
										   - u*np.exp(-q)*self.eos_new(q,r)/r**2 if r>0 else 0])
		# Define the initial values
		r = r0
		q0 = np.log(p0)
		q = q0
		u = u0

		# Initialise the result arrays (for total and inside)
		results = np.array([[r0, u0, np.exp(q0), self.eos_new(q0,r0)]])
		results_small = np.array([[r0, u0, np.exp(q0), self.eos_new(q0,r0)]])
		
		# Make sure that errors are raised
		np.seterr(all = "raise")
		 
		# Actually integrate the equations
		while r <= rend:
			if r <= R:
				try:
					# Define the new values of u,p via the increment calculated by the function specified above
					u, q = [u + du(r, u, q, dr)[0], q + du(r, u, q, dr)[1]]
					# Increase the radial coordinate
					r = r + dr
					# Append the calculated values to the results array
					results = np.concatenate((results, np.array([[r, u, np.exp(q), self.eos_new(q,r)]])))
					results_small = np.concatenate((results_small, np.array([[r, u, np.exp(q), self.eos_new(q,r)]])))
				except:
					break
				if np.exp(q) <= 0:
					break
			else:
				u, q = [u, 0]
				r = r + dr
				results = np.concatenate((results, np.array([[r, u, np.exp(q), self.eos_new(q,r)]])))
		# Return the results and a True value for success of the method
		if r == 0:
			return results, results_small, False, r
		elif r < R:
			print("Warning: Solving only possible up to r = " + str(results[-1][0]) + " with m = " + str(results[-1][1]) + " and p = " + str(results[-1][2]))
			return results, results_small, True, r
		elif r > R:
			return results, results_small, True, r