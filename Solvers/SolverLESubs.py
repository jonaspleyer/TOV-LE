#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
class DiffEqSolverLESubs(DiffEqSolver):
	# Solves the Tov equation (also for different additional terms present in equation)
	def solveTOV(self, r0, u0, p0, R, rend, dr, terms=0, exponent=None, suppressWarning=False, noconvert=False):
		# Set the exponent n = 1/(gamma-1) if chose different to init
		if exponent == None:
			gamma = self.gamma
		else:
			gamma = 1+1/exponent
		
		g = gamma
		n = 1/(g-1)
		A = self.factor
		K = A**(-g)
		
		rho0 = A*p0**(1/g)
		alpha =  np.sqrt((n+1)*rho0**(1/n-1)*K/(4*np.pi))
		xi0 = r0/alpha
		xi_end = rend/alpha
		xi_R = R/alpha
		dxi = dr/alpha
		
		# We have to adjust the EOS to accomodate for the change of variables.
		self.eos_new = lambda T, xi: rho0*T**n
		
		# Define the function that returns the increment
		# The parameter "terms" defines how many additional terms are present when comparing
		# the LE and TOV equation
		if terms==0:
			du = self.RK4(lambda xi, u, T: [4 * np.pi * alpha**3 * self.eos_new(T,xi) * xi ** 2,
										   -1/((n+1)*K*rho0**(1/n))*1/(alpha*xi**2)*(1 + K*rho0**(1/n)*T) * (4*np.pi*xi**3*alpha**3*K*rho0**(1+1/n)*T**(n+1) + u) / (1 - 2 * u /(alpha*xi)) if xi>0 else 0])
		elif terms==1:
			du = self.RK4(lambda xi, u, T: [4 * np.pi * alpha**3 * self.eos_new(T,xi) * xi ** 2,
										   -1/((n+1)*K*rho0**(1/n))*1/(alpha*xi**2)*(1 + K*rho0**(1/n)*T) * (4*np.pi*xi**3*alpha**3*K*rho0**(1+1/n)*T**(n+1) + u) if xi>0 else 0])
		elif terms==2:
			du = self.RK4(lambda xi, u, T: [4 * np.pi * alpha**3 * self.eos_new(T,xi) * xi ** 2,
										   -1/((n+1)*K*rho0**(1/n))*u/(alpha*xi**2)*(1 + K*rho0**(1/n)*T) if xi>0 else 0])
		elif terms==3:
			du = self.RK4(lambda xi, u, T: [4 * np.pi * alpha**3 * self.eos_new(T,xi) * xi ** 2,
										   -1/((n+1)*K*rho0**(1/n))*u/(alpha*xi**2) if xi>0 else 0])
		# Define the initial values
		xi = xi0
		T0 = 1
		T = T0
		u = u0

		# Initialise the result arrays (for total and inside)
		if noconvert==False:
			results = np.array([[xi0*alpha, u0, K*rho0**(1+1/n)*T**(n+1), self.eos_new(T0,xi0)]])
			results_small = np.array([[xi0*alpha, u0, K*rho0**(1+1/n)*T**(n+1), self.eos_new(T0,xi0)]])
		else:
			results = np.array([[xi0, u0, T, self.eos_new(T0,xi0)]])
			results_small = np.array([[xi0, u0, T, self.eos_new(T0,xi0)]])
		# Make sure that errors are raised
		np.seterr(all = "raise")
		 
		# Actually integrate the equations
		while xi <= xi_end:
			if xi <= xi_R:
				try:
					# Define the new values of u,p via the increment calculated by the function specified above
					u, T = [u + du(xi, u, T, dxi)[0], T + du(xi, u, T, dxi)[1]]
					# Increase the radial coordinate
					xi = xi + dxi
					# Append the calculated values to the results array
					if noconvert==False:
						results = np.concatenate((results, np.array([[xi*alpha, u, K*rho0**(1+1/n)*T**(n+1), self.eos_new(T,xi)]])))
						results_small = np.concatenate((results_small, np.array([[xi*alpha, u, K*rho0**(1+1/n)*T**(n+1), self.eos_new(T,xi)]])))
					else:
						results = np.concatenate((results, np.array([[xi, u, T, self.eos_new(T,xi)]])))
						results_small = np.concatenate((results_small, np.array([[xi, u, T, self.eos_new(T,xi)]])))
				except:
					break
				if T <= 0:
					break
			else:
				u, T = [u, 0]
				xi = xi + dxi
				if noconvert==False:
					results = np.concatenate((results, np.array([[xi*alpha, u, K*rho0**(1+1/n)*T**(n+1), self.eos_new(T,xi)]])))
				else:
					results = np.concatenate((results, np.array([[xi, u, T, self.eos_new(T,xi)]])))
		# Return the results and a True value for success of the method
		if xi == 0:
			if noconvert==False:
				return results, results_small, False, xi*alpha
			else:
				return results, results_small, False, xi
		elif xi < xi_R:
			if suppressWarning==False:
				print("Warning: Solving only possible up to r = " + str(results[-1][0]) + " with m = " + str(results[-1][1]) + " and p = " + str(results[-1][2]))
			if noconvert==False:
				return results, results_small, True, xi*alpha
			else:
				return results, results_small, True, xi
		elif xi > xi_R:
			if noconvert==False:
				return results, results_small, True, xi*alpha
			else:
				return results, results_small, True, xi