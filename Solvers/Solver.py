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
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

class DiffEqSolver:
	def __init__(self, gamma, A):
		self.gamma = gamma
		self.factor = A
		# Define the equation of state rho(p,r) = A*p^g
		self.eos = lambda p,r: A*p**(1/gamma) if p > 0 else 0

	# Implementation of the Runge-Kutta method
	def RK4(self, f):
		# This returns a function that upon provided with r, u, dr gives dr*f(r,u,p) via
		# the Runge-Kutta-4 method.
		return lambda r, u, p, dr: (
			lambda du1, dp1: (
				lambda du2, dp2: (
					lambda du3, dp3: (
						lambda du4, dp4: [(du1 + 2 * du2 + 2 * du3 + du4) / 6, (dp1 + 2 * dp2 + 2 * dp3 + dp4) / 6]
					)(dr * f(r + dr	 , u + du3	   , p + dp3)[0]	   , dr * f(r + dr	 , u + du3	   , p + dp3)[1])
				)(dr * f(r + dr / 2 , u + du2 / 2   , p + dp2 / 2)[0]   , dr * f(r + dr / 2 , u + du2 / 2   , p + dp2 / 2)[1])
			)(dr * f(r + dr / 2 , u + du1 / 2   , p + dp1 / 2)[0]   , dr * f(r + dr / 2 , u + du1 / 2   , p + dp1 / 2)[1])
		)(dr * f(r		  , u			 , p)[0]			 , dr * f(r		  , u			 , p)[1])

	# Solves the Tov equation (also for different additional terms present in equation)
	def solveTOV(self, r0, u0, p0, R, rend, dr, terms=0, exponent=None, suppressWarning=False):
		# Set the exponent n = 1/(gamma-1) if chose different to init
		if not exponent == None:
			self.gamma = 1+1/exponent
		
		# Define the function that returns the increment
		# The parameter "terms" defines how many additional terms are present when comparing
		# the LE and TOV equation
		if terms==0:
			du = self.RK4(lambda r, u, p: [4 * np.pi * self.eos(p,r) * r ** 2,
										   -(p + self.eos(p,r)) * (4 * np.pi * r ** 3 * p + u) / (r * (r - 2 * u)) if r>0 else 0])
		elif terms==1:
			du = self.RK4(lambda r, u, p: [4 * np.pi * self.eos(p,r) * r ** 2,
										   - 1/r**2 *(self.eos(p,r) + p) * (4 * np.pi * r ** 3 * p + u) if r>0 else 0])
		elif terms==2:
			du = self.RK4(lambda r, u, p: [4 * np.pi * self.eos(p,r) * r ** 2,
										   - u/r**2 *(p + self.eos(p,r)) if r>0 else 0])
		elif terms==3:
			du = self.RK4(lambda r, u, p: [4 * np.pi * self.eos(p,r) * r ** 2,
										   - u*self.eos(p,r)/r**2 if r>0 else 0])
		# Define the initial values
		r = r0
		p = p0
		u = u0

		# Initialise the result arrays (for total and inside)
		results = np.array([[r0, u0, p0, self.eos(p0,r0)]])
		results_small = np.array([[r0, u0, p0, self.eos(p0,r0)]])

		# Actually integrate the equations
		while r <= rend:
			if r <= R:
				try:
					# Define the new values of u,p via the increment calculated by the function specified above
					u, p = [u + du(r, u, p, dr)[0], p + du(r, u, p, dr)[1]]
					# Increase the radial coordinate
					r = r + dr
					# Append the calculated values to the results array
					results = np.concatenate((results, np.array([[r, u, p, self.eos(p,r)]])))
					results_small = np.concatenate((results_small, np.array([[r, u, p, self.eos(p,r)]])))
				except:
					break
				if p <= 0:
					break
			else:
				u, p = [u, 0]
				r = r + dr
				results = np.concatenate((results, np.array([[r, u, p, self.eos(p,r)]])))
		# Return the results and a True value for success of the method
		if r == 0:
			return results, results_small, False, r
		elif r < R:
			if suppressWarning==False:
				print("Warning: Solving only possible up to r = " + str(results[-1][0]) + " with m = " + str(results[-1][1]) + " and p = " + str(results[-1][2]))
			return results, results_small, True, r
		elif r > R:
			return results, results_small, True, r

	# Solves the Lane-Emden euation	
	def solveLE(self, xi0, T0, dT0, xi_max, dxi, exponent=None, suppressWarning=False):
		if exponent == None:
			exponent = 1/(self.gamma-1)
		dy = self.RK4(lambda x, y1, y2:[y2, 0 if x <= 0 else -2 / x * y2 - y1**exponent if y1 > 0 else 0])
		
		# Define initial values
		x = xi0
		y1 = T0
		y2 = dT0
		
		# Initialise the result array
		results = np.array([[x,y1,y2]])
		
		# Actually integrate the equations
		while x <= xi_max:
			try:
				y1, y2 = [y1 + dy(x,y1,y2,dxi)[0], y2 + dy(x,y1,y2,dxi)[1]]
				if y1 >= 0:
					x = x + dxi
					results = np.concatenate((results, np.array(([[x,y1,y2]]))))
				else:
					break
			except:
				break
  
		if x == 0:
			return results, False, x
		elif x <= xi_max:
			if suppressWarning == False:
				print("Warning: Solving only possible up to xi = " + str(round(results[-1][0],2)) +"/" + str(round(xi_max,2)) + " with T = " + str(results[-1][1]) + " and dT = " + str(results[-1][2]))
			return results, True, x
		elif x > xi_max:
			return results, True, x
	
	# Calculate the density from the mass (if one does not calculate it via EOS)
	def getDensity(self,results):
		# Calculate the density inside the star
		# Make an exception if interpolation fails due to not enough results
		try:
			u_interpolate = interp1d(results[i][:, 0], results[i][:, 1], kind='cubic')
			# Calculate the density from u(r)
			rho = lambda s: 1/(4*np.pi) * derivative(u_interpolate, s, dx=dr/100) / s ** 2 if s <= R else 0
			# Define a reduced range since the interpolation is not valid at the endpoints
			r_reduced = results[i][:, 0][1:results[i].shape[0] - 2]
			
			return r_reduced, rho, True
		except:
			print("Warning: Interpolation was not successful!")
			return [], False, False
	
	# Converts results from LE equation to compare with TOV results
	def convertSolveLE(self, r0, u0, p0, R, rend, dr, exponent=None, suppressWarning=False, suppressOutput=False, noconvert=False):
		# Gather all important variables
		g = self.gamma
		if exponent == None:
			A = self.factor
		else:
			A = exponent
		rho0 = A*p0**(1/g)
		alpha =  np.sqrt(g/(g-1)*rho0**(g-2)/(4*np.pi)/A**(g))
		xi_end = R/alpha
		dxi = dr/alpha
		
		# Solve the LE equation
		results_LE, succ, xi_max = self.solveLE(0,1,0,xi_end, dxi, exponent=exponent, suppressWarning=suppressWarning)
		if noconvert==True:
			return results_LE, succ, xi_max
		r_max = xi_max*alpha
		
		# If this has no success, we do not need to further evaluate
		if succ == False:
			return [], False
		if suppressOutput == False:
			print('alpha     = ' + str(round(alpha,2)))
			print('xi_end[r] = ' + str(round(xi_max*alpha,2)))
		
		# Transform the LE result into the TOV Parameters to compare them
		# this computes the pressure from the LaneEmden function Theta
		conversion_func = lambda x: 1/A**g*rho0**(g)*x**(g/(g-1)) if x > 0 else 0
		p_results_LE_transf = np.array([[res[0]*alpha, conversion_func(res[1]) ] for res in results_LE ])
		rho_results_LE_transf = np.array([[r,self.eos(p,r)] for r,p in p_results_LE_transf])
		
		# Calculate the LE-Mass in dependece of r
		# Use an interpolation for this. If it fails communicate that
		try:
			rho_LE_interpolate = interp1d(rho_results_LE_transf[:,0],rho_results_LE_transf[:,1], kind='cubic')
			m_LE_calculated = np.array([[r,integrate.quad(lambda x: 4*np.pi*x**2*rho_LE_interpolate(x),0,r)[0]] for r in rho_results_LE_transf[:,0]])
			succ = True
		except:
			m_LE_calculated = np.array([[r,0] for r in rho_results_LE_transf[:,0]])
			succ = False
		
		# Return an array with all the results. The array is built [[r,m,p,rho],...,]
		results = np.array([[
			results_LE[i][0]*alpha,
			m_LE_calculated[i][1],
			p_results_LE_transf[i][1],
			rho_results_LE_transf[i][1]
			] for i in range(0,len(results_LE))])
		
		# Return the values of the array
		return results, succ, r_max
		
	# Calculate the total TOV-mass for a given configuration
	def getMass(self, r0, u0, p0, R, rend, dr):
		# A function to quickly calculate the total mass of the star
		# First we need to solve the equations
		results, results_small, succ, r_max = self.solveTOV(r0, u0, p0, R, rend, dr)
		# Check if the solving was successful
		if succ == True :
			# We obtain the mass via the maximum of u(r) since u(r) only increases and is constant for r>=R
			# Also return value True if solving was succesful
			return results[:, 1][-1], results[:,0][-1], True
		else:
			return False, False, False