#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 18:07:16 2020

@author: jonas
"""
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from Solvers.Solver import DiffEqSolver
import Standards.PlottingStandards as standards

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

class Plotter(DiffEqSolver):
	def calculateThetaWithInitialB(self, b0, steps_max):
		b = [b0]*steps_max
		for k in range(0,steps_max-1):
			try:
				b[k+1] = -sum([b[i]*b[k-i] for i in range(0,k+1)])/((2*k+2)*(2*k+3))
			except:
				print("Stopped calculating")
				break
			
		# Determine the radius of convergence
		conv_radius = [np.float(1/(abs(d)**(1/(2*n+2)))) for n, d in enumerate(b[1:])]
		# Create a function that calculates the value of theta at value xi
		T = lambda xi:sum([c*xi**(2*i) for i, c in enumerate(b)])
		
		return b, conv_radius, T
	
	def compareResults(self, b0, steps_max, xi_max, dxi):
		#  Calculate the coefficients stored in b, the radius of convergence and the function theta(xi)
		b, conv_radius, T = self.calculateThetaWithInitialB(b0, steps_max)
		
		# Calculate the LE results
		results, succ, xi_max = self.solveLE(0, 1, 0, xi_max, dxi, exponent=2)
		
		# Use the same 
		x_vals = results[:,0][results[:,0]<=conv_radius[-1]]
		y_vals = [T(xi) for xi in x_vals]
		
		# Initialise plot with right size
		cm = 1/2.54
		plt.figure(figsize=[16*cm,12*cm])
		
		ax1 = plt.subplot(2,2,1)
		# Set the plot title
		ax1.set_title(r'Series $\theta_{ser}=\sum b_n\xi^{2n}$')
		# Plot the LE results from series expansion
		ax1.plot(x_vals, y_vals, label=r'$\theta_{ser}$', c='k', linestyle=standards.linestyles[0])
		# Plot the solution of the differential equation as well
		
		# Set the ticks to show a tick where the last radius value was calculated
		ticks = np.arange(0,conv_radius[-1],1)
		tickslabels = [str(round(ticks[i],0)) for i in range(len(ticks))]
		ticks = np.append(ticks,conv_radius[-1])
		tickslabels.append(str(round(conv_radius[-1],3)))
		ax1.set_xticks(ticks)
		ax1.set_xticklabels(tickslabels)
		# Also show a legend
		ax1.legend()
		
		# Plot the radius of convergence
		ax2 = plt.subplot(2,2,2)
		ax2.set_title("Radius of convergence")
		ax2.plot(range(2,2*len(b),2), conv_radius, c='k', label='$R_n$')
		ax2.legend()
		
		ax3 = plt.subplot(2,2,3)
		mask = x_vals<=conv_radius[-1]*0.93
		diff = np.array([results[:,1][i]-y_vals[i] for i in range(len(y_vals))])[mask]
		ax3.plot(x_vals[mask], diff, label=r'$\Delta$', c='k', linestyle=standards.linestyles[1])
		ax3.legend()
		ax3.set_title(r'Difference $\Delta=\theta_{calc}-\theta_{ser}$')
		ax3.set_xticks(ticks)
		ax3.set_xticklabels(tickslabels)
		
		ax4 = plt.subplot(2,2,4)
		diff_rel = np.array([(results[:,1][i]-y_vals[i])/results[:,1][i] for i in range(len(y_vals))])[mask]
		ax4.plot(x_vals[mask], diff_rel, label=r'$\Delta/\theta_{calc}$', c='k', linestyle=standards.linestyles[0])
		ax4.set_title(r'Relative Difference $\Delta/\theta_{calc}$')
		ax4.legend()
		ax4.set_xticks(ticks)
		ax4.set_xticklabels(tickslabels)
		
		# Show the final plot
		plt.tight_layout()
		plt.subplots_adjust(hspace=0.3)
		plt.savefig('pictures/LE-ExactN2.svg')
		plt.show()
		matplotlib.use("pgf")
		matplotlib.rcParams.update({
		    'font.family': 'serif',
		    'text.usetex': True,
		    'pgf.rcfonts': False,
		})
		plt.savefig("pictures/LE-ExactN2.pgf", dpi=1000, bbox_inches='tight')

# Initialise Plotter class
Solver = Plotter()

# Describes the number of steps to take to calculate the exponents
steps_max = 100
# This is the initial value of Theta_0 and also the first coefficient of the series expansion of Theta.
b0 = 1
xi_max = 10
dxi = 0.03

# Main method
Solver.compareResults(b0, steps_max, xi_max, dxi)