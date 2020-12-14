#!/bin/spyder

import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import matplotlib.pyplot as plt
from scipy.special import kv
import numpy as np
import matplotlib
from scipy.interpolate import interp1d

import Standards.PlottingStandards

ns = [0.01,0.5,2.5]
n = 0.1
gamma = 1+1/n
p0 = 1
factors = [1,50,2500]
linestyles = ['-', '--', '-.', ':']
curve_colours = ["black","teal","darkorange","lime"]

plt.figure(figsize=[8,6])

for i, factor in enumerate(factors):
	# Define pressure of alpha function
	# First we need to define the function that we want to invert
	f_alpha = lambda a: factor/(kv(2,a)*a**2)*np.exp(-a*(kv(1,a)+kv(3,a))/(2*kv(2,a)))
	a0 = 0.000000001
	
	# Define the x- and y-vals for interpolation
	y_vals = np.linspace(a0,100)
	x_vals = f_alpha(y_vals)
	# interpolate the inverse by switching y and x 
	alpha = interp1d(y_vals, x_vals)
	# Define the EOS
	eos_new = lambda p:p*(1+alpha(p)*(kv(1,alpha(p))+kv(3,alpha(p)))/(2*kv(2,alpha(p))))
	
	# Choose the multiplication factor such that polytropic and rel eos results
	# can be compared later on
	A = eos_new(p0)/p0**(1/gamma)
	eos = lambda p: A*p**(1/gamma)
	
	# Define plot Range
	factor_range = 1
	# Set plotting values
	x_vals = np.arange(0.000001, factor_range*p0,factor_range*p0/10000)/p0
	y_rel = eos_new(x_vals)/eos_new(p0)
	
	plt.plot(x_vals, y_rel, label=r'$\rho(p)_{rel}$ for $A=$'+str(factor), c=curve_colours[i], linestyle=linestyles[0])

for j, n in enumerate(ns):
	gamma = 1+1/n
	eos = lambda p: A*p**(1/gamma)
	x_vals = np.linspace(0.01, factor_range*p0)/p0
	y_cla = eos(x_vals)/eos_new(p0)
	plt.plot(x_vals, y_cla, label=r'$\rho(p)_{cla}$ for $n=$'+str(n), c=curve_colours[j], linestyle=linestyles[1])
	
# for j, n in enumerate(ns):
	
	
# Set ticks for better display (first x-axis)
x_N_ticks = 4 # >= 2
x_tickstep = (max(x_vals)-min(x_vals))/(x_N_ticks-1)
x_ticks = np.arange(min(x_vals), max(x_vals) + x_tickstep/2, x_tickstep)
x_labels = [str(round(tick,1))+"$p_0$" for tick in x_ticks]

# Do the same for y axis
y_N_ticks = 6 # >= 2
y_up = max(max(y_cla),max(y_rel))
y_down = min(min(y_cla),min(y_rel))
y_tickstep = (y_up-y_down)/(y_N_ticks-1)
y_ticks = np.arange(y_down, y_up + y_tickstep/2, y_tickstep)
y_labels = [str(round(tick,1))+r"$\rho_{0,rel}$" for tick in y_ticks]
	
plt.xticks(x_ticks, x_labels)
plt.yticks(y_ticks, y_labels)
plt.legend()
plt.savefig("RelEOS.svg", dpi=1000, bbox_inches='tight')
plt.show()