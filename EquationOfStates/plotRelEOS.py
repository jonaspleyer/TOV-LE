#!/bin/spyder

# import os, sys
# currentdir = os.path.dirname(os.path.realpath(__file__))
# parentdir = os.path.dirname(currentdir)
# sys.path.append(parentdir)

import matplotlib.pyplot as plt
from scipy.special import kv
import numpy as np
import matplotlib
from scipy.interpolate import interp1d

import Standards.PlottingStandards as standards

ns = [0.2,1,2]
n = 0.5
gamma = 1+1/n
p0 = 5
#Factor needs to be larger than 17
factors = [20,200]

linestyles_group1 = [(0,(4,1)),(0,(4,4)),(0,(4,8))]
linestyles_group2 = [(0,(0.5,1)),(0,(0.5,3)),(0,(0.5,6))]

curve_colours = ["black","gray","lightgrey"]

cm = 1/2.54
plt.figure(figsize=[16*cm,12*cm])

for i, factor in enumerate(factors):
	# Enlarge the factor to make inversion via interpolation successfull
# 	factor = p0*factor*(np.exp(1)**2*2)*1.1
	
	# Define pressure of alpha function
	# First we need to define the function that we want to invert
	p_alpha = lambda a: p0*factor/(kv(2,a)*a**2)*np.exp(-a*(kv(1,a)+kv(3,a))/(2*kv(2,a))) if a > 0 else factor/(np.exp(1)**2*2)
	
	# Calculate the highest alpha value until which we need to interpolate
	a_max = 1
	while p_alpha(a_max) > 0.01*p0:
		a_max += 1
	a_max += -1
	
	# Define the x- and y-vals for interpolation
	a_vals = np.arange(0,a_max,0.1)
	p_vals = [p_alpha(x) for x in a_vals]
	# interpolate the inverse by switching y and x 
	alpha = interp1d(p_vals, a_vals)
	# Define the EOS
	eos_new = lambda p:p*(1+alpha(p)*(kv(1,alpha(p))+kv(3,alpha(p)))/(2*kv(2,alpha(p))))
	
	# Define plot Range
	factor_range = 1
	# Set plotting values
	x_vals = np.linspace(p_vals[-2], p0)
	y_rel = eos_new(x_vals)/eos_new(p0)
	x_vals = x_vals/p0
	
	plt.plot(x_vals, y_rel, label=r'$\rho(p)_{rel}$ for $A='+str(factor)+"$", c=curve_colours[0], linestyle=linestyles_group1[i])

for j, n in enumerate(ns):
	gamma = 1+1/n
	eos = lambda p: p**(1/gamma)
	x_vals = np.linspace(0.01*p0, factor_range*p0)
	y_cla = eos(x_vals)/p0**(1/gamma)
	x_vals = x_vals/p0
	plt.plot(x_vals, y_cla, label=r'$\rho(p)_{cla}$ for $n='+str(n)+"$", c=curve_colours[0], linestyle=linestyles_group2[j])	
	
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
y_labels = [str(round(tick,1))+r"$\rho_{0}$" for tick in y_ticks]
	
plt.xticks(x_ticks, x_labels)
plt.yticks(y_ticks, y_labels)
plt.legend()
plt.savefig("RelEOS.svg", dpi=1000, bbox_inches='tight')
plt.show()