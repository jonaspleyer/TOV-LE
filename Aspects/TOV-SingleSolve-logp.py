#!/bin/spyder

# Only use this to import the package DiffEqSolver 
# from ../Solvers/ correctly
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from Solvers.SolverLogp import DiffEqSolverLogp
import matplotlib.pyplot as plt

class Plotter(DiffEqSolverLogp):
	def solveAndPlotResults(self, r0, u0, p0, R, rend, dr):
		# First we need to solve the equations
		results_TOV, results_TOV_small, succ_TOV, r_max_TOV = self.solveTOV(r0, u0, p0, R, rend, dr)
		results_LE, succ_LE, r_max_LE = self.convertSolveLE(r0, u0, p0, R, rend, dr)
		# Check if the solving was succ_TOVessful
		if succ_TOV == True:
			# Create Plot with different subplots
			plt.figure(figsize=[8,6])
			plt.tight_layout()
			plt.subplots_adjust(hspace=0.3)
			
			# Plot pressure
			plt.subplot(2,2,1)
			plt.title("Pressure")
			plt.plot(results_TOV[:, 0], results_TOV[:, 2], label=r'$p_{TOV}(r)$', linestyle='-', c='black')
			plt.plot(results_LE[:, 0], results_LE[:, 2], label=r'$p_{LE}(r)$', linestyle='-.', c='black')
			plt.legend()
			
			# Plot u(r)
			plt.subplot(2,2,2)
			plt.title("Mass")
			plt.plot(results_TOV[:, 0], results_TOV[:, 1], label=r'$m_{TOV}(r)$', linestyle='-', c='black')
			plt.plot(results_LE[:, 0], results_LE[:, 1], label=r'$m_{LE}(r)$', linestyle='-.', c='black')
			plt.legend()
			
			# Plot density
			plt.subplot(2,2,3)
			plt.title("Density")
			plt.plot(results_TOV[:,0], results_TOV[:,3], label=r'$\rho_{TOV} (r)$', linestyle='-', c='black')
			plt.plot(results_LE[:,0], results_LE[:,3], label=r'$\rho_{LE} (r)$', linestyle='-.', c='black')
			plt.legend()
			
			# Create table with all values
			plt.subplot(2,2,4)
			# plt.title(r'Values of the different plots')
			plt.axis('off')
			plt_vals = [
				["TOV", r'$r_{max}$', "%.4g" % r_max_TOV],
				["TOV", r'$p(r_{max})$', "%.4g" % results_TOV_small[:,2][-1]],
				["TOV", r'$\rho(r_ {max})$', "%.4g" % results_TOV_small[:,3][-1]],
				["LE,", r'$r_{max}$', "%.4g" % r_max_LE],
				["LE",  r'$p(r_{max})$', "%.4g" % results_LE[:,2][-1]],
				["LE",  r'$\rho(r_ {max})$', "%.4g" % results_LE[:,3][-1]]
				]
			polyTable = plt.table(cellText=plt_vals,
						 colLabels=['Equation','Parameter', 'Value'],
						 colColours=['#dedede','#dedede','#dedede'],
						 loc='center',
						 cellLoc='center'
						 )
			polyTable.auto_set_font_size(False)
			polyTable.set_fontsize(11)
			polyTable.scale(1.2,1.675)
			
			# Save the figure
			plt.savefig('pictures/TOV-SingleSolve-OverR.svg')
			plt.show()
		else:
			print("Solving was not possible.")

# Define initial values
r0 = 0
u0 = 0.0
p0 = 1
R  = 10
rend = R
dr = 0.01

# Create an instance of the Solver with polytropic EOS
n = 2
gamma = 1+1/n
A = 5

# Initialise Solver
Solver = Plotter(gamma, A)

Solver.solveAndPlotResults(r0, u0, p0, R, rend, dr)