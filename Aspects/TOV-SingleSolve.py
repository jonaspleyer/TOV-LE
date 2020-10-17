# Only use this to import the package DiffEqSolver 
# from ../Solvers/ correctly
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from Solvers.Solver import DiffEqSolver
import matplotlib.pyplot as plt

class Plotter(DiffEqSolver):
	def solveAndPlotResults(self, r0, u0, p0, R, rend, dr):
		# First we need to solve the equations
		results, results_small, succ, r_max = self.solveTOV(r0, u0, p0, R, rend, dr)
		# Check if the solving was successful
		if succ == True:
			# Plot pressure
			plt.plot(results[:, 0], results[:, 2], label=r'Pressure $p(r)$', linestyle='-', c='black')
			# Plot u(r)
			plt.plot(results[:, 0], results[:, 1], label=r'$m(r)$', linestyle='--', c='black')
			# Plot density
			plt.plot(results[:,0], results[:,3], label=r'Density $\rho (r)$', linestyle='-.', c='black')
			plt.legend()
			plt.savefig('pictures/TOV-SingleSolve.svg')
			plt.show()
		else:
			print("Solving was not possible.")

# Define initial values
r0 = 0
u0 = 0
p0 = 0.4
R  = 100
rend = R
dr = 0.01

# Create an instance of the Solver with polytropic EOS
n = 2
gamma = 1+1/n
A = 1

# Initialise Solver
Solver = Plotter(gamma, A)

Solver.solveAndPlotResults(r0, u0, p0, R, rend, dr)