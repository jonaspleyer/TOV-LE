# TOV-LE

## Structure:

Aspects - Contains routines that depend on Solvers and yield plots

- LE-Exponents		
	-:arrow_right:	Solves the LE equation for different exponents n and plots the value of xi0 at which theta=0
- TOV-LE-Combi		
	-:arrow_right:	Plots TOV and LE results
- TOV-OverR		
	-:arrow_right: 	Supposed to plot TOV results for 1/r solving method (currently not working)
- TOV-SingleSolve	
	-:arrow_right:	Solve and Plot TOV equation
- TOV-Terms		
	-:arrow_right:	Solve the TOV equation for different amounts of terms present

Solvers - Contains differnt ways to solving the LE and TOV equations

- Solver  		
	-:arrow_right:	Standard Solver. Solves TOV and LE as written down
- SolverOverR  		
  -:arrow_right:	Supposed to solve with (p*r) instead of p. Currently not working

## Usage

Usage for all Aspects is the same. Enter the directory and execute the corresponding python file in a python shell.

Example:

```bash
cd Aspects

python TOV-SingleSolve.py
```
