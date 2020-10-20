# TOV-LE

## Structure:

Aspects - Contains routines that depend on Solvers and yield plots

* LE-Exponents
  * :arrow_right:	Plots the value of xi_0 (theta(xi_0)=0) for the LE equation for different exponents n
* TOV-LE-Combi
  * :arrow_right:	Plots TOV and LE results
* TOV-SingleSolve-Logp
  * :arrow_right: 	Plots TOV results for 1/r solving method (currently not working)
* TOV-SingleSolve
  * :arrow_right:	Plots TOV results
* TOV-Terms
  * :arrow_right:	Plots TOV results for different amounts of terms present

Solvers - Contains differnt ways to solving the LE and TOV equations

* Solver
  * :arrow_right:	Standard Solver. Solves TOV and LE as written down
* SolverLogp
  * :arrow_right:	Solve with q=log(p) instead of p

## Usage

Usage for all Aspects is the same. Enter the directory and execute the corresponding python file in a python shell.

Example:

```bash
cd Aspects

python TOV-SingleSolve.py
```
