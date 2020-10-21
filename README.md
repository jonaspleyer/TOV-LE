# TOV-LE

## Structure:

### Aspects - Contains routines that depend on Solvers and yield plots

* :arrow_right: LE-Exponents
  * Plots the value of xi_0 (theta(xi_0)=0) for the LE equation for different exponents n
  * Uses Solver
* :arrow_right: TOV-LE-Combi
  * Plots TOV and LE results
  * Uses Solver
* :arrow_right: TOV-Exponents-Logp
  * Same as LE-Exponents but also for TOV results with q=log(p) substitution
  * Uses SolverLogp
* :arrow_right: TOV-Exponents-LESubs
  * Same as LE-Exponents but also for TOV results with same substitution as LE equation
  * Uses SolverLESubs
* :arrow_right: TOV-SingleSolve
  * Plots TOV results
  * Uses Solver
* :arrow_right: TOV-SingleSolve-Logp
  * Plots TOV results for q=log(p) solving method
  * Uses Solverlogp
* :arrow_right: TOV-Terms
  * Plots TOV results for different amounts of terms present
  * Uses Solver

### Solvers - Contains differnt ways to solving the LE and TOV equations

* :arrow_right: Solver
  * Standard Solver. Solves TOV and LE as written down
* :arrow_right: SolverLogp
  * Solve with q=log(p) instead of p
* :arrow_right: SolverLESubs
  * Solve with same substitution as the LE equation


## Usage

Usage for all Aspects is the same. Enter the directory and execute the corresponding python file in a python shell.

Example:

```bash
cd Aspects

python TOV-SingleSolve.py
```
