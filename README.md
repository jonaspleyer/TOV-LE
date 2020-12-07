# TOV-LE

## Structure:

### Aspects
Contains routines that depend on Solvers and yield plots

Aspect		|Solver		|Description				
--------------------|---------------|---------------------------------------
LE-Exponents                |Solver 		|Plots the value of xi_0 (theta(xi_0)=0) for the LE equation for different exponents n
LE-SingleSolve              |Solver         |Plots LE results (for multiple exponents)
TOV-LE-Combi                |Solver         |Plots TOV and LE results
TOV-RelEOS-Comparison       |Solver         |Plots the TOV equation for a polytropic and relativistic EOS
TOV-Exponents-Logp          |SolverLogp     |Same as LE-Exponents but also for TOV results with q=log(p) substitution
TOV-Exponents-LESubs|SolverLESubs   |Same as LE-Exponents but also for TOV results with same substitution as LE equation
TOV-SingleSolve             |Solver         |Plots TOV results
TOV-LE-SingleSolve-LESubs   |SolverLESubs   |Plots TOV and LE Results with theta substitution
TOV-LE-SingleSolve-LESubs-InitialVals   |SolverLESubs   |Plots TOV and LE Results with theta substitution and different initial Values for A and p0
TOV-Exponents-LESubs-InitialVals-Database-CalculateOnly |SolverLESubs   |Calculates TOV and LE Results and stores them in a mongo database. Requires running mongo server!
TOV-Exponents-LESubs-InitialVals-Database-PlotResults   |-              |Only for plotting results calculated by above method.
TOV-LE-SingleSolve-Logp        |SolverLogp     |Plots TOV results for q=log(p) solving method
TOV-Terms                   |Solver         |Plots TOV results for different amounts of terms present

### Solvers
Contains differnt ways to solving the LE and TOV equations

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
