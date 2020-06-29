# Joint-Chance-Constraints
We consider a n-player game where each player's strategy set contains some stochastic linear constraints. The existence of a **Nash Equilibrium** under certain conditions has been proved earlier *(General sum games with joint chance constraints, 2018)*. We further analyse this problem via a constructed example of a Cournot Competition among electricity firms over a network. We run simulations of an Iterative scheme to find the Nash Equilibrium in this context and describe the results.

### Part 1
We have considered different instances (20) in each we have taken a different random matrix (Multivariate Normal Distribution) (Mean and Covariance matrices are different) to ensure existence of Nash Equilibrium in mostly all cases, for similicity we have ignored the changes in other variables which can be incorporated as well. Data obtained from the simulation is also attached. 

>**cournot_latest.m** - main file. Code is quite general, No variables are hard coded, even chance constraint probability distribution can be modified easily

>**converging.m** - plots the payoffs vs iterations depicting the convergence for all different instances.

>**plot_T_I.m** - plots time taken and total iterations done for all different instances, it also provides some information in output.

### Part 2
For a particular instance obtained from above simulation, we simulate for different initial points. This ensures we obtain solution irrespective of the starting point in almost all cases.

>**cournot_latest_changing_x.m** - main file.

>**converging.m** - plots the payoffs vs iterations depicting the convergence for all different instances.

>**plot_T_I.m** - plots time taken and total iterations done for all different instances, it also provides some information in output.


**opRes.pdf** - General sum games with joint chance constraints, 2018
**jcc.pdf** - Joint Chance Constraints - Solving for Nash Equilibrium, 2020
