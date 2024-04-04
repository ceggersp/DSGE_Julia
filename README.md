# DSGE_Julia
 `Julia` functions to solve log-linearized DSGE models

## `BK_solve.jl`

`BK_solve.jl` implements the Blanchard & Kahn (1980) method to solve dynamic stochastic linear systems:

$$
E_t X_{t+1}=AE_t X_{t+1}+BE_t X_{t} 
$$

After doing log-linearization, `BK_solve.jl` needs $A$, and $B$ as an input. Also, specifying the indices for endogenous state variables, exogenous state variables, dynamic control variables, and static variables.

After the inputs are given, a matrix $W$ is built such that the system is written only in dynamic terms:

$$
E_t X^D_{t+1}=W X^D_{t} 
$$

If the number of stable eigen-values of $W$ is the same as the number of state variables, then the system is solved and a unique path is found.
