# Intro 

R function I made to solve a PDE with fixed boundary conditions via Liebmann's method

$T_{i,j} = \frac{T_{i+1,j} + T_{i-1,j} + T_{i,j+1} + T_{i,j-1} }{4}$

and overrelaxation applied via formula:

$T_{i,j}^{new}=\lambda T_{i,j}^{new} + (1-\lambda) T_{i,j}^{old}$

# Usage

`libm_fixed_bound <- function(A, bounds, relax_fac, tol = 0.01, maxit = 500)`

Where `A` is a numeric matrix, `bounds` is a vector of the boundary conditions of the order `[bottom, top, left, right]`, 
`relax_fac` is a relaxation factor between 1 and 2 to accelerate convergence, `tol` is the accepted tolerance for error with default value being 1% error, 
and `maxit` is the maximum number of iterations to perform before stopping incase of divergence.

## Examples


# Liscence 

Use however you want no need for credit lmao
