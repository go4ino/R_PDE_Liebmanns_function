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

Using an example from Numerical Methods for Engineers 7th Edition (Steven Chapra, Raymond Canale):

![](https://github.com/go4ino/R_PDE_Liebmanns_function/blob/main/29.4%20copy.png)

29.2 Use Liebmann’s method to solve for the temperature of the square heated plate in Fig. 29.4, but with the upper boundary condition increased to 120C and the left boundary decreased to 60C. Use a relaxation factor of 1.2 and iterate to epsilon = 1%

```R
matty <- matrix(c(0,0,0,
              0,0,0,
              0,0,0), byrow = TRUE, ncol = 3 , nrow = 3)

bcon <- c(0, 120, 60, 50)


libm_fixed_bound(A = matty, bcon, relax_fac =  1.2)
```

Output:

```R
> libm_fixed_bound(A = matty, bcon, relax_fac =  1.2)
[[1]]
     [,1]      [,2]      [,3]      [,4] [,5]
[1,]   60 120.00000 120.00000 120.00000   50
[2,]   60  80.74282  83.87091  77.24861   50
[3,]   60  59.03858  57.52364  54.72136   50
[4,]   60  37.86663  32.42950  34.29923   50
[5,]   60   0.00000   0.00000   0.00000   50

[[2]]
[1] 0.003730279

[[3]]
[1] 5
```
Hella now you can avoid paying for this answer on chegg/quizlet

# Liscence 

Use however you want no need for credit lmao
