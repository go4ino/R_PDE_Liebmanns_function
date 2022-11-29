libm_fixed_bound <- function(A, bounds, relax_fac, tol = 0.01, maxit = 500) {
					# func to calc liebman's method based solution
					# for hotplate problems in Chapra Canale TB
					#
					# inputs
					# A: matrix for inner points (i,j), initially should be all zeros
					#   but does support non zero initial conditions just in case
					# bounds: vector of boundary conditions. Format is [bottom, top, left, right] 
					# relax_fac: relaxation factor
					# tol: error tolerance, default is 1%
					# maxit: max iterations incase we diverge
					
					
					# get length and width for looping var
					n <- dim(A)[1] #rows, aka j values
					m <- dim(A)[2] #columns, aka i values
					
					# add 2 rows and 2 columns for the boundary conditions
					# to avoid having to use a ton of if statements each time a new value is calculated
					
					# bind bottom boundary
					A <- rbind(A, t(rep(bounds[1] , m)))
					# bind top boundary
					A <- rbind(t(rep(bounds[2] , m)) , A)
					# bind right boundary
					A <- cbind(A, rep(bounds[4] , n+2)) # +2 to account for new rows
					# bind left boundary
					A <- cbind(rep(bounds[3] , n+2) , A) # +2 to account for new rows
					
					# past value matrix
					A_prev = A
					
					
					# stopping var
					eps_a <- 4 #epsilon / local error
					n_iter <- 0
					
					# main loop to iterate val
					# 2 dimensional loop
					# outer loop is by row
					# inner loop is by column
					# encased in a while statement with the 2 stopping conditions
					
					while ( (n_iter < maxit) & (eps_a > tol)) {
						# calc new matrix
						# since the entirety of the surrounding rows and cols are boundary conditions
						# we only modify the inside
						
						
						for (j in seq(n+1, 2, by = -1)) {# rows cycle
							
							for (i in 1:m+1) {# column cycle
								# eg: like (i,j) = (1,1) , (2,1) , (3,1) etc
								
								A[i,j] = (A[i+1, j] + A[i-1, j] + A[i, j+1] + A[i, j-1]) / 4
								
								# apply relaxing factor
								A[i,j] = relax_fac * A[i,j] + (1 - relax_fac) * A_prev[i,j]
								
							}
							
						}
						
						
						# adjust stopping paramters
						eps_a <- max(abs(A[2:m+1, 2:n+1] - A_prev[2:m+1, 2:n+1]) / A[2:m+1, 2:n+1]) 
						# we stop when the max error meets our desired tolerance^
						n_iter <- n_iter + 1
						# update prev entry matrix
						A_prev <- A
						
					}
					# return matrix solution, and stopping parameters
					# TODO: remove the added rows and columns of boundary conditions on A
					return(list(A, eps_a, n_iter))
					
					
				}
