######################
##
## Numerical integration with trapezoid rule
##
## Int( f,a,b ) = h/2 ( f(x0) + f(xn) ) + h sum( f(xi) ); xi = a + ih; h = ( b - a )/n
##
## Given p processes, each process can work on n/p segments (assuming n/p an integer)
##
## process 0 -> [a, a+nh/p]
## process 1 -> [a+nh/p, a+2nh/p]
##
## process p-1 -> [ a + (p-1)nh/p, b]

"""
test_integration.py -- Parallel trapezoidal rule

input: none
output: Estimate of the integral from a to b of f(x) using the trapezoidal rule and n trapezoids

algorith:
	1. Each process calculates "its" interval of integration
	2. Each process estimates the integral of f(x) over its interval using the trapezoidal rule
	3a. Each process !=0 sends its integral to process 0
	3b. Process 0 sums the calculations received from the individual processes and prints the results

The number of processes (p) should evenly divide the number of trapezoids (n = 1024)

parameters:
	int my_rank	 		- My process rank
	int p 			 		- The number of processes
	float a = 0.0 	- Left endpoint
	float b = 1.0 	- Right endpoint
	int n = 1024		- Number of trapezoids
	float h 				- Trapezoid base length
	float local_a		- left endpoint my process
	float local_b 	- right endpoint my process
	int local_n			- number of trapezoids for my calculation
	float integral	- integral over my interval
	float total=-1 	- Total integral
	int source			- process sending integral
	int dest = 0 		- all messages go o 0
"""

import numpy
from test_func import f, Trap, Get_data
from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
p = comm.Get_size()

if my_rank == 0:
	print('Start computation')
comm.Barrier()

a, b, n = Get_data( my_rank, p, comm )
#a = 0.0; b = 1.0; n = 1024
dest = 0; total = -1.0

h = ( b - a )/n; local_n = n/p 
local_a = a + my_rank * local_n*h 
local_b = local_a + local_n*h 
integral = Trap( local_a, local_b, local_n, h )
print("Local integral computed by process: ",my_rank)

total = comm.reduce( integral, root = 0)
# if my_rank == 0:
# 	total = integral
# 	for source in range(1,p):
# 		integral = comm.recv( source = source )
# 		print("PE ", my_rank," <- ", source, ",", integral)
# 		total = total + integral
# else:
# 	print("PE ", my_rank," -> ", dest, ",", integral)
# 	comm.send( integral, dest = 0 )

print("PE ", my_rank," -> ", dest, ",", total)

if my_rank == 0:
	print("With n= ", n,", trapezoids, \n")
	print("integral from ", a, " to ",b, " = ", total, "\n")

comm.Barrier()
MPI.Finalize

print( my_rank, '--> hi')
