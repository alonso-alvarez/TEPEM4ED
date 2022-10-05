import numpy as np

def f(x):
	return x*x

def Trap( a, b, n, h ):
	integral = ( f(a) + f(b) )/2.0
	x = a 
	for i in range( 1,int(n) ):
		x = x + h
		integral = integral + f(x)

	return integral*h

def Get_data( my_rank, p, comm ):
	a = None
	b = None
	n = None
	if my_rank == 0:
		print("Rank ",my_rank,": Enter a, b and n")
		a = float( input("Enter a: "))
		b = float( input("Enter b: "))
		n = float( input("Enter n: "))
		print("ready to broadcast\n")

	a = comm.bcast( a, root = 0 )
	b = comm.bcast( b, root = 0 )
	n = comm.bcast( n, root = 0 )

	return a, b, n