#########################################################################
#Author: Cricket Bergner
#Date: 02/10/2026
#########################################################################

# Adaptive Integration

#########################################################################

# import libraries
import numpy as np
import math as m
from prettytable import PrettyTable as pt

# initialize variables
counter = itl = e = 0
N = 1
ans = 1.005702542825726 #from https://www.integral-calculator.com
t = pt(["Number of Slices", "Integral Estimate", "Estimate of Error"])

# trapezoid rule
def trap(f, a, b, N):
  dx = (b-a)/N
  s = 0 # sum  
  for i in range(N):
    xi = a + (i*dx)
    xip = a + ((i+1)*dx)
    s += ((f(xi) + f(xip))* (dx / 2))
  return s

# while loop to iterate until conditions of problem are met
while counter < 13:
  N *= 2
  counter += 1
  itl = trap(lambda x: np.sin(np.sqrt(100*x))**2, 0, 2, N) # calculate integral
  e =  np.abs(itl - ans) # calculate error
  t.add_row([N, round(itl, 7), round(e, 7)]) # add to table
  
print(t) # table
print("The correct output to the integral is ", ans, ".")
print("It took 8192 intervals for the trapezoid rule to approximate this with an accuracy of 10^-6.")

#########################################################################

# Gaussian Quadrature

#########################################################################





#########################################################################

# Subplots

#########################################################################

# import libraries
import scipy as sp
import matplotlib as plt

# initialize variables
roots, weights = sp.special.roots_legendre(N)
t2 = pt(["", "p1", "p2", "p3", "p4"]) # initializing table with headers

for i in range(4): # rows
    for j in range(4): # columns
       


       plt.subplot(4, 4, j+1)

    plt.subplot(4, 4, i+1)



#
#
#

