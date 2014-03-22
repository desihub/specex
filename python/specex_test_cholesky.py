#!/usr/bin/env python

import numpy,scipy,scipy.linalg
from specex_cholesky import *


A=numpy.array([[1.,1.,1.],[1.,2.,1.],[1.,1.,4.]])
B=numpy.array([3.,3.,3.])
print "A=",A
print "B=",B
print 
print

X,Ai=cholesky_solve_and_invert(A,B)
print "A=",A
print "Ai=",Ai
print "A.Ai=",A.dot(Ai)
print "X=",X
print "A.X=",A.dot(X)


