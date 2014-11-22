import ROOT
import numpy
from numpy import *
import scipy
from scipy import integrate
from scipy.stats import poisson, norm
import matplotlib.pyplot as plt

def function(liste):
	ff=0
	for z in liste:
		ff+=z/2
	return ff

def another_function(**kwargs):
	ff=0
	for key, z in kwargs.items():
		ff+=z/2
	return ff

def main():
	n=6
	liste=[q for q in xrange(1,n+1)]
	
	func= lambda x1,x2,x3: function([x1,x2,x3])
	print integrate.nquad(func,[[0,1],[0,1],[0,1]], opts=[{"epsabs": 1.49e-20},{"epsabs": 1.49e-20},{"epsabs": 1.49e-20}])
	
	
	
	
	func2= lambda y1,y2,y3: another_function(**{"eins": y1,"zwei": y2,"drei": y3})
	print integrate.nquad(func2,[[0,1],[0,1],[0,1]], opts=[{"epsabs": 1.49e-20},{"epsabs": 1.49e-20},{"epsabs": 1.49e-20}])
main()
