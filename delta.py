import numpy
from numpy import *
import scipy, math
from scipy import integrate
from scipy.stats import poisson, norm, lognorm
import matplotlib.pyplot as plt

def Delta():
	sigma=1
	func1= lambda xx: norm.pdf(xx,0,sigma)
	#func2= lambda xx: norm.pdf(xx,0,0.1)
	func3= lambda z: lognorm.pdf(z,1)
	#print 1 - integrate.quad(func1,-4*sigma,4*sigma)[0]
	print 1 - integrate.quad(func3,0,10)[0]
	
	
	x=arange(-1,10,0.1)
	#y1=map(func1,x)
	#y2=map(func2,x)
	y3=map(func3,x)
	#plt.plot(x,y1)
	#plt.plot(x,y2)
	#plt.plot(x,y3)
	#plt.show()

	jeffreys= lambda z: 1/math.sqrt(z+1.5)
	jeffreys2= lambda z: 1/math.sqrt(z+15)
	xxx=arange(0,2.6,0.1)
	y=map(jeffreys,xxx)
	yy=map(jeffreys2,xxx)
	plot1, =plt.plot(xxx,y)
	plot2, =plt.plot(xxx,yy,'r')
	plt.xlabel('Signal strength')
	plt.ylabel('Jeffreys prior')
	plt.legend([plot1,plot2],['s=1.0, b=1.5','s=1.0, b=15.0'])
	#plt.show()
	
	print jeffreys(0)/jeffreys(1)
Delta()	
