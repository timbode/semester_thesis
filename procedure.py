import ROOT
import numpy
from numpy import *
import scipy
from scipy import integrate
from scipy.stats import poisson, norm, lognorm
import matplotlib.pyplot as plt

#ZH 15 0.445458561182 0.0197227299213
#Zj2b 15 0.742596149445 0.212676465511

def ExtendedPoisson(N,z):
	if z > 0:
		return poisson(z).pmf(N)
	else:
		return 0

# ------------------------------------------------------------------------------

N=1
likelihood=lambda z1,z2: ExtendedPoisson(N,z1+z2)
gauss=lambda z: norm.pdf(z,0.742596149445,0.212676465511)
bayes_bg=lambda z1: ExtendedPoisson(N,z1)*norm.pdf(z1,0.742596149445,0.212676465511)
bayes_signal= lambda z1,z2: ExtendedPoisson(N,z1+z2)*norm.pdf(z1,0.742596149445,0.212676465511)*norm.pdf(z2,0.445458561182,0.0197227299213)
implicit=lambda z: norm.pdf(z,1,0.1)

'''
x=arange(0,10,0.1)
y1=map(likelihood,x)
y2=map(gauss,x)
y3=map(bayes_bg,x)
y4=map(implicit,x)
y5=map(bayes_signal,x)
plt.plot(x,y1)
#plt.plot(x,y2)
plt.plot(x,y3)
#plt.plot(x,y4)
plt.plot(x,y5)
plt.show()
'''

print integrate.quad(bayes_bg,-10,10)
print integrate.nquad(bayes_signal,[[-10,10],[-10,10]])
