import itertools
import numpy as np
import scipy
from scipy import integrate
import random
import math
from scipy.stats import poisson, norm, lognorm
from time import gmtime, strftime

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#    Author: Tristan Snowsill
#    License: Creative Commons Attribution-Noncommercial-Share Alike license
#    Package Index Owner: tristan.snowsill
#    DOAP record: mcint-0.1dev5.xml

def integrate(integrand, sampler, measure=1.0, n=1000):
    # Sum elements and elements squared
    total = 0.0
    total_sq = 0.0
    for x in itertools.islice(sampler, n):
        f = integrand(x)
        total += f
        total_sq += (f**2)
    # Return answer
    sample_mean = total/n
    sample_var = (total_sq - ((total/n)**2)/n)/(n-1.0)
    return (measure*sample_mean, measure*math.sqrt(sample_var/n))
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def integrand1(z):
	return lognorm.pdf(z,1)
	
def integrand2(z):
	return norm.pdf(z,0,1)

def sampler1(num):
	while True:
		z=random.uniform(0,10)
                yield z
                
def sampler2(num):
	while True:
		z=random.uniform(-4,4)
                yield z

err1=0
err2=0
for i in xrange(0,100):
	err1+=integrate(integrand1,sampler1(1),measure=10.0, n=1000)[1]/integrate(integrand1,sampler1(1),measure=10.0, n=1000)[0]
	err2+=integrate(integrand2,sampler2(1),measure=40.0, n=1000)[1]/integrate(integrand2,sampler2(1),measure=40.0, n=1000)[0]
print err1/100, err2/100


#print (int(strftime("%H", gmtime())[0:2])+2), int(strftime("%M", gmtime())[0:2])

#func1= lambda *args: function(**Assign(args,**IntVars))
#func2= lambda *args: sum(args)
#print func2(1,1,1,1,1)
#print integrate.nquad(func1,[[0,1],[0,1],[0,1],[0,1]])
