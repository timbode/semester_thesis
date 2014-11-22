from numpy import *
from scipy.stats import norm, poisson
from scipy import integrate
import matplotlib.pyplot as plt
from time import gmtime, strftime

# Gaussian for s times Gaussian for b times Poissonian alla (mu*s+b)

start_time=strftime("%Y-%m-%d %H:%M:%S", gmtime())
'''
N=100.0 # ACHTUNG! muss wegen des poissonians ganzzahlig sein...
mean_s=5.5
sigma_s=0.5
mean_b=100.1
sigma_b=0.5
'''
# dieser Satz von Werten produziert die disparaten Resultate: Daten leicht unter dem background...
#N=6.0 # ACHTUNG! muss wegen des poissonians ganzzahlig sein...
#mean_s=0.5
#sigma_s=0.1
#mean_b=6.1
#sigma_b=0.5

N=6.0 # ACHTUNG! muss wegen des poissonians ganzzahlig sein...
mean_s=0.5
sigma_s=0.1
mean_b=6.1
sigma_b=0.5

def ExtendedPoisson(M,z):
    if z > 0.0:
        return poisson(z).pmf(M)
    else:
        return 0

# Bayesian method: marginalisation

def Posterior():
	func= lambda s,b,mumu: norm.pdf(s,mean_s,sigma_s)*norm.pdf(b,mean_b,sigma_b)*ExtendedPoisson(N,mumu*s+b)
	func2= lambda mumumu: integrate.nquad(func,[[mean_s - 5*sigma_s,mean_s + 5*sigma_s],[mean_b - 5*sigma_b,mean_b + 5*sigma_b]],[mumumu])[0]#/integrate.nquad(func,[[mean_s - 5*sigma_s,mean_s + 5*sigma_s],[mean_b - 5*sigma_b,mean_b + 5*sigma_b],[1e-15,10]])[0]

	x=arange(0,5.0,0.2)
	y=map(func2,x)
	plt.plot(x,y)
	plt.show()

 
# Frequentist method:

# finding the maximum analytically:
def pqFormel(p,q):
	return [-0.5*p+sqrt(0.25*p*p-q), -0.5*p-sqrt(0.25*p*p-q)]

def AnaMax(mumu,M):
	if mumu > 0.0:
		# log der Distribution nach s und b ableiten, Hilfsvariable z=mu*s+b
		p=mumu*mumu*sigma_s*sigma_s+sigma_b*sigma_b-(mumu*mean_s+mean_b)
		q=-M*(mumu*mumu*sigma_s*sigma_s + sigma_b*sigma_b)
		z=max(pqFormel(p,q)) # max?
		b=mean_b - sigma_b*sigma_b + M*sigma_b*sigma_b/z
		s=(z-b)/mumu
		return s,b

	elif mumu == 0.0:
		return mean_s,max(mean_b/2+sqrt(0.25*mean_b*mean_b+sigma_b*sigma_b*(M-1)),mean_b/2+sqrt(0.25*mean_b*mean_b+sigma_b*sigma_b*(M-1)))
	else:
		return 'error'
		

# global max. is given by s=mean_s, b=mean_b, mu=(N-b)/s

def Log(s,b,mumu,M):
	liste=[q for q in xrange(1,int(M+1))]
	return -log(2*pi)-log(sigma_s)-log(sigma_b)-((s-mean_s)*(s-mean_s)/(0.5*sigma_s*sigma_s))-((b-mean_b)*(b-mean_b)/(0.5*sigma_b*sigma_b))+M*log(mumu*s+b)-sum(liste)-(mumu*s+b)

def Log2(b,M):
	liste=[q for q in xrange(1,int(M+1))]
	return -0.5*log(2*pi)-log(sigma_b)-((b-mean_b)*(b-mean_b)/(0.5*sigma_b*sigma_b))+M*log(b)-sum(liste)-b

def Profile(mumu):
	if mumu > 0.0:
		s0,b0 = AnaMax(0.0,N)
		s1,b1 = AnaMax(1.0,N)

		param_s0=random.normal(s0,sigma_s)
		param_b0=random.normal(b0,sigma_b)
		param_s1=random.normal(s1,sigma_s)
		param_b1=random.normal(b1,sigma_b)

		if mumu*param_s0+param_b0 < 0: # ?
			toy0=0.0
		else:
			toy0=random.poisson(mumu*param_s0+param_b0)

		if mumu*param_s1+param_b1 < 0: # ?
			toy1=0.0
		else:
			toy1=random.poisson(mumu*param_s1+param_b1)

		# constraints from the procedure paper:
		if (toy0-mean_b)/mean_s < 0.0:
			muMax0=0.0
		elif (toy0-mean_b)/mean_s > 1.0:
			muMax0=1.0
		else:
			muMax0=(toy0-mean_b)/mean_s

		# constraints from the procedure paper:
		if (toy1-mean_b)/mean_s < 0.0:
			muMax1=0.0
		elif (toy1-mean_b)/mean_s > 1.0:
			muMax1=1.0
		else:
			muMax1=(toy1-mean_b)/mean_s

		ss0,bb0 = AnaMax(mumu,toy0)
		ss1,bb1 = AnaMax(mumu,toy1) # Ist das so jetzt richtig?

		return [-2*Log(ss0,bb0,mumu,toy0)+2*Log(mean_s,mean_b,muMax0,toy0),-2*Log(ss1,bb1,mumu,toy1)+2*Log(mean_s,mean_b,muMax1,toy1)]

# observed value of test statistics
def q_obs(mumu):
	ss,bb=AnaMax(mumu,N)

	# constraints from the procedure paper:
	if (N-mean_b)/mean_s < 0.0:
		muMax=0.0
	elif (N-mean_b)/mean_s > 1.0:
		muMax=1.0
	else:
		muMax=(N-mean_b)/mean_s

	return -2*Log(ss,bb,mumu,N)+2*Log(mean_s,mean_b,muMax,N)

# create toy sample
sample0=[]
for q in xrange(0,100000):
	sample0.append(Profile(1.0)[0])

sample1=[]
for q in xrange(0,100000):
	sample1.append(Profile(1.0)[1])

plt.hist(sample0,50,histtype='step')
plt.hist(sample1,50,histtype='step')
plt.axvline(q_obs(0.0), color='blue', lw=2)
plt.axvline(q_obs(1.0), color='green', lw=2)
plt.show()
print '_______'
print q_obs(1.0),q_obs(0.0)
print '_______'
print Posterior()


end_time=strftime("%Y-%m-%d %H:%M:%S", gmtime())
print start_time
print end_time
