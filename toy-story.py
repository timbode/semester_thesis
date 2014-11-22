import ROOT
import numpy
from numpy import *
import scipy
from scipy import integrate
from scipy.stats import poisson, norm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def readDC(datacard):
	# path to the data card folder
	path="/home/lappet/Desktop/BDTBlindFinalBoxOpening/"
	
	with open(path+datacard) as f:
		
		processes=[] # declare 'processes' as empty (needed below)
		for line in f:
			if "shapes" in line:
				rootFile=line.strip().split()
			if "bin" in line:
				channel=line.strip().split()
			if "process" in line:
				if not processes: # use emptiness of 'processes' to avoid overwriting on appearance of 2nd string "process"
					processes=line.strip().split()
					processes=processes[1:]
		f.close()


	
	stats=[path,rootFile[3],channel[1],processes]
	return stats

# works only BELOW tfile.cd(...)!
def loadHisto(name):
	return ROOT.gDirectory.Get(name).Clone() # Klonen!

def addHistos(arg,mu):
	h=loadHisto(arg[0])
	h.Scale(mu)
	#h.Sumw2() # run time warning
	for i in xrange(len(arg)-1):
		h.Add(loadHisto(arg[i+1]))
	return h

def yieldBinContent(datacard):
	stats=readDC(datacard)
	tfile=ROOT.TFile(stats[0]+stats[1])
	tfile.cd((stats[2]))
	processes=stats[3]
	
	k=0
	while k < loadHisto(processes[0]).GetNbinsX():
		sk=loadHisto(processes[0]).GetBinContent(k+1)
		b=[]
		i=0
		while i < len(processes)-1:
			b.append(loadHisto(processes[i+1]).GetBinContent(k+1))
			i+=1
		bk=sum(b)
		yield [sk,bk]
		k+=1

# ------------------------------------
# probabilistic functions

# f(Nk|mu,sk,bk)
def extended_poisson(Nk,mu,sk,bk):
	if (mu*sk+bk) > 0:
		return poisson(mu*sk+bk).pmf(Nk)
	else:
		return 0

# I(mu) for bin k:
def fisher(mu,sk,bk):
	return power(sk,2)/(mu*sk+bk)

# L(Nk,mu|s1b1 ... snbn)=Pi(mu/s1b1 ... snbn)*f(N1|mu,s1b1)* ... *f(Nn|mu,snbn)
def L(N,mu,datacard):
	var=yieldBinContent(datacard)
	I=0
	k=0
	ff=1
	for x in var:
		I+=fisher(mu,x[0],x[1])
		ff*=extended_poisson(N[k],mu,x[0],x[1]) # k runs from 0 to 14, while the bins run from 1 to 15
		k+=1
	# returns the prior, the product of f-likelihoods AND the probuct of these
	return [sqrt(I),ff,sqrt(I)*ff]

# marginalization
def marginal(N,datacard):
	func= lambda mumu: L(N,mumu,datacard)[2]
	Lm=integrate.quad(func,0,1)
	return Lm[0] # integrate.quad() returns [value,error]

# P(mu|N=N1 ... Nn, s1b1 ... snbn)
def posterior(mu,N,datacard):
	return L(N,mu,datacard)[2]/marginal(N,datacard)

# ------------------------------------

def MakeToyData(mu_given,plot_counter,datacard):
	stats=readDC(datacard)
	tfile=ROOT.TFile(stats[0]+stats[1])
	tfile.cd((stats[2]))
	processes=stats[3]

	s_test=addHistos(processes,mu_given)
	
	# generate toy data and fill ROOT histo with toy data
	hist=ROOT.TH1F("toyData","testing mu",15,-1,1)
	toy_data=zeros(s_test.GetNbinsX()+1)
	for i in xrange(1,s_test.GetNbinsX()+1):
		toy_data[i]=random.poisson(s_test.GetBinContent(i))
		hist.SetBinContent(i,toy_data[i])
	'''
	# Plots
	c=ROOT.TCanvas("c","toy data")
	c.SetLogy()
	hist.Draw("e1 x0")
	hist.SetMarkerStyle(11)
	s_test.Draw("same hist")
	toy_string="toy"
	c.Print("toy_mu"+str(mu_given)+"_"+str(plot_counter))
	'''
	hist.SetDirectory(0)
	return toy_data[1:]


def test():
	datacard="vhbb_DC_TH_BDT_M110_ZeeHighPt_8TeV.txt"
	plot_counter=0
	
	
	var=yieldBinContent(datacard)
	k=0
	N1=0
	N0=0
	maxi=1000
	for i in xrange(0,maxi):
		N1+=MakeToyData(10,plot_counter,datacard)
		N0+=MakeToyData(9,plot_counter,datacard)

	N1=int_(N1/maxi)
	N0=int_(N0/maxi)
	print N1
	print N0
	#N1=MakeToyData(1,plot_counter,datacard)
	#N0=MakeToyData(0,plot_counter,datacard)
	
	for x in var:
		#print N[k], x[0],x[1]
		#print poisson(x[0]+x[1]).pmf(N[k]),poisson(x[1]).pmf(N[k])
		xx=arange(0,x[1]+max(x[1],30),1)
		func1= lambda Nk1: poisson(10*x[0]+x[1]).pmf(Nk1)
		func2= lambda Nk2: poisson(9*x[0]+x[1]).pmf(Nk2)
		#func= lambda mumu: extended_poisson(N[k],mumu,x[0],x[1])
		y1=map(func1,xx)
		y2=map(func2,xx)
		#print y
		plt.plot(xx,y1, label="mu=10")
		plt.plot(xx,y2, color='r', label="mu=9")
		plt.axvline(x=N1[k], ymin=0, ymax=1)
		plt.axvline(x=N0[k], ymin=0, ymax=1, color='r')
		plt.title("Bin "+str(k+1))
		plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
		plt.show()
		
		k+=1
	
	'''
	q=10
	xxx=arange(0,2.2*q,0.1)
	Y=zeros(len(xxx))
	for i in xrange(0,10000):
		plot_counter+=1
		N=MakeToyData(q,plot_counter,datacard)
		func= lambda mumu: L(N,mumu,datacard)[1]
		y=map(func,xxx)
		Y+=y
	print xxx[argmax(Y)]
	plt.plot(xxx,Y)
	plt.title("mu="+str(q))
	plt.axvline(x=xxx[argmax(Y)], ymin=0, ymax=1)
	plt.axvline(x=q, ymin=0, ymax=1, color='r')
	plt.savefig("likelihood_mu_"+str(q)+".png")
	plt.close()
	'''
test()
