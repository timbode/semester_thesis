import ROOT
import numpy
from numpy import *
import scipy
from scipy import integrate
from scipy.stats import poisson, norm
import matplotlib.pyplot as plt

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

def stack(name,histos):
	colorDict=[2,401,41,5,922,920,596,840] #{'DYlight':401,'DYc':402,'DY2b':5,'DY1b':41,'TT':596,'ST':840,'VVlight':922,'VVb':920,'ZZb':18,'WZb':18,'ZH':2,'DYlc':401,'':1}
	hs=ROOT.THStack("hs",name)
	for i in xrange(len(histos)-1,-1,-1): #ACHTUNG: looping backwards
		#print histos[i], "in colour:", colorDict[i]
		h1=loadHisto(histos[i])
		h1.SetFillColor(colorDict[i])
		h1.SetLineColor(1)
		hs.Add(h1)
		
	return hs

def yieldBinContent(datacard):
	stats=readDC(datacard)
	tfile=ROOT.TFile(stats[0]+stats[1])
	tfile.cd((stats[2]))
	processes=stats[3]
	
	k=0
	while k < loadHisto(processes[0]).GetNbinsX():
		sk=loadHisto(processes[0]).GetBinContent(k+1)
		sk_e=loadHisto(processes[0]).GetBinError(k+1)
		bk=[]
		bk_e=[]
		i=0
		while i < len(processes)-1:
			bk.append(loadHisto(processes[i+1]).GetBinContent(k+1))
			bk_e.append(loadHisto(processes[i+1]).GetBinError(k+1))
			i+=1
		yield [sk,bk,sk_e,bk_e]
		k+=1

def CombineDCs(datacards):
	k=0
	while k < len(datacards):
		gen=yieldBinContent(datacards[k])
		for z in gen:
			yield z		
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

# L(N,mu|s1b1 ... snbn)=Pi(mu/s1b1 ... snbn)*f(N1|mu,s1b1)* ... *f(Nn|mu,snbn)
def L(N,mu,datacards):
	gen=CombineDCs(datacards)
	I=0
	k=0
	ff=1
	for x in gen:
		I+=fisher(mu,x[0],sum(x[1]))
		ff*=extended_poisson(N[k],mu,x[0],sum(x[1])) # k runs from 0 to 14, while the bins run from 1 to 15
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

# priors for sk, bk
def prior(mean,sigma):
	
# ------------------------------------

def PlotHisto(datacard):

	stats=readDC(datacard)
	tfile=ROOT.TFile(stats[0]+stats[1])
	tfile.cd((stats[2]))
	processes=stats[3]	

	# data
	data=loadHisto("data_obs")
	data.SetMarkerStyle(11)
	data.SetLineColor(1)

	# signal
	s=addHistos(processes,1)
	# background
	b=addHistos(processes,0)

	# required for pdf labelling...
	mass=stats[1].split('_')
	mass=mass[3]

	
	# signal + background
	c1=ROOT.TCanvas("c1","signal + background")
	c1.SetLogy()
	sb=stack("Signal, background & data",processes) # ACHTUNG: y range fixen --- bis 10^(-2)
	
	sb.Draw("hist")
	data.Draw("same e1 x0")
	c1.Print(mass+"_"+stats[2])

def test():
	datacard1="vhbb_DC_TH_BDT_M110_ZeeHighPt_8TeV.txt"
	datacard2="vhbb_DC_TH_BDT_M110_ZeeLowPt_8TeV.txt"
	datacards=[datacard1,datacard2]
	test=CombineDCs(datacards)
	for x in test:
		print x

	N=[2,46,150,300,400,350,280,170,100,70,50,20,10,3,0,20,100,220,380,550,572,630,500,390,270,180,120,56,20,15]
	print "-------"
	print L(N,0,datacards)
	print marginal(N,datacards)
	print posterior(0,N,datacards)

test()
