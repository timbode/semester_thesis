import ROOT
import numpy
from numpy import *
import scipy
from scipy import integrate
from scipy.stats import poisson, norm
import matplotlib.pyplot as plt

import pickle
import math


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
	for i in xrange(len(histos)):
		print histos[i]
		h1=loadHisto(histos[i])
		h1.SetFillColor(colorDict[i])
		h1.SetLineColor(1)
		hs.Add(h1)
		
	return hs

def extended_poisson(d,Nk):
	if d > 0:
		return poisson(d).pmf(Nk)
	else:
		return 0

# get marginalised likelihood in bin k --------------- arguments: k: bin --- dk: histogram (b or s) --- Nk: data counting rate --- intLimits: mean plus/minus intLimits times sigma
def Lk(k,dk,Nk,intLimits,prior):
	# "prior" times "likelihood conditioned on d"

	#print dk.GetBinContent(k)

	if prior==0:
		# Gaussian
		pLdk= lambda d: norm.pdf(d,dk.GetBinContent(k),dk.GetBinError(k))*extended_poisson(d,Nk) # Der Gaussian kann NEGATIVE WERTE annehmen... In den bins, auf denen er das tut, ist die Wahrscheinlichkeit aber verschwindend... Also okay, erstmal weiterzumachen
	else:
		# uniform
		pLdk= lambda d: extended_poisson(d,Nk)
		# Log-normal
		#pLdk= lambda d: (1/d)*norm.pdf(log(d),log(dk.GetBinContent(k)),(dk.GetBinError(k))/(dk.GetBinContent(k)))*poisson(d).pmf(Nk)
	'''
	x=arange(dk.GetBinContent(k)-intLimits*dk.GetBinError(k),dk.GetBinContent(k)+intLimits*dk.GetBinError(k),0.5)
	y=map(pLdk,x)
	for i in xrange(len(x)):
		if math.isnan(y[i])==1:
			print x[i], y[i]#FEHLER STECKT IN DER POISSON-FKT
	'''
	# marginalised likelihood
	Lk=integrate.quad(pLdk,dk.GetBinContent(k)-intLimits*dk.GetBinError(k),dk.GetBinContent(k)+intLimits*dk.GetBinError(k))
	#if (dk.GetBinContent(k) < 0):
		#print dk.GetBinContent(k), Lk[0]
	return Lk[0] # integrate.quad returns a list=[value,error]

# get the product of the bin likelihoods
def L(dk,data,intLimits,prior):
	L=1.0
	for i in xrange(1,(dk.GetNbinsX())+1):
		L*=Lk(i,dk,data.GetBinContent(i),intLimits,prior)# ACHTUNG: Lk vertraegt als drittes Argument nur integers - wegen der bloeden Poisson-Fkt.
	return L

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

def MakeToyData(hist,mu_given,processes,plot_counter,datacard):

	s_test=addHistos(processes,mu_given)

	# generate toy data and fill ROOT histo with toy data
	toy_data=zeros(s_test.GetNbinsX()+1)
	#toyData=ROOT.TH1F("toyData","testing mu",15,-1,1)
	for i in xrange(1,s_test.GetNbinsX()+1):
		toy_data[i]=random.poisson(s_test.GetBinContent(i))
		hist.SetBinContent(i,toy_data[i])
		#print toy_data[i], s_test.GetBinContent(i)
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
	return hist
	

def GetFromDC(mu,processes,toyData,datacard):

	s=addHistos(processes,mu)
	
	intLimits=3
	return L(s,toyData,intLimits,0)#ACHTUNG: PRIOR

def toy(mu_given,processes,plot_counter,datacard):

	hist=ROOT.TH1F("toyData","testing mu",15,-1,1)
	toyData=MakeToyData(hist,mu_given,processes,plot_counter,datacard)

	func= lambda mumu: GetFromDC(mumu,processes,toyData,"vhbb_DC_TH_BDT_M110_ZeeHighPt_8TeV.txt") #GRRRRRRRRRRRRRRRRR: toyData muessen integers sein... ----> diskrete Poisson-Fkt.
	y=map(func,arange(-2.5,5.5,0.05))
	return y


def test(datacard):

	stats=readDC(datacard)
	tfile=ROOT.TFile(stats[0]+stats[1])
	tfile.cd((stats[2]))
	processes=stats[3]
	
	x=arange(-2.5,5.5,0.05)
	plot_counter=0
	Y=zeros(len(x))
	for i in xrange(0,1000):
		plot_counter+=1
		Y+=toy(1,processes,plot_counter,"vhbb_DC_TH_BDT_M110_ZeeHighPt_8TeV.txt")
		#print Y
	pickle.dump(Y,open("1000mus.py","wb"))
	itemlist=pickle.load(open("1000mus.py","rb"))
	print itemlist
	print x[argmax(Y)]
	plt.plot(x,Y)
	plt.axvline(x=1, ymin=0, ymax=1)
	plt.show()
	
test("vhbb_DC_TH_BDT_M110_ZeeHighPt_8TeV.txt")
