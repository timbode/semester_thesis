import ROOT
import numpy
from numpy import *
import scipy
from scipy import integrate
from scipy.stats import poisson, norm
import matplotlib.pyplot as plt



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

def extended_poisson(d,Nk):
	if d > 0:
		return poisson(d).pmf(Nk)
	else:
		return 0
	

# get marginalised likelihood in bin k --------------- arguments: k: bin --- dk: histogram (b or s) --- Nk: data counting rate --- intLimits: mean plus/minus intLimits times sigma
def Lk(k,dk,Nk,intLimits,prior):
	# "prior" times "likelihood conditioned on d"

	if prior==0:
		# Gaussian
		pLdk= lambda d: norm.pdf(d,dk.GetBinContent(k),dk.GetBinError(k))*extended_poisson(d,Nk) # Der Gaussian kann NEGATIVE WERTE annehmen... In den bins, auf denen er das tut, ist die Wahrscheinlichkeit aber verschwindend... Also okay, erstmal weiterzumachen
	else:
		# uniform
		pLdk= lambda d: extended_poisson(d,Nk)
		# Log-normal
		#pLdk= lambda d: (1/d)*norm.pdf(log(d),log(dk.GetBinContent(k)),(dk.GetBinError(k))/(dk.GetBinContent(k)))*poisson(d).pmf(Nk)
	
	# marginalised likelihood
	Lk=integrate.quad(pLdk,dk.GetBinContent(k)-intLimits*dk.GetBinError(k),dk.GetBinContent(k)+intLimits*dk.GetBinError(k))
	return Lk[0] # integrate.quad returns a list=[value,error]

# get the product of the bin likelihoods
def L(dk,data,intLimits,prior):
	L=1.0
	for i in xrange(1,(dk.GetNbinsX())+1):
		L*=Lk(i,dk,data.GetBinContent(i),intLimits,prior)# ACHTUNG: Lk vertraegt als drittes Argument nur integers - wegen der bloeden Poisson-Fkt.
	return L

def readDC(datacard):
	# path to the data card folder
	path="/home/lappet/Desktop/Semesterarbeit/datacards/125/"
	
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


def GetFromDC(mu,datacard):

	stats=readDC(datacard)
	tfile=ROOT.TFile(stats[0]+stats[1])
	tfile.cd((stats[2]))
	processes=stats[3]	
	
	print processes
	histo=loadHisto("Zj2b")
	print histo.GetBinContent(15)
	signal=loadHisto("ZH")
	print signal.GetBinContent(15)
	# data
	data=loadHisto("data_obs")
	data.SetMarkerStyle(11)
	data.SetLineColor(1)

	s=addHistos(processes,mu)
	b=addHistos(processes,0)


	# required for pdf labelling...
	mass=stats[1].split('_')
	mass=mass[3]

	# signal + background
	c1=ROOT.TCanvas("c1","signal + background")
	c1.SetLogy()
	c1.SetFillStyle(4000)
	sb=stack("signal + background",processes) # ACHTUNG: y range fixen --- bis 10^(-2)
	
	sb.Draw("hist")
	data.Draw("same e1 x0")
	c1.Print(mass+"_"+stats[2]+"_s-b-data")

	'''
	# background
	c2=ROOT.TCanvas("c2","background")

	b.SetFillColor(31)
	b.SetLineColor(31)	
	b.Draw("hist")
	c2.Print("b")

	# signal + background
	c3=ROOT.TCanvas("c2","background")
	
	s.SetFillColor(31)
	s.SetLineColor(31)	
	s.Draw("hist")
	c3.Print("s-b")
	'''

	intLimits=3
	return L(s,data,intLimits,0)
	
def test():
	GetFromDC(1,"vhbb_DC_TH_BDT_M125_ZeeHighPt_8TeV.txt")

	'''
	# M110_Zee (HighPt and LowPt) * M110_Zmm (HighPt and LowPt)
	func= lambda mumu: GetFromDC(mumu,"vhbb_DC_TH_BDT_M110_ZeeHighPt_8TeV.txt")*GetFromDC(mumu,"vhbb_DC_TH_BDT_M110_ZeeLowPt_8TeV.txt")*GetFromDC(mumu,"vhbb_DC_TH_BDT_M110_ZmmHighPt_8TeV.txt")*GetFromDC(mumu,"vhbb_DC_TH_BDT_M110_ZmmLowPt_8TeV.txt")
	x=arange(-2.5,5.5,0.05)
	y=map(func,arange(-2.5,5.5,0.05))
	#print x,y
	plt.plot(x,y)
	plt.show()
	
	
	mu=0
	for i in xrange(125,130,5):
		a=GetFromDC(mu,"vhbb_DC_TH_BDT_M"+str(i)+"_ZeeHighPt_8TeV.txt")*GetFromDC(mu,"vhbb_DC_TH_BDT_M"+str(i)+"_ZeeLowPt_8TeV.txt")
		b=GetFromDC(mu,"vhbb_DC_TH_BDT_M"+str(i)+"_ZmmHighPt_8TeV.txt")*GetFromDC(mu,"vhbb_DC_TH_BDT_M"+str(i)+"_ZmmLowPt_8TeV.txt")
		print "M"+str(i)+"_8TeV:", a*b
	'''
test()
