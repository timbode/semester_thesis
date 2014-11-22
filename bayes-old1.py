import ROOT
import numpy
from numpy import *
import scipy
from scipy import integrate
from scipy.stats import poisson, norm
import matplotlib.pyplot as plt

def read2(datacard):
	# path to the data card folder
	path="/home/lappet/Desktop/BDTBlindFinalBoxOpening/"

	with open(path+datacard) as f:
		processes=[] # declare 'processes' as empty (needed below)
		lnNs=[]
		lnNs_dict={}
		stats=[]
		stats_dict={}
		shapes=[]
		shapes_dict={}
		i=0
		k=0
		kk=0
		for line in f:
			if "shapes" in line:
				rootFile=line.strip().split()
			if "bin" in line:
				channel=line.strip().split()
			if "process" in line:
				if not processes: # use emptiness of 'processes' to avoid overwriting on appearance of 2nd string "process"
					processes=line.strip().split()
					processes=processes[1:]
			
			if "lnN" in line:
				lnNs.append(line.strip().split())
				lnNs_dict[lnNs[i][0]]=lnNs[i][2:]
				lnNs[i]=lnNs[i][:1]
				i+=1

			if "shape" in line:
				if "stats" in line:
					stats.append(line.strip().split())
					stats_dict[stats[k][0]]=stats[k][2:]
					stats[k]=stats[k][:1]
					k+=1
				else:
					if not "shapes" in line:
						shapes.append(line.strip().split())
						shapes_dict[shapes[kk][0]]=shapes[kk][2:]
						shapes[kk]=shapes[kk][:1]
						kk+=1
				
		f.close()
	
	lnNs=[q[0] for q in lnNs]

	stats=[q[0] for q in stats]
	shapes=[q[0] for q in shapes]
	
	return [path,rootFile[3],channel[1],processes,lnNs,lnNs_dict,stats,stats_dict,shapes,shapes_dict]

def dataa(datacard):
	read=read2(datacard)
	tfile=ROOT.TFile(read[0]+read[1])
	tfile.cd((read[2]))
	processes=read[3]
	stats=read[6]
	stats_dict=read[7]
	shapes=read[8]
	shapes_dict=read[9]
	
	
	# obsolete
	BinContents={}
	for name in processes: #for i, string in enumerate(shapes):
		h=loadHisto(name)
		liste=[]
		for k in xrange(0,h.GetNbinsX()):
			liste.append(h.GetBinContent(k+1))
		BinContents[name+"_BinContents"]=liste
	# obsolete
	BinErrors={}
	for name in processes: #for i, string in enumerate(shapes):
		h=loadHisto(name)
		liste=[]
		for k in xrange(0,h.GetNbinsX()):
			liste.append(h.GetBinError(k+1))
		BinErrors[name+"_BinErrors"]=liste

	'''
	test=10
	print loadHisto("ZHCMS_vhbb_stats_ZH_ZeeHighPt_8TeVUp").GetBinContent(test)
	print loadHisto("ZHCMS_vhbb_stats_ZH_ZeeHighPt_8TeVDown").GetBinContent(test)
	print loadHisto(processes[0]).GetBinError(test)
	print loadHisto(processes[0]).GetBinContent(test)
	print loadHisto(processes[0]).GetBinContent(test)-loadHisto("ZHCMS_vhbb_stats_ZH_ZeeHighPt_8TeVUp").GetBinContent(test)
	print loadHisto(processes[0]).GetBinContent(test)-loadHisto("ZHCMS_vhbb_stats_ZH_ZeeHighPt_8TeVDown").GetBinContent(test)
	'''
	
	vhbb_stats={}
	for i, string in enumerate(stats):
		h_Up=loadHisto(processes[i]+string+"Up")
		h_Down=loadHisto(processes[i]+string+"Down")
		liste=[]
		for k in xrange(0,h_Up.GetNbinsX()):
			nochneliste=[]
			nochneliste.append(h_Up.GetBinContent(k+1)) # +1-shift nicht vergessen!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			nochneliste.append(h_Down.GetBinContent(k+1))
			liste.append(nochneliste)
		vhbb_stats[processes[i]]=liste #ACHTUNG: naming dict entries after the process, NOT the name in the datacard
	
	shape_systematics={}
	for string in shapes:
		Dict={}
		for jj, strang in enumerate(processes):
			h_Up=loadHisto(strang+string+"Up")
			h_Down=loadHisto(strang+string+"Down")
			liste=[]
			for k in xrange(0,h_Up.GetNbinsX()):
				nochneliste=[]
				nochneliste.append(h_Up.GetBinContent(k+1)) # +1-shift nicht vergessen!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				nochneliste.append(h_Down.GetBinContent(k+1))
				liste.append(nochneliste)
			Dict[strang]=liste
		shape_systematics[string]=Dict
	return [vhbb_stats,shape_systematics]

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

# ------------------------------------
# probabilistic functions

# statistical prior for sk, bk
def prior(var,mean,sigma,shape):
	if mean==0:
		return 1
	if shape==1:
		return norm.pdf(var,mean,sigma)
	else:
		return 1/(2*sigma) #1/((mean+sigma)-(mean-sigma))

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
def L(N,mu,s,b,datacards): # s and b must have len(gen)
	gen=CombineDCs(datacards)
	k=0
	P=1
	I=0
	ff=1
	for z in gen:
		P*=prior(s[k],z[0],z[2],1)
		for i in xrange(len(b[k])):
			P*=prior(b[k][i],z[1][i],z[3][i],1) # completes the multiplication of statistical priors
		I+=fisher(mu,s[k],sum(b[k])) # ACHTUNG: b[k] must be a list in itself, containing all backgrounds on one bin
		ff*=extended_poisson(N[k],mu,s[k],sum(b[k])) 
		k+=1
	# returns the priors, the product of f-likelihoods AND the probuct of these
	return [[P,sqrt(I),ff],P*sqrt(I)*ff]

def what(liste):
	return power(liste[1],2)

def funcy(x,y):
	return y*x*x

def help(k,length,var):
	liste=zeros(length)
	liste[k]=var
	return liste

def marginalize(f,k,length):
	func= lambda var: f(help(k,length,var))
	Int=integrate.quad(func,0,1)
	return Int
	
def whoot(liste):
	ff=1
	for z in liste:
		ff*=z
	return ff

	
# ------------------------------------# ------------------------------------# ------------------------------------
# ------------------------------------# ------------------------------------# ------------------------------------
# new functions

# "liste" ist eine Liste der Systematiken...
def specialGaussian(z,liste,sigma):
	return norm.pdf(z,sum(liste),sigma)

# gibt das Produkt aus Prozess-prior und systematischen priors auf bin k+1 --- zk ist die Prozess-Variable, die statistisch zu integrieren ist
def systPriors(zk,liste,process,k,datacard):
	# get the statistical sigma
	sigma_k=abs(dataa(datacard)[0][process][k][0]-dataa(datacard)[0][process][k][1])/2
	print specialGaussian(zk,liste,sigma_k)
	
	# aus diesen up und down values muessen der Mittelwert und das Sigma fuer die Systematik-priors berechnet werden
	for string in dataa(datacard)[1]: #for i, j in dataa(datacard)[1].items():
		print string, dataa(datacard)[1][string][process][k]
	return 0
	
# ------------------------------------

def test():
	datacard1="vhbb_DC_TH_BDT_M110_ZeeHighPt_8TeV.txt"
	datacard2="vhbb_DC_TH_BDT_M110_ZeeLowPt_8TeV.txt"
	datacards=[datacard1,datacard2]

	'''
	N=[2,46,150,300,400,350,280,170,100,70,50,20,10,3,1,20,100,220,380,550,572,630,500,390,270,180,120,56,20,15]
	#print L(N,0,datacards)
	#print marginal(N,datacards)
	#print posterior(0,N,datacards)
	#print Pi(datacards)
	print marginalize(what,1,3)
	print integrate.quad(funcy,0,1,1)
	func= lambda x1,x2,x3: whoot([x1,x2,x3])
	print integrate.nquad(func,[[0,1],[0,1],[0,1]])
	'''
	num=0
	length=0
	for string, vals in read2(datacard1)[9].items():
		if vals[num]=="1.0":
			length+=1
	# Liste hat die Laenge, die der Anzahl der Integrationsvariablen entspricht (je nach dem, ob alle Systematiken einen Einfluss auf den Prozess "num" haben)
	liste=zeros(length)
	#print systPriors(1,liste,"ZH",1,datacard1)
	print dataa(datacard1)[1]

test()
