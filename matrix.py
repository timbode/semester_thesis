# Author: Tim Lappe
import ROOT
import numpy
from numpy import *
import scipy
import math
from scipy import integrate
from scipy.stats import poisson, norm, lognorm
import matplotlib.pyplot as plt
import itertools
from time import gmtime, strftime
from progbar import progbar
import sys

start_time=strftime("%Y-%m-%d %H:%M:%S", gmtime())

# works only BELOW tfile.cd(...)!
def loadHisto(name):
	return ROOT.gDirectory.Get(name).Clone() # Klonen!


# Checklist for new datacards:
# 1: "process" appears twice (below, odd elements of "processes" are removed for this reason)
# 2: "shape" and "shapes" appear as they do... 
# 3: "bin" appears twice

def read(datacards):
	# path to the datacard folder
	path="datacards/125/"
	rootFile=[]
	channel=[]
	processes=[]
	lnNs=[]
	shapes=[]
	# the Wln are normalized...
	rates=[]
	# in the Znn datacards, "processes" appears thrice...
	for k, datacard in enumerate(datacards):
		c3=0
		with open(path+datacard) as f:
			for line in f:
				if line.startswith("#"):
					continue
			
				if "shapes" in line:
					rootFile.append(line.strip().split()[3])
					
				if "bin" in line:
					# the stat_bin diagrams...
					if not "shape" in line:
						channel.append(line.strip().split()[1])
			
				if "process" in line:
					# in the Znn datacards, "processes" appears thrice...
					if "Znunu" in datacard:
						if c3 == 0:
							c3+=1
							continue
							
					x=line.strip().split()
					# special case for "Wln" where no histogram "QCD" exists: remove the entire QCD process for "Wln" and "Znn"...
					# criterion to find out the datacards from which to remove "QCD":
					if len(x) > 10:
						x=x[:-1]
					processes.append(x)
					
				if "lnN" in line:
					lnNs.append(line.strip().split())
					
				if "shape" in line:
					if not "shapes" in line:
						if not "stat" in line:
							shapes.append(line.strip().split())
							
				if "rate" in line:
					rates.append(line.strip().split())
		f.close()
		
	# remove the second appearances of "bin"
	channel=[q for k, q in enumerate(channel) if k % 2 == 0] # AFTER the datacard loop
	
	# remove the second appearances of "process"
	processes=[q[1:] for k, q in enumerate(processes) if k % 2 == 0]
	
	# die Wmn-datacards sind fehlerhaft: die Zeile "bin" liest sich "Wmu" statt "Wmn" --- Korrektur allenfalls in "schneiderei.py"
	for k, chan in enumerate(channel):
		if "Wmu" in chan:
			channel[k]="Wmn"+chan[3:]
	# ------------------------------------------------------------------------------
	stats=[]
	for i, Bin in enumerate(channel):
		liste=[]
		for j in xrange(len(processes[i])):
			if ("Zee" in Bin or "Zmm" in Bin):
				liste.append("CMS_vhbb_stats_"+processes[i][j]+"_"+Bin)
			else:
				if ("Wen" in Bin or "Wmn" in Bin):
					if processes[i][j] == "s_Top":
						liste.append("CMS_vhbb_stat"+"sTop"+"_"+Bin)
					else:
						liste.append("CMS_vhbb_stat"+processes[i][j]+"_"+Bin)
				else:
					liste.append("CMS_vhbb_stat"+processes[i][j]+"_"+Bin)
		stats.append(liste)
	# ------------------------------------------------------------------------------
	
	lnNs=[q[0] for q in lnNs] # requires that the "lnN name" be the first element in the sublist/datacard line
	lnNs=sorted(list(set(lnNs)))
	
	shapes=[q[0] for q in shapes] # requires that the "shape name" be the first element in the sublist/datacard line
	shapes=sorted(list(set(shapes)))
	
	# get rid of most of the "unimportant" shape to minimize the marginalisation sum
	approx_shapes=[]
	for shape in shapes:
		if "eff_b" in shape:
			approx_shapes.append(shape)
		'''
		if "fake_b" in shape:
			approx_shapes.append(shape)
		
		if "res_j" in shape:
			approx_shapes.append(shape)
		if "scale_j" in shape:
			approx_shapes.append(shape)
		'''
	shapes=approx_shapes
	
	rates=[q[1:] for q in rates]
	
	# Integration
	IntVars={}
	
	# make lnN_matirx
	lnN_matrix={}
	for string in lnNs:
	
		# Integration
		IntVars[string]=1
		
		lnN_matrix[string]=[]
		for num, datacard in enumerate(datacards):
			with open(path+datacard) as f:
				c=0
				for line in f:
					if string in line:
						lnN_matrix[string].append(line.strip().split()[2:])
						c=1
			f.close()
			if c==0:
				lnN_matrix[string].append(["-" for q in xrange(len(processes[num]))])
				
	# make shape_matirx
	shape_matrix={}
	for string in shapes:
	
		# Integration
		IntVars[string]=0
		
		shape_matrix[string]=[]
		for num, datacard in enumerate(datacards):
			with open(path+datacard) as f:
				c=0
				for line in f:
					if string in line:
						shape_matrix[string].append(line.strip().split()[2:])
						c=1
			f.close()
			if c==0:
				shape_matrix[string].append(["-" for q in xrange(len(processes[num]))])
# ------------------------------------------------------------------------------	

	datacards_dict={q: k for k, q in enumerate(datacards)}
	
	PrettyHugeDict={}
	data_obsDict={}
	processes_dict={}
	for datacard in datacards:
		tfile=ROOT.TFile(path+rootFile[datacards_dict[datacard]])
		tfile.cd((channel[datacards_dict[datacard]]))
		
		# get the observed data
		data_obs=loadHisto("data_obs")
		data_obsDict[datacard]=[data_obs.GetBinContent(q) for q in xrange(1,data_obs.GetNbinsX()+1)]

		processes_dict[datacard]={q: k for k, q in enumerate(processes[datacards_dict[datacard]])}
	
		BigDict={}
		for process in processes[datacards_dict[datacard]]:
			# get the (process) histo
			hist=loadHisto(process)
			if ("Wen" in datacard or "Wmn" in datacard):
				hist.Scale(float(rates[datacards_dict[datacard]][processes_dict[datacard][process]]))
	
			# get the statistical histos
			if ("Zee" in stats[datacards_dict[datacard]][processes_dict[datacard][process]] or "Zmm" in stats[datacards_dict[datacard]][processes_dict[datacard][process]]):
				h_Down=loadHisto(process+stats[datacards_dict[datacard]][processes_dict[datacard][process]]+"Down") 
				h_Up=loadHisto(process+stats[datacards_dict[datacard]][processes_dict[datacard][process]]+"Up")
			else:
				h_Down=loadHisto(process+"_"+stats[datacards_dict[datacard]][processes_dict[datacard][process]]+"Down") 
				h_Up=loadHisto(process+"_"+stats[datacards_dict[datacard]][processes_dict[datacard][process]]+"Up")
				if ("Wen" in datacard or "Wmn" in datacard):
					h_Down.Scale(float(rates[datacards_dict[datacard]][processes_dict[datacard][process]]))
					h_Up.Scale(float(rates[datacards_dict[datacard]][processes_dict[datacard][process]]))
	
			# beim ersten Durchlauf von "string in shapes" die statistischen Werte (fuer die statistischen sigmas) mit rausholen -> c
			#
			# ACHTUNG: Ordnung der sublists: Down, "middle" (hist), Up
			c=0
			Dict={}
			Dict[stats[datacards_dict[datacard]][processes_dict[datacard][process]]]=[]
			for string in shapes:
				# QCD hat z. B. bei allen shapes "-" stehen - der loop fuer die stats muss aber auf jeden Fall ausgefuehrt werden, darf also nicht innerhalb der folgenden if-Klausel stehen!
				# Koennte sein, dass bei der dadurch noetig gewordenen Umstellung ein Fehler hineingekommen ist...
				shape_entry=shape_matrix[string][datacards_dict[datacard]][processes_dict[datacard][process]]
				if (shape_entry == "1.0" or shape_entry == "1.00"):
					if ("Zee" in stats[datacards_dict[datacard]][processes_dict[datacard][process]] or "Zmm" in stats[datacards_dict[datacard]][processes_dict[datacard][process]]):
						H_Down=loadHisto(process+string+"Down")
						H_Up=loadHisto(process+string+"Up")
					else:
						H_Down=loadHisto(process+"_"+string+"Down")
						H_Up=loadHisto(process+"_"+string+"Up")
						if ("Wen" in datacard or "Wmn" in datacard):
							H_Down.Scale(float(rates[datacards_dict[datacard]][processes_dict[datacard][process]]))
							H_Up.Scale(float(rates[datacards_dict[datacard]][processes_dict[datacard][process]]))
			
					liste=[]
				for k in xrange(0,data_obs.GetNbinsX()):
					if c < data_obs.GetNbinsX():
						Dict[stats[datacards_dict[datacard]][processes_dict[datacard][process]]].append([h_Down.GetBinContent(k+1),hist.GetBinContent(k+1),h_Up.GetBinContent(k+1)])
						# Integration:
						# only for bins with non-zero bin content
						# bei Wln sind die Ups und Downs auf den Bins durcheinandergeworfen! Der Abstand scheint aber wieder zu stimmen --- abs()
						if abs(hist.GetBinContent(k+1) - h_Down.GetBinContent(k+1)) > 1:#0.1:
							# Set the perliminary value to the upper integration boundary (4 sigma) --- now obsolete...
							IntVars[datacard+"-"+process+"-"+str(k+1)]=hist.GetBinContent(k+1)+ 4*(hist.GetBinContent(k+1) - h_Down.GetBinContent(k+1))
						else:
							IntVars[datacard+"-"+process+"-"+str(k+1)]="-"
								
								# Beachte folgenden Spezialfall:
								# bin content ist null, aber die Systematiken sind es nicht!
								# In diesem Fall gibt es keine statistische Integration,
								# sondern nur eine Summation ueber die drei Moeglichkeiten...
								# zk muss dann "direkt" jeweils einen der drei Werte annehmen.
								
					if (shape_entry == "1.0" or shape_entry == "1.00"):
						nochneliste=[]
						nochneliste.append(H_Down.GetBinContent(k+1))
						nochneliste.append(hist.GetBinContent(k+1))
						nochneliste.append(H_Up.GetBinContent(k+1))		
						liste.append(nochneliste)
						
					c+=1
					
				if (shape_entry == "1.0" or shape_entry == "1.00"):
					Dict[string]=liste
			BigDict[process]=Dict
		PrettyHugeDict[datacard]=BigDict
		
	return [PrettyHugeDict,lnN_matrix,IntVars,datacards_dict,processes_dict,stats,data_obsDict,shapes]
	
	
# new functions --------------------------------------------
# f(Nk|mu,sk,bk)
def ExtendedPoisson(Nk,mu,sk,bk):
	if (mu*sk+bk) > 0:
		return poisson(mu*sk+bk).pmf(Nk)
	else:
		return 0
		
def ExtendedLog(x):
	if x <= 0:
		print bla
		
# I(mu) for bin k:
def Fisher(mu,sk,bk):
	return sk*sk/(mu*sk+bk)

# **kwargs contains the to-be-integrated variables
# Structure of **kwargs:
# see IntVars
def P(mu,data_obsDict,PrettyHugeDict,lnN_matrix,datacards_dict,processes_dict,stats,**kwargs):
	ff=0
	I=0
	c=0
	plus=0
	for datacard in PrettyHugeDict.iterkeys():
		# ACHTUNG: "for k" and "for processes" have been interchanged...
		for k in xrange(0,len(data_obsDict[datacard])):# is there always a process "data_obs" which gives the representative bin number?
			sk=0
			bk=0
			for process in PrettyHugeDict[datacard].iterkeys():
				bin_k=PrettyHugeDict[datacard][process][stats[datacards_dict[datacard]][processes_dict[datacard][process]]][k]
				# get the nominal value:
				mu_k=bin_k[1]
				
				# get the statistical sigma:
				# are the Ups and Downs symmetric around the nominal value? - First look: yes...
				
				# bei Wln sind die Ups und Downs auf den Bins durcheinandergeworfen! Der Abstand scheint aber wieder zu stimmen --- abs()
				sigma_k=abs(bin_k[1]-bin_k[0])
				
				# beim ersten Durchlauf die lnN-priors mit rausholen -> c
				scale_factors=1
				for scale in lnN_matrix.iterkeys():
					# lnN-priors dranmultiplizieren
					if c < len(lnN_matrix):
						ff+=log(lognorm.pdf(kwargs[scale],1))
						if lognorm.pdf(kwargs[scale],1) == 0.0:
							print lognorm.pdf(kwargs[scale],1)
						c+=1
					# scale_factors generieren
					xxx=lnN_matrix[scale][datacards_dict[datacard]][processes_dict[datacard][process]]
					if xxx != "-":
						scale_factors*=power(kwargs[scale],(float(xxx) - 1.0))
				# ------------------------------------------------------------
				GBC_sum=0
				for syst in PrettyHugeDict[datacard][process]:
					if not "stat" in syst:
						# GBC(IntVars[syst])
						GBC_sum+=(PrettyHugeDict[datacard][process][syst][k][kwargs[syst]] - mu_k)
						
				# sum of the backgrounds (required for the likelihood)
				kwarg=kwargs[datacard+"-"+process+"-"+str(k+1)]
				if kwarg == "-":
					if process != ("ZH" or "WH"):
						bk+=scale_factors*(mu_k + GBC_sum)
					else:
						sk+=scale_factors*(mu_k + GBC_sum)
				else: # if kwargs[datacard+"-"+process+"-"+str(k+1)] != "-": --- # if sigma_k > 0.1:
					ff+=log(norm.pdf(kwarg,scale_factors*(mu_k + GBC_sum),sigma_k))
				
					if process != ("ZH" or "WH"):
						bk+=kwarg
					else:
						sk+=kwarg
				# ------------------------------------------------------------
			#sk und bk verarbeiten
			ff+=log(ExtendedPoisson(data_obsDict[datacard][k],mu,sk,bk))
			
			# Jeffreys-Prior vorbereiten
			I+=Fisher(mu,sk,bk)
			
	# Jeffreys-Prior
	# I kann stellenweise negativ werden - sollte aber keine Rolle spielen, weil der Poissonian dort eh verschwindet...			
	I=sqrt(abs(I))
	
	# Comment out for flat prior:
	#ff+=log(I)
	
	if ff == -float("Inf"):
		return 0.0
	else:
		return -ff#exp(ff+4300)
# ------------------------------------------------------------------------------

# lnNs und stats durch Integration marginalisieren:
# Hilfsfunktionen:
def Assign(liste,Keys,shapes,**IntVars):
	k=0
	for key in Keys:
		if not key in shapes:
			if IntVars[key] != "-":
				#if IntVars[key] != 1: # removes the lnN variables... Must be undone later!
				IntVars[key]=liste[k]
				k+=1
	return IntVars

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#    Author: Tristan Snowsill
#    License: Creative Commons Attribution-Noncommercial-Share Alike license
#    Package Index Owner: tristan.snowsill
#    DOAP record: mcint-0.1dev5.xml

def integrate(x,shapes,mu,data_obsDict,PrettyHugeDict,lnN_matrix,datacards_dict,processes_dict,stats,integrand,Keys,sampler,measure=1.0,n=10,**IntVars):# ACHTUNG: Unordnung in der Argumentliste (wegen Syntax)
    # Sum elements and elements squared
    total = 0.0
    total_sq = 0.0
    for y in itertools.islice(sampler, n):
        f = integrand(x,y,shapes,mu,data_obsDict,PrettyHugeDict,lnN_matrix,datacards_dict,processes_dict,stats,**IntVars)
        total += f
        total_sq += (f**2)
    # Return answer
    sample_mean = total/n
    sample_var = (total_sq - ((total/n)**2)/n)/(n-1.0)
    return (measure*sample_mean, measure*math.sqrt(sample_var/n))
    
def OuterIntegration(i,shapes,mu,data_obsDict,PrettyHugeDict,lnN_matrix,datacards_dict,processes_dict,stats,integrand,Keys,sampler2,measure=1.0,n=10,**IntVars):
    # Sum elements and elements squared
    total = 0.0
    total_sq = 0.0
    for x in itertools.islice(sampler2, n):
    	num,a,b=Boundaries(i,x,shapes,Keys,lnN_keys,lnN_matrix,datacards_dict,processes_dict,PrettyHugeDict)
        f = integrate(x,shapes,mu,data_obsDict,PrettyHugeDict,lnN_matrix,datacards_dict,processes_dict,stats,integrand,Keys,sampler(num,a,b),measure=1.0,n=100,**IntVars)[0]
        total += f
        total_sq += (f**2)
    # Return answer
    sample_mean = total/n
    sample_var = (total_sq - ((total/n)**2)/n)/(n-1.0)
    return (measure*sample_mean, measure*math.sqrt(sample_var/n))

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def integrand(x,z,shapes,mu,data_obsDict,PrettyHugeDict,lnN_matrix,datacards_dict,processes_dict,stats,**IntVars):
	return P(mu,data_obsDict,PrettyHugeDict,lnN_matrix,datacards_dict,processes_dict,stats,**Assign(append(x,z),Keys,shapes,**IntVars))

def sampler(num,a,b):
	while True:
		z=zeros(num)
		for k in xrange(0,num):
			z[k]=random.uniform(a[k],b[k])
                yield z
                
def Boundaries(i,x,shapes,Keys,lnN_keys,lnN_matrix,datacards_dict,processes_dict,PrettyHugeDict):
	num=0
	a=[]
	b=[]
	Measure=1
	cc=0
	for key in Keys:
		if not key in shapes:
			if IntVars[key] != "-":
				if IntVars[key] != 1:
					scale_factors=1
					for k, lnN_key in enumerate(lnN_keys):
						scale_factor=lnN_matrix[lnN_key][datacards_dict[key.split(".")[0]+".txt"]][processes_dict[key.split(".")[0]+".txt"][key.split(".")[1].split("-")[1]]]
						if scale_factor != "-":
						# spaeter IntVars[lnN_key] durch x ersetzen!!! 
							scale_factors*=power(x[k],(float(scale_factor) - 1.0))
				
					num+=1
					shift=0.0
					for string in PrettyHugeDict[key.split(".")[0]+".txt"][key.split(".")[1].split("-")[1]]:
						if "stat" in string:
							# Beachte Indexverschiebung um -1: die keys sind wie die bins von 1 bis 15 nummeriert
							xxx=PrettyHugeDict[key.split(".")[0]+".txt"][key.split(".")[1].split("-")[1]][string][int(key.split(".")[1].split("-")[2]) - 1]
							mean=xxx[1]
							sigma=xxx[1] - xxx[0]
							sigma=abs(sigma)
						else:
							# The mean must be given by "GBC_sum" rather than mu_k!
							xxx=PrettyHugeDict[key.split(".")[0]+".txt"][key.split(".")[1].split("-")[1]][string][int(key.split(".")[1].split("-")[2]) - 1]
							shift+=xxx[i[shapes.index(string)]] - xxx[1]
							
							# Set the integration variable:
							IntVars[string]=i[shapes.index(string)]
							
					upper=(scale_factors*(mean+shift)) + 4*sigma
					lower=(scale_factors*(mean+shift)) - 4*sigma
					a.append(lower)
					b.append(upper)
					if Measure > 1e+100:
						Measure*=1e-100
						cc+=1
					Measure*=upper - lower
	#print num, Measure, cc
	return [num,a,b]					
	
# ------------------------------------------------------------------------------

mH="M125"
#datacards=["vhbb_DC_TH_BDT_"+mH+"_ZnunuMedPt_8TeV.txt"]
#datacards=["vhbb_DC_TH_BDT_"+mH+"_ZeeHighPt_8TeV.txt"]
#datacards=["vhbb_DC_TH_BDT_"+mH+"_ZnunuHighPt_8TeV.txt"]
datacards=["vhbb_DC_BDT_"+mH+"_WmnHighPt_8TeV.txt","vhbb_DC_TH_BDT_"+mH+"_ZnunuHighPt_8TeV.txt","vhbb_DC_TH_BDT_"+mH+"_ZmmLowPt_8TeV.txt","vhbb_DC_TH_BDT_"+mH+"_ZnunuLowPt_8TeV.txt","vhbb_DC_BDT_"+mH+"_WenHighPt_8TeV.txt","vhbb_DC_BDT_"+mH+"_WenLowPt_8TeV.txt","vhbb_DC_BDT_"+mH+"_WenMidPt_8TeV.txt","vhbb_DC_BDT_"+mH+"_WmnLowPt_8TeV.txt","vhbb_DC_BDT_"+mH+"_WmnMidPt_8TeV.txt","vhbb_DC_TH_BDT_"+mH+"_ZeeHighPt_8TeV.txt","vhbb_DC_TH_BDT_"+mH+"_ZeeLowPt_8TeV.txt","vhbb_DC_TH_BDT_"+mH+"_ZmmHighPt_8TeV.txt","vhbb_DC_TH_BDT_"+mH+"_ZnunuMedPt_8TeV.txt"]
readout=read(datacards)
PrettyHugeDict=readout[0]
lnN_matrix=readout[1]
IntVars=readout[2]
datacards_dict=readout[3]
processes_dict=readout[4]
stats=readout[5]
data_obsDict=readout[6]
shapes=readout[7]

#print PrettyHugeDict["vhbb_DC_TH_BDT_"+mH+"_ZeeLowPt_8TeV.txt"]["Zj1b"]["CMS_vhbb_stats_Zj1b_ZeeLowPt_8TeV"][7]

# Liste der keys (um fixe Reihenfolge festzulegen):
Keys=[]
lnN_keys=[]
for key in IntVars.iterkeys():
	for string in lnN_matrix.iterkeys():
		if key == string:
			lnN_keys.append(key)
		else:
			Keys.append(key)
Keys=sorted(set(Keys))

for key in IntVars.iterkeys():
	for string in lnN_matrix.iterkeys():
		if key == string:
			Keys.remove(key)

for k, key in enumerate(lnN_keys):
	Keys.insert(k,key)


# Test:
mu=float(sys.argv[1])
print "mu=", mu

# Integration boundaries of the lognormals:
A=zeros(len(lnN_keys))
B=A+10
Num=len(A)

pb=progbar(3)

S=0
for i in itertools.product(xrange(3),repeat=len(shapes)):
	S+=OuterIntegration(i,shapes,mu,data_obsDict,PrettyHugeDict,lnN_matrix,datacards_dict,processes_dict,stats,integrand,Keys,sampler(Num,A,B),measure=1.0,n=100,**IntVars)[0]
	pb.move()
	
print "mu=", mu
print "Flat prior!"
print mH
print "N=", 100
print "n=", 100
print S

print start_time
print strftime("%Y-%m-%d %H:%M:%S", gmtime())
# ------------------------------------------------------------------------------
