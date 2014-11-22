def Schneiderei(datacard):
	# path to the data card folder
	path="/home/lappet/Desktop/Semesterarbeit/datacards/Wln/"

	head=[]
	shapes=[]
	bins=[]
	obs=[]
	cc=0
	Lines=[]
	with open(path+datacard,"r+") as f:
		for line in f:
			if "max" in line:
				head.append(line.strip().split())
			
			if "shapes *" in line:
				shapes.append(line.strip().split())
				
			if "bin" in line:
				bins.append(line.strip().split())
				cc+=1
				
			if "observation" in line:
				obs.append(line.strip().split())
				
			if cc >= 2:
				if not line.startswith("-"):
					Lines.append(line.strip().split())
	f.close()
	
	bins=bins[:2]
	bins[1]=bins[1][1:]
	obs=obs[0]

	liste=[[],[],[]]

	
	for k in xrange(0,4):
		prepend=Lines[k][0]
		Lines[k]=Lines[k][1:]
		for i in xrange(1,4):
			xx=Lines[k][bins[1].index(bins[0][i]):(len(Lines[k]) - 1) - bins[1][::-1].index(bins[0][i]) + 1]
			xx.insert(0,prepend)
			liste[i-1].append(xx)
	
	for k in xrange(4,len(Lines)):
		prepend0=Lines[k][0]
		prepend1=Lines[k][1]
		Lines[k]=Lines[k][2:]
		for i in xrange(1,4):
			xx=Lines[k][bins[1].index(bins[0][i]):(len(Lines[k]) - 1) - bins[1][::-1].index(bins[0][i]) + 1]
			xx.insert(0,prepend0)
			xx.insert(1,prepend1)
			liste[i-1].append(xx)
			
	Pts=["LowPt","MidPt","HighPt"]
	for k in xrange(0,3):
		with open(path+datacard[0:12]+"M"+datacard[13:20]+Pts[k]+datacard[20:29],"w") as f:
			for line in head:
				for string in line:
					f.write(string+"\t")
				f.write("\n")
			f.write("\n")
		
			for string in shapes[k]:
				if "vhbb_WS" in string:
					f.write(string[:5]+"TH"+string[7:12]+"M"+string[13:]+"\t")
				else:
					f.write(string+"\t")
			f.write("\n")
		
			f.write("\n")
		
			f.write(bins[0][0]+"\t")
			f.write(bins[0][k+1][0:3]+Pts[k]+datacard[20:25]+"\t")
			f.write("\n")
		
			f.write("\n")
		
			f.write(obs[0]+"\t")
			f.write(obs[k+1]+"\n")
		
			f.write("\n")
		
			for line in liste[k]:
				for string in line:
					f.write(string+"\t")
				f.write("\n")

def main():
	masses=[110,115,120,125,130,135,140,145,150]
	for mass in masses:
		Schneiderei("vhbb_DC_BDT_H"+str(mass)+"_Wmn_8TeV.txt")

main()
