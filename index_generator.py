import os,glob,re

notebooks=[]
directories = glob.glob(os.getcwd()+'/**'+'/*.ipynb',recursive=True)
names = []
areaofphys = []
descriptions = []

for i in range(0,len(directories)):
	redName=None
	redDesc=None
	with open(directories[i]) as f:
		content=f.readlines()
	for j in range(0,len(content)):
		if(content[j].find("#NAME:")!=-1):
			redName=content[j][content[j].find("#NAME:")+len("#NAME:"):]
			names.append(redName[0:redName.find("\"")].strip())
			break
	if redName==None:
		names.append("")
	for j in range(0,len(content)):
		if(content[j].find("#DESCRIPTION:")!=-1):
			redDesc=content[j][content[j].find("#DESCRIPTION:")+len("#DESCRIPTION:"):]
			descriptions.append(redDesc[0:redDesc.find("\"")].strip())
			break
	if redDesc==None:
		descriptions.append("")
	directories[i]="."+directories[i][len(os.getcwd()):]
	if directories[i] ==  "./index.ipynb":
		iIndex=i
	#fails if ipynb isn't in a sub directory
	areaofphys.append(directories[i][2:])
	areaofphys[i]=areaofphys[i][0:areaofphys[i].find("/")]
	areaofphys[i]=re.sub(r"(\w)([A-Z])", r"\1 \2", areaofphys[i])

for i in range(0,len(directories)):
	if names[i]=='':
		pass
	else:
		notebooks.append([names[i],areaofphys[i],descriptions[i],directories[i]])

notebooks.sort(key=lambda x: x[1])

indexNotebook = open("indexgen.ipynb",'w')
with open("./.indexraw.txt","r") as p:
	lines=p.readlines()

insertFrom=None
for k in range(0,len(lines)):
	if(lines[k].find("# PyCav Demo Index")!=-1):
		insertFrom=k+1
		break
lineset1=lines[0:k+1]
lineset2=lines[k+1:]
prevSection=None
for nb in notebooks:
	if(prevSection!=str(nb[1])):
		lineset1.append("\"\\n\",\n")
		lineset1.append("\"## "+str(nb[1])+"\\n\",\n")
		lineset1.append("\"\\n\",\n")
	lineset1.append("\"\\n\",\n")
	lineset1.append("\"["+str(nb[0])+"]("+str(nb[3])+"): "+str(nb[2])+"\\n\",\n")
	lineset1.append("\"\\n\",\n")
	prevSection=str(nb[1])
lineset1.append("\"\\n\"\n")
lines=lineset1+lineset2
for l in lines:
	indexNotebook.write(l)
indexNotebook.close()
