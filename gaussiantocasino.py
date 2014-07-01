#!/usr/bin/env python

import numpy, os, sys, argparse
from copy import deepcopy

#### parse the command line arguments

parser = argparse.ArgumentParser(description='Convert Gaussian09 output to CASINO input')
parser.add_argument('-fchk', default='gaussian.fchk',help='Input formatted Gaussian checkpoint file')
parser.add_argument('-gout', default='gaussian.out', help='File containing the stdout from a Gaussian run')
parser.add_argument('-gwfn', default='gwfn.data',help='name for the gwfn.data file, if different')
parser.add_argument('-corr', default='correlation.data',help='name for the correlation.data file, if different' )
parser.add_argument('-cut', default=1.0,type=float,help='Sum of squares of coefficients for truncating CIS/CASSCF expansion')
parser.add_argument('-maxdets',default = 99999,type=int,help='Maximum number of determinants to include')
parser.add_argument('-state',default = 1, type=int,help='Which CIS excited state to use') 
args = parser.parse_args()

## OPEN all of the relevant files ##

try:
	fchk = open(vars(args)['fchk'],'r')
except:
	print "File " + vars(args)['fchk'] + " does not exist!\n"
	quit()

try:
	output = open(vars(args)['gout'],'r')
except:
	print "File " + vars(args)['gout'] + " does not exist!\n"
	quit()

gwfn = open(vars(args)['gwfn'],'w')

## define useful functions

def normCoeff(coeff, expt, angmom):
	if ((angmom == 1) or (angmom == 2)): 
		normed = numpy.float(coeff) * ((2. * numpy.float(expt) / numpy.pi) ** (3./4.))
	elif (angmom == 3):
		normed = numpy.float(coeff) * ((2. * numpy.float(expt) / numpy.pi) ** (3./4.)) * 2. * (numpy.float(expt) ** 0.5) 
	elif (angmom == 4):
		normed = numpy.float(coeff) * ((2. * numpy.float(expt) / numpy.pi) ** (3./4.)) * 2. * numpy.float(expt)
	elif (angmom == 5):
		normed = numpy.float(coeff) * ((2. * numpy.float(expt) / numpy.pi) ** (3./4.)) * ((4 * numpy.float(expt)) ** 1.5) / numpy.sqrt(8.)
	return normed
	
def nnrepulsion(coordsinRows, atomnumbers):
	sum = 0.
	for i in range(len(atomnumbers) - 1):
		for j in range(i + 1,len(atomnumbers)):
			distance = numpy.sqrt(numpy.sum((numpy.asarray(coordsinRows[i]) - numpy.asarray(coordsinRows[j])) ** 2)) #/ 0.5291772109
			nn = numpy.float(atomnumbers[i]) * numpy.float(atomnumbers[j])
			sum += (nn / distance)
	sum /= len(atomnumbers)
	return sum

def parseCAS09(config, core, ground, detno, writeout,signlist):
	groundlist = list(ground)
	configlist = list(config)
	missingalpha = []
	addedalpha = []
	missingbeta = []
	addedbeta = []
	groundOrbUp = []
	groundOrbDown = []
	exOrbUp = []
	exOrbDown = []
	for i in range(len(configlist)):
		if (groundlist[i] != configlist[i]):
			if ((groundlist[i] == '0') and (configlist[i] == 'a')):
				addedalpha.append(i) 
				exOrbUp.append(i)
			elif ((groundlist[i] == '0') and (configlist[i] == 'b')):
				addedbeta.append(i)
				exOrbDown.append(i)
			elif ((groundlist[i] == '0') and (configlist[i] == '1')):
				addedalpha.append(i)
				addedbeta.append(i)
				exOrbUp.append(i)
				exOrbDown.append(i)
			elif ((groundlist[i] == 'a') and (configlist[i] == '0')):
                                missingalpha.append(i)
                                groundOrbUp.append(i)
                        elif ((groundlist[i] == 'a') and (configlist[i] == '1')):
                                addedbeta.append(i)
                                groundOrbUp.append(i)
                                exOrbUp.append(i)
                                exOrbDown.append(i)
			elif ((groundlist[i] == 'a') and (configlist[i] == 'b')):
                                addedbeta.append(i)
                                groundOrbUp.append(i)
				missingalpha.append(i)
				exOrbDown.append(i)
			elif ((groundlist[i] == 'b') and (configlist[i] == '0')):
                                missingbeta.append(i)
                                groundOrbDown.append(i)
                        elif ((groundlist[i] == 'b') and (configlist[i] == '1')):
                                addedalpha.append(i)
                                groundOrbDown.append(i)
                                exOrbUp.append(i)
                                exOrbDown.append(i)
                        elif ((groundlist[i] == 'b') and (configlist[i] == 'a')):
                                addedalpha.append(i)
				missingbeta.append(i)
				groundOrbDown.append(i)
				exOrbUp.append(i)
			elif ((groundlist[i] == '1') and (configlist[i] == '0')):
                                missingbeta.append(i)
				missingalpha.append(i)
				groundOrbUp.append(i)
				groundOrbDown.append(i)
				exOrbDown.append(i)
                        elif ((groundlist[i] == '1') and (configlist[i] == 'a')):
                                missingbeta.append(i)
                                groundOrbUp.append(i)
                                groundOrbDown.append(i)
                                exOrbUp.append(i)
			elif ((groundlist[i] == '1') and (configlist[i] == 'b')):
                                missingalpha.append(i)
                                groundOrbUp.append(i)
                                groundOrbDown.append(i)
                                exOrbDown.append(i)
			elif ((groundlist[1] == '1') and (configlist[i] == '1')):
				groundOrbUp.append(i)
				groundOrbDown.append(i)
                                exOrbUp.append(i)
                                exOrbDown.append(i)
                      	elif ((groundlist[1] == 'a') and (configlist[i] == 'a')):
				groundOrbUp.append(i)
                                exOrbUp.append(i)
			elif ((groundlist[1] == 'b') and (configlist[i] == 'b')):
				groundOrbDown.append(i)
                                exOrbDown.append(i)                           
	for i in range(len(missingalpha)):
		orbout = core + missingalpha[i] + 1
		orbin = core + addedalpha[i] + 1 
		writeout.append('DET ' + str(detno) + ' 1 PR ' + str(orbout) + ' 1 ' + str(orbin) + ' 1\n')
	for i in range(len(missingbeta)):
		orbout = core +	missingbeta[i]	+ 1
                orbin =	core + addedbeta[i] + 1 
                writeout.append('DET ' + str(detno) + ' 2 PR ' + str(orbout) + ' 1 ' + str(orbin) + ' 1\n')  
	for i in range(len(groundOrbUp):
		for j in range(len(exOrbUp)):
			if ((groundOrbUp(i) == exOrbUp(j)) and (i != j)):
				orb = ExOrbUp.pop(j)
				ExOrbUp.insert(i,orb)
				signlist[len(signlist) - 1] *= -1
		

def CAS09(gout):
	line = ' '
	while (line != 'BOTTOM WEIGHT'):
		newline = gout.readline()
		try:
			line = newline.split('=')[0].strip()
		except:
			pass
	newline = gout.readline()
	line = newline.split()[0].strip()
	configs = []
	while (line == 'Configuration'):
		config = newline.split()[4].strip()
		configs.append(config)		
		newline = gout.readline()
		line = newline.split()[0].strip()
	while (line != "NO. OF ORBITALS"):
		newline = gout.readline()
		line = newline.split("=")[0].strip()
	numorbitals = newline.split("=")[1].strip()
	orbitals = numpy.int(numorbitals)
	newline = gout.readline()
	numelectrons = newline.split('=')[1].strip()
	electrons = numpy.int(numelectrons)
	while (line != "Enter MCSCF program."):
		newline = gout.readline()
		line = newline.strip()
	line = gout.readline()
	numcore = line.split('=')[2].strip()
	core = numpy.int(numcore)
	return (orbitals, electrons, core, configs)	
	
def writeCAS09(gwfn,output,cutoff,maxdets):
	detsout = []
	energies = []
	detcount = 0
        (orbitals, electrons, core, configs) = CAS09(output)
	line = ''
	while (line != 'EIGENVALUE'):
		newline = output.readline()
		try:
			line = newline.split(')')[1].split()[0].strip()
		except:	
			pass
	newline = output.readline()
	eiglist = []
	indexlist = []
	while (line != 'Final'):
		line = newline.replace('(',')').split(')')
		line.pop(0)				
		while (len(line) > 1):
			index = line.pop(0)
			index = numpy.int(index)
			indexlist.append(index)
			eigenvalue = line.pop(0)
			eiglist.append(eigenvalue.strip())
		newline = output.readline()
		line = newline.split()[0].strip()
	for i in range(len(indexlist)):
		if (indexlist[i] == 1):
			groundindex = i
			break
		else:
			print i, indexlist[i]
	sumForCutoff = 0.
	signList = []
	for i in range(len(indexlist)):
		if ((numpy.float(eiglist[i]) != 0.) and (sumForCutoff <= cutoff) and (detcount < maxdets)):
			detcount += 1
			signList.append(1)
			parseCAS09(configs[indexlist[i] - 1], core, configs[groundindex], detcount, detsout,signList)
			energies.append(eiglist[i])	
			sumForCutoff += numpy.float(eiglist[i]) ** 2.
	output.close()
	gwfn.write('MD\n')
	gwfn.write(str(detcount) + '\n')
	f = open(vars(args)['corr'],'a')
	f.write('START MDET\n')
	f.write('Title\n')
	f.write('MDET Example\n')
	f.write('Multideterminant/excitation specification (see manual)\n')
	f.write('MD\n')
	f.write(str(detcount) +'\n')
	counter = 1
	for i in range(len(energies)): #these aren't actually energies. they're coefficients. i am dumb.
		item = energies[i] * signList[i]
		gwfn.write(item + '\n')
		if (counter == 1):
			f.write(item + ' 1 0\n')
			counter += 1
		else:
			f.write(item + ' ' + str(counter) + ' 1 \n')
			counter += 1
	for item in detsout:
		gwfn.write(item)
		f.write(item)
	f.write('END MDET \n')
	f.close()

def CISparse(gwfn, output,cutoff,maxdets):
	detone = []
	dettwo = []
	coeffs = []
	state = str(vars(args)['state'])
	line = ' '
	targetline = ' Excited State   ' + state
	while (line != targetline):
		newline = output.readline()
		line = newline.split(':')[0]
	while ((line != 'SavETr:') and (line != 'This')):
		newline = output.readline()
		if (line == ''):
			newline = 'This'
		line = newline.split()[0].strip()
		if ((line != 'SavETr:') and (line != 'This')):
			line = newline.split('->')
			detone.append(line[0].strip())
			dettwo.append(line[1].split()[0].strip())
			coeffs.append(line[1].split()[1].strip())
	abslist = []
	for entry in coeffs:
		abslist.append(numpy.abs(numpy.float(entry)))
	a = numpy.vstack((detone, dettwo,coeffs,abslist)).T
	totalExcitations = sorted(a, key=lambda a_entry: a_entry[3], reverse=True)
	writeDet = []
	writeCoeff = []
	detcount = 0
	sumCoeff = 0
	while ((sumCoeff <= cutoff) and (detcount < maxdets) and (len(totalExcitations) > 0)):
		entry = totalExcitations.pop(0)
		sumCoeff += numpy.float(entry[3]) ** 2.
		detcount += 1			
		writeDet.append('DET ' + str(detcount) + ' 1 PR ' + str(entry[0]) +' 1 ' + str(entry[1] + ' 1\n'))
		writeCoeff.append(entry[2])
	print "Squared sum of coefficients is " + str(sumCoeff) + '\n'
	gwfn.write('MD\n')
	gwfn.write(str(detcount) + '\n')
	f = open(vars(args)['corr'],'a')
	f.write('START MDET\n')
	f.write('Title\n')
	f.write('MDET Example\n')
	f.write('Multideterminant/excitation specification (see manual)\n')
	f.write('MD\n')
	f.write(str(detcount) +'\n')
	for i in range(len(writeCoeff)):
		gwfn.write(writeCoeff[i] + ' \n')
		if (i == 0):
			f.write(writeCoeff[i] + ' ' + str(i + 1) + ' 0\n')
		else:
			f.write(writeCoeff[i] + ' ' + str(i + 1) + ' 1\n')
	for item in writeDet:
		gwfn.write(item)
		f.write(item)
	f.write('END MDET \n')
	f.close()
	output.close()
## Read the Gaussian output to get the basic information portion##

line = ' '
while (line != 'This is part of the Gaussian'):
	newline = output.readline()
	line = newline.split(r'(R)')[0].lstrip()
code = newline.split(r'(R)')[1].split()[0].strip()

# Read the fchk file for the key information ##

title = fchk.readline().strip()
method = fchk.readline().split()[1].strip()
numberatoms = fchk.readline().split('I')[1].strip()

while (line != "Number of electrons"):
	newline = fchk.readline()
	line = newline.split('I')[0].strip()
numberelectrons = newline.split('I')[1].strip()

while (line != "Number of basis functions"):
	newline = fchk.readline()
	line = newline.split('I')[0].strip()
numbasisfunctions = newline.split('I')[1].strip()

atomnumbers = []
while (line != "Atomic numbers"):
	newline = fchk.readline()
	line = newline.split('I')[0].strip()

while (line != "Nuclear charges"):
	newline = fchk.readline()
	line = newline.split('R')[0].strip()
	if (line != "Nuclear charges"):
		atoms = newline.split()
		for item in atoms:
			atomnumbers.append(item)

valencecharges = []
while (line != "Current cartesian coordinates"):
	newline = fchk.readline()
	line = newline.split('R')[0].strip()
	if (line != "Current cartesian coordinates"):
		charges = newline.split()
		for charge in charges:
			valencecharges.append(charge.strip())
	
coords = []
while (line != "Force Field"):
	newline = fchk.readline()
	line = newline.split('I')[0].strip()
	if (line != "Force Field"):
		coordline = newline.split()
		for item in coordline:
			coords.append(item.strip())
for i in range(len(coords)):
	coords[i] = '%.13E' % numpy.float(coords[i])

coordsinRows = [[coords[i],coords[i + 1], coords[i + 2]] for i in range(0,len(coords),3)]
coordsinRowsfloat = [[numpy.float(coords[i]),numpy.float(coords[i+1]),numpy.float(coords[i+2])] for i in range(0,len(coords),3)]

while (line != "Number of contracted shells"):
	newline = fchk.readline()
	line = newline.split('I')[0].strip()
contractshells = newline.split('I')[1].strip()

while (line != 'Number of primitive shells'):
	newline = fchk.readline()
	line = newline.split('I')[0].strip()
primshells = newline.split('I')[1].strip()

while (line != 'Highest angular momentum'):
	newline = fchk.readline()
	line = newline.split('I')[0].strip()
highestangular = newline.split('I')[1].strip()
highestangular = str(numpy.int(highestangular) + 1)

angulardict = {'0': 1, '1': 3, '-1': 2, '2': -4, '-2':4, '3':-5, '-3':5, '4': -6, '-4':6, '5': -7, '-5':7, '6': -8, '-6': 8}

while (line != 'Shell types'):
	newline = fchk.readline()
	line = newline.split('I')[0].strip()

shelltypes = []
haveSP = 0
while (line != 'Number of primitives per shell'):
	newline = fchk.readline()
	line = newline.split('I')[0].strip()
	if (line != 'Number of primitives per shell'):
		shells = newline.split()
		for shell in shells:
			shelltypes.append(angulardict[shell.strip()])
			if (angulardict[shell.strip()] == 2):
				haveSP = 1

shellmap = [] # number of primitives per shell stored here
while (line != "Shell to atom map"):
	newline = fchk.readline()
	line = newline.split('I')[0].strip()
	if (line != "Shell to atom map"):
		contractedprims = newline.split()
		for contracted in contractedprims:
			shellmap.append(contracted.strip())

atommap = []
while (line != "Primitive exponents"):
	newline = fchk.readline()
	line = newline.split('R')[0].strip()
	if (line != "Primitive exponents"):
		sequences = newline.split()
		for sequence in sequences:
			atommap.append(sequence.strip())	

sequence = [1]
counter = 0
positionofeachshell = [coordsinRows[counter]]
for i in range(1,len(atommap)):
	if (atommap[i] != atommap[i - 1]):
		sequence.append(i + 1)
		counter += 1
	positionofeachshell.append(coordsinRows[counter])
sequence.append(len(atommap) + 1)


primitiveexpts = []
while (line != 'Contraction coefficients'):
	newline = fchk.readline()
	line = newline.split('R')[0].strip()
	if (line != 'Contraction coefficients'):
		expline = newline.split()
		for exponent in expline:
			primitiveexpts.append(exponent.strip()) 
for i in range(len(primitiveexpts)):
	primitiveexpts[i] = '%.13E' % numpy.float(primitiveexpts[i])

contractioncoeffs = []
if (haveSP == 0):
	stoppingPoint = 'Coordinates of each shell'
else:
	stoppingPoint = 'P(S=P) Contraction coefficients'

while (line != stoppingPoint):
	newline = fchk.readline()
 	line = newline.split('R')[0].strip()
	if (line != stoppingPoint):
		coeffline = newline.split()
		for coeff in coeffline:
			contractioncoeffs.append(coeff.strip())
for i in range(len(contractioncoeffs)):
	contractioncoeffs[i] = '%.13E' % numpy.float(contractioncoeffs[i])


spcontractions = []
if (haveSP == 1):
	while (line != 'Coordinates of each shell'):
		newline = fchk.readline()
		line = newline.split('R')[0].strip()
		if (line != 'Coordinates of each shell'):
			spconline = newline.split()
			for spcons in spconline:
					spcontractions.append(spcons.strip())
	for i in range(len(spcontractions)):
		spcontractions[i] = '%.13E' % numpy.float(spcontractions[i])

shellcoords = []
while (line != "Constraint Structure"):
	newline = fchk.readline()
	line = newline.split('R')[0].strip()
	if (line != "Constraint Structure"):
		shellcoordline = newline.split()
		for shellcoord in shellcoordline:
			shellcoords.append(shellcoord.strip())

while (line != "Alpha MO coefficients"):
	newline = fchk.readline()
	line = newline.split('R')[0].strip()

alphaMOs = []
while ((line != "Total") and (line != "Beta")):
	newline = fchk.readline()
	line = newline.split('R')[0].strip().split()[0].strip()
	if ((line != "Total") and (line != "Beta")):
		MOline = newline.split()
		for MO in MOline:
			alphaMOs.append(MO.strip())
for i in range(len(alphaMOs)):
	alphaMOs[i] = '%.13E' % numpy.float(alphaMOs[i])

betaMOs = []
unrestricted = 0
if (line == "Beta MO coefficients"):
	unrestricted = 1
	while (line != "Total SCF Density"):
		newline = fchk.readline()
		line = newline.split('R')[0].strip()
		if (line != "Total SCF Density"):
			MOline = newline.split()
			for MO in MOline:
				betaMOs.append(MO.strip())
for i in range(len(betaMOs)):
	betaMOs[i] = '%.13E' % numpy.float(betaMOs[i])
fchk.close()

##############################################################

# process angmom
listedAngMom = []
for i in range(len(shellmap)):
	total = numpy.int(shellmap[i])
	for j in range(total):
		listedAngMom.append(numpy.int(shelltypes[i]))

# process d- and f- coefficients
MOcounter = 0
while (MOcounter < len(alphaMOs)):
	for shell in shelltypes:
		if (shell == 1):
			MOcounter += 1
		elif (shell == 2):
			MOcounter += 4
		elif (shell == 3):
			MOcounter += 3
                elif ((shell == 4) and (unrestricted == 0)):
                        alphaMOs[MOcounter] = '%.13E' % (numpy.float(alphaMOs[MOcounter]) / numpy.sqrt(3.))
                        alphaMOs[MOcounter + 1] = '%.13E' % (numpy.float(alphaMOs[MOcounter + 1]) * 2.)
                        alphaMOs[MOcounter + 2] = '%.13E' % (numpy.float(alphaMOs[MOcounter + 2]) * 2.)
                        alphaMOs[MOcounter + 4] = '%.13E' % (numpy.float(alphaMOs[MOcounter + 4]) * 2.)
                        MOcounter += 5
                elif ((shell == 4) and (unrestricted == 1)):
                        alphaMOs[MOcounter] = '%.13E' % (numpy.float(alphaMOs[MOcounter]) / numpy.sqrt(3.))
                        betaMOs[MOcounter] = '%.13E' % (numpy.float(betaMOs[MOcounter]) / numpy.sqrt(3.))
                        alphaMOs[MOcounter + 1] = '%.13E' % (numpy.float(alphaMOs[MOcounter + 1]) * 2.)
                        betaMOs[MOcounter + 1] = '%.13E' % (numpy.float(betaMOs[MOcounter + 1]) * 2.)
                        alphaMOs[MOcounter + 2] = '%.13E' % (numpy.float(alphaMOs[MOcounter + 2]) * 2.)
                        betaMOs[MOcounter + 2] = '%.13E' % (numpy.float(betaMOs[MOcounter + 2]) * 2.)
                        alphaMOs[MOcounter + 4] = '%.13E' % (numpy.float(alphaMOs[MOcounter + 4]) * 2.)
                        betaMOs[MOcounter + 4] = '%.13E' % (numpy.float(betaMOs[MOcounter + 4]) * 2.)
                        MOcounter += 5
                elif ((shell == 5) and (unrestricted == 0)):
			alphaMOs[MOcounter] = '%.13E' % (numpy.float(alphaMOs[MOcounter]) / numpy.sqrt(30) * 4.)
			alphaMOs[MOcounter + 1] = '%.13E' % (numpy.float(alphaMOs[MOcounter+ 1]) / numpy.sqrt(5) * 2. /3.)
			alphaMOs[MOcounter + 2] = '%.13E' % (numpy.float(alphaMOs[MOcounter+ 2]) / numpy.sqrt(5) * 2. /3.)
			alphaMOs[MOcounter + 3] = '%.13E' % (numpy.float(alphaMOs[MOcounter+ 3]) / numpy.sqrt(2) * 2. /15.)
			alphaMOs[MOcounter + 4] = '%.13E' % (numpy.float(alphaMOs[MOcounter+ 4]) / numpy.sqrt(2) * 2. /15.)
			alphaMOs[MOcounter + 5] = '%.13E' % (numpy.float(alphaMOs[MOcounter+ 5]) / numpy.sqrt(3) / 15.)
			alphaMOs[MOcounter + 6] = '%.13E' % (numpy.float(alphaMOs[MOcounter+ 6]) / numpy.sqrt(3) / 15.)
                        MOcounter += 7
                elif ((shell == 5) and (unrestricted == 1)):
			alphaMOs[MOcounter] = '%.13E' % (numpy.float(alphaMOs[MOcounter]) /  numpy.sqrt(30) * 4.)
			betaMOs[MOcounter] = '%.13E' % (numpy.float(betaMOs[MOcounter]) /  numpy.sqrt(30) * 4.)
			alphaMOs[MOcounter + 1] = '%.13E' % (numpy.float(alphaMOs[MOcounter+ 1]) / numpy.sqrt(5) * 2. /3.)
			betaMOs[MOcounter + 1] = '%.13E' % (numpy.float(betaMOs[MOcounter+ 1]) / numpy.sqrt(5) * 2. /3.)
			alphaMOs[MOcounter + 2] = '%.13E' % (numpy.float(alphaMOs[MOcounter+ 2]) / numpy.sqrt(5) * 2. /3.)
			betaMOs[MOcounter + 2] = '%.13E' % (numpy.float(betaMOs[MOcounter+ 2]) / numpy.sqrt(5) * 2. /3.)
			alphaMOs[MOcounter + 3] = '%.13E' % (numpy.float(alphaMOs[MOcounter+ 3]) / numpy.sqrt(2) * 2. /15.)
			betaMOs[MOcounter + 3] = '%.13E' % (numpy.float(betaMOs[MOcounter+ 3]) / numpy.sqrt(2) * 2. /15.)
			alphaMOs[MOcounter + 4] = '%.13E' % (numpy.float(alphaMOs[MOcounter+ 4]) / numpy.sqrt(2) * 2. /15.)
			betaMOs[MOcounter + 4] = '%.13E' % (numpy.float(betaMOs[MOcounter+ 4]) / numpy.sqrt(2) * 2. /15.)
			alphaMOs[MOcounter + 5] = '%.13E' % (numpy.float(alphaMOs[MOcounter+ 5]) / numpy.sqrt(3) / 15.)
			betaMOs[MOcounter + 5] = '%.13E' % (numpy.float(betaMOs[MOcounter+ 5]) / numpy.sqrt(3) / 15.)
			alphaMOs[MOcounter + 6] = '%.13E' % (numpy.float(alphaMOs[MOcounter+ 6]) / numpy.sqrt(3) / 15.)
			betaMOs[MOcounter + 6] = '%.13E' % (numpy.float(betaMOs[MOcounter+ 6]) / numpy.sqrt(3) / 15.)
                        MOcounter += 7


#########################################################################

## Actually write the goddamn gwfn.data file based on our reading ##

gwfn = open('gwfn.data','w')

gwfn.write(title + '\n')
gwfn.write('\n')
gwfn.write('BASIC_INFO\n')
gwfn.write('----------\n')
gwfn.write('Generated by:\n')
gwfn.write('Gaussian ' + code + '\n')
gwfn.write('Method:\n')
gwfn.write(method + '\n')
gwfn.write('DFT Functional:\n')
gwfn.write('none\n')
gwfn.write('Periodicity:\n')
gwfn.write('0\n')
gwfn.write('Spin unrestricted:\n')
if (unrestricted == 0):
	gwfn.write('.false.\n')
else:
	gwfn.write('.true.\n') 

nnrep = nnrepulsion(coordsinRowsfloat, valencecharges)
gwfn.write('nuclear-nuclear repulsion energy (au/atom):\n')
gwfn.write('%.14f \n' % nnrep)
gwfn.write('Number of electrons per primitive cell:\n')
gwfn.write(numberelectrons +'\n')
gwfn.write('\n')
gwfn.write('GEOMETRY\n')
gwfn.write('--------\n')
gwfn.write('Number of atoms:\n')
gwfn.write(str(len(atomnumbers)) + '\n')
gwfn.write('Atomic positions (au):\n')
for coord in coordsinRows:
	gwfn.write(str(coord[0]).rjust(20) + str(coord[1]).rjust(20) + str(coord[2]).rjust(20) + ' \n')
gwfn.write('Atomic numbers for each atom:\n')
counter = 0
while (len(atomnumbers) > 0):
	string = atomnumbers.pop(0)
	gwfn.write(string.rjust(10))
	counter += 1
	if (counter == 8):
		gwfn.write('\n')
		counter = 0
if (counter != 0):
	gwfn.write('\n')
gwfn.write('Valence charges for each atom:\n')
counter = 0
while (len(valencecharges) > 0):
        string = valencecharges.pop(0)
        gwfn.write(string.rjust(20))
        counter += 1
        if (counter == 4):    
                gwfn.write('\n')
                counter = 0
if (counter != 0):
	gwfn.write('\n')
gwfn.write('\n')
gwfn.write('BASIS SET\n')
gwfn.write('---------\n')
gwfn.write('Number of Gaussian centres\n')
gwfn.write(numberatoms.rjust(10) + '\n')
gwfn.write('Number of shells per primitive cell\n')
gwfn.write(contractshells.rjust(10) +'\n')
gwfn.write("Number of basis functions ('AO') per primitive cell\n")
gwfn.write(numbasisfunctions.rjust(10) + '\n')
gwfn.write('Number of Gaussian primitives per primitive cell\n')
gwfn.write(primshells.rjust(10) +'\n')
gwfn.write('Highest shell angular momentum (s/p/d/f... 1/2/3/4...)\n')
gwfn.write(highestangular.rjust(10) + '\n')
gwfn.write('Code for shell types (s/sp/p/d/f... 1/2/3/4/5...)\n')
counter = 0
while (len(shelltypes) > 0):
	string = str(shelltypes.pop(0))
	gwfn.write(string.rjust(10))
	counter += 1
	if (counter == 8):
		gwfn.write('\n')
		counter = 0
if (counter != 0):
	gwfn.write('\n') 
gwfn.write('Number of primitive Gaussians in each shell\n')
counter = 0
while (len(shellmap) > 0):
	string = shellmap.pop(0)
	gwfn.write(string.rjust(10))
	counter += 1
	if (counter == 8):
		gwfn.write('\n')
		counter = 0
if (counter != 0):
	gwfn.write('\n')
gwfn.write('Sequence number of first shell on each centre\n')
counter = 0
while (len(sequence) > 0):
	string = str(sequence.pop(0))
	gwfn.write(string.rjust(10))
	counter += 1
	if (counter == 8):
		counter = 0
		gwfn.write('\n')
if (counter != 0):
	gwfn.write('\n')
gwfn.write('Exponents of Gaussian primitives\n')
counter = 0
normfix = deepcopy(primitiveexpts)
normfixsecond = deepcopy(primitiveexpts)
while (len(primitiveexpts) > 0):
	string = primitiveexpts.pop(0)
	gwfn.write(string.rjust(20))
	counter += 1
	if (counter == 4):
		counter = 0
		gwfn.write('\n')
if (counter != 0):
	gwfn.write('\n')
gwfn.write('Normalized contraction coefficients\n')
counter	= 0
while (len(contractioncoeffs) > 0):
        coeff = contractioncoeffs.pop(0)
	expt = normfix.pop(0)
	angmom = listedAngMom.pop(0)
	normed = normCoeff(coeff, expt, angmom)
        string = '%.13E' % normed
	gwfn.write(string.rjust(20))
        counter	+= 1
        if (counter == 4):
                counter	= 0
                gwfn.write('\n')
if (counter != 0):
	gwfn.write('\n')
if (haveSP == 1):
	gwfn.write('2nd contraction coefficients (p coeff. for sp shells, 0 otherwise)\n')
	counter = 0
	while (len(spcontractions) > 0):
		coeff = spcontractions.pop(0)
		expt = normfixsecond.pop(0)
		normed = normCoeff(coeff, expt, 3)
		string = '%.13E' % normed
		gwfn.write(string.rjust(20))
		counter += 1
		if (counter == 4):
			counter = 0
			gwfn.write('\n')
	if (counter != 0):
		gwfn.write('\n')
gwfn.write('Position of each shell (au)\n')
for row in positionofeachshell:
	gwfn.write(row[0].rjust(20) + row[1].rjust(20) + row[2].rjust(20) + '\n')
gwfn.write(' \n')
gwfn.write('MULTIDETERMINANT INFORMATION\n')
gwfn.write('----------------------------\n')
cutoff = numpy.float(vars(args)['cut'])
maxdets = vars(args)['maxdets']
if ((method == 'CASSCF') and (code == '09')):
	writeCAS09(gwfn,output,cutoff,maxdets)
elif (method == 'RCIS-FC'):
	CISparse(gwfn, output,cutoff,maxdets)
	output.close()
else:
	gwfn.write('GS\n')
	output.close()
gwfn.write(' \n')
gwfn.write('ORBITAL COEFFICIENTS\n')
gwfn.write('------------------------\n')
for i in range(len(alphaMOs)):
	string = alphaMOs[i]
	gwfn.write(string.rjust(20))
	if (i % 4 == 3):
		gwfn.write('\n')
if (unrestricted == 0):
	gwfn.write('\n\n\n\n')
else:
	for i in range(len(betaMOs)):
		string = betaMOs[i]
		gwfn.write(string.rjust(20))
		if (i % 4 == 3):
			gwfn.write('\n')
gwfn.write('\n\n\n\n\n')
gwfn.close()	
	
