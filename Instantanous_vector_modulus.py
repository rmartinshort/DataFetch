#!/usr/bin/env python
#Calculate the instantanous vector modulus (p(t)) of a seismogram 
#from the E, N and Z components. This follows Lomax & Michelini (1987)

import obspy as op
import numpy as np 
import os
import glob


def to_IVM(Ecomp,Ncomp,Zcomp):

	'''Convert to IVM'''

	print 'Working on %s %s %s' %(Ecomp,Ncomp,Zcomp)

	e = op.read(Ecomp,format='mseed')
	n = op.read(Ncomp,format='mseed')
	z = op.read(Zcomp,format='mseed')

	p = z.copy()
	theta = z.copy()
	phi = z.copy()

	#Demean and detrend the data (for production code this should be done
	#before calculation of p)

	e.detrend(type='demean')
	e.detrend(type='linear')

	n.detrend(type='demean')
	n.detrend(type='linear')

	z.detrend(type='demean')
	z.detrend(type='linear')

	#Write the vector modulus and other components to file

	p[0].data = np.sqrt(e[0].data**2 + n[0].data**2 + z[0].data**2)
	theta[0].data = np.arcsin(z[0].data/p[0].data)
	phi[0].data = np.arctan(n[0].data/e[0].data)

	pname = Zcomp.replace("BHZ","BHP")
	thetaname = Zcomp.replace("BHZ","BHTheta")
	phiname = Zcomp.replace("BHZ","BHPhi")

	p.write(pname,format='mseed')
	theta.write(thetaname,format='mseed')
	phi.write(phiname,format='mseed')



def ConvertLoop(indir):

	pwd = os.getcwd()

	os.chdir(indir)

	#Group all the E, N and Z file
	#Make a list of all the events, by their start times

	Zfiles = glob.glob('*BHZ*.mseed')
	ev_comps = []

	for mseedname in Zfiles:
		starttime = mseedname.split('_')[2]
		station = mseedname.split('..')[0]
		compchar = '%s*%s*.mseed' %(station,starttime)
		if compchar not in ev_comps:
			ev_comps.append(compchar)

	for compchar in ev_comps:

		components = glob.glob(compchar)

		if len(components) == 3:
			Ecomp = components[0]
			Ncomp = components[1]
			Zcomp = components[2]

			to_IVM(Ecomp,Ncomp,Zcomp)



		

	os.chdir(pwd)







if __name__ == '__main__':

	ConvertLoop('waveforms')