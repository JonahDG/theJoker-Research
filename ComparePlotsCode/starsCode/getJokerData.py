import astropy as ap 
import astropy.units as u 
import astropy.table as at 
from astropy.io import fits
import numpy as np 
import os
import glob
def main():
	samplesFiles='/scratch/jdg577/theJoker/Data/samples'
	idListFile='/scratch/jdg577/theJoker/Data/starsData/apogee_id-ticid-filtered.fits'
	idHDUL=fits.open(idListFile)
	starIDs=idHDUL[1].data
	for i in range(len(starIDs['TICID'])):
		ticid=starIDs[i]['TICID']
		apogeeid=starIDs[i]['APOGEE_ID']
		jokerSamples=getJokerSampling(samplesFiles,apogeeid)
		makeFits(apogeeid,ticid,jokerSamples)
	print('Done!')

def getJokerSampling(samplesFolder,apogeeid):
	genericPath=samplesFolder+'/**/'+str(apogeeid)+'.fits'
	perList=[]
	eList=[]
	omegaList=[]
	M0List=[]
	sList=[]
	KList=[]
	v0List=[]
	ln_priorList=[]
	ln_likelihoodList=[]
	samplesPath=glob.glob(genericPath)
	jokerSampling=at.Table()
	sampling=at.Table.read(samplesPath[0])
	perList.append(sampling['P'])
	eList.append(sampling['e'])
	omegaList.append(sampling['omega'])
	M0List.append(sampling['M0'])
	sList.append(sampling['s'])
	KList.append(sampling['K'])
	v0List.append(sampling['v0'])
	ln_priorList.append(sampling['ln_prior'])
	ln_likelihoodList.append(sampling['ln_likelihood'])

	perCol=at.Column(perList,dtype=float,unit=u.d)
	eCol=at.Column(eList,dtype=float)
	omegaCol=at.Column(omegaList,dtype=float,unit=u.rad)
	M0Col=at.Column(M0List,dtype=float,unit=u.rad)
	sCol=at.Column(sList,dtype=float,unit=u.km/u.s)
	KCol=at.Column(KList,dtype=float,unit=u.km/u.s)
	v0Col=at.Column(v0List,dtype=float,unit=u.km/u.s)
	ln_priorCol=at.Column(ln_priorList,dtype=float)
	ln_likelihoodCol=at.Column(ln_likelihoodList,dtype=float)

	jokerSampling['tjP']=perList
	jokerSampling['e']=eCol
	jokerSampling['omega']=omegaCol
	jokerSampling['M0']=M0Col
	jokerSampling['s']=sCol
	jokerSampling['tjK']=KCol
	jokerSampling['v0']=v0Col
	jokerSampling['ln_prior']=ln_priorCol
	jokerSampling['ln_likelihood']=ln_likelihoodCol

	return jokerSampling

def makeFits(APOGEE_ID,TICID,JOKERSAMPLES):
	file='/scratch/jdg577/theJoker/Data/starsData/'+str(APOGEE_ID)+'_'+str(TICID)+'/'+str(APOGEE_ID)+'_'+str(TICID)+'-Joker.fits'
	if os.path.exists('/scratch/jdg577/theJoker/Data/starsData/'+str(APOGEE_ID)+'_'+str(TICID)+'/')==False:
		os.mkdir('/scratch/jdg577/theJoker/Data/starsData/'+str(APOGEE_ID)+'_'+str(TICID)+'/')
	else:
		pass
	hdr=fits.Header()
	hdr['APOGEEID']=str(APOGEE_ID)
	hdr['TICID']=str(TICID)
	dataHDU=fits.BinTableHDU(data=JOKERSAMPLES,header=hdr)
	dataHDU.writeto(file,overwrite=True)

main()

