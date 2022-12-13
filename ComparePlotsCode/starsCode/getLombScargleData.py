import astropy as ap 
import astropy.units as u 
import astropy.table as at 
import astropy.coordinates as coords 
from astropy.io import ascii,fits
from astropy.time import Time
from astropy.timeseries import LombScargle, BoxLeastSquares
import numpy as np 
import glob
import os

def main():
	idListFile='/scratch/jdg577/theJoker/Data/starsData/apogee_id-ticid-filtered.fits'
	idHDUL=fits.open(idListFile)
	starIDs=idHDUL[1].data
	allStarLite='/scratch/jdg577/theJoker/Data/metadataDR17.fits'
	metadata=at.Table.read(allStarLite)
	starFolder='/scratch/jdg577/theJoker/Data/starsData/'
	for i in range(len(starIDs)):
		ticid=starIDs[i]['TICID']
		apogeeid=starIDs[i]['APOGEE_ID']
		star=metadata[metadata['APOGEE_ID']==apogeeid]
		map_p=star['MAP_P']
		LC=pullLC(ticid,apogeeid,starFolder)
		getBoxLeastSquares(LC, map_p, ticid, apogeeid)
		getLombScargle(LC, map_p, ticid, apogeeid)

def pullLC(TIC,APOGEE,STARS):
	path=STARS+str(APOGEE)+'_'+str(TIC)+'/'+str(APOGEE)+'_'+str(TIC)+'-LightCurve.fits'
	LCHDUL=fits.open(path)
	LCData=LCHDUL[1].data 
	return LCData

def getLombScargle(LCData,MAPP,TIC,APOGEE):
	time=LCData['Time']
	flux=LCData['Flux']
	nTerms=[1,4,16]
	for n in nTerms:
		periodogramType=str(n)+'TermLombScargle'
		data=at.Table() 
		frequency,power=LombScargle(time,flux,nterms=n).autopower(minimum_frequency=(0.05/MAPP),maximum_frequency=(20./MAPP))
		print('Done '+str(n))
		period=1./frequency
		data['Frequency']=frequency
		data['Period']=period
		data['Power']=power 
		makeFits(TIC, APOGEE, periodogramType, data)