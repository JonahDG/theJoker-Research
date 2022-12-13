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
import multiprocessing as mp

def main():
	idListFile='/scratch/jdg577/theJoker/Data/starsData/apogee_id-ticid-filtered.fits'
	idHDUL=fits.open(idListFile)
	starIDs=idHDUL[1].data 
	allStarLite='/scratch/jdg577/theJoker/Data/metadataDR17.fits'
	metadata=at.Table.read(allStarLite)
	starFolder='/scratch/jdg577/theJoker/Data/starsData/'
	processes=[]
	pList=[]
	for i in range(len(starIDs)):
		ticid=starIDs[i]['TICID']
		apogeeid=starIDs[i]['APOGEE_ID']
		star=metadata[metadata['APOGEE_ID']==apogeeid]
		mapp=star['MAP_P']
		LC=pullLC(ticid,apogeeid,starFolder)
		low=2
		high=4
		p=mp.Process(target=getLowHighBLS,args=(LC,mapp,ticid,apogeeid,h,l,starFolder))
		p.start()
		processes.append(p)
		pList.append([apogeeid,ticid])
		print('Process Added')
	for process, star in zip(processes, pList):
		print(str(star)+'Joined')
		process.join()
def pullLC(TIC,APOGEE,STARS):
	path=STARS+str(APOGEE)+'_'+str(TIC)+'/'+str(APOGEE)+'_'+str(TIC)+'-LightCurve.fits'
	LCHDUL=fits.open(path)
	LCData=LCHDUL[1].data 
	return LCData
def getLowHighBLS(LCDATA,MAPP,TIC,APOGEE,HIGH,LOW,STARS):
	time=LCDATA['Time']*u.d 
	flux=LCDATA['FLux']
	perioType='BLS_'+str(LOW)+str(HIGH)
	lowPow=LOW/2.-1.
	highPow=HIGH/2.-1.
	minP=20**lowPow
	maxP=20**highPow
	data=at.Table()
	results=BoxLeastSquares(time,flux).autopower(duration=minP/2.,minimum_period=minP,maximum_period=maxP)
	period=results.period 
	power=results.power 
	frequency=1./period 
	data['Frequency']=frequency
	data['Period']=period 
	data['Power']=power 
	makeFits(TIC,APOGEE,perioType,data,STARS)
def makeFits(TIC,APOGEE,PERIODOGRAM_TYPE,PERIODOGRAM_DATA,STARS):
	file=STARS+str(APOGEE)+'_'+str(TIC)+'/'+str(APOGEE)+'_'+str(TIC)+'-'+PERIODOGRAM_TYPE+'.fits'
	if os.path.exists(STARS+str(APOGEE)+'_'+str(TIC)+'/')==False:
		os.mkdir(STARS+str(APOGEE)+'_'+str(TIC)+'/')
	else:
		pass
	hdr=fits.Header()
	hdr['APOGEEID']=str(APOGEE)
	hdr['TICID']=str(TIC)
	hdr['TYPE']=str(PERIODOGRAM_TYPE)
	hdr['DR']='17'
	dataHDU=fits.BinTableHDU(data=PERIODOGRAM_DATA,header=hdr)
	dataHDU.writeto(file,overwrite=True)

main()

