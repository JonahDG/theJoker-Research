import astropy as ap 
import astropy.units as u 
import astropy.table as at 
from astropy.io import fits
import lightkurve as lk
import numpy as np 
import os


def main():
	idListFile='/scratch/jdg577/theJoker/Data/starsData/apogee_id-ticid-filtered.fits'
	idHDUL=fits.open(idListFile)
	starIDs=idHDUL[1].data
	for i in range(len(starIDs['TICID'])):
		ticid=starIDs[i]['TICID']
		apogeeid=starIDs[i]['APOGEE_ID']
		lcdata=getLCData(ticid)
		makeFits(apogeeid,ticid,lcdata)
	print('Done!')

def getLCData(TICID):
	ticStr='TIC'+str(TICID)
	LCData=at.Table()
	lcCollection=lk.search_lightcurve(target=ticStr,mission='TESS').download_all()
	lcStitch=lcCollection.stitch().remove_nans().remove_outliers()
	time=np.ascontiguousarray(lcStitch.time.value)
	flux=np.ascontiguousarray(lcStitch.flux)
	fluxErr=lcStitch.flux_err
	timeCol=at.Column(time,dtype=float)
	fluxCol=at.Column(flux,dtype=float)
	fluxErrCol=at.Column(fluxErr,dtype=float)
	LCData['Time']=timeCol
	LCData['Flux']=fluxCol
	LCData['FluxErr']=fluxErrCol
	return LCData 

def makeFits(APOGEE_ID,TICID,LCData):
	file='/scratch/jdg577/theJoker/Data/starsData/'+str(APOGEE_ID)+'_'+str(TICID)+'/'+str(APOGEE_ID)+'_'+str(TICID)+'-LightCurve.fits'
	if os.path.exists('/scratch/jdg577/theJoker/Data/starsData/'+str(APOGEE_ID)+'_'+str(TICID)+'/')==False:
		os.mkdir('/scratch/jdg577/theJoker/Data/starsData/'+str(APOGEE_ID)+'_'+str(TICID)+'/')
	else:
		pass
	hdr=fits.Header()
	hdr['APOGEEID']=str(APOGEE_ID)
	hdr['TICID']=str(TICID)
	dataHDU=fits.BinTableHDU(data=LCData,header=hdr)
	dataHDU.writeto(file,overwrite=True)

main()



