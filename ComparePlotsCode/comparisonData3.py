import astropy as ap 
import astropy.units as u 
import astropy.table as at 
import astropy.coordinates as coords 
from astropy.io import ascii,fits
from astropy.time import Time
from astropy.timeseries import LombScargle, BoxLeastSquares
import matplotlib.pyplot as plt 
import numpy as np 
import lightkurve as lk 
from PyPDF2 import PdfFileMerger,PdfFileReader
import glob

def main():
    samplesFiles='/scratch/jdg577/theJoker/Data/samples'
    separationFile='/scratch/jdg577/theJoker/Data/allStarLite-r12-l33-tess_2min-max_20arcsec-xm.fits'
    allStarLite='/scratch/jdg577/theJoker/Data/metadataDR17.fits'
    apogeeWTICID=matchApogeeDataTICID(allStarLite,separationFile)
    print('APOGEE+TIC')
    print(len(apogeeWTICID['APOGEE_ID']))
    apogeeWTICID=apogeeWTICID[apogeeWTICID['MAP_P']<10*u.d]
    apogeeWTICID=apogeeWTICID[apogeeWTICID['separation']<2*u.arcsec]
    print('APOGEE+TIC+Filter')
    print(len(apogeeWTICID['APOGEE_ID']))
    jokerPercent=getK16Percentile(samplesFiles)
    print('Joker')
    print(len(jokerPercent))
    jokerPercent=jokerPercent[jokerPercent['K16percentile']>1]
    print('Joker+Filter')
    print(len(jokerPercent))
    combine1=at.join(apogeeWTICID,jokerPercent,keys='APOGEE_ID')
    print('COmbine1')
    print(len(combine1))
    flagLC=getLCFlag(combine1)
    print('flag')
    print(len(flagLC))
    flagLC=flagLC[flagLC['Flag']==True]
    print('flag+filter')
    print(len(flagLC))
    combined=at.join(combine1,flagLC,keys='TICID')
    print(combined)
    stars=combined['APOGEE_ID','TICID']
    print(stars)

def matchApogeeDataTICID(starFile,separationFile):
	stars=at.QTable.read(starFile)
	seps=at.QTable.read(separationFile)
	apogeeData_TIC=at.join(stars,seps,keys='APOGEE_ID')

	return apogeeData_TIC['APOGEE_ID','TICID','MAP_P','separation']

def getK16Percentile(samplesFolder):
	genericPath=samplesFolder+'/**/2M*.fits' # This ignores any apogee ids that aren't 2Mass but that is an negligible amount
	'''I did a quick check of how many non-2mass stars made it through the <10 days MAP P filter and the <2 arcsec separation
	filter and it was just one star, so I assume (at least for now) that it is prettty safe to ignore'''
	samplesPaths=glob.glob(genericPath)
	sampleIDList=[]
	K16percentileList=[]
	jokerPercentile=at.Table()
	for samplesPath in samplesPaths:
		idx=samplesPath.find('2M')
		sampleID=samplesPath[idx:-5]
		sampling=at.Table.read(samplesPath)
		sampleIDList.append(sampleID)
		K16percentile=np.percentile(sampling['K'],16)
		K16percentileList.append(K16percentile)
	sampleIDCol=at.Column(sampleIDList,dtype=str)
	K16percentileCol=at.Column(K16percentileList,dtype=float)
	jokerPercentile['APOGEE_ID']=sampleIDCol
	jokerPercentile['K16percentile']=K16percentileCol
	return jokerPercentile

def getLCFlag(apogeeTic):
	flagList=[]
	ticList=[]
	LCFlag=at.Table()
	for i in range(len(apogeeTic['TICID'])):
		print(i)
		tic=apogeeTic[i]['TICID'] 
		ticStr='TIC'+str(tic)
		ticList.append(tic)
		try:
			lcCollection=lk.search_lightcurve(target=ticStr,mission='TESS').download_all()
			lcStitch=lcCollection.stitch().remove_nans().remove_outliers()
			flagList.append(True)
			print('LIGHT CURVE SUCCESSFUL ON '+ticStr+' | '+str(apogeeTic[i]['APOGEE_ID']))
		except:
			flagList.append(False)
			print('LIGHT CURVE Failed ON '+ticStr+' | '+str(apogeeTic[i]['APOGEE_ID']))
			pass
	ticCol=at.Column(ticList,dtype=int)
	flagCol=at.Column(flagList,dtype=bool)
	LCFlag['TICID']=ticCol
	LCFlag['Flag']=flagList
	return LCFlag

def combine(apogeeTic,percentile,flags):
	combinedData=at.join(apogeeTic,percentile,keys='APOGEE_ID')
	combinedData=at.join(combinedData,flags,keys='TICID')
	finalData=combinedData['APOGEE_ID','TICID','MAP_P','separation','K16percentile','Flag']
	return finalData
main()
