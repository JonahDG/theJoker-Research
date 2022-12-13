# Imports 
import astropy as ap
import astropy.units as u
import astropy.table as at
import astropy.coordinates as coords
from astropy.io import ascii,fits
from astropy.time import Time
from astropy.timeseries import LombScargle
import matplotlib.pyplot as plt
import numpy as np
import lightkurve as lk
import thejoker as tj
from PyPDF2 import PdfFileMerger,PdfFileReader
import glob

def getStarData(apogeeVisistsFile,starFile,jokerSampleFile,SeparationFile):
    apogeeVisists=at.Table.read(apogeeVisistsFile)
    stars=at.QTable.read(starFile)
    # print(stars.colnames)
    separations=at.QTable.read(separationFile,format='fits')
    # print(separations.colnames)
    genericPath=jokerSampleFile+'/**/2M*.fits'
    sampleFitsPaths=glob.glob(genericPath)
    sampleIDList=[]
    for sampleFitsPath in sampleFitsPaths:
        index=sampleFitsPath.find('2M')
        sampleID=sampleFitsPath[index:-5]
        sampleIDList.append(sampleID)
    samplesTable=at.Table([sampleIDList,sampleFitsPaths],names=('APOGEE_ID','FITS_FILE'),dtype=('str','str'))
    sources=at.join(stars,separations,keys='APOGEE_ID')
    sources=at.join(sources,samplesTable,keys='APOGEE_ID')
    sources=sources[sources['MAP_P']<10*u.d]
    sources=sources[sources['separation']<2*u.arcsec]
    sources=sources[2]
    return apogeeVisists,sources

def getLCData(sources):
    ticStr='TIC'+str(sources['TICID'])
    lcCollection=lk.search_lightcurve(target=ticStr,mission='TESS').download_all()
    lcStitch=lcCollection.stitch().remove_nans().remove_outliers()
    time=np.ascontiguousarray(lcStitch.time.value)
    flux=np.ascontiguousarray(lcStitch.flux)
    fluxErr=lcStitch.flux_err
    lcDataTable=at.QTable([sources['TICID'],time,flux,fluxErr],names=('TICID','Time','Flux','Flux_err'))
    return lcDataTable

def getLSPeriodogram(sources):
    lcData=getLCData(sources)
    for i in range(len(sources['MAP_P'])):
        mapP=float(sources[i]['MAP_P']/u.d)
        ticID=sources[i]['TICID']
        timeViaLC=lcData[ticID]['Time']
        fluxViaLC=lcData[ticID]['Flux']
        freq,power=LombScargle(timeViaLC,fluxViaLC).autopower(minimum_frequency=(.1/mapP),maximum_frequency=(10./mapP))
        period=1./freq
        periodogramDataTable=at.QTable([per,power,mapP],names=('Period','Power','Map_P'))
    return periodogramDataTable

def getPlots(sources):
    periodogramData=getLSPeriodogram(sources)
    supTitle='TESTING' + str(    
allStarLite='/scratch/jdg577/theJoker/Data/allStarLite-metadata.fits'
apogeeVisitsFile='/scratch/jdg577/theJoker/Data/allVisit-r12-l33.fits'
separationFile='/scratch/jdg577/theJoker/Data/allStarLite-r12-l33-tess_2min-max_20arcsec-xm.fits'
samplesFile='/scratch/jdg577/theJoker/Data/samples'
apogeeData,sourceData=getStarData(apogeeVisitsFile,allStarLite,samplesFile,separationFile)
print(sourceData['APOGEE_ID','TICID','FITS_FILE'])
