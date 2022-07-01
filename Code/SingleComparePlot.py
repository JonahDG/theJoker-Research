# Imports
import astropy as ap
import astropy.units as u
import astropy.table as at
import astropy.coordinates as coords
from astropy.io import ascii,fits
from astropy.time import Time
from astropy.timeseries import LombScargle
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import lightkurve as lk
import thejoker as tj
from PyPDF2 import PdfFileMerger,PdfFileReader
import glob

def getStarData(apogeeFile,starFile,separationFile,jokerSampleFile):
    apogeeVisits=at.Table.read(apogeeFile)
    stars=at.QTable.read(starFile)
    separations=at.QTable.read(separationFile)

    sampleIDList=[]
    GenPath=jokerSampleFile+'/**/2M*.fits'
    sampleFitsPaths=glob.glob(GenPath)
    for sampleFitsPath in sampleFitsPaths:
        index=sampleFitsPath.find('2M')
        sampleID=sampleFitsPath[index:-5]
        # print(sampleID)
        sampleIDList.append(sampleID)
    samplesTable=at.Table([sampleIDList,sampleFitsPaths],names=('APOGEE_ID','FITS_FILE'))

    #Filters
    sources=at.join(stars,separations,keys='APOGEE_ID')
    sources=at.join(sources,samplesTable,keys='APOGEE_ID')
    sourceLCs=sources[sources['MAP_P']<10*u.d]
    sources=sources[sources['separation']<2*u.arcsec]

    return apogeeVisits,sources

def getLCData(sources):
    ticList=[]
    timeList=[]
    fluxList=[]
    fluxErrList=[]
    for tic in sources['TICID']:
        ticStr='TIC'+str(tic)
        ticList.append(ticStr)
        print('Attempting Light Curve on ' +ticStr)
        try:
            lcCollection=lk.search_lightcurve(target=ticStr,mission='TESS').download_all()
            lcStitch=lcCollection.stitch().remove_nans().remove_outliers()
            time=lcStitch.time.value
            flux=lcStitch.flux
            fluxErr=lcStitch.flux_err
            tiqme=np.ascontiguousarray(time)
            flux=np.ascontiguousarray(flux)
            timeList.append(time)
            fluxList.append(flux)
            fluxErrList.append(fluxErr)
            print('Light Curve Successful on '+ticStr)
        except:
            print('Light Curve Failed on '+ticStr)
            timeList.append([])
            fluxList.append([])
            fluxErrList.append([])
            pass
    lcDataTable=at.QTable([ticList,timeList,fluxList,fluxErrList],
    names=('TICID','Time','Flux','Flux_err'))
    return lcDataTable

allStarFile='/scratch/jdg577/theJoker/Data/allStarLite-metadata.fits'
apogeeVisitsFile='/scratch/jdg577/theJoker/Data/allVisit-r12-l33.fits'
# goldSampleFile='/scratch/jdg577/theJoker/Data/gold_sample.fits'
separationFile='/scratch/jdg577/theJoker/Data/allStarLite-r12-l33-tess_2min-max_20arcsec-xm.fits'
samplesFile='/scratch/jdg577/theJoker/Data/samples'

apogeeData,allStarSources=getStarData(apogeeVisitsFile,allStarFile,separationFile,samplesFile)
# apogeeData,goldSources=getStarData(apogeeVisitsFile,goldSampleFile,separationFile,samplesFile)
print(len(allStarSources))
lcData=getLCData(allStarSources)
print(len(lcData))
