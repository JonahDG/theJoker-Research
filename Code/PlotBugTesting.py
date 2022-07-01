
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

def getStarData(apogeeFilePath,sampleFilePath,separationFilePath,jokerSamplesFilePath):
    ###### FUTRE: MIGHT NEED TO PULL APOGEE ID FROM GOLD TBD
    #APOGEE VISITS
    apogeeVisits=at.Table.read(apogeeFilePath)

    #Gold Sample and separation File
    goldSample=at.QTable.read(sampleFilePath)
    separationSample=at.QTable.read(separationFilePath)
    # Gold Sample and separation joining
    separationSample=separationSample[separationSample['separation']<2.*u.arcsec]
    sources=at.join(goldSample,separationSample,keys='APOGEE_ID')
    # samples joining
    # sampleFitsListFull=[]
    sampleIDListFull=[]
    samplesFitsListPath=glob.glob('/scratch/jdg577/theJoker/Data/samples/**/2M*.fits')
    for sampleFitsPath in samplesFitsListPath:
        index=sampleFitsPath.find('2M')
        sampleFits=sampleFitsPath[index:]
        sampleID=sampleFits[:-5]
        sampleIDListFull.append(sampleID)
    samplesTableFull=at.Table([sampleIDListFull,samplesFitsListPath],names=('APOGEE_ID','FITS_FILE'),dtype=('str','str'))
    sources=at.join(sources,samplesTableFull,keys='APOGEE_ID')
    # sources=sources[sources['MAP_P']<10.*u.d]
    #period filtering
    print(len(sources))

    print('Data Tables Formed')
    return apogeeVisits,sources

def findNearest(periodArray,MAP_P):
    periodArray=np.asarray(periodArray)
    idx=(np.abs(periodArray-MAP_P)).argmin()
    return periodArray[idx]

allStarLite='/scratch/jdg577/theJoker/Data/allStarLite-metadata.fits'
apogeeVisitsFile='/scratch/jdg577/theJoker/Data/allVisit-r12-l33.fits'
separationFile='/scratch/jdg577/theJoker/Data/allStarLite-r12-l33-tess_2min-max_20arcsec-xm.fits'
samplesFile='/scratch/jdg577/theJoker/Data/samples'
goldSampleFile='/scratch/jdg577/theJoker/Data/gold_sample.fits'
# apogeeVisits,sources=getStarData(apogeeVisitsFile,allStarLite,separationFile,samplesFile)
apogeeVisits,sources=getStarData(apogeeVisitsFile,allStarLite,separationFile,samplesFile)
print(sources['APOGEE_ID','FITS_FILE'])
TESTID='2M09010280+0237517'
testSource=sources[sources['APOGEE_ID']==TESTID]
sources=sources[sources['MAP_P']<10*u.d]
testSource2=sources[sources['APOGEE_ID']==TESTID]

print(testSource)
print(testSource2)
