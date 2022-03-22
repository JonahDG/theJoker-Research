# Imports
#region astropy imports
import astropy as ap
import astropy.units as u
import astropy.table as at
import astropy.coordinates as coords
from astropy.io import ascii
from astropy.io import fits
from astropy.time import Time
from astropy.timeseries import LombScargle
#endregion
#region other imports
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import lightkurve as lk
import thejoker as tj
from fpdf import FPDF
#endregion

# getStarData returns apogeeVisits and sources
def getStarData(apogeeFileName,tessFileName,binariesFile):
    # APOGEE Vists Data
    apogeeVisits=at.Table.read(apogeeFileName)
    # Binary Catalog
    binaries=at.QTable.read(binariesFile)
    # TESS Data
    tessData=at.QTable.read(tessFileName)
    tessData=tessData[tessData['separation']<2.*u.arcsec] # Separation Filter
    # Join Tables
    sources=at.join(binaries,tessData,keys='APOGEE_ID')
    return apogeeVisits, sources

# getJokerRVData returns RV Data from Joker
def getJokerRVData(allVisits,sources):
    data=[]
    for row in sources:
        visits=allVisits[allVisits['APOGEE_ID']==row['APOGEE_ID']] # Filters to individual APOGEE Visits
        datum=tj.RVData(Time(visits['JD'],format='jd'),
            rv=visits['VHELIO']*u.km/u.s,
            rv_err=visits['VRELERR']*u.km/u.s)
        data.append(datum)
    return data

# getLightCurveData returns Light Curve Data Table from lightkurve as a Data Table
def getLightCurveData(sources):
    ticList=[]
    timeList=[]
    fluxList=[]
    fluxErrList=[]
    for tic in sources['TICID']:
        ticStr='TIC'+str(tic)
        ticList.append(ticStr)
        print(ticStr)
        try:
            lcCollection=lk.search_lightcurve(target=ticStr,mission='TESS').download_all()
            lcStitch=lcCollection.stitch().remove_nans().remove_outliers()
            time=lcStitch.time.value
            flux=1e3*(lcStitch.flux-1)
            fluxErr=1e3*(lcStitch.flux_err)
            time=np.ascontiguousarray(time)
            flux=np.ascontiguousarray(flux)
            timeList.append(time)
            fluxList.append(flux)
            fluxErrList.append(fluxErr)
            print('Success')
        except:
            print('Fail')
            print('Probably MergeConflictError and/or TableMergeError')
            timeList.append([])
            fluxList.append([])
            fluxErrList.append([])
            pass
        lcDataTable=at.QTable([ticList,timeList,fluxList,fluxErrList],names=('TICID',
                                                                            'Time',
                                                                            'Flux',
                                                                            'Flux_err'))
    return lcDataTable

# getLsPeriodogram returns Lomb Scargle Periodogram from Astropy as a Data Table
def getLsPeriodogram(sources):
    freqList=[]
    powList=[]
    lcData=getLightCurveData(sources)
    for i in range(len(sources['MAP_P'])):
        period=float(sources[i]['MAP_P']/u.d)
        timeViaLC=lcData[i]['Time']
        fluxViaLC=lcData[i]['Flux']
        freq,pow=LombScargle(timeViaLC,fluxViaLC).autopower(minimum_frequency=(.1/period),maximum_frequency=(10./period))
        freqList.append(freq)
        powList.append(pow)
    periodogramDataTable=at.QTable([freqList,powList],names=('Frequency,Power'))
    return periodogramDataTable

# getPlots returns nothing but saves individual plots as pngs and one file as pdf
def getPlots(osurces,jData,periodogram,joker,prior_sample):
    pdf=FPDF(unit='in',format=[8.5,11])

# INCOMPLETE