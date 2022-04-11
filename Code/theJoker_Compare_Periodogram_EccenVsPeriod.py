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
from PyPDF2 import PdfFileMerger,PdfFileReader
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
    sources=sources[sources['MAP_P']<10*u.d]
    sources=sources[:250]
    print('Data Tables Formed')
    return apogeeVisits, sources

# getJokerRVData returns RV Data from Joker
def getJokerRVData(allVisits,sources):
    data=[]
    for row in sources:
        visits=allVisits[allVisits['APOGEE_ID']==row['APOGEE_ID']] # Filters to individual APOGEE Visits
        print('Joker Running on: ', visits['APOGEE_ID'])
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
        print('Attempting Light Curve on '+ticStr)
        try:
            lcCollection=lk.search_lightcurve(target=ticStr,mission='TESS').download_all()
            lcStitch=lcCollection.stitch().remove_nans().remove_outliers()
            time=lcStitch.time.value
            flux=(lcStitch.flux)
            fluxErr=(lcStitch.flux_err)
            time=np.ascontiguousarray(time)
            flux=np.ascontiguousarray(flux)
            timeList.append(time)
            fluxList.append(flux)
            fluxErrList.append(fluxErr)
            print('Light Curve Success: '+ticStr)
        except:
            print('Light Curve Fail: '+ticStr)
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
    powList=[]
    perList=[]
    mapPList=[]
    lcData=getLightCurveData(sources)
    for i in range(len(sources['MAP_P'])):
        period=float(sources[i]['MAP_P']/u.d)
        mapPList.append(period)
        timeViaLC=lcData[i]['Time']
        fluxViaLC=lcData[i]['Flux']
        try:
            freq,pow=LombScargle(timeViaLC,fluxViaLC).autopower(minimum_frequency=(.1/period),maximum_frequency=(10./period))
            per=1./freq
            perList.append(per)
            powList.append(pow)
            print('Lomb Scargle Success: ',sources[i]['TICID'])
        except:
            perList.append([])
            powList.append([])
            print('Lomb Scargle Fail: ',sources[i]['TICID'])
            pass
    periodogramDataTable=at.QTable([perList,powList,mapPList],names=('Period','Power','MAP_P'))
    return periodogramDataTable

# getPlots returns nothing but saves individual plots as pngs and one file as pdf
def getPlots(sources,jData,periodogram):
    indivPlotArray=[]
    for i in range(len(sources['APOGEE_ID'])):
        # Titles
        supTitle='APOGEE ID: '+str(sources[i]['APOGEE_ID']+' | TIC ID: TIC'+str(sources[i]['TICID']))
        plot1Title='Eccentricity vs Period'
        plot2Title='Lomb Scargle Periodogram'
        # Figure Creation
        fig,(ax1,ax2)=plt.subplots(nrows=2,ncols=1,figsize=(15,10),facecolor='w',sharex='all')
        fig.suptitle(supTitle)
        # plot 1
        prior=tj.JokerPrior.default(P_min=periodogram[i]['MAP_P']/10.*u.day,
        P_max=periodogram[i]['MAP_P']/.1*u.day,
        sigma_K0=300*u.km/u.s, sigma_v=100*u.km/u.s)
        joker=tj.TheJoker(prior)
        prior_sample=prior.sample(100_00)
        sample=joker.rejection_sample(jData[i],prior_sample)
        ax1.plot(sample['P'],sample['e'],'o',color='black',rasterized=True)
        ax1.set_xlabel('Period (d)')
        ax1.set_ylabel('Eccentricity')
        ax1.set_title(plot1Title)
        ax1.set_xscale('log')
        # ax1.axvline(x=periodList[i]/0.1,color='red')
        # ax1.axvline(x=periodList[i]/10,color='red')
        # Plot 2
        freq=periodogram[i]['Period']
        pow=periodogram[i]['Power']
        ax2.plot(freq,pow,rasterized=True)
        ax2.set_xlabel('Period (d)')
        ax2.set_ylabel('Power')
        ax2.set_title(plot2Title)
        # ax2.axvline(x=periodList[i]/0.1,color='red')
        # ax2.axvline(x=periodList[i]/10,color='red')
        # save figure
        # local path
        # pngPath='/Users/jonahgoldfine/Desktop/theJoker-Research/Plots/Compare_Periodogram_PeriodVsEccentricity/PNGs/Periodogram_EccenPeriod_'+str(sources[i]['APOGEE_ID'])+'.png'
        indivPdfPath='/scratch/jdg577/theJoker/Plots/PNGs/PerGram_EccenPer_'+str(sources[i]['APOGEE_ID'])+'.pdf'
        fig.savefig(indivPdfPath)
        indivPlotArray.append(indivPdfPath)
        print(supTitle+' DONE')
    merger=PdfFileMerger()
    for subPDF in indivPlotArray:
        merger.append(PdfFileReader(subPDF,'rb'))
    merger.write('/scratch/jdg577/theJoker/Plots/PDFs/PerGram_EccenPer_Full.pdf')
apogeeFile='/scratch/jdg577/theJoker/Data/allVisit-r12-l33.fits'
tessFile='/scratch/jdg577/theJoker/Data/allStarLite-r12-l33-tess_2min-max_20arcsec-xm.fits'
binaryMetadataFile='/scratch/jdg577/theJoker/Data/allStarLite-metadata.fits'
apogeeData,sourceData=getStarData(apogeeFile,tessFile,binaryMetadataFile)
jokerRVData=getJokerRVData(apogeeData,sourceData)
lsPeriodogramData=getLsPeriodogram(sourceData)
getPlots(sourceData,jokerRVData,lsPeriodogramData)


