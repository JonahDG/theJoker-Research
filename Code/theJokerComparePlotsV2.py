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

# Star Data
def getStarData(apogeeFileName,tessFileName,tessMetaFileName):
    # APOGEE Visits
    apogeeVisits=at.Table.read(apogeeFileName)
    # Tess Data and metadata
    tessData=at.QTable.read(tessFileName)
    tessData=tessData[tessData['separation']<2.*u.arcsec] # filter out large separationts
    tessMetaData=at.QTable.read(tessMetaFileName)
    # Join Tess Files and  Filter Data
    sources=at.join(tessMetaData,tessData,key='APOGEE_ID')
    sources=sources[sources['MAP_P']<10*u.d] # Filter out long periods
    sources=sources[:250] #first 250 stars
    print('Data Tables Formed')
    return apogeeVisits, sources

# Get Joker Data
def getJokerRVData(allVisits,sources):
    data=[]
    for row in sources:
        visists=allVisit[allVisit['APOGEE_ID']==row['APOGEE_ID']]
        datum=tj.RVData(Time(visists['JDG'],format='jd'),
        rv=visists['VHELIO']*u.km/u.s,
        rv_err=visists['VRELERR']*u.km/u.s)
        data.append(datum)
    return data

# get LC time, lc flux, lc flux err
def getLCData(sources):
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
            flux=lcStitch.flux
            fluxErr=lcStitch.flux_err
            time=np.ascontiguousarray(time)
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

# Periodogram
def getLSPeriodogram(sources):
    powList=[]
    perList=[]
    mapPList=[]
    lcData=getLCData(sources)
    for i in range(len(sources['MAP_P'])):
        mapP=float(sources[i]['MAP_P']/u.d)
        mapPList.append(mapP)
        timeViaLC=lcData[i]['Time']
        fluxViaLC=lcData[i]['Flux']
        print('Attempting Periodogram on TIC'+str(sources[i]['TICID']))
        try:
            frequency,power=LombScargle(timeViaLC,fluxViaLC).autopower(minimum_frequency=(.1/mapP),
            maximum_frequency=(10./mapP))
            period=1./frequency
            perList.append(period)
            powList.append(power)
            print('Lomb Scargle Successful on TIC'+str(sources[i]['TICID']))
        except:
            perList.append([])
            powList.append([])
            print('Lomb Scargle Failed on TIC'+str(sources[i]['TICID']))
            pass
    periodogramDataTable=at.QTable([perList,powList,mapPList],
    names=('Period','Power','MAP_P'))
    return periodogramDataTable

# Create and save plots
def getPlots(sources,jData,periodogramData):
    indivPlotArray=[]
    for i in range(len(sources['APOGEE_ID'])):
        # Titles
        supTitle='APOGEE ID: '+str(sources[i]['APOGEE_ID'])+\
        ' | TIC ID: '+ str(sources[i]['TICID'])
        plot1Title='Eccentricity vs Log(Period)'
        plot2Title='Lomb Scargle Periodogram'
        #figure
        fig,(ax1,ax2)=plt.subplots(nrows=2,ncols=1,
        figsze=(15,10),facecolor='w',sharex='all')
        fig.supTitle(supTitle)
        # Plot 1
        prior=tj.JokerPrior.default(P_min=periodogramData[i]['MAP_P']/10.*u.d,
        P_max=periodogramData[i]['MAP_P']/.1*u.d,
        sigma_K0=300*u.km/u.s,sigma_v=100*u.km/u.s)
        joker=tj.TheJoker(prior)
        prior_sample=prior.sampe(100_00)
        sampe=joker.rejection_sample(jData[i],prior_sample)
        ax1.plot(sample['P'],sample['e'],'o',color='Black',rasterized=True)
        ax1.set_xlabel('Log(Period) (d)')
        ax1.set_ylabel('Eccentricity')
        ax1.set_title('plot1Title')
        ax1.set_xscale('log')
        ax1.axvline(x=periodogramData[i]['MAP_P'],color='red',alpha=.75)
        # Plot 2
        period=periodogramData[i]['Period']
        power=periodogramData[i]['Power']
        ax2.plt(period,power,rasterized=True)
        ax2.set_xlabel('Log(Period) (d)')
        ax2.set_ylabel('Power')
        ax2.axvline(x=periodogramData[i]['MAP_P'],color='red',alpha=.75)

        #saving files
        indivPDF='/scratch/jdg577/theJoker/Plots/PNGs/ComparePlots_'+str(sources[i]['APOGEE_ID'])+'.pdf'
        fig.savefig(indivPDF)
        indivPlotArray.append(indivPDF)
        print(supTitle+' Plotted')
    merger=PdfFileMerger()
    for subPDF in indivPlotArray:
        merger.append(PdfFileMerger(subPDF,'rb'))
    merger.write('/scratch/jdg577/theJoker/Plots/PDFs/ComparePlotsFull.pdf')
apogeeFile='/scratch/jdg577/theJoker/Data/allVisit-r12-l33.fits'
tessFile='/scratch/jdg577/theJoker/Data/allStarLite-r12-l33-tess_2min-max_20arcsec-xm.fits'
tessMetaDataFile='/scratch/jdg577/theJoker/Data/allStarLite-metadata.fits'
apogeeData,sourceData=getStarData(apogeeFile,tessFile,tessMetaDataFile)
jokerRVData=getJokerRVData(apogeeData,sourceData)
lsPeriodogramData=getLSPeriodogram(sourceData)
getPlots(sourceData,jokerRVData,lsPeriodogramData)