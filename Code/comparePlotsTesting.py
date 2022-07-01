# Imports
import astropy as ap
import astropy.units as u
import astropy.table as at
import astropy.coordinates as coords
from astropy.io import ascii,fits
from astropy.time import Time
from astropy.timeseries import LombScargle as LS
import matplotlib.pyplot as plt
import numpy as np
import lightkurve as lk
import thejoker as tj
from PyPDF2 import PdfFileMerger,PdfFileReader
import glob

def getStarData(apogeeFile,starFile,separationFile,jokerSampleFile):
    # Read Files
    apogeeVisits=at.Table.read(apogeeFile)
    stars=at.QTable.read(starFile)
    separations=at.QTable.read(separationFile)
    # Generate Samples File Paths
    sampleIDList=[]
    genPath=jokerSampleFile+'/**/2M*.fits' # Generic path for glob
    sampleFitsPaths=glob.glob(genPath) # Generates paths for each sample 
    # Generate list of APOGEE IDs & adds to table with File path
    for sampleFitsPath in sampleFitsPaths:
        index=sampleFitsPath.find('2M') # finds where apogee id starts in file path
        sampleID=sampleFitsPath[index:-5] # pulls apogee id from file path
        sampleIDList.append(sampleID)
    sampleTable=at.Table([sampleIDList,sampleFitsPaths],names=('APOGEE_ID','FITS_FILE'))
    # Filters
    sources=at.join(stars,separations,keys='APOGEE_ID')
    sources=at.join(sources,sampleTable,keys='APOGEE_ID')
    sources=sources[sources['MAP_P']<10*u.d]
    sources=sources[sources['separation']<2*u.arcsec]
    sources=sources[:10]
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

def getPlots(sources,periodogramData):
    indivPlotArray=[]
    for i in range(len(sources['APOGEE_ID'])):
        # Titles
        supTitle='APOGEE ID: '+str(sources[i]['APOGEE_ID'])+\
        ' | TIC ID: '+ str(sources[i]['TICID'])
        plot1Title='Eccentricity vs Log(Period)'
        plot2Title='Lomb Scargle Periodogram'
        #figure
        fig,(ax1,ax2)=plt.subplots(nrows=2,ncols=1,
        figsize=(15,10),facecolor='w',sharex='all')
        fig.suptitle(supTitle)
        # Plot 1
        samplesFits=sources[i]['FITS_FILE']
        sample=at.Table.read(samplesFits)
        ax1.plot(sample['P'],sample['e'],'o',color='Black',rasterized=True)
        ax1.set_xlabel('Log(Period) (d)')
        ax1.plot(sample['P'],sample['e'],'o',color='Black',rasterized=True)
        ax1.set_xlabel('Log(Period) (d)')
        ax1.set_ylabel('Eccentricity')
        ax1.set_title(plot1Title)
        ax1.set_xscale('log')
        ax1.set_ylim(0,1)
        ax1.axvline(x=periodogramData[i]['MAP_P'],color='red',alpha=.75)
        # Plot 2
        period=periodogramData[i]['Period']
        power=periodogramData[i]['Power']
        ax2.plot(period,power,rasterized=True)
        ax2.set_xlabel('Log(Period) (d)')
        ax2.set_ylabel('Power')
        ax2.axvline(x=periodogramData[i]['MAP_P'],color='red',alpha=.75)

        #saving files
        indivPDF='/scratch/jdg577/theJoker/Plots/PNGs/ComparePlots_Testing_'+str(sources[i]['APOGEE_ID'])+'.pdf'
        fig.savefig(indivPDF)
        indivPlotArray.append(indivPDF)
        print(supTitle+' Plotted')
    merger=PdfFileMerger()
    for subPDF in indivPlotArray:
        merger.append(PdfFileReader(open(subPDF,'rb'),strict=False))
    merger.write('/scratch/jdg577/theJoker/Plots/PDFs/Compare_Plots_Testing.pdf')

allStarLite='/scratch/jdg577/theJoker/Data/allStarLite-metadata.fits'
apogeeVisitsFile='/scratch/jdg577/theJoker/Data/allVisit-r12-l33.fits'
separationFile='/scratch/jdg577/theJoker/Data/allStarLite-r12-l33-tess_2min-max_20arcsec-xm.fits'
samplesFile='/scratch/jdg577/theJoker/Data/samples'
apogeeData,sourceData=getStarData(apogeeVisitsFile,allStarLite,separationFile,samplesFile)
lsPeriodogramData=getLSPeriodogram(sourceData)
getPlots(sourceData,lsPeriodogramData)


