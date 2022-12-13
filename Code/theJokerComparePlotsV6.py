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
def main():
    allStarLite='/scratch/jdg577/theJoker/Data/metadataDR17.fits'
    apogeeVisitsFile='/scratch/jdg577/theJoker/Data/allVisit-r12-l33.fits'
    separationFile='/scratch/jdg577/theJoker/Data/allStarLite-r12-l33-tess_2min-max_20arcsec-xm.fits'
    samplesFile='/scratch/jdg577/theJoker/Data/samples'

    apogeeVisits,sources=getStarData(apogeeVisitsFile,allStarLite,samplesFile,separationFile)
    combData=combineData(sources)
    print(combData['APOGEE_ID','TICID'])
    getPlots(combData)

def getStarData(apogeeVisitsFile,starFile,jokerSampleFile,separationFile):
    apogeeVisits=at.Table.read(apogeeVisitsFile)
    stars=at.QTable.read(starFile)
    separations=at.QTable.read(separationFile,format='fits')
    genericPath=jokerSampleFile+'/**/2M*.fits'
    sampleFitsPaths=glob.glob(genericPath)
    sampleIDList=[]
    for sampleFitsPath in sampleFitsPaths:
        idx=sampleFitsPath.find('2M')
        sampleID=sampleFitsPath[idx:-5]
        sampleIDList.append(sampleID)
    samplesTable=at.Table([sampleIDList,sampleFitsPaths],names=('APOGEE_ID','FITS_FILE'),dtype=('str','str'))
    sources=at.join(stars,separations,keys='APOGEE_ID')
    sources=at.join(sources,samplesTable,keys='APOGEE_ID')
    sources=sources[sources['MAP_P']<10*u.d]
    sources=sources[sources['separation']<2*u.arcsec]
    #sources=sources[:10]
    return apogeeVisits,sources

def combineData(sources):
    ticList=[]
    timeList=[]
    fluxList=[]
    fluxErrList=[]
    powList=[]
    perList=[]
    printList=[]
    fullTable=at.Table()
    for i in range(len(sources['TICID'])):
        tic=sources[i]['TICID']
        ticStr='TIC'+str(tic)
        #fullTableTest[i]['TICID']=tic
        try:
            # Light Curve
            lcCollection=lk.search_lightcurve(target=ticStr,mission='TESS').download_all()
            lcStitch=lcCollection.stitch().remove_nans().remove_outliers()
            time=np.ascontiguousarray(lcStitch.time.value)
            flux=np.ascontiguousarray(lcStitch.flux)
            fluxErr=lcStitch.flux_err

            # Periodogram
            mapP=float(sources[i]['MAP_P']/u.d)
            freq,power=LombScargle(time,flux).autopower(minimum_frequency=(.1/mapP),maximum_frequency=(10./mapP))
            period=1./freq
            
            ticList.append(tic)
            timeList.append(time)
            fluxList.append(flux)
            fluxErrList.append(fluxErr)
            perList.append(period)
            powList.append(power)
            # SUCCESS
            print('Light Curve and Periodogram Successful on ' +ticStr)
        except:
            print('Light Curve and/or Periodogram Failed on '+ticStr)
            pass 
    fullTable['TICID']=ticList
    fullTable['Time']=timeList
    fullTable['Flux']=fluxList
    fullTable['Flux_Err']=fluxErrList
    fullTable['Period']=perList
    fullTable['Power']=powList
    #fullTable['MAP_P']=mapPList

    fullTable=at.join(fullTable,sources,keys='TICID')
    fullTable=at.unique(fullTable,keys='APOGEE_ID')
    return fullTable
def getPlots(allData):
    indivPlotArray=[]
    for i in range(len(allData['APOGEE_ID'])):
        # MAPP
        map_p=allData[i]['MAP_P']/u.d
        half_map_p=0.5*map_p
        twice_map_p=2.*map_p
        tenth_map_p=.25*map_p
        tence_map_p=4.*map_p
        #Titles
        supTitle='APOGEE ID: ' +str(allData[i]['APOGEE_ID'])+' | TIC ID: ' +str(allData[i]['TICID'])
        plot1Title='Eccentricity vs Log of Period'
        plot2Title='Lomb Scargle Periodogram'
        #FIGURE
        fig,(ax1,ax2)=plt.subplots(nrows=2,ncols=1,figsize=(15,10),facecolor='w',sharex='all')
        fig.suptitle(supTitle)
        # plot 1
        sampleFits=allData[i]['FITS_FILE']
        sample=at.Table.read(sampleFits)
        ax1.plot(sample['P'],sample['e'],'o',color='Black',rasterized=True)
        ax1.set_xlabel('Period - Log Scale (d)')
        ax1.set_ylabel('Eccentricity')
        ax1.set_title(plot1Title)
        ax1.set_xscale('log')
        ax1.set_ylim(0,1)
        ax1.axvline(x=map_p,color='red',alpha=.75,label='MAP_P')
        ax1.axvline(x=half_map_p,color='purple',linestyle='-.',alpha=0.75,label=r'$\frac{1}{2}$ MAP_P')
        ax1.axvline(x=twice_map_p,color='darkolivegreen',linestyle='-.',alpha=0.75,label='2 MAP_P') 
        ax1.axvline(x=tenth_map_p,color='purple',linestyle=':',alpha=0.75,label=r'$\frac{1}{10}$ MAP_P')
        ax1.axvline(x=tence_map_p,color='darkolivegreen',linestyle=':',alpha=0.75,label='10 MAP_P')
        ax1.legend()   
        #plot2
        ax2.plot(allData[i]['Period'],allData[i]['Power'],rasterized=True)
        ax2.set_xlabel('Period - Log Scale (d)')
        ax2.set_ylabel('Power')
        ax2.axvline(x=allData[i]['MAP_P']/u.d,color='red',alpha=.75)
        ax2.axvline(x=half_map_p,color='purple',linestyle='-.',alpha=0.75,label=r'$\frac{1}{2}$ MAP_P')
        ax2.axvline(x=twice_map_p,color='darkolivegreen',linestyle='-.',alpha=0.75,label='2 MAP_P')
        ax2.axvline(x=tenth_map_p,color='purple',linestyle=':',alpha=0.75,label=r'$\frac{1}{10}$ MAP_P')
        ax2.axvline(x=tence_map_p,color='darkolivegreen',linestyle=':',alpha=0.75,label='10 MAP_P')
        ax2.legend()
        # Save Files
        indivPDF='/scratch/jdg577/theJoker/Plots/indivPdfs/ComparePlotsV6_'+str(allData[i]['APOGEE_ID'])+'.pdf'
        print(supTitle+' Plotted')
        fig.savefig(indivPDF)
        indivPlotArray.append(indivPDF)
        print(supTitle+' Saved')
    merger=PdfFileMerger()
    for subPDF in indivPlotArray:
        merger.append(PdfFileReader(open(subPDF,'rb'),strict=False))
    merger.write('/scratch/jdg577/theJoker/Plots/PDFs/ComparePlotsV6_test.pdf')

main()
