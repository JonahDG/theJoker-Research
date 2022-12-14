from astropy.io import fits
import astropy.table as at 
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import glob
from PyPDF2 import PdfFileMerger,PdfFileReader
from fpdf import FPDF

def main():
	metaListFile='/scratch/jdg577/theJoker/Data/starsData/apogee_id-ticid-filtered.fits'
	metaList=at.Table.read(metaListFile)
	plotPathList=[]
	for star in metaList:
		apogeeid=star['APOGEE_ID']
		ticid=star['TICID']
		mapp=star['MAP_P']
		jokerFile=f'/scratch/jdg577/theJoker/Data/starsData/%s_%s/%s_%s-Joker.fits'\
		%(apogeeid,ticid,apogeeid,ticid)
		LS1File=f'/scratch/jdg577/theJoker/Data/starsData/%s_%s/%s_%s-1TermLombScargle.fits'\
		%(apogeeid,ticid,apogeeid,ticid)
		LS4File=f'/scratch/jdg577/theJoker/Data/starsData/%s_%s/%s_%s-4TermLombScargle.fits'\
		%(apogeeid,ticid,apogeeid,ticid)
		LS16File=f'/scratch/jdg577/theJoker/Data/starsData/%s_%s/%s_%s-16TermLombScargle.fits'\
		%(apogeeid,ticid,apogeeid,ticid)
		# BLS02File=f'/scratch/jdg577/theJoker/Data/starsData/%s_%s/%s_%s-BLS_0_2.fits'\
		# %(apogeeid,ticid,apogeeid,ticid)
		BLS13File=f'/scratch/jdg577/theJoker/Data/starsData/%s_%s/%s_%s-BLS_1_3.fits'\
		%(apogeeid,ticid,apogeeid,ticid)
		BLS24File=f'/scratch/jdg577/theJoker/Data/starsData/%s_%s/%s_%s-BLS_2_4.fits'\
		%(apogeeid,ticid,apogeeid,ticid)
		# print(jokerFile)
		# print(LS1File)
		# print(LS4File)
		# print(LS16File)
		# print(BLS13File)
		# print(BLS24File)
		joker=at.Table.read(jokerFile)
		LS1=at.Table.read(LS1File)
		LS4=at.Table.read(LS4File)
		LS16=at.Table.read(LS16File)
		# BLS02=at.Table.read(BLS02FIle)
		BLS13=at.Table.read(BLS13File)
		BLS24=at.Table.read(BLS24File)
		plotPath=getPlots(joker,LS1,LS4,LS16,BLS13,BLS24,apogeeid,ticid,mapp)
		plotPathList.append(plotPath)
	savePNGsToSinglePDF(plotPathList)
	print('Done!')

def getPlotsUnfinished(jokerSamples,\
	LombScargle1Term,LombScargle4Term,LombScargle16Term,\
	BoxLeastSquares02,BoxLeastSquares13,BoxLeastSquares24,\
	APOGEE_ID,TICID,MAP_P):
	NotImplemented

def getPlots(jokerSamples,\
	LombScargle1Term,LombScargle4Term,LombScargle16Term,\
	BoxLeastSquares13,BoxLeastSquares24,\
	APOGEE_ID,TICID,MAP_P):
	fig,axes=plt.subplots(nrows=5,ncols=1,figsize=(25,35),\
		facecolor='w',constrained_layout=True,dpi=75.0,sharex=True)
	fig.suptitle(f'Comparison of Joker Samples with Periodograms\n\
		APOGEE ID: %s\n\
		TIC ID: %s'%(APOGEE_ID,TICID))
	subTitles=['Joker Samples',\
	'1-Term Lomb Scargle Periodogram',\
	'4-Term Lomb Scargle Periodogram',\
	'16-Term Lomb Scargle Periodogram',\
	'Box Least Squares Periodogram']
	for i,ax in enumerate(axes.flat):
		if i==0:
			ax.plot(jokerSamples['tjP'],\
				jokerSamples['e'],'o',color='black',rasterized=True)
			ax.axvline(MAP_P/3.,color='black',alpha=0.25)
			ax.axvline(MAP_P/2.,color='black',alpha=0.5)
			ax.axvline(2.*MAP_P/3.,color='black',alpha=0.75)
			ax.axvline(MAP_P,color='black',alpha=1,label=['Modes of MAP_P'])
			ax.axvline(3.*MAP_P/2.,color='black',alpha=0.75)
			ax.axvline(2.*MAP_P,color='black',alpha=0.5)
			ax.axvline(3.*MAP_P,color='black',alpha=0.25)
			ax.set_title(subTitles[i])
			ax.set_xlabel('Period (days)')
			ax.set_ylabel('Eccentricity')
			ax.set_xscale('log')
			ax.set_ylim(0,1)
		if i==1:
			ax.plot(LombScargle1Term['Period'],\
				LombScargle1Term['Power'],rasterized=True)
			ax.axvline(MAP_P/3.,color='black',alpha=0.25)
			ax.axvline(MAP_P/2.,color='black',alpha=0.5)
			ax.axvline(2.*MAP_P/3.,color='black',alpha=0.75)
			ax.axvline(MAP_P,color='black',alpha=1,label=['Modes of MAP_P'])
			ax.axvline(3.*MAP_P/2.,color='black',alpha=0.75)
			ax.axvline(2.*MAP_P,color='black',alpha=0.5)
			ax.axvline(3.*MAP_P,color='black',alpha=0.25)
			ax.set_title(subTitles[i])
			ax.set_xlabel('Period (days)')
			ax.set_ylabel('Power')
		if i==2:
			ax.plot(LombScargle4Term['Period'],\
				LombScargle4Term['Power'],rasterized=True)
			ax.axvline(MAP_P/3.,color='black',alpha=0.25)
			ax.axvline(MAP_P/2.,color='black',alpha=0.5)
			ax.axvline(2.*MAP_P/3.,color='black',alpha=0.75)
			ax.axvline(MAP_P,color='black',alpha=1,label=['Modes of MAP_P'])
			ax.axvline(3.*MAP_P/2.,color='black',alpha=0.75)
			ax.axvline(2.*MAP_P,color='black',alpha=0.5)
			ax.axvline(3.*MAP_P,color='black',alpha=0.25)
			ax.set_title(subTitles[i])
			ax.set_xlabel('Period (days)')
			ax.set_ylabel('Power')
		if i==3:
			ax.plot(LombScargle16Term['Period'],\
				LombScargle16Term['Power'],rasterized=True)
			ax.axvline(MAP_P/3.,color='black',alpha=0.25)
			ax.axvline(MAP_P/2.,color='black',alpha=0.5)
			ax.axvline(2.*MAP_P/3.,color='black',alpha=0.75)
			ax.axvline(MAP_P,color='black',alpha=1,label=['Modes of MAP_P'])
			ax.axvline(3.*MAP_P/2.,color='black',alpha=0.75)
			ax.axvline(2.*MAP_P,color='black',alpha=0.5)
			ax.axvline(3.*MAP_P,color='black',alpha=0.25)
			ax.set_title(subTitles[i])
			ax.set_xlabel('Period (days)')
			ax.set_ylabel('Power')
		if i==4:
			ax.plot(BoxLeastSquares13['Period'],\
				BoxLeastSquares13['Power'],\
				color='#984ea3',rasterized=True)
			ax.plot(BoxLeastSquares24['Period'],\
				BoxLeastSquares24['Power'],\
				color='#a65628',rasterized=True)
			ax.axvline(MAP_P/3.,color='black',alpha=0.25)
			ax.axvline(MAP_P/2.,color='black',alpha=0.5)
			ax.axvline(2.*MAP_P/3.,color='black',alpha=0.75)
			ax.axvline(MAP_P,color='black',alpha=1,label=['Modes of MAP_P'])
			ax.axvline(3.*MAP_P/2.,color='black',alpha=0.75)
			ax.axvline(2.*MAP_P,color='black',alpha=0.5)
			ax.axvline(3.*MAP_P,color='black',alpha=0.25)
			ax.set_title(subTitles[i])
			ax.set_xlabel('Period (days)')
			ax.set_ylabel('Power')
	pngFile=f'/scratch/jdg577/theJoker/Plots/starPlots/%s_%s/%s_%s-Compare_Periodograms_Joker.png'%(str(APOGEE_ID),str(TICID),str(APOGEE_ID),str(TICID))
	if os.path.exists(f'/scratch/jdg577/theJoker/Plots/starPlots/%s_%s/'%(str(APOGEE_ID),str(TICID)))==False:
 		os.mkdir(f'/scratch/jdg577/theJoker/Plots/starPlots/%s_%s/'%(str(APOGEE_ID),str(TICID)))
	else:
		pass
	fig.savefig(pngFile)
	print('Figure Done!')
	return pngFile

def savePNGsToSinglePDF(plots):
	pdf=FPDF(unit='in',format=[8.5,11])
	for subPNG in plots:
		pdf.add_page()
		pdf.image(subPNG,w=6.5,h=10)
		print(f'%s Saved to PDF'%subPNG)
	pdf.output('/scratch/jdg577/theJoker/Plots/starPlots/PlotTypePDFs/Compare_Periododgrams_Joker.pdf')

main()




