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
	starsFile='/scratch/jdg577/theJoker/Data/apogee_id-ticid-filtered.fits'
	stars=at.Table.read(starsFile)
	plotPathList=[]
	plotPathList02=[]
	plotPathList13=[]
	plotPathList24=[]
	plotPathListComp=[]
	'''Anything related to the low end frequencies has been removed 
	until HPC Greene Script working again'''
	print('test')
	for star in stars:
		apogeeid=star['APOGEE_ID']
		ticid=star['TICID']

		# file02=f'/scratch/jdg577/theJoker/Data/starsData/%s_%s/%s_%s-BLS_0_2.fits'\
		# %(apogeeid,ticid,apogeeid,ticid)
		# _02=at.Table.read(file02)
		# freq02=_02['Frequency']
		# pow02=_02['Power']
		# bestFreq02=getBestFreq(freq02,pow02)
		# bestPer02=1./bestFreq02

		file13=f'/scratch/jdg577/theJoker/Data/starsData/%s_%s/%s_%s-BLS_1_3.fits'\
		%(apogeeid,ticid,apogeeid,ticid)
		_13=at.Table.read(file13)
		freq13=_13['Frequency']
		pow13=_13['Power']
		bestFreq13=getBestFreq(freq13,pow13)
		bestPer13=1./bestFreq13

		file24=f'/scratch/jdg577/theJoker/Data/starsData/%s_%s/%s_%s-BLS_2_4.fits'\
		%(apogeeid,ticid,apogeeid,ticid)
		_24=at.Table.read(file24)
		freq24=_24['Frequency']
		pow24=_24['Power']
		bestFreq24=getBestFreq(freq24,pow24)
		bestPer24=1./bestFreq24

		# composite=at.vstack([_02,_13,_24])
		composite=at.vstack([_13,_24])
		freqComp=composite['Frequency']
		powComp=composite['Power']
		bestFreqComp=getBestFreq(freqComp,powComp)
		bestPerComp=1./bestFreqComp

		fileLC=f'/scratch/jdg577/theJoker/Data/starsData/%s_%s/%s_%s-LightCurve.fits'\
		%(apogeeid,ticid,apogeeid,ticid)
		LightCurve=at.Table.read(fileLC)

		# BLSFracList=['02','13','24','COMP']
		BLSFracList=['13','24','COMP']
		for frac in BLSFracList:
			if frac=='02':
				NotImplemented
			if frac=='13':
				plotPath=getPlots(bestPer13,LightCurve,star,frac)
				plotPathList13.append(plotPath)
				plotPathList.append(plotPath)
			if frac=='24':
				plotPath=getPlots(bestPer24,LightCurve,star,frac)
				plotPathList24.append(plotPath)
				plotPathList.append(plotPath)
			if frac=='COMP':
				plotPath=getPlots(bestPerComp,LightCurve,star,frac)
				plotPathListComp.append(plotPath)
				plotPathList.append(plotPath)
	savePNGsToSinglePDF(plotPathList13,'BLS_13')
	savePNGsToSinglePDF(plotPathList24,'BLS_24')
	savePNGsToSinglePDF(plotPathListComp,'BLS_COMP')
	savePNGsToSinglePDF(plotPathList,'BLS_Full')



def getBestFreq(frequency,power):
	idx=np.argmax(power)-1
	if idx<0:
		return frequency[0]
	if idx+2>=len(frequency):
		return frequency[-1]
	negb=0.5*(power[idx]-power[idx+2])
	c=power[idx]-2.*power[idx+1]+power[idx+2]
	return frequency[idx+1]+0.5*(frequency[idx+2]-frequency[idx])*negb/c

def getPlots(period,LCData,metaRow,BLSType):
	apogeeid=metaRow['APOGEE_ID']
	ticid=metaRow['TICID']
	modes=period*np.array([0,1./3,1./2,2./3,1.,3./2,2.,3.])
	#Super Title
	supTitle=f'Light Curves Folded @ Modes of sections of BLS Periodogram\n\
	APOGEE ID: %s\n TIC ID: %s'%(str(apogeeid),str(ticid))
	#Figure
	fig,axes=plt.subplots(nrows=4,ncols=2,figsize=(32,10),\
		facecolor='w',constrained_layout=True,dpi=75)
	fig.suptitle(supTitle)
	for i,ax in enumerate(axes.flat):
		print(i,ax)
		if i==0:
			print('Raw Data')
			plot=ax.scatter(LCData['Time'],LCData['Flux'],\
				marker='o',s=3,edgecolor='None',alpha=1,\
				c=LCData['Time'],cmap='plasma')
			ax.set_xlabel('Period (d)')
			ax.set_ylabel('Flux')
			ax.set_title('Raw Data')
			fig.colorbar(plot,ax=axes.ravel().tolist(),\
				label='TESS Time (d)')
		else:
			print('else')
			foldFrac=modes[i]/period
			subTitle=f'Folded @ %.2f (%.2f Times Best Period)'%(modes[i],foldFrac)
			print(subTitle)
			plot=ax.scatter((LCData['Time']%modes[i]),\
				LCData['Flux'], marker='o',s=3,\
				edgecolor='None',alpha=0.5,c=LCData['Time'],cmap='plasma')
			ax.set_xlabel(f'Phase (d %% %.2f)'%(modes[i]))
			ax.set_ylabel('Flux')
			ax.set_title(subTitle)
	pngFile=f'/scratch/jdg577/theJoker/Plots/starPlots/%s_%s/%s_%s-%s_Fold.png'%(str(apogeeid),str(ticid),str(apogeeid),str(ticid),BLSType)
	if os.path.exists(f'/scratch/jdg577/theJoker/Plots/starPlots/%s_%s/'%(str(apogeeid),str(ticid)))==False:
		os.mkdir(f'/scratch/jdg577/theJoker/Plots/starPlots/%s_%s/'%(str(apogeeid),str(ticid)))
	else:
		pass
	fig.savefig(pngFile)
	print(supTitle+'Done!')
	return pngFile

def savePNGsToSinglePDF(pathList,BLSType):
	pdf=FPDF(unit='in',format=[11,8.5])
	for subPng in pathList:
		pdf.add_page()
		pdf.image(subPng,w=10,h=3.125)
		print(f'%s Saved to PDF'%subPng)
	pdf.output(f'/scratch/jdg577/theJoker/Plots/starPlots/PlotTypePDFs/%s_Fold_From_PNGs.pdf'%BLSType)

main()

