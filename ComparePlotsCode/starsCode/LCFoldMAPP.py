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
	starsFile='/scratch/jdg577/theJoker/Data/starsData/apogee_id-ticid-filtered.fits'
	stars=at.Table.read(starsFile)
	plotPathList=[]
	for star in stars:
		apogeeid=star['APOGEE_ID']
		ticid=star['TICID']

		lightCurveFile=f'/scratch/jdg577/theJoker/Data/starsData/%s_%s/%s_%s-LightCurve.fits'\
		%(apogeeid,ticid,apogeeid,ticid)
		lightCurve=at.Table.read(lightCurveFile)
		plotPath=getPlots(star,lightCurve)
		plotPathList.append(plotPath)
	
	savePNGsToSinglePDF(plotPathList)

def getPlots(starData,LCData):
	# Data
	apogeeid=starData['APOGEE_ID']
	ticid=starData['TICID']
	MAPP=starData['MAP_P']
	modes=MAPP*np.array([0,1./3.,1./2.,2./3.,1.,3./2.,2.,3.])
	print(modes)
	#Titles
	supTitle=f'Light Curves Folded @ Modes of MAP_P\n APOGEE \
	ID: %s\n TIC ID: %s' % (str(apogeeid), str(ticid))
	

	#Make Fig
	fig,axes=plt.subplots(nrows=4,ncols=2,figsize=(32,10),\
		facecolor='w',constrained_layout=True,dpi=75.0) #DPI 25 for regular 100+ for poster
	fig.suptitle(supTitle)

	#make sub plots
	
	for i,ax in enumerate(axes.flat):
		print(i,ax)
		if i==0:
			print('RAW DATA')
			plot=ax.scatter(LCData['Time'],LCData['Flux'],\
				marker='o',s=3,edgecolor='None',alpha=1,\
				c=LCData['Time'],cmap='plasma')
			ax.set_xlabel('Period (d)')
			ax.set_ylabel('Flux')
			ax.set_title('Raw Data')
			fig.colorbar(plot,ax=axes.ravel().tolist(),label='TESS Time (d)')
		else:
			subTitle=f'Folded @ %.2f (d)'%(modes[i])
			print(subTitle)
			plot=ax.scatter((LCData['Time']%modes[i]),\
				LCData['Flux'],marker='o',s=3,edgecolor='None',\
				alpha=0.5,c=LCData['Time'],cmap='plasma')
			ax.set_xlabel(f'Phase (d %% %.2f)'%(modes[i]))
			ax.set_ylabel('Flux')
			foldFrac=modes[i]/MAPP
			ax.set_title(f'Folded @ %.2f (%.2f Times MAP_P)'%(modes[i],foldFrac))
		pngFile='/scratch/jdg577/theJoker/Plots/starPlots/'\
			+str(apogeeid)+'_'+str(ticid)+'/'\
		 	+str(apogeeid)+'_'+str(ticid)+'-MAPP_Fold.png'
		if os.path.exists('/scratch/jdg577/theJoker/Plots/starPlots/'+str(apogeeid)+'_'+str(ticid)+'/')==False:
			os.mkdir('/scratch/jdg577/theJoker/Plots/starPlots/'+str(apogeeid)+'_'+str(ticid)+'/')
		else:
			pass
	fig.savefig(pngFile)
	print(supTitle+' DONE')
	return pngFile

def savePDFsToSinglePDF(plots):
	# Deprecated b/c file too big 
	# Maybe for use as for paper/poster b/c images are high quality
	fullPDFFIle='/scratch/jdg577/theJoker/Plots/starPlots/PlotTypePDFs/MAPP_Fold.pdf'
	merger=PdfFileMerger()
	for subPDF in plots:
		merger.append(PdfFileReader(subPDF,'rb'))
	merger.write(fullPDFFIle)

def savePNGsToSinglePDF(plots):
	pdf=FPDF(unit='in',format=[11,8.5])
	for subPNG in plots:
		pdf.add_page()
		pdf.image(subPNG,w=10,h=3.125)
		print(f'%s Saved to PDF'%subPNG)
	pdf.output('/scratch/jdg577/theJoker/Plots/starPlots/PlotTypePDFs/MAPP_Fold_From_PNGs.pdf','F')
	print('Done!')

main()

