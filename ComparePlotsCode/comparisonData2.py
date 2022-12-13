import astropy as ap 
import astropy.units as u 
import astropy.table as at 
import astropy.coordinates as coords 
from astropy.io import ascii,fits
from astropy.time import Time
from astropy.timeseries import LombScargle, BoxLeastSquares
import matplotlib.pyplot as plt 
import numpy as np 
import lightkurve as lk 
from PyPDF2 import PdfFileMerger,PdfFileReader
import glob

def main():
    samplesFile='/scratch/jdg577/theJoker/Data/samples'
    separationFile='/scratch/jdg577/theJoker/Data/allStarLite-r12-l33-tess_2min-max_20arcsec-xm.fits'
    allStarLite='/scratch/jdg577/theJoker/Data/metadataDR17.fits'
    apogeeWTICID=matchApogeeDataTICID(allStarLite,separationFile)
    jokerSampling=getJokerSampling(samplesFile)
    print(jokerSampling.info)
    allData=combineFilterData(apogeeWTICID,jokerSampling)
    print(allData.info)
    makeFits(allData)

def matchApogeeDataTICID(starFile,separationFile):
	stars=at.QTable.read(starFile)
	seps=at.QTable.read(separationFile)
	apogeeData_TIC=at.join(stars,seps,keys='APOGEE_ID')
	return apogeeData_TIC

def getLCData(apogeeTic): # Input Data set subject to change based on filtering needs 
	ticList=[]
	timeList=[]
	fluxList=[]
	fluxErrList=[]

	LCData=at.Table()
	for i in range(len(apogeeTic['TICID'])):
		tic=apogeeTic[i]['TICID']
		ticStr='TIC'+str(tic)
		try:
			# Light curve
			lcCollection=lk.search_lightcurve(target=ticStr,mission='TESS').download_all()
			lcStitch=lcCollection.stitch().remove_nans().remove_outliers()
			time=np.ascontiguousarray(lcStitch.time.value)
			print(time.shape)
			flux=np.ascontiguousarray(lcStitch.flux)
			print(flux.shape)
			fluxErr=lcStitch.flux_err
			ticList.append(tic)
			print(ticList)
			timeList.append(time)
			fluxList.append(flux)
			fluxErrList.append(fluxErr)
			print('LIGHT CURVE SUCCESSFUL ON '+ticStr+' | '+str(apogeeTic[i]['APOGEE_ID']))
		except:
			print('LIGHT CURVE Failed ON '+ticStr+' | '+str(apogeeTic[i]['APOGEE_ID']))
			pass
	print(ticList[:4])
	print(timeList[:4])
	LCData['TICID']=ticList
	LCData['Time']=timeList
	LCData['Flux']=fluxList
	LCData['FluxErr']=fluxErrList
	print(LCData.info)
	return LCData

def getJokerSampling(samplesFolder):
	genericPath=samplesFolder+'/**/2M*.fits' # This ignores any apogee ids that aren't 2Mass but that is an negligible amount
	'''I did a quick check of how many non-2mass stars made it through the <10 days MAP P filter and the <2 arcsec separation
	filter and it was just one star, so I assume (at least for now) that it is prettty safe to ignore'''
	samplesPaths=glob.glob(genericPath)
	sampleIDList=[]
	perList=[]
	eList=[]
	omegaList=[]
	M0List=[]
	sList=[]
	KList=[]
	K16percentileList=[]
	v0List=[]
	ln_priorList=[]
	ln_likelihoodList=[]

	jokerSampling=at.Table()
	for samplesPath in samplesPaths:
		idx=samplesPath.find('2M')
		sampleID=samplesPath[idx:-5]
		sampling=at.Table.read(samplesPath)
		sampleIDList.append(sampleID)
		K16percentile=np.percentile(sampling['K'],16)
		perList.append(sampling['P'][0])
		eList.append(sampling['e'][0])
		omegaList.append(sampling['omega'][0])
		M0List.append(sampling['M0'][0])
		sList.append(sampling['s'][0])
		KList.append(sampling['K'][0])
		K16percentileList.append(K16percentile)
		v0List.append(sampling['v0'][0])
		ln_priorList.append(sampling['ln_prior'][0])
		ln_likelihoodList.append(sampling['ln_likelihood'][0])
	# print([perList])
	sampleIDCol=at.Column(sampleIDList,dtype=str)
	perCol=at.Column(perList,dtype=float,unit=u.d)
	eCol=at.Column(eList,dtype=float)
	omegaCol=at.Column(omegaList,dtype=float,unit=u.rad)
	M0Col=at.Column(M0List,dtype=float,unit=u.rad)
	sCol=at.Column(sList,dtype=float,unit=u.km/u.s)
	KCol=at.Column(KList,dtype=float,unit=u.km/u.s)
	K16percentileCol=at.Column(K16percentileList,dtype=float)
	v0Col=at.Column(v0List,dtype=float,unit=u.km/u.s)
	ln_priorCol=at.Column(ln_priorList,dtype=float)
	ln_likelihoodCol=at.Column(ln_likelihoodList,dtype=float)

	jokerSampling['APOGEE_ID']=sampleIDCol
	jokerSampling['tjP']=perList
	jokerSampling['e']=eCol
	jokerSampling['omega']=omegaCol
	jokerSampling['M0']=M0Col
	jokerSampling['s']=sCol
	jokerSampling['tjK']=KCol
	jokerSampling['K16percentile']=K16percentileCol
	jokerSampling['v0']=v0Col
	jokerSampling['ln_prior']=ln_priorCol
	jokerSampling['ln_likelihood']=ln_likelihoodCol
	return jokerSampling

def combineFilterData(apogeeDataWTic,jokerData):
	apogeeJoker=at.join(apogeeDataWTic,jokerData,keys='APOGEE_ID')
	apogeeJoker=apogeeJoker[apogeeJoker['MAP_P']<10*u.d]
	apogeeJoker=apogeeJoker[apogeeJoker['separation']<2*u.arcsec]
	apogeeJoker=apogeeJoker[apogeeJoker['K16percentile']>1]
	lcData=getLCData(apogeeJoker)
	unsortedCombinedData=at.join(apogeeJoker,lcData,keys='TICID')
	combinedData=unsortedCombinedData['APOGEE_ID','TICID','separation','n_visits',\
	'MAP_P','MAP_P_err','MAP_e','MAP_e_err','MAP_omega','MAP_omega_err','MAP_M0','MAP_M0_err',\
	'MAP_K','MAP_K_err','MAP_v0','MAP_v0_err','MAP_s','MAP_s_err','MAP_t0_bmjd','t_ref_bmjd','baseline','MAP_ln_likelihood',\
	'MAP_ln_prior','max_unmarginalized_ln_likelihood','max_phase_gap','periods_spanned','phase_coverage',\
	'phase_coverage_per_period','unimodal','joker_completed','mcmc_completed','mcmc_status','gelman_rubin_max',\
	'constant_ln_likelihood','robust_constant_ln_likelihood','robust_constant_mean','robust_constant_scatter',\
	'robust_constant_success','robust_linear_ln_likelihood','robust_linear_a','robust_linear_b','robust_linear_scatter',\
	'robust_linear_success','tjP','e','omega','M0','s','tjK','K16percentile','v0','ln_prior','ln_likelihood',\
	'Time','Flux','FluxErr']
	return combinedData

# def periodogramData(fullData):
# 	necessaryData=fullData['APOGEE_ID','TICID','MAP_P','MAP_P_err','Time','Flux','FluxErr']
# 	periodogramData=at.Table()
# 	idArray=[]
# 	BLSarray=[]
# 	LS1array=[]
# 	LS4array=[]
# 	LS16array=[]
# 	for apogeeID in necessaryData['APOGEE_ID']:
		

def makeFits(fullData):
	# check if finite

	fitsFile='/scratch/jdg577/theJoker/Data/crossmatch-apogee-tess_LC_joker.fits' #need better name
	fullData.write(fitsFile,overwrite=True)

main()
