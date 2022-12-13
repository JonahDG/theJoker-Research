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
from PyPDF2 import PdfFileMerger,PdfFileReader
import glob

def main():
    samplesFile='/scratch/jdg577/theJoker/Data/samples'
    separationFile='/scratch/jdg577/theJoker/Data/allStarLite-r12-l33-tess_2min-max_20arcsec-xm.fits'
    apogeeVisitsFile='/scratch/jdg577/theJoker/Data/allVisit-r12-l33.fits'
    allStarLite='/scratch/jdg577/theJoker/Data/metadataDR17.fits'
    apogeeWTICID=matchApogeeDataTICID(apogeeVisitsFile,allStarLite,separationFile)
    jokerSampling=getJokerSampling(samplesFile)
    allData=combineFilterData(apogeeWTICID,jokerSampling)
    print(allData.columns)
    print(allData.dtypes)
    #makeFits(allData)


# adds both the separation data (for filtering later) and TICID to apoogee visits
def matchApogeeDataTICID(visitsFile, starFile,separationFile):
	visits=at.Table.read(visitsFile) # I don't think this is necessary
	stars=at.QTable.read(starFile)
	seps=at.QTable.read(separationFile)
	apogeeData_TIC=at.join(stars,seps,keys='APOGEE_ID')
	return apogeeData_TIC

# Uses table of TICIDS matched to apogee visits to produce light curve data of every star in list
# Then using LC data power and period from LS periodogram is computer
def getLCData(apogeeTic): # Input Data set subject to change based on filtering needs 
	ticList=[]
	timeList=[]
	fluxList=[]
	fluxErrList=[]
	perList=[]
	powList=[]
	LCData=at.Table()
	for i in range(len(apogeeTic['TICID'])):
		tic=apogeeTic[i]['TICID']
		ticStr='TIC'+str(tic)
		try:
			# Light curve
			lcCollection=lk.search_lightcurve(target=ticStr,mission='TESS').download_all()
			lcStitch=lcCollection.stitch().remove_nans().remove_outliers()
			time=np.ascontiguousarray(lcStitch.time.value)
			flux=np.ascontiguousarray(lcStitch.flux)
			fluxErr=lcStitch.flux_err
			# lomb scargle periodogram
			mapP=float(apogeeTic[i]['MAP_P']/u.d)
			freq,power=LombScargle(time,flux).autopower(minimum_frequency=(0.05/mapP),maximum_frequency=(20./mapP))
			period=1./freq
			# add to lists
			ticList.append(tic)
			timeList.append(time)
			fluxList.append(flux)
			fluxErrList.append(fluxErr)
			perList.append(period)
			powList.append(power)


			print('LIGHT CURVE AND PERIODOGRAM SUCCESSFUL ON '+ticStr)
		except:
			print('LIGHT CURVE AND/OR PERIDOGRAM FAILED ON '+ticStr)
			pass
	LCData['TICID']=ticList
	LCData['Time']=timeList
	LCData['Flux']=fluxList
	LCData['FluxErr']=fluxErrList
	LCData['LSPeriod']=perList
	LCData['Power']=powList
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
		sampleIDList.append(sampleID)
		sampling=at.Table.read(samplesPath)
		K16percentile=np.percentile(sampling['K'],16)
		perList.append(sampling['P'])
		eList.append(sampling['e'])
		omegaList.append(sampling['omega'])
		M0List.append(sampling['M0'])
		sList.append(sampling['s'])
		KList.append(sampling['K'])
		K16percentileList.append(K16percentile)
		v0List.append(sampling['v0'])
		ln_priorList.append(sampling['ln_prior'])
		ln_likelihoodList.append(sampling['ln_likelihood'])
	jokerSampling['APOGEE_ID']=sampleIDList
	jokerSampling['tjP']=perList
	jokerSampling['e']=eList
	jokerSampling['omega']=omegaList
	jokerSampling['M0']=M0List
	jokerSampling['s']=sList
	jokerSampling['tjK']=KList
	jokerSampling['K16percentile']=K16percentileList
	jokerSampling['v0']=v0List
	jokerSampling['ln_prior']=ln_priorList
	jokerSampling['ln_likelihood']=ln_likelihoodList
	return jokerSampling

def combineFilterData(apogeeDataWTic,jokerData):
	apogeeJoker=at.join(apogeeDataWTic,jokerData,keys='APOGEE_ID')
	apogeeJoker=apogeeJoker[apogeeJoker['MAP_P']<10*u.d]
	apogeeJoker=apogeeJoker[apogeeJoker['separation']<2*u.arcsec]
	apogeeJoker=apogeeJoker[apogeeJoker['K16percentile']>1]
	lcData=getLCData(apogeeJoker)
	unsortedCombinedData=at.join(apogeeJoker,lcData,keys='TICID')
	combinedData=unsortedCombinedData['APOGEE_ID','TICID','n_visits','separation',
	'MAP_P','MAP_P_err','MAP_e','MAP_e_err','MAP_omega','MAP_omega_err','MAP_M0','MAP_M0_err',
	'MAP_K','MAP_K_err','MAP_v0','MAP_v0_err','MAP_s','MAP_s_err','MAP_t0_bmjd','t_ref_bmjd','baseline','MAP_ln_likelihood',
	'MAP_ln_prior','max_unmarginalized_ln_likelihood','max_phase_gap','periods_spanned','phase_coverage',
	'phase_coverage_per_period','unimodal','joker_completed','mcmc_completed','mcmc_status','gelman_rubin_max',
	'constant_ln_likelihood','robust_constant_ln_likelihood','robust_constant_mean','robust_constant_scatter',
	'robust_constant_success','robust_linear_ln_likelihood','robust_linear_a','robust_linear_b','robust_linear_scatter',
	'robust_linear_success','tjP','e','omega','M0','s','tjK','K16percentile','v0','ln_prior',
	'ln_likelihood','Time','Flux','FluxErr','LSPeriod','Power']
	return combinedData

def makeFits(fullData):
	fitsFile='/scratch/jdg577/theJoker/Data/crossmatch-apogee-tess_LC_joker.fits' #need better name
	fullData.write(fitsFile,format='fits')


main()


	
