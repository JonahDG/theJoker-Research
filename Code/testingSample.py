import astropy as ap
import astropy.units as u
import astropy.table as at
import astropy.coordinates as coords
from astropy.io import ascii,fits
from astropy.time import Time
from astropy.timeseries import LombScargle
import glob


allStarLite='/scratch/jdg577/theJoker/Data/metadataDR17.fits'
stars=at.QTable.read(allStarLite)
genericPath='/scratch/jdg577/theJoker/Data/samples/**/2M*.fits'
samplesFitsPaths=glob.glob(genericPath)
checkNans=stars['MAP_P_err'].isnull()
nansNum=0
valsNum=0
for check in checkNans:
    if check==True:
       nansNum+=1
    else:
        valsNum+=1
print(nansNum)
print(valsNum)

























