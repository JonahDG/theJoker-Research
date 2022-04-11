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
from fpdf import FPDF
#endregion
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
    sources=sources[0]
    print('Data Tables Formed')
    return apogeeVisits, sources

apogeeFile='/scratch/jdg577/theJoker/Data/allVisit-r12-l33.fits'
tessFile='/scratch/jdg577/theJoker/Data/allStarLite-r12-l33-tess_2min-max_20arcsec-xm.fits'
binaryMetadataFile='/scratch/jdg577/theJoker/Data/allStarLite-metadata.fits'

apogeeData,sourceData=getStarData(apogeeFile,tessFile,binaryMetadataFile)
print(apogeeData)
print(sourceData)
for row in sourceData:
    visits=apogeeData[apogeeData['APOGEE_ID']==row['APOGEE_ID']]
    print(visits)