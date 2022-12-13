from fpdf import FPDF
from astropy.io import fits
import astropy.table as at 
starsFile='/scratch/jdg577/theJoker/Data/starsData/apogee_id-ticid-filtered.fits'
stars=at.Table.read(starsFile)
for star in stars:
    apogeeid=star['APOGEE_ID']
    ticid=star['TICID']
    