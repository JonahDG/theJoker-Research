import astropy.table as at
import astropy.units as u
import glob
separationFile='/scratch/jdg577/theJoker/Data/allStarLite-r12-l33-tess_2min-max_20arcsec-xm.fits'
separations=at.QTable.read(separationFile)
apogeeVisitsFile='/scratch/jdg577/theJoker/Data/allVisit-r12-l33.fits'
visits=at.QTable.read(apogeeVisitsFile)
allStarLite='/scratch/jdg577/theJoker/Data/metadataDR17.fits'
stars=at.QTable.read(allStarLite)
print(separations.columns)
print(visits.columns)
print(stars.columns)
sources=at.join(stars,separations,keys='APOGEE_ID')
sources=sources[sources['MAP_P']<10*u.d]
sources=sources[sources['separation']<2*u.arcsec]
massCounter=0
nonMassCounter=0
samplesList=glob.glob('/scratch/jdg577/theJoker/Data/samples/**/2M*.fits')
test=at.Table.read(samplesList[0])
#print(test.columns)
