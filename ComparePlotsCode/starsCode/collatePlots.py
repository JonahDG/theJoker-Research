from fpdf import FPDF
from astropy.io import fits
import astropy.table as at 
def main():
    plotFileList=[]
    starsFile='/scratch/jdg577/theJoker/Data/starsData/apogee_id-ticid-filtered.fits'
    stars=at.Table.read(starsFile)
    for star in stars:
        apogeeid=star['APOGEE_ID']
        ticid=star['TICID']
        print(apogeeid,ticid)
        comparisonPlotFile=f'/scratch/jdg577/theJoker/Plots/starPlots/%s_%s/%s_%s-Compare_Periodograms_Joker.png'\
        %(apogeeid,ticid,apogeeid,ticid)
        plotFileList.append(comparisonPlotFile)
        mappPlotFile=f'/scratch/jdg577/theJoker/Plots/starPlots/%s_%s/%s_%s-MAPP_Fold.png'\
        %(apogeeid,ticid,apogeeid,ticid)
        plotFileList.append(mappPlotFile)
        LC1FoldPlotFile=f'/scratch/jdg577/theJoker/Plots/starPlots/%s_%s/%s_%s-LC1_Fold.png'\
        %(apogeeid,ticid,apogeeid,ticid)
        plotFileList.append(LC1FoldPlotFile)
        LC4FoldPlotFile=f'/scratch/jdg577/theJoker/Plots/starPlots/%s_%s/%s_%s-LC4_Fold.png'\
        %(apogeeid,ticid,apogeeid,ticid)
        plotFileList.append(LC4FoldPlotFile)
        LC16FoldPlotFile=f'/scratch/jdg577/theJoker/Plots/starPlots/%s_%s/%s_%s-LC16_Fold.png'\
        %(apogeeid,ticid,apogeeid,ticid)
        plotFileList.append(LC16FoldPlotFile)
        # BLS02FoldFile=f'/scratch/jdg577/theJoker/Plots/starPlots/%s_%s/%s_%s-'\
        # %(apogeeid,ticid,apogeeid,ticid)
        # plotFileList.append(BLS02FoldFile)
        BLS13FoldFile=f'/scratch/jdg577/theJoker/Plots/starPlots/%s_%s/%s_%s-13_Fold.png'\
        %(apogeeid,ticid,apogeeid,ticid)
        plotFileList.append(BLS13FoldFile)
        BLS24FoldFile=f'/scratch/jdg577/theJoker/Plots/starPlots/%s_%s/%s_%s-24_Fold.png'\
        %(apogeeid,ticid,apogeeid,ticid)
        plotFileList.append(BLS24FoldFile)
        BLSCOMPFoldFile=f'/scratch/jdg577/theJoker/Plots/starPlots/%s_%s/%s_%s-COMP_Fold.png'\
        %(apogeeid,ticid,apogeeid,ticid)
        plotFileList.append(BLSCOMPFoldFile)
    savePNGsToPDF(plotFileList)

def savePNGsToPDF(fileList):
    pdf=FPDF(unit='in',format=[8.5,11])
    for subPNG in fileList:
        print(subPNG,' Creating')
        pdf.add_page()
        pdf.image(subPNG,w=6.5,h=10)
        print(subPNG, 'Added')
    pdf.output('/scratch/jdg577/theJoker/Plots/starPlots/PlotTypePDFs/collatedPlots.pdf')
    print('ALL DONE! :)')

main()  