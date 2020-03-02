"""

Module:         Testing Anisotropic Reconstructions
Author:         James R. Beattie
First Created:  29 / 01 / 2020

File:           the data processing file for the anisotropic reconstructions
Description:    this file processes the .dat files generated from main.py

"""

########################################################################################################################

import PlottingFunctions as PF
reload(PF)
import EllipseFunctions as EF
reload(EF)
import PowerSpectrumFunctions as PSF
reload(PSF)
import AnisoReconFunctions as AR
reload(AR)

from PickleData import *
from header import *
from scipy.stats import pearsonr as cor

########################################################################################################################

ap 			= argparse.ArgumentParser(description = 'Just a bunch of input arguments')
ap.add_argument('-viz','--viz',default=False,help='an argument for turning the visualisation plots on: True / False',type=bool)
ap.add_argument('-print','--print',default=True,help='an argument for printing out the variance recon. output: True / False',type=bool)
args 		= vars(ap.parse_args())

########################################################################################################################


# Functions
########################################################################################################################

def lineReader(file,type=None):
    # open a file
    g = open(file,'r')
    lineCounter = 0

    # define the new variables for each file
    lineName        = []
    lineFileCount   = []
    line3DVar1Sigma = []
    line3DVar       = []
    lineReconVar    = []
    lineReconVarSD  = []
    lineAnisoVar    = []
    lineProlateVar  = []
    lineOblateVar   = []

    line3DSkew             = []
    lineReconSkew          = []
    lineReconSkewSD        = []
    lineReconSkew1Sigma    = []
    lineAnisoSkew          = []

    line3DKurt             = []
    lineReconKurt          = []
    lineReconKurtSD        = []
    lineReconKurt1Sigma    = []
    lineAnisoKurt          = []

    for line in g:
        if lineCounter == 0:
            lineCounter +=1
            continue

        line  = line.strip('\n').split(',')

        # File information and anisotropy
        lineName.append(line[0])
        lineFileCount.append(float(line[1]))
        # Variance
        lineAnisoVar.append(float(line[11]))
        line3DVar.append(float(line[5]))
        lineReconVar.append(float(line[8]))
        lineReconVarSD.append(float(line[9]))
        lineProlateVar.append(float(line[6]))
        lineOblateVar.append(float(line[7]))
        # Skewness
        lineAnisoSkew.append(float(line[19]))
        line3DSkew.append(float(line[14]))
        lineReconSkew.append(float(line[17]))
        lineReconSkewSD.append(float(line[18]))
        # Kurtosis
        lineAnisoKurt.append(float(line[27]))
        line3DKurt.append(float(line[22]))
        lineReconKurt.append(float(line[25]))
        lineReconKurtSD.append(float(line[26]))

        # add here for skewness and kurtosis

        lineCounter+=1

    if type != None:
        g.close()
        return lineName, lineFileCount, line3DVar, line3DKurt, line3DSkew, lineAnisoKurt, lineAnisoSkew, lineAnisoVar, lineReconSkew, lineReconKurt, lineReconVar, lineProlateVar, lineOblateVar

    # Correlations are taken before averaging
    fileCorVar.append(cor(np.array(line3DVar),np.array(lineReconVar)))

    # take means and 1sigma fluctations of variance
    line3DVar1Sigma     = np.std(np.array(line3DVar))
    line3DVar           = np.mean(np.array(line3DVar))
    lineReconVar1Sigma  = np.std(np.array(lineReconVar))
    lineReconVar        = np.mean(np.array(lineReconVar))
    lineReconVarSD      = np.mean(np.array(lineReconVarSD))
    lineAnisoVar        = np.mean(np.sqrt(1 - np.array(lineAnisoVar)))

    # take means and 1sigma fluctations of skewness
    line3DSkew1Sigma     = np.std(np.array(line3DSkew))
    line3DSkew           = np.mean(np.array(line3DSkew))
    lineReconSkew1Sigma  = np.std(np.array(lineReconSkew))
    lineReconSkew        = np.mean(np.array(lineReconSkew))
    lineReconSkewSD      = np.mean(np.array(lineReconSkewSD))
    lineAnisoSkew        = np.mean(np.sqrt(1 - np.array(lineAnisoSkew)))

    # take means and 1sigma fluctations of kurtosis
    line3DKurt1Sigma     = np.std(np.array(line3DKurt))
    line3DKurt           = np.mean(np.array(line3DKurt))
    lineReconKurt1Sigma  = np.std(np.array(lineReconKurt))
    lineReconKurt        = np.mean(np.array(lineReconKurt))
    lineReconKurtSD      = np.mean(np.array(lineReconKurtSD))
    lineAnisoKurt        = np.mean(np.sqrt(1 - np.array(lineAnisoKurt)))

    # append to the global dataset
    fileMach.append(machData[lineName[0]]['M'])
    fileMach1Sigma.append(machData[lineName[0]]['MStd'])
    fileMA.append(machData[lineName[0]]['MA'])
    fileName.append(lineName[0])

    # Variance
    file3DVar.append(line3DVar)
    file3DVar1Sigma.append(line3DVar1Sigma)
    fileReconVar.append(lineReconVar)
    fileReconVarSD.append(lineReconVarSD)
    fileReconVar1Sigma.append(lineReconVar1Sigma)
    fileAnisoVar.append(lineAnisoVar)

    # Skewness
    file3DSkew         .append(line3DSkew)
    file3DSkew1Sigma   .append(line3DSkew1Sigma)
    fileReconSkew      .append(lineReconSkew)
    fileReconSkewSD    .append(lineReconSkewSD)
    fileReconSkew1Sigma.append(lineReconSkew1Sigma)
    fileAnisoSkew      .append(lineAnisoSkew)

    # Kurtosis
    file3DKurt         .append(line3DKurt)
    file3DKurt1Sigma   .append(line3DKurt1Sigma)
    fileReconKurt      .append(lineReconKurt)
    fileReconKurtSD    .append(lineReconKurtSD)
    fileReconKurt1Sigma.append(lineReconKurt1Sigma)
    fileAnisoKurt       .append(lineAnisoKurt)

    g.close()

def oneToOnePlot(measured,measure1Sigma,recon,recon1Sigma,output,type):
    recon       = np.array(recon)
    measured    = np.array(measured)

    x       = np.linspace(0,25,100)
    f, ax   = plt.subplots(1,1,dpi=200)
    ax.plot(x,x,ls='--',color='black',zorder=1)
    ax.plot(x,x/2.,ls='--',color='blue',zorder=1)
    ax.plot(x,2*x,ls='--',color='blue',zorder=1)
    ax.errorbar(y=recon, x=measured, xerr=measure1Sigma, yerr=recon1Sigma,linestyle="None",color='black',elinewidth=1,capsize=2,alpha=0.4,zorder=2)
    p1      = ax.scatter(y=recon, x=measured,c=fileMA,cmap=plt.cm.plasma,zorder=3,norm=colors.LogNorm(vmin=min(fileMA), vmax=max(fileMA)))
    cbar    = plt.colorbar(p1)
    cbar.set_label(r"$\mathcal{M}_{A0}$",fontsize=16,x=-0.5)
    if type == "var":
        ax.set_xlabel(r"$\sigma^2_{\text{true}}$",fontsize=fs)
        ax.set_ylabel(r"$\sigma^2_{\text{recon}}$",fontsize=fs)
    elif type == "skew":
        ax.set_xlabel(r"$\mathscr{S}_{\text{true}}$",fontsize=fs)
        ax.set_ylabel(r"$\mathscr{S}_{\text{recon}}$",fontsize=fs)
    elif type == "kurt":
        ax.set_xlabel(r"$\mathscr{K}_{\text{true}}$",fontsize=fs)
        ax.set_ylabel(r"$\mathscr{K}_{\text{recon}}$",fontsize=fs)

    #ax.set_xlim(0.2,20)
    #ax.set_ylim(0.2,20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.savefig(output,dpi=200)
    plt.close()

def timeDomainPlot(fileName):
    lineName, lineFileCount, line3DVar, line3DKurt, line3DSkew, lineAnisoKurt, lineAnisoSkew, lineAnisoVar, lineReconSkew, lineReconKurt, lineReconVar, lineProlateVar, lineOblateVar = lineReader(fileName,type=1)

    lineFileCount   = np.array(lineFileCount)
    line3DVar       = np.array(line3DVar)
    lineProlateVar  = np.array(lineProlateVar)
    lineOblateVar   = np.array(lineOblateVar)
    lineReconVar    = np.array(lineReconVar)

    f, ax   = plt.subplots(1,1,dpi=200,sharex=True)
    f.subplots_adjust(wspace=0.32,hspace=0.05,top=0.98,right=0.99,left=0.1,bottom=0.12)
    # ax.plot(x,x,ls='--',color='black',zorder=1)
    # ax.plot(x,x/2.,ls='--',color='blue',zorder=1)
    # ax.plot(x,2*x,ls='--',color='blue',zorder=1)
    # ax[0].plot(lineFileCount*0.1,lineAnisoVar)
    # ax[0].set_ylabel(r"$e$",fontsize=fs)
    # ax[0].set_xlabel(r"$t/T$",fontsize=fs)

    ax.plot(lineFileCount*0.1,line3DVar,label=r"true $\sigma_{3D}^2$",color="k")
    ax.plot(lineFileCount*0.1,lineProlateVar,label=r"$\sigma_{3D,p}^2$",linestyle="--",color="r")
    ax.plot(lineFileCount*0.1,lineOblateVar,label=r"$\sigma_{3D,o}^2$",linestyle="--",color="b")
    ax.set_ylabel(r"$\sigma^2$",fontsize=fs)
    ax.set_xlabel(r"$t/T$",fontsize=fs)
    plt.legend()
    plt.savefig("sigmaAsTurnover_{}.png".format(fileName),dpi=200)
    plt.close()

########################################################################################################################


machData = {"M2MA0.1":{"M":2.6,"MStd":0.2,"MA":0.133},"M4MA0.1":{"M":5.2,"MStd":0.4,"MA":0.13},
              "M10MA0.1":{"M":12,"MStd":1,"MA":0.125}, "M20MA0.1":{"M":24,"MStd":1,"MA":0.119},
              "M2MA0.5":{"M":2.2,"MStd":0.2,"MA":0.54},"M4MA0.5":{"M":4.4,"MStd":0.4,"MA":0.54},
              "M10MA0.5":{"M":10.5,"MStd":0.5,"MA":0.52}, "M20MA0.5":{"M":21,"MStd":1,"MA":0.53},
              "M2MA1":{"M":2.0,"MStd":0.1,"MA":0.98},"M4MA1":{"M":3.8,"MStd":0.3,"MA":0.95},
              "M10MA1":{"M":9.3,"MStd":0.5,"MA":0.93}, "M20MA1":{"M":19,"MStd":1,"MA":0.93},
              "M2MA2":{"M":1.66,"MStd":0.05,"MA":1.7},"M4MA2":{"M":3.5,"MStd":0.1,"MA":1.73},
              "M10MA2":{"M":9.0,"MStd":0.4,"MA":1.8}, "M20MA2":{"M":18,"MStd":1,"MA":1.8},
              "M2MA10":{"M":1.80,"MStd":0.08,"MA":9.0},"M4MA10":{"M":3.7,"MStd":0.1,"MA":9.2},
              "M10MA10":{"M":9.2,"MStd":0.4,"MA":9.2}, "M20MA10":{"M":19,"MStd":1,"MA":9.3}}


# Global variable dataset
########################################################################################################################

# File Information
fileMach                = []
fileMach1Sigma          = []
fileMA                  = []
fileName                = []

# Variance
file3DVar               = []
file3DVar1Sigma         = []
fileReconVar            = []
fileReconVarSD          = []
fileReconVar1Sigma      = []
fileAnisoVar            = []

# Skewness
file3DSkew              = []
file3DSkew1Sigma        = []
fileReconSkew           = []
fileReconSkewSD         = []
fileReconSkew1Sigma     = []
fileAnisoSkew           = []

# Kurtosis
file3DKurt              = []
file3DKurt1Sigma        = []
fileReconKurt           = []
fileReconKurtSD         = []
fileReconKurt1Sigma     = []
fileAnisoKurt           = []

# Correlations
fileCorVar              = []

os.system("ls output*.dat > TempDatFiles.txt")
f = open("TempDatFiles.txt",'r')
fileCounter     = 0
for file in f:
    file  = file.strip('\n')
    fileCounter +=1
    # line reader fills the global dataset with time averaged quantities
    lineReader(file)
f.close()

oneToOnePlot(file3DVar,file3DVar1Sigma,fileReconVar,fileReconVar1Sigma,"var.png","var")
oneToOnePlot(file3DSkew,file3DSkew1Sigma,fileReconSkew,fileReconSkew1Sigma,"skew.png","skew")
oneToOnePlot(file3DKurt,file3DKurt1Sigma,fileReconKurt,fileReconKurt1Sigma,"kurt.png","kurt")


f = open("TempDatFiles.txt",'r')
fileCounter     = 0
for file in f:
    file  = file.strip('\n')
    fileCounter +=1
    # line reader fills the global dataset with time averaged quantities
    timeDomainPlot(file)
f.close()
