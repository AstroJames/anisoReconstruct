"""

Module:         Testing Anisotropic Reconstructions
Author:         James R. Beattie
First Created:  22 / 10 / 2019

File:           main file for anisotropic (and isotropic) reconstructions
Description:    this is the main file

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

########################################################################################################################

ap 			= argparse.ArgumentParser(description = 'Just a bunch of input arguments')
ap.add_argument('-viz','--viz',default=False,help='an argument for turning the visualisation plots on: True / False',type=bool)
ap.add_argument('-print','--print',default=True,help='an argument for printing out the variance recon. output: True / False',type=bool)
args 		= vars(ap.parse_args())

########################################################################################################################

if __name__ == "__main__":

    # Load in the data
    MachData    = "M2MA0.1"
    iter        = "50"
    file3DDir   = "./testData/"+ MachData +"/3DDens/"
    fileProjDir = "./testData/"+ MachData +"/Proj/"

    # Calculate the 3D variance and delete the 3D field to save memory
    dens3D      = load_obj(file3DDir + "Turb_hdf5_plt_cnt_00{}_{}_dens".format(iter,MachData))
    var3D       = (dens3D / dens3D.mean() - 1).var()
    del dens3D

    # Load the density projections and power spectra
    densProj    = load_obj(fileProjDir + "Turb_hdf5_plt_cnt_00{}_{}_proj".format(iter,MachData))
    varProj     = (densProj / densProj.mean() - 1).var()

    # Calculate the 2D power spectrum
    densPS, k   = PSF.PowerSpectrum(densProj)
    logDensPS   = np.log10(densPS)

    # Define some visualisation parameters
    sigma       = 10                        # for Guassian kernel for contour plots
    isobars     = np.array([0,-20,40])      # define the domain for the isobars on contour plots in log10 [to, from, by]

    # Plot the column density and the power spectrum
    if args['viz'] == True:
        PF.kAndDensityFieldPlot(densProj,densPS,logDensPS,sigma,isobars,fileName)

    # Fit and Calculate the ellipses to all of the k-space power spectra
    aniso, anisoStd, prinAxisKeep, center, kperpAxis, kparAxis  = EF.calculateEllipses(logDensPS,sigma,isobars,"kperp")

    # Calculate the anisotropic variances for the prolate and oblate ellipse.
    varProlate, varOblate,  rr, cc, relError, error             = AR.calculateAnisoVar(densPS,center,aniso,viz=args['viz'])

    # Calculate the isotropic variance for comparison.
    R, var3DIso                                                 = PSF.calculateIsoVar(densPS,k,varProj)

    # Print off the values for the reconstructions
    if args['print'] == True:
        print("Relative residual power: {}".format(relError))
        print("Absolute residual power correction: {}".format(error))
        print("The 2D projection variance is: {}".format(varProj))
        print("The 3D variance is: {}".format(var3D))
        print("The variance from the 3D prolate: {}".format(varProlate))
        print("The variance from the 3D oblate: {}".format(varOblate))
        print("The variance from the isotropic recon: {}".format(var3DIso))
