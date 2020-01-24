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

    # Read directory
    readDir     = "/Volumes/JamesBe/MHD/"
    writeDir    = "./reconStructData/"
    outData     = {} # dictionary for all of the output data




    # Directory names
    dirNames   = ["M2MA1", "M2MA10", "M2MA2",
                  "M4MA0.1","M4MA0.5", "M4MA1", "M4MA10", "M4MA2",
                  "M10MA0.1","M10MA0.5", "M10MA1", "M10MA10", "M10MA2",
                  "M20MA0.1","M20MA0.5", "M20MA1", "M20MA10", "M20MA2"]

    for dir in dirNames:

        print("Overwriting output_{}.dat file".format(dir))
        os.system("touch output_{}.dat".format(dir))
        f = open("output_{}.dat".format(dir),"w+")

        print("Starting on directory: {}".format(dir))

        # File iteration number
        fileNumbers         = xrange(50,101)
        # the file count
        fileCount           = 0

        # for each file
        for fileNumber in fileNumbers:
            fileCount += 1
            try:

                print("Starting on file number: {} in directory: {}".format(fileNumber,dir))

                # Load in the data
                file3DDir   = readDir + dir +"/3DDens/"
                fileProjDir = readDir + dir +"/Proj/"
                try:
                    # Calculate the 3D variance and delete the 3D field to save memory
                    print("Reading in 3D density.")
                    dens3D      = loadObj(file3DDir + "Turb_hdf5_plt_cnt_00{}_{}_dens".format(fileNumber,dir))
                    var3D       = (dens3D / dens3D.mean() - 1).var()
                    del dens3D

                    print("Reading in 2D projection.")
                    # Load the density projections and power spectra
                    densProj    = loadObj(fileProjDir + "Turb_hdf5_plt_cnt_00{}_{}_proj".format(fileNumber,dir))
                    varProj     = (densProj / densProj.mean() - 1).var()
                except:
                    print("Skipping file number: {}".format(fileNumber))

                    # Check if we are at the end of the files (there are 100)
                    if fileNumber > 98:
                        print("Saving object because you have reached then end of the files.")
                        print("------------------------------------------------- \n")
                        f.close()
                        break

                    continue

                print("Calculating 2D power spectra.")
                # Calculate the 2D power spectrum
                densPS, k   = PSF.PowerSpectrum(densProj)
                logDensPS   = np.log10(densPS)

                # Define some visualisation parameters
                sigma       = 10                        # for Guassian kernel for contour plots
                isobars     = np.array([0,-20,40])      # define the domain for the isobars on contour plots in log10 [to, from, by]

                # Plot the column density and the power spectrum
                if args['viz'] == True:
                    print("Creating a denisty and power spectrum plot.")
                    PF.kAndDensityFieldPlot(densProj,densPS,logDensPS,sigma,isobars,fileName)

                print("Fitting ellipses.")
                # Fit and Calculate the ellipses to all of the k-space power spectra
                aniso, anisoStd, prinAxisKeep, center, kperpAxis, kparAxis  = EF.calculateEllipses(logDensPS,sigma,isobars,"kperp")

                print("Reconstructing anisotropic dispersio.n")
                # Calculate the anisotropic variances for the prolate and oblate ellipse.
                varProlate, varOblate,  rr, cc, relError, error, aniso      = AR.calculateAnisoVar(densPS,center,aniso,viz=args['viz'])

                print("Reconstructing isotropic dispersion.")
                # Calculate the isotropic variance for comparison.
                R, var3DIso                                                 = PSF.calculateIsoVar(densPS,k,varProj)

                print("------------------------------------------------- \n")

                # Print off the values for the reconstructions
                if args['print'] == True:
                    print("Relative residual power: {}".format(relError))
                    print("Absolute residual power correction: {}".format(error))
                    print("The 2D projection variance is: {}".format(varProj))
                    print("The 3D variance is: {}".format(var3D))
                    print("The variance from the 3D prolate: {}".format(varProlate))
                    print("The variance from the 3D oblate: {}".format(varOblate))
                    print("The variance from the isotropic recon: {}".format(var3DIso))
                    print("------------------------------------------------- \n")

                if fileCount == 1:
                    f.write("directory, fileNumber, relError, error, varProj, var3D, varProlate, varOblate, var3DIso \n".format(fileNumber))

                f.write("{} {}, {}, {}, {}, {}, {}, {}, {} \n".format(dir,fileNumber,relError,error,varProj,var3D,varProlate,varOblate,var3DIso))
            except:
                print("Skipping number.")
                print("------------------------------------------------- \n")
                continue
