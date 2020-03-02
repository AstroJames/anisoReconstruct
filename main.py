"""

Module:         Testing Anisotropic Reconstructions
Author:         James R. Beattie
First Created:  22 / 10 / 2019

File:           main file for anisotropic (and isotropic) reconstructions
Description:    this is the main file

"""

########################################################################################################################

import PlottingFunctions as PF
import EllipseFunctions as EF
import PowerSpectrumFunctions as PSF
import AnisoReconFunctions as AR

from PickleData import *
from header import *

#########################################################################################################################################

ap 			= argparse.ArgumentParser(description = 'Just a bunch of input arguments')
ap.add_argument('-viz','--viz',default=False,help='an argument for turning the visualisation plots on: True / False',type=bool)
ap.add_argument('-print','--print',default=True,help='an argument for printing out the variance recon. output: True / False',type=bool)
args 		= vars(ap.parse_args())

#########################################################################################################################################

if __name__ == "__main__":

    os.system("touch shell_out.dat")
    g = open("shell_out.dat","w+")
    # Read directory
    readDir     = "/Volumes/JamesBe/MHD/"
    writeDir    = "./reconStructData/"
    outData     = {} # dictionary for all of the output data

    # Directory names
    dirNames   = ["M2MA0.5","M2MA1", "M2MA2", "M2MA10",
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
                    var3D       = (dens3D / dens3D.mean() - 1).var()                    # calculate the true 3D variance
                    skew3D      = skew(dens3D / dens3D.mean() - 1,axis=None)                      # calculate the true 3D skewness
                    kurt3D      = kurtosis(dens3D / dens3D.mean() - 1,fisher=False,axis=None)     # calculate the true 3D kurtosis
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

                ################################################################
                # Dispersion Statistics
                ################################################################

                print("Dispersion: Calculating 2D power spectra")
                # Calculate the 2D power spectrum
                densPS, k   = PSF.PowerSpectrum(densProj)
                logDensPS   = np.log10(densPS)

                # Define some visualisation parameters
                sigma       = 10                        # for Guassian kernel for contour plots
                isobars     = np.array([0,-20,40])      # define the domain for the isobars on contour plots in log10 [to, from, by]

                # Plot the column density and the power spectrum
                if args['viz'] == True:
                    print("Dispersion: Creating a denisty and power spectrum plot.")
                    PF.kAndDensityFieldPlot(densProj,densPS,logDensPS,sigma,isobars,fileName)

                print("Dispersion: Fitting ellipses.")
                # Fit and Calculate the ellipses to all of the k-space power spectra
                aniso, anisoStd, prinAxisKeep, center, kperpAxis, kparAxis  = EF.calculateEllipses(logDensPS,sigma,isobars,"kperp")

                print("Dispersion: Reconstructing anisotropic dispersion.")
                # Calculate the anisotropic variances for the prolate and oblate ellipse.
                varProlate, varOblate,  rr, cc, relError, error, aniso      = AR.calculateAnisoVar(densPS,center,aniso,viz=args['viz'])
                averVar = np.mean([varProlate, varOblate])
                stdVar  = abs(varProlate - varOblate) /2.

                print("Dispersion: Reconstructing isotropic dispersion.")
                # Calculate the isotropic variance for comparison.
                R, var3DIso                                                 = PSF.calculateIsoVar(densPS,k,varProj)

                # Write to a shell_out file
                g.write("Dispersion statistics\n")
                g.write("-----------------------------------------------------\n")
                g.write("Dispersion: Relative residual power: {}\n".format(relError))
                g.write("Dispersion: Absolute residual power correction: {}\n".format(error))
                g.write("Dispersion: The anisotropy factor is: {}\n".format(aniso))
                g.write("Dispersion: The 2D projection variance is: {}\n".format(varProj))
                g.write("Dispersion: The 3D variance is: {}\n".format(var3D))
                g.write("Dispersion: The variance from the 3D prolate transform: {}\n".format(varProlate))
                g.write("Dispersion: The variance from the 3D oblate transform: {}\n".format(varOblate))
                g.write("Dispersion: The variance from the isotropic recon: {}\n".format(var3DIso))
                g.write("------------------------------------------------- \n")

                ################################################################
                # Skewness Statistics: < X^3 > / < X^2 > ^ (3 / 2)
                ################################################################

                print("Skewness: Calculating 2D power spectra")
                # Construct the < X^3 > term
                densPS_skew, k    = PSF.PowerSpectrum( densProj ** (3./2.) )
                logDensPS_skew    = np.log10(densPS_skew)

                # Define some visualisation parameters
                sigma       = 10                        # for Guassian kernel for contour plots
                isobars     = np.array([0,-20,40])      # define the domain for the isobars on contour plots in log10 [to, from, by]

                # Plot the column density and the power spectrum
                if args['viz'] == True:
                    print("Skewness: Creating a denisty and power spectrum plot.")
                    PF.kAndDensityFieldPlot(densProj,densPS_skew,logDensPS_skew,sigma,isobars,fileName)

                print("Skewness: Fitting ellipses.")
                # Fit and Calculate the ellipses to all of the k-space power spectra
                aniso_skew, anisoStd_skew, prinAxisKeep, center, kperpAxis, kparAxis  = EF.calculateEllipses(logDensPS_skew,sigma,isobars,"kperp")

                print("Skewness: Reconstructing anisotropic dispersion.")
                # Calculate the anisotropic variances for the prolate and oblate ellipse.
                skewProlate, skewOblate,  rr, cc, relError_skew, error_skew, aniso_skew    = AR.calculateAnisoVar(densPS_skew,center,aniso_skew,viz=args['viz'])

                print("Skewness: Calculating skewness")
                # constructing the skewness for each of the transforms
                skewProlate = skewProlate / varProlate**(3./2.)
                skewOblate  = skewOblate / varOblate**(3./2.)

                # the average skewnness
                averSkew = np.mean([skewProlate, skewOblate])
                stdSkew  = abs(skewProlate - skewOblate) /2.

                # Write to a shell_out file
                g.write("Skewness statistics\n")
                g.write("-----------------------------------------------------\n")
                g.write("Skewness: Relative residual power: {}\n".format(relError_skew))
                g.write("Skewness: Absolute residual power correction: {}\n".format(error_skew))
                g.write("Skewness: The anisotropy factor is: {}\n".format(aniso_skew))
                g.write("Skewness: The 3D skewness is: {}\n".format(skew3D))
                g.write("Skewness: The skewness from the 3D prolate transform: {}\n".format(skewProlate))
                g.write("Skewness: The skewness from the 3D oblate transform: {}\n".format(skewOblate))
                g.write("------------------------------------------------- \n")

                ################################################################
                # Kurtosis Statistics: < X^4 > / < X^2 > ^ 2
                ################################################################

                print("Kurtosis: Calculating 2D power spectra")
                # Construct the < X^4 > term
                densPS_kurt, k     = PSF.PowerSpectrum( densProj ** 2. )
                logDensPS_kurt     = np.log10(densPS_kurt)

                # Define some visualisation parameters
                sigma       = 10                        # for Guassian kernel for contour plots
                isobars     = np.array([0,-20,40])      # define the domain for the isobars on contour plots in log10 [to, from, by]

                # Plot the column density and the power spectrum
                if args['viz'] == True:
                    print("Kurtosis: Creating a denisty and power spectrum plot.")
                    PF.kAndDensityFieldPlot(densProj,densPS_kurt,logDensPS_kurt,sigma,isobars,fileName)

                print("Kurtosis: Fitting ellipses.")
                # Fit and Calculate the ellipses to all of the k-space power spectra
                aniso_kurt, anisoStd_kurt, prinAxisKeep, center, kperpAxis, kparAxis  = EF.calculateEllipses(logDensPS_kurt,sigma,isobars,"kperp")

                print("Kurtosis: Reconstructing anisotropic dispersion.")
                # Calculate the anisotropic variances for the prolate and oblate ellipse.
                kurtProlate, kurtOblate,  rr, cc, relError_kurt, error_kurt, aniso_kurt    = AR.calculateAnisoVar(densPS_kurt,center,aniso_kurt,viz=args['viz'])

                print("Kurtosis: Calculating skewness")
                #
                kurtProlate = kurtProlate / varProlate ** 2.
                kurtOblate  = kurtOblate / varOblate ** 2.

                averKurt = np.mean([kurtProlate, kurtOblate])
                stdKurt  = abs(kurtProlate - kurtOblate) / 2.

                # Write to a shell_out file
                g.write("Kurtosis statistics\n")
                g.write("-----------------------------------------------------\n")
                g.write("Kurtosis: Relative residual power: {}\n".format(relError_kurt))
                g.write("Kurtosis: Absolute residual power correction: {}\n".format(error_kurt))
                g.write("Kurtosis: The anisotropy factor is: {}\n".format(aniso_kurt))
                g.write("Kurtosis: The 3D kurtosis is: {}\n".format(kurt3D))
                g.write("Kurtosis: The kurtosis from the 3D prolate: {}\n".format(kurtProlate))
                g.write("Kurtosis: The kurtosis from the 3D oblate: {}\n".format(kurtOblate))
                g.write("------------------------------------------------- \n")

                ################################################################
                # Printing Statistics
                ################################################################

                # Print off the values for the reconstructions
                if args['print'] == True:
                    print("-------------------------------------------------")
                    print("Dispersion statistics")
                    print("------------------------------------------------- \n")
                    print("Relative residual power: {}".format(relError))
                    print("Absolute residual power correction: {}".format(error))
                    print("The anisotropy factor is: {}".format(aniso))
                    print("The 2D projection variance is: {}".format(varProj))
                    print("The 3D variance is: {}".format(var3D))
                    print("The variance from the 3D prolate: {}".format(varProlate))
                    print("The variance from the 3D oblate: {}".format(varOblate))
                    print("The variance from the isotropic recon: {}".format(var3DIso))
                    print("------------------------------------------------- \n\n")

                    print("-------------------------------------------------")
                    print("Skewness statistics")
                    print("------------------------------------------------- \n")
                    print("Relative residual power: {}".format(relError_skew))
                    print("Absolute residual power correction: {}".format(error_skew))
                    print("The anisotropy factor is: {}".format(aniso_skew))
                    print("The 3D skewness is: {}".format(skew3D))
                    print("The skewness from the 3D prolate: {}".format(skewProlate))
                    print("The skewness from the 3D oblate: {}".format(skewOblate))
                    print("------------------------------------------------- \n\n")

                    print("-------------------------------------------------")
                    print("Kurtosis statistics")
                    print("------------------------------------------------- \n")
                    print("Relative residual power: {}".format(relError_kurt))
                    print("Absolute residual power correction: {}".format(error_kurt))
                    print("The anisotropy factor is: {}".format(aniso_kurt))
                    print("The 3D kurtosis is: {}".format(kurt3D))
                    print("The kurtosis from the 3D prolate: {}".format(kurtProlate))
                    print("The kurtosis from the 3D oblate: {}".format(kurtOblate))
                    print("------------------------------------------------- \n")

                if fileCount == 1:
                    f.write("dir_01, fileNumber_02, relError_03, error_04, "+
                            "varProj_05, var3D_06, varProlate_07, varOblate_08, averVar_09, stdVar_10, var3DIso_11, varAniso_12, "+
                            "relErrorSkew_13, errorSkew_14, skew3D_15, skewProlate_16, skewOblate_17, averSkew_18, stdSkew_19, skewAniso_20, "+
                            "relErrorKurt_21, errorKurt_22, kurt3D_23, kurtProlate_24, kurtOblate_25, averKurt_26, stdKurt_27, kurtAniso_28 \n")

                f.write("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, ".format(dir,fileNumber,relError,error,varProj,var3D,varProlate,varOblate,averVar,stdVar,var3DIso,aniso) +
                        "{}, {}, {}, {}, {}, {}, {}, {}, ".format(relError_skew,error_skew,skew3D,skewProlate,skewOblate,averSkew,stdSkew,aniso_skew) +
                        "{}, {}, {}, {}, {}, {}, {}, {} \n".format(relError_kurt,error_kurt,kurt3D,kurtProlate,kurtOblate,averKurt,stdKurt,aniso_kurt))
            except:
                print("Skipping number: {}".format(fileNumber))
                print("------------------------------------------------- \n")
                continue

    g.close()
