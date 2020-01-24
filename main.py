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
from scipy.integrate import quad
from PickleData import *

########################################################################################################################

ap 			= argparse.ArgumentParser(description = 'Just a bunch of input arguments')
ap.add_argument('-type','--type',default="iso",help='aniso, iso',type=str)
args 		= vars(ap.parse_args())

def extractEllipseCircum(kperp,kpar,aniso):
    """
    Calculate the circumference of the ellipse in k-space.

    INPUTS:
    ----------
    kperp   - the modulus of the principle axis for the ellipse perpendicular to the B-field
    kpar    - the modulus of the principle axis for the ellipse parallel to the B-field
    aniso   - the anisotropy factor: kpar / kperp


    OUTPUTS:
    ----------
    circ    - the circumference of the ellipse

    """

    if aniso > 1.:
        print("Swapping axis for circumference")
        aniso = 1. / aniso

    # Define the eccentricity of the ellipse
    e = np.sqrt( 1 - aniso**2 )

    # Circumference function
    f       = lambda theta, e: np.sqrt( 1 - e**2 * np.sin(theta)**2  )

    # Using quadrature integration
    F, err  = quad(f,0,np.pi/2.,args=(e,))

    if err > 1e-5:
        raise Exception("Circumference integration has failed.")

    # Circumference = 4 * a * int_0^pi/2 f d\theta
    circ    = 4 * kperp * F

    return circ

def extractProlateEllipse(kperp,kpar,aniso):
    """
    Calculate the surface area of a prolate ellipse in k-space

    INPUTS:
    ----------
    kperp   - the modulus of the principle axis for the ellipse perpendicular to the B-field
    kpar    - the modulus of the principle axis for the ellipse parallel to the B-field
    aniso   - the anisotropy factor: kpar / kperp


    OUTPUTS:
    ----------
    surface    - the prolate surface area of an ellipse

    """

    if aniso > 1.:
        print("Swapping axis for prolate ellipse")
        aniso = 1. / aniso
        # kperpDummy  = kperp
        # kperp       = kpar
        # kpar        = kperpDummy


    # Define the eccentricity of the ellipse
    e = np.sqrt( 1. - aniso**2 )


    # the prolate surface area
    surface = 2. * np.pi * kpar**2. * ( 1. + (1. / ( e * np.sqrt(  1. - e**2. ) ) ) * np.arcsin(e) )

    return surface


def extractOblateEllipse(kperp,kpar,aniso):
    """
    Calculate the oblate surface area of the ellipse in k-space.

    INPUTS:
    ----------
    kperp   - the modulus of the principle axis for the ellipse perpendicular to the B-field
    kpar    - the modulus of the principle axis for the ellipse parallel to the B-field
    aniso   - the anisotropy factor: kpar / kperp


    OUTPUTS:
    ----------
    surface     - the oblate surface area of an ellipse

    """

    if aniso > 1.:
        print("Swapping axis for oblate ellipse")
        aniso = 1. / aniso
        # kperpDummy  = kperp
        # kperp       = kpar
        # kpar        = kperpDummy

    # Define the eccentricity of the ellipse
    e = np.sqrt( 1. - aniso**2. )


    # the oblate surface area
    surface  = 2. * np.pi * kperp**2. * ( 1. + ( (1. - e**2.) / e ) * np.arctanh(e) )

    return surface

def createAnisotropicK(powerSpectrum,center,aniso):
    """
    Create the anisotropic k space vectors based on the average anisotropy

    INPUTS:
    ----------
    powerSpectrum   - the 2D power spectrum ( untransformed )
    aniso           - the averaged anistropy factor


    OUTPUTS:
    ----------
    kperp       -
    kpar        -

    """

    if aniso > 1.:
        anisoNew = 1. / aniso
        padDim      = int( np.round( powerSpectrum.shape[0] / ( anisoNew ) ) )
    else:
        padDim      = int( np.round( powerSpectrum.shape[0] / ( aniso) ) )

    # the amount to pad the power spectrum

    padAmount   = padDim - powerSpectrum.shape[0]/2.
    paddedPS    = np.pad(powerSpectrum, (padAmount, padAmount), 'constant', constant_values=(0, 0))
    center      = np.array(center) + padAmount

    # the kperp and kpar components of the ellipses
    if aniso < 1.:
        kperp       = np.arange(1,padDim-1)
        kpar        = np.floor(aniso * kperp).astype(int)
    else:
        kpar      = np.arange(1,padDim-1)
        kperp     = np.floor(anisoNew * kpar).astype(int)


    return paddedPS, kpar, kperp, center, padAmount

def calculateAnisoVar(powerSpectrum,center,aniso,viz):
    """
    Calculate the 3D variance from the anisotropic 2D-projected power spectrum

    INPUTS:
    ----------
    powerSpectrum   - the 2D power specutrm (untransformed)
    center          - the center of the 2D power spectrum for reference
    ansio           - the anisotropy array
    viz             - a parameter for visualisation

    OUTPUTS:
    ----------
    varProlate  - the variance contructed from rotating the 2D power spectrum into a
                - 3D prolate ellipse
    varOblate   - the vairance constructed from rotating the 2D power spectrum into a
                - 3D oblate ellipse

    """

    # Take the anisotropy factor and big the most isotropic anisotropy.
    # the reason is because the most isotropic modes hold the most power
    # especially in the high Mach regime where anistropic shocks are on smaller
    # k scales than the driving.
    print("Picked the anistropy for the k-modes with the most power.")
    if np.max(aniso) > 1.0:
        aniso = np.min(aniso)
    else:
        aniso = np.max(aniso)


    # Note, this is what will change as a function of scale in a later implementation
    paddedPS, kpars, kperps, center, padAmount = createAnisotropicK(powerSpectrum,center,anisoMean)
    # Plot for checking the elliptic fits


    counter             = 0
    prolateVolume       = []
    oblateVolume        = []
    totalProlatePower   = []
    totalOblatePower    = []

    # make a copy of the power spectrum for visualisations
    if viz:
        modifiedPS = paddedPS.copy()

    for kpar, kperp in zip(kpars, kperps):
        print("fitting an ellipse with kperp scale: {}".format(kperp))

        # get the ellipse coordinates on the power spectrum
        rr, cc = skimage.draw.ellipse_perimeter(int(center[0]), int(center[1]), kpar, kperp)

        if(kperp == 0 or kpar == 0):
            # make sure not to include a 0 in either of the coordinates
            # which could happen because of the floor function in createAnisotropicK
            continue

        # Create a new function that compares the new and old rr and cc
        if counter > 0:
            # Make sure that for the first 50 k vectors that they don't contain the
            # same coordinates (this is very expensive for large k, but also there is very
            # little power in large k)
            if(kperp < 50):
                rr, cc, toDelete    = compareIndexes(rr,cc,rrOld,ccOld)

        powerAtKparProj         = np.sum(paddedPS[rr,cc])

        if viz:
            if np.mod(counter,1) == 0:
                modifiedPS[rr,cc] = 1

        # Calculate the power in each of the rotations of the power spectra
        prolateVolumeFactor     = extractProlateEllipse(kperp,kpar,anisoMean) / extractEllipseCircum(kperp,kpar,aniso)
        powerAtProlateKpar3D    = powerAtKparProj * prolateVolumeFactor
        oblateVolumeFactor      = extractOblateEllipse(kperp,kpar,anisoMean) / extractEllipseCircum(kperp,kpar,aniso)
        powerAtOblateKpar3D     = powerAtKparProj * oblateVolumeFactor

        # Append an arrays of powers and volume factors
        prolateVolume.append(prolateVolumeFactor)
        oblateVolume.append(oblateVolumeFactor)
        totalProlatePower.append(powerAtProlateKpar3D)
        totalOblatePower.append(powerAtOblateKpar3D)
        counter += 1

        # Store the old indices for the Ellipse
        # because these will need to be omitted
        rrOld = rr
        ccOld = cc

    error      = sum(sum((modifiedPS!=1)*paddedPS))
    relError   =  sum(sum((modifiedPS!=1)*paddedPS)) / sum(sum(paddedPS))

    varProlate  = np.sum(np.array(totalProlatePower)) + relError*np.sum(np.array(totalProlatePower))
    varOblate   = np.sum(np.array(totalOblatePower)) + relError*np.sum(np.array(totalProlatePower))

    print("relative residual power: {} absolute residual power: {}".format(relError,error))

    if viz:
        regionPS = modifiedPS[padAmount:(modifiedPS.shape[0]-padAmount),padAmount:(modifiedPS.shape[1]-padAmount)]
        f, ax = plt.subplots(1,1,dpi=200)
        ax.imshow(regionPS,extent=[-256.5, 256.5, -256.5, 256.5],vmin=10**-15,norm=mpl.colors.LogNorm(),cmap=plt.cm.plasma)
        ax.set_ylabel(r"$k_{\parallel}$",fontsize=fs)
        ax.set_xlabel(r"$k_{\perp}$",fontsize=fs)
        plt.show()

    return varProlate, varOblate, rr, cc


def compareIndexes(rrs,ccs,rrOlds,ccOlds):
    """
    A function that removes power spectrum coordiantes that are the
    same between successive elliptic fits.

    """

    # create a better way of doing this
    #newArr = np.vstack([rr,cc])
    #oldArr = np.vstacl([rrOld,ccOld])

    iterCount   = 0
    toDelete    = []

    for rr, cc in zip(rrs,ccs):
        for rrOld, ccOld in zip(rrOlds,ccOlds):
            if rr == rrOld and cc == ccOld:
                toDelete.append(iterCount)
                break

        iterCount += 1


    rrMod = np.delete(rrs,toDelete)
    ccMod = np.delete(ccs,toDelete)

    return rrMod, ccMod, toDelete

def rotateByPixel():

    pass

########################################################################################################################

if __name__ == "__main__":

    if args['type'] == "aniso":

        MachData    = "M20MA0.1"
        iter        = "50"
        file3DDir   = "./testData/"+ MachData +"/3DDens/"
        fileProjDir = "./testData/"+ MachData +"/Proj/"

        # 3D Variance
        dens3D      = load_obj(file3DDir + "Turb_hdf5_plt_cnt_00{}_{}_dens".format(iter,MachData))
        var3D       = (dens3D / dens3D.mean() - 1).var()
        del dens3D

        # Projections and Power Spectrum
        densProj    = load_obj(fileProjDir + "Turb_hdf5_plt_cnt_00{}_{}_proj".format(iter,MachData))
        varProj     = (densProj / densProj.mean() - 1).var()
        densPS, k   = PSF.PowerSpectrum(densProj)
        logDensPS   = np.log10(densPS)


        sigma       = 10                       # for the Gaussin kernel
        fs          = 14                       # global fontsize
        isobars     = np.array([0,-20,40])     # define the domain for the isobars


        f, ax = plt.subplots(2,1,dpi=300)
        plt.subplots_adjust(left=0.01, bottom=0.05, right=0.9, top=0.98, wspace=0.05, hspace=0.03)
        xi = densProj/densProj.mean() -1
        plot1 = ax[0].imshow(xi,cmap=plt.cm.plasma,vmin=xi.min(), vmax=xi.max())
        PF.annotatePlot(ax[0],2,10,xi.max(),xi.min())
        cbar1 = plt.colorbar(plot1,pad=0.01,ax=ax[0])
        cbar1.set_label(r"$\Sigma/\Sigma_0 - 1$",fontsize=fs)
        ax[0].set_xticks([])
        ax[0].set_yticks([])
        plot2 = ax[1].imshow(densPS,cmap=plt.cm.plasma,extent=[-256.5, 256.5, -256.5, 256.5],vmin=10**-15,norm=mpl.colors.LogNorm())
        ax[1].set_ylabel(r"$k_{\parallel}$",fontsize=fs)
        ax[1].set_xlabel(r"$k_{\perp}$",fontsize=fs)
        cbar2 = plt.colorbar(plot2,pad=0.01,ax=ax[1])
        cbar2.set_label(r"$\mathscr{P}_{\Sigma/\Sigma_0 - 1}$",fontsize=fs)
        aniso, anisoStd, prinAxisKeep, center, kperpAxis, kparAxis = EF.EllipsePlotter(ax[1],logDensPS,sigma,isobars,"kperp")
        EF.ContourPlot(ax[1],logDensPS,sigma,isobars,center)
        #plt.tight_layout()
        plt.savefig("Figure1_M2Ma01.png",bbox_inches='tight',dpi=300)
        plt.close()


        varProlate, varOblate,  rr, cc,   = calculateAnisoVar(densPS,center,aniso,viz=True)
        R, var3DIso                       = PSF.calculateIsoVar(densPS,k,varProj)
        print("The 2D projection variance is: {}".format(varProj))
        print("The 3D variance is: {}".format(var3D))
        print("The variance from the 3D prolate: {}".format(varProlate))
        print("The variance from the 3D oblate: {}".format(varOblate))
        print("The variance from the isotropic recon: {}".format(var3DIso))

    if args['type'] == "iso":

        def powerSpec(dens):
            densPS, k   = PSF.PowerSpectrum(dens)
            logDensPS   = np.log10(densPS)
            return logDensPS

        isoDir      = "./testData/isotropic/"
        M2iso       = powerSpec(load_obj(isoDir + "Turb_hdf5_plt_cnt_0050_M2MA10_proj"))
        M4iso       = powerSpec(load_obj(isoDir + "Turb_hdf5_plt_cnt_0050_M4MA10_proj"))
        M10iso      = powerSpec(load_obj(isoDir + "Turb_hdf5_plt_cnt_0050_M10MA10_proj"))
        M20iso      = powerSpec(load_obj(isoDir + "Turb_hdf5_plt_cnt_0050_M20MA10_proj"))

        sigma       = 10                       # for the Gaussin kernel
        fs          = 16                       # global fontsize
        isobars     = np.array([0,-20,40])     # define the domain for the isobars

        f, ax = plt.subplots(2,2,dpi=200,sharex=True,sharey=True)
        ax[0,0].imshow(M2iso,cmap=plt.cm.plasma,extent=[-256.5, 256.5, -256.5, 256.5],vmin=-15,vmax=0)
        aniso, anisoStd, prinAxisKeep, center, kperpAxis, kparAxis = EF.EllipsePlotter(ax[0,0],M2iso,sigma,isobars,"kperp")
        EF.ContourPlot(ax[0,0],M2iso,sigma,isobars,center)

        ax[0,1].imshow(M4iso,cmap=plt.cm.plasma,extent=[-256.5, 256.5, -256.5, 256.5],vmin=-15,vmax=0)
        aniso, anisoStd, prinAxisKeep, center, kperpAxis, kparAxis = EF.EllipsePlotter(ax[0,1],M4iso,sigma,isobars,"kperp")
        EF.ContourPlot(ax[0,1],M4iso,sigma,isobars,center)

        ax[1,0].imshow(M10iso,cmap=plt.cm.plasma,extent=[-256.5, 256.5, -256.5, 256.5],vmin=-15,vmax=0)
        aniso, anisoStd, prinAxisKeep, center, kperpAxis, kparAxis = EF.EllipsePlotter(ax[1,0],M10iso,sigma,isobars,"kperp")
        EF.ContourPlot(ax[1,0],M10iso,sigma,isobars,center)

        ax[1,1].imshow(M20iso,cmap=plt.cm.plasma,extent=[-256.5, 256.5, -256.5, 256.5],vmin=-15,vmax=0)
        aniso, anisoStd, prinAxisKeep, center, kperpAxis, kparAxis = EF.EllipsePlotter(ax[1,1],M20iso,sigma,isobars,"kperp")
        EF.ContourPlot(ax[1,1],M20iso,sigma,isobars,center)




        #ax.set_ylabel(r"$k_{\parallel}$",fontsize=fs)
        #ax.set_xlabel(r"$k_{\perp}$",fontsize=fs)
        #cbar = plt.colorbar(plot,pad=0.01)
        #cbar.set_label(r"$\mathscr{P}(k_{\perp},k_{\parallel})$",fontsize=fs)
        plt.show()

    if args['type'] == "viz":
        import matplotlib.path as mpath
        def colorBarPos(ax,loc):
            axIn = inset_axes(ax,
                        width="50%",  # width = 50% of parent_bbox width
                        height="5%",  # height : 5%
                        loc=loc)
            return axIn

        circle = mpath.Path.unit_circle()
        verts = np.copy(circle.vertices)
        verts[:, 0] *= 2
        ellipMarker2 = mpath.Path(verts, circle.codes)

        verts = np.copy(circle.vertices)
        verts[:, 1] *= 2
        ellipMarker1 = mpath.Path(verts, circle.codes)


        Mach    = np.array([2.6,2.2,2.0,1.66,1.8,5.2,4.4,3.8,3.5,3.7])
        MachA   = np.array([0.133,0.54,0.98,1.7,9.2,0.13,0.54,0.95,1.73,9.2])

        realVar = np.array([0.543,0.574,0.512,0.674,0.603,1.575,1.749,1.997,2.260,2.346])
        proVar  = np.array([0.295,0.318,0.369,0.494,0.558,0.871,1.263,1.514,1.802,1.589])
        oblVar  = np.array([0.548,0.626,0.568,0.743,0.687,1.633,1.820,2.068,2.483,2.234])
        isoVar  = np.array([0.379,0.464,0.373,0.585,0.634,1.076,1.450,1.850,2.094,1.956])

        linear = lambda x: x
        xdomain = np.linspace(0,10,1000)

        fig, ax = plt.subplots(1,2,dpi=200)
        # ax.scatter(Mach,realVar,marker='X',c=MachA,norm=colors.LogNorm(vmin=min(MachA),vmax=max(MachA)),
        #            cmap=plt.cm.plasma,edgecolors="black",linewidths=0.5,zorder=10,s=100,label=r"True $\sigma_s^2$")
        ax[0].scatter(realVar[0:4],isoVar[0:4],marker='o',c=MachA[0:4],norm=colors.LogNorm(vmin=min(MachA),vmax=max(MachA)),
                   cmap=plt.cm.plasma,edgecolors="black",linewidths=0.5,zorder=10,s=100,label="Brunt+2010 estimate")
        ax[0].scatter(realVar[0:4],proVar[0:4],marker=ellipMarker1,c=MachA[0:4],norm=colors.LogNorm(vmin=min(MachA),vmax=max(MachA)),
                   cmap=plt.cm.plasma,edgecolors="black",linewidths=0.5,zorder=10,s=100,label="prolate estimate")
        ax[0].scatter(realVar[0:4],oblVar[0:4],marker=ellipMarker2,c=MachA[0:4],norm=colors.LogNorm(vmin=min(MachA),vmax=max(MachA)),
                   cmap=plt.cm.plasma,edgecolors="black",linewidths=0.5,zorder=10,s=100,label="oblate estimate")
        ax[0].plot(xdomain,linear(xdomain),color="red",ls='--')
        ax[0].set_xlim(0.2,0.8)
        ax[0].set_ylim(0.2,0.8)
        ax[0].set_xlabel(r"$\sigma^2_{\text{true}}$",fontsize=14)
        ax[0].set_ylabel(r"$\sigma^2_{\text{estimated}}$",fontsize=14)

        p = ax[1].scatter(realVar[5:9],isoVar[5:9],marker='o',c=MachA[5:9],norm=colors.LogNorm(vmin=min(MachA),vmax=max(MachA)),
                   cmap=plt.cm.plasma,edgecolors="black",linewidths=0.5,zorder=10,s=100,label="Brunt+2010 estimate")
        ax[1].scatter(realVar[5:9],proVar[5:9],marker=ellipMarker1,c=MachA[5:9],norm=colors.LogNorm(vmin=min(MachA),vmax=max(MachA)),
                   cmap=plt.cm.plasma,edgecolors="black",linewidths=0.5,zorder=10,s=100,label="prolate estimate")
        ax[1].scatter(realVar[5:9],oblVar[5:9],marker=ellipMarker2,c=MachA[5:9],norm=colors.LogNorm(vmin=min(MachA),vmax=max(MachA)),
                   cmap=plt.cm.plasma,edgecolors="black",linewidths=0.5,zorder=10,s=100,label="oblate estimate")
        ax[1].plot(xdomain,linear(xdomain),color="red",ls='--')
        ax[1].set_xlim(0.8,3)
        ax[1].set_ylim(0.8,3)
        ax[1].set_xlabel(r"$\sigma^2_{\text{true}}$",fontsize=14)
        ax[1].set_ylabel(r"$\sigma^2_{\text{estimated}}$",fontsize=14)
        cax     = colorBarPos(ax[1],'upper right')
        cbar2   = plt.colorbar(p, orientation='horizontal', cax=cax)
        cbar2.set_label(r"$\mathcal{M}_{\text{A}0}$", labelpad=-28,fontsize=16,x=-0.2)

        plt.show()