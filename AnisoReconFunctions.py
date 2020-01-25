from header import *

import PlottingFunctions as PF
reload(PF)
import EllipseFunctions as EF
reload(EF)
import PowerSpectrumFunctions as PSF
reload(PSF)
from PickleData import *

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
        #print("Swapping axis for circumference")
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
        #print("Swapping axis for prolate ellipse")
        aniso = 1. / aniso

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
        #print("Swapping axis for oblate ellipse")
        aniso = 1. / aniso

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
        padDim      = int( np.round( powerSpectrum.shape[0] / ( aniso ) ) )

    # the amount to pad the power spectrum
    padAmount   = padDim - powerSpectrum.shape[0]/2
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
    #print("Picked the anistropy for the k-modes with the most power.")
    if np.max(aniso) > 1.:
        aniso = np.min(aniso)
    else:
        aniso = np.max(aniso)

    if aniso > 1.:
        aniso = 1./aniso

    if np.isreal(aniso) == False:
        aniso = np.abs(aniso)


    #print("aniso: {}".format(aniso))
    # Note, this is what will change as a function of scale in a later implementation
    paddedPS, kpars, kperps, center, padAmount = createAnisotropicK(powerSpectrum,center,aniso)
    # Plot for checking the elliptic fits


    counter             = 0
    prolateVolume       = []
    oblateVolume        = []
    totalProlatePower   = []
    totalOblatePower    = []

    # make a copy of the power spectrum for visualisations
    modifiedPS = paddedPS.copy()

    for kpar, kperp in zip(kpars, kperps):
        #print("fitting an ellipse with kperp scale: {}".format(kperp))

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

        # Keeping track of all of the modes that are hit by the ellipses in
        # 2d
        if np.mod(counter,1) == 0:
            modifiedPS[rr,cc] = 1

        # Calculate the power in each of the rotations of the power spectra
        prolateVolumeFactor     = extractProlateEllipse(kperp,kpar,aniso) / extractEllipseCircum(kperp,kpar,aniso)
        powerAtProlateKpar3D    = powerAtKparProj * prolateVolumeFactor
        oblateVolumeFactor      = extractOblateEllipse(kperp,kpar,aniso) / extractEllipseCircum(kperp,kpar,aniso)
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

    #print("Relative residual power: {} Absolute residual power: {}".format(relError,error))

    if viz:
        regionPS = modifiedPS[padAmount:(modifiedPS.shape[0]-padAmount),padAmount:(modifiedPS.shape[1]-padAmount)]
        f, ax = plt.subplots(1,1,dpi=200)
        ax.imshow(regionPS,extent=[-256.5, 256.5, -256.5, 256.5],vmin=10**-15,norm=mpl.colors.LogNorm(),cmap=plt.cm.plasma)
        ax.set_ylabel(r"$k_{\parallel}$",fontsize=fs)
        ax.set_xlabel(r"$k_{\perp}$",fontsize=fs)
        plt.show()

    return varProlate, varOblate, rr, cc, relError, error, aniso


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
