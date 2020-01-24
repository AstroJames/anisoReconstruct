from header import *

def PowerSpectrum(data,type=None,variance=None):
    """
    Calculate the power spectrum.

    INPUTS:
    ----------
    data        - The 2D image data
    n           - the coef. of the fourier transform = grid size.
    type        - either '2D' or 'aziAverage' for 2D or averaged power spectrum calculations.
    variance    - passed to the azimuthal averaging if you want to calculate the variance of
                the average.

    OUTPUTS:
    ----------
    Pspec       - the 2D power spectrum
    k           - the 2D k vector array
    """

    data    = data.astype(float)                # make sure the data is float type

    # This removes the k = 0 wavevalue and makes integrating easier.
    data    = data/data.mean() - 1             # Pre-processing following Federrath et al. 2016

    ft      = (1./ (data.shape[0]*data.shape[1] ) )*fftpack.fft2(data)                # 2D fourier transform
    ft_c    = fftpack.fftshift(ft)              # center the transform so k = (0,0) is in the center
    PSpec   = np.abs(ft_c*np.conjugate(ft_c))   # take the power spectrum

    # Take the azimuthal average of the powr spectrum, if required.
    if type == 'aziAverage':

        if not variance:

            PSpec  = azimuthalAverage(PSpec)

            return data
        else:

            PSpec, var  = azimuthalAverage(PSpec,variance=True)

            return PSpec, var

    # Create the kx and ky vector components.
    kx      = np.round(np.linspace( -( PSpec.shape[0] + 1 )/2, ( PSpec.shape[0] + 1 )/2, PSpec.shape[0]))
    ky      = np.round(np.linspace( -( PSpec.shape[1] + 1 )/2, ( PSpec.shape[1] + 1 )/2, PSpec.shape[1]))
    kx, ky  = np.meshgrid(kx,ky,indexing="xy")
    k       = np.hypot(kx,ky)

    return PSpec, k

def PowerSpectrumAveraging(files,densOrderd,run):
    """
    this functions averages over the power spectrums and returns a dictionary with the averaged
    power spectrums in it.

    INPUTS:
    ----------
    files           - all of the file names for each simulations
    densOrdered     - the density dictionary ordered by plot order
    run             - if the function needs to be rerun for recompiling of the density plots

    OUTPUTS:
    ----------
    PSpecAver       - the average power spectrum as a dictionary, for each of the simulations
    PSpecVar        - the variance power spectrum as a dictionary, for each of the simulations

    """

    if run == 0: # if the power spectrum need to be recompiled.
        PSpecAver   = {}
        PSpecVar    = {}
        fileCounter = 0
        # Average the power spectrum, from 5T to 10T
        for iter in xrange(50,101):
            print("Currently on iteration {}".format(iter))
            # Load the density files
            try:
                density     = LoadPickles(files,iter)
            except IndexError:
                print "Index error, I'm going to break the loop."
                break
            plotCounter = 0
            for i in xrange(0,5):
                for j in xrange(0,4):
                    if fileCounter == 0:
                        dens        = density[densOrderd[plotCounter]]          # column density
                        PSpec, k    = PowerSpectrum(dens)                       # the power spectrum and wavevector
                        PSpecAver[densOrderd[plotCounter]] = PSpec              # add a new power spectrum to the dictionary
                        PSpecVar[densOrderd[plotCounter]]  = PSpec**2           # for constructing the variance
                    else:
                        dens        = density[densOrderd[plotCounter]]          # column density
                        PSpec, k    = PowerSpectrum(dens)                       # the power spectrum and wavevector
                        PSpecAver[densOrderd[plotCounter]]  += PSpec            # add the power spectrum together
                        PSpecVar[densOrderd[plotCounter]]   += PSpec**2         # for constructing the variance
                    plotCounter +=1 #update the plot
            fileCounter +=1 #update the file
        # Average the power spectrum and take the log10 transform
        for key in PSpecAver.keys():
            PSpecAver[key]  = PSpecAver[key]/fileCounter
            PSpecVar[key]   = (PSpecVar[key]/fileCounter - PSpecAver[key]**2)**0.5
        save_obj(PSpecAver,"AveragePSpec")
        save_obj(PSpecVar,"StdPSpec")
    else:
        PSpecAver   = load_obj("AveragePSpec.pkl")
        PSpecVar    = load_obj("StdPSpec.pkl")
    return PSpecAver, PSpecVar

def calculateIsoVar(PowerSpectrum,k,var2D):
    """
    Assuming isotropy of the kz, this function calculates R = sigma^2_2 / sigma^2_3

    INPUTS:
    ------------------------------------------------------------------------------------------
    PowerSpectrum   - the 2D power spectrum.
    k               - the k wavevector as a 2D grid.
    var2D           - the varianace of the 2D column density.

    OUTPUTS:
    ------------------------------------------------------------------------------------------
    R       - the ratio between the 2D and 3D variance
    var3D   - the estimated 3D variance

    """

    # Differentials for integration
    dkx = k[0,0]-k[0,1]
    dky = k[0,0]-k[1,0]
    dk  = np.hypot(dkx,dky)

    # Calculate the integrals over the 2D and 3D power spectrum, assuming isotropy
    P2D = 2* np.pi* sum( sum( PowerSpectrum ) ) * dk
    P3D = 4* np.pi* sum( sum( PowerSpectrum * k ) ) * dk

    # Calculate R from Brunt et al. 2010, and the 3D variance.
    R       = P2D / P3D
    var3D   = var2D / R

    return R, var3D
