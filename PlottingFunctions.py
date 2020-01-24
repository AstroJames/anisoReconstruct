from header import *
import EllipseFunctions as EF
reload(EF)
import PowerSpectrumFunctions as PSF
reload(PSF)

def annotatePlot(ax,M,Ma,sMax,sMin,type=None,logType=None):
    """


    """

    # Domain for the Mach label
    x_coord_M = 0.078
    y_coord_M = 0.9

    # Domain for the Mach Alven label
    x_coord_Ma = 0.65
    y_coord_Ma = 0.895

    # Domain for the turnover time label
    x_coord_T = 0.04
    y_coord_T = 0.8

    # Domain for the scale label
    x_coord_D = 0.078
    y_coord_D = 0.05

    # the shift in the shadow
    eps = 0.002

    sMax = np.round(sMax,1)
    sMin = np.round(sMin,1)

    fs = 10; # font size

    if type is not None:
        if type == "proj":
            character = "\Sigma"
        elif type == "slice":
            character = r"\rho"

    if logType is not None:
        if logType == "natural":
            logChar = r"\ln"
        elif logType == "base10":
            logChar = r"\log_{10}"

    # Turbulent Mach number
    ax.annotate(r'$\mathcal{M} \approx $' + ' {}'.format(M),
                     xy=(x_coord_M-eps, y_coord_M+eps),xycoords='axes fraction',
                     xytext=(x_coord_M-eps, y_coord_M+eps),textcoords='axes fraction',
                     fontsize=fs,color='black')
    ax.annotate(r'$\mathcal{M} \approx $' + ' {}'.format(M),
                     xy=(x_coord_M, y_coord_M),xycoords='axes fraction',
                     xytext=(x_coord_M, y_coord_M),textcoords='axes fraction',
                     fontsize=fs,color='w')

    # Turbulent Mach Alven number
    ax.annotate(r'$\mathcal{M}_{\text{A}0} \approx $' + ' {}'.format(Ma),
                     xy=(x_coord_Ma-eps, y_coord_Ma+eps),xycoords='axes fraction',
                     xytext=(x_coord_Ma-eps, y_coord_Ma+eps),textcoords='axes fraction',
                     fontsize=fs,color='black')
    ax.annotate(r'$\mathcal{M}_{\text{A}0} \approx $' + ' {}'.format(Ma),
                     xy=(x_coord_Ma, y_coord_Ma),xycoords='axes fraction',
                     xytext=(x_coord_Ma, y_coord_Ma),textcoords='axes fraction',
                     fontsize=fs,color='w')

    # Density Domain
    if type is not None and logType is not None:
        ax.annotate(r'${}$'.format(logChar) + r'$\left({}/{}_0\right)$:'.format(character,character) +r' $\left[{}, {}\right]$'.format(sMin,sMax),
                         xy=(x_coord_D-eps, y_coord_D+eps),xycoords='axes fraction',
                         xytext=(x_coord_D-eps, y_coord_D+eps),textcoords='axes fraction',
                         fontsize=fs,color='black')
        ax.annotate(r'${}$'.format(logChar) + r'$\left({}/{}_0\right)$:'.format(character,character) + r' $\left[{}, {}\right]$'.format(sMin,sMax),
                         xy=(x_coord_D, y_coord_D),xycoords='axes fraction',
                         xytext=(x_coord_D, y_coord_D),textcoords='axes fraction',
                         fontsize=fs,color='w')


# This is a plot for main.py that should be run if the
def kAndDensityFieldPlot(densProj,densPS,logDensPS,sigma,isobars,fileName):
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
    aniso, anisoStd, prinAxisKeep, center, kperpAxis, kparAxis = EF.calculateEllipses(logDensPS,sigma,isobars,"kperp",ax[1])
    EF.ContourPlot(ax[1],logDensPS,sigma,isobars,center)
    #plt.tight_layout()
    plt.savefig("{}.png".format(fileName),bbox_inches='tight',dpi=300)
    plt.close()
