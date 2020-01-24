from header import *

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
