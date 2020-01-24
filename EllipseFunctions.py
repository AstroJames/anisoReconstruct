from header import *

def ContourPlot(ax,data,sigma,isobars,center):
    """
    Creates contours along the isobars of 2D data.

    INPTS:

    ----------
    ax      - the plot that the contours will be added to.
    data    - the 2D data that the contours will be created on.
    sigma   - the standard deviations of the Gaussin kernel for smoothing.
    isobars - the domain of the isobars.

    OUTPUTS:
    ----------
    contours added to a plot, on an axis, ax

    """

    # Error handling
    if len(isobars) != 3:
        raise Exception("Incorrect dimensions for the length of the isobar range.")

    # Read in the isobar domain
    bottom = isobars[0]
    top    = isobars[1]
    number = isobars[2]

    # A counter for making sure the plots only get labelled once.
    counter = 0

    # Smooth the 2D data, and create the contours on the plot
    for parms in np.linspace(bottom,top,number): #the contour domain, i.e. the power of each isobar
        contours = measure.find_contours( filters.gaussian(data, sigma), parms)

        for n, contour in enumerate(contours):
            if counter == 0:
                ax.plot(contour[:, 1] - center[0], contour[:, 0] - center[1], linewidth=1, color='red',label='Data',zorder=1)
                counter += 1
            else:
                ax.plot(contour[:, 1] - center[0], contour[:, 0] - center[1], linewidth=1, color='red',zorder=1)


def fitEllipse(x,y):
    """
    Fit an ellipse of the form D(x^2, xy, y^2, x, y, 1) * a(a_xx, a_yy, a_x, a_y, a_1) = 0
    to a set of data (x,y), using constrained least squares.

    INPUTS:
    ----------
    x - the x data
    y - the y data

    OUTPUTS:
    ----------
    a - the coefficient for the quadratic (elliptic) form

    """

    # Define the x and y data for the fitting
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]

    # Define the quadratic form of the ellipse D(x^2,xy,y^2,x,y,1)
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))

    # Define the S matrix, D^TD
    S = np.dot(D.T,D)

    # Define the constraint matrix, C
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1

    # Calculate the eigen values of S{^-1}C
    E, V =  eig(np.dot(inv(S), C))

    # Find the largest eigen value
    n = np.argmax(E)

    # a is the eigen vector correspdonding to the largest eigen value
    a = V[:,n]

    # return a(a_xx, a_xy, a_yy, a_x, a_y, a_1)
    return a

def findEllipseCenter(a):
    """
    This function takes in the coefficient vector, fit using constrained least squares of the
    ellipse and calculates the center coordinates.

    INPUT:
    ----------
    a - the coefficients for the ellipse quadratic form:
    D(x^2, xy, y^2, x, y, 1) * a(a_xx, a_yy, a_x, a_y, a_1) = 0

    OUTPUTS:
    ----------
    x0  - the centroid for the x coordinate
    y0  - the centroid for the y coordinate

    """

    b, c, d, f, g, a = a[1]/2., a[2], a[3]/2., a[4]/2., a[5], a[0]

    num = b*b-a*c
    x0  = (c*d-b*f)/num
    y0  = (a*f-b*d)/num

    return x0,y0

def findEllipseRotation(a):
    """
    This function takes in the coefficient vector, fit using contrained least squares of the
    ellipse and calculates the degree of rotation.

    INPUT:
    ----------
    a - the coefficients for the ellipse quadratic form:
    D(x^2, xy, y^2, x, y, 1) * a(a_xx, a_yy, a_x, a_y, a_1) = 0

    OUTPUTS:
    ----------
    the rotation angle of the ellipse

    """

    b, c, d, f, g, a = a[1]/2., a[2], a[3]/2., a[4]/2., a[5], a[0]

    return 0.5*np.arctan(2*b/(a-c))



def findEllipseAxis(a):
    """
    This function takes in the coefficient vector, fit using contrained least squares of the
    ellipse and calculates the major and minor axis.

    INPUT:
    ----------
    a - the coefficients for the ellipse quadratic form:
    D(x^2, xy, y^2, x, y, 1) * a(a_xx, a_yy, a_x, a_y, a_1) = 0

    OUTPUTS:
    ----------
    a - major axis length
    b - minor axis length

    """

    b, c, d, f, g, a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]


    up      = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1   = (b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2   = (b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))

    # a: the axis of the ellipse along the x axis,
    # b: the axis of the ellipse along the y axis.
    a       = np.sqrt(up/down1)
    b       = np.sqrt(up/down2)

    return a, b


def EllipseFitter(data,sigma,isobar):
    """
    This function takes in the coefficient vector, fit using contrained least squares of the
    ellipse and calculates the center coordinates.

    INPUT:
    ----------
    data    - the 2D power spectrum data.
    sigma   - the variance of the Gaussian kernel.
    isobar  - the value of the isobar.

    OUTPUTS:
    ----------
    x           - the x equation of the ellipse
    y           - the y equation of the ellipse
    a           - the length of the major axis
    b           - the length of the minor axis
    [x0, y0]    - the centroid coordinates

    """

    eps     = 5.0 # threshold for where the center is

    x, y    = ContourCoords(data,sigma,isobar)
    C       = fitEllipse(x,y)
    x0, y0  = findEllipseCenter(C)
    phi     = findEllipseRotation(C)
    a, b    = findEllipseAxis(C)

    # Check if the center has not fit properly.
    if np.isnan(np.array([a, b])).any() == True:
        return

    # Check if the center is too far away from the center of the power spectrum
    # x coordinate
    if abs(x0) > data.shape[0]/2.0 + eps or abs(x0) < data.shape[0]/2.0 - eps:
        return

    # y coordinate
    if abs(y0) > data.shape[1]/2.0 + eps or abs(y0) < data.shape[1]/2.0 - eps:
        return

    print("Ellipse fitting: center = [{},{}].".format(x0,y0))
    print("Ellipse fitting: angle of rotation = {}.".format(phi))
    print("Ellipse fitting: axes = [{},{}].".format(a,b))

    # the angle to evaluate the elliptic equation
    theta = np.arange(0,2*np.pi, 0.01)

    # parametric elliptic equations
    x = x0 + a*np.cos(theta)*np.cos(phi) - b*np.sin(theta)*np.sin(phi)
    y = y0 + a*np.cos(theta)*np.sin(phi) + b*np.sin(theta)*np.cos(phi)

    # Extract a coordinate for the major axis
    theta = 0
    xMax = x0 + a*np.cos(theta)*np.cos(phi) - b*np.sin(theta)*np.sin(phi)
    yMax = y0 + a*np.cos(theta)*np.sin(phi) + b*np.sin(theta)*np.cos(phi)

    # Extract a coordinate for the minor axis
    theta = np.pi/2.
    xMin = x0 + a*np.cos(theta)*np.cos(phi) - b*np.sin(theta)*np.sin(phi)
    yMin = y0 + a*np.cos(theta)*np.sin(phi) + b*np.sin(theta)*np.cos(phi)

    return x, y, a, b, [x0,y0], [xMax,yMax,xMin,yMin]

def EllipsePlotter(ax,PowerSpectrum,sigma,isobars,anisotropy):
    """
    This function takes in the coefficient vector, fit using contrained least squares of the
    ellipse and calculates the center coordinates.

    INPUTS:
    ----------
    ax              - the axis for the visualisation.
    PowerSpectrum   - the 2D power spectrum data.
    sigma           - the variance of the Gaussian kernel.
    isobars         - all of the isobars

    OUTPUTS:
    ----------
    an ellipse plotted fit ontop of the power spectrum, to each of the isobars provided
    also the anisotropic ratio is calculated per fit

    """

    # Error handling
    if len(isobars) != 3:
        raise Exception("Incorrect dimensions for the length of the isobar range.")

    # Make an ansiotropic list for averaging later
    aniso   = []

    # Make a counter so that we only label it once
    counter         = 0     # label counter
    prinAxisKeep    = []    # the principle axis of the first ellipse
    centerKeep      = []    # the center of the ellipse
    aAxis           = []
    bAxis           = []
    block           = False # parameter for extracting only the block of ellipses

    # Read in the isobar domain
    bottom = isobars[0]
    top    = isobars[1]
    number = isobars[2]

    for isobar in np.linspace(top,bottom,number):

        try:
            xEllipse, yEllipse, a, b, center, prinAxis  = EllipseFitter(PowerSpectrum,sigma,isobar)
            block = True
        except:
            # Test if we are fitting ellipses in a block
            if not block:
                print "Skipping an isobar at small k."
                continue
            else:
                print "Skipping an isobar at large k."
                prinAxisKeep    = prinAxis   # the principle axis of the first ellipse
                centerKeep      = center
                break
        # if the principle axis is too large (i.e. there isn't a close contour)
        if max(a,b) > 0.5*PowerSpectrum.shape[0]:
            print "Skipping an isobar at large k."
            break
        else:
            ax.plot(xEllipse-center[0],yEllipse-center[1],'blue',alpha=0.8,zorder=2)
            prinAxisKeep    = prinAxis   # the principle axis of the first ellipse
            centerKeep      = center


        # append the ansistropic ratio (make sure the principle axis is always on the bottom)
        if anisotropy == "principle":
            if a < b:
                aniso.append(a/b)
            elif a >= b:
                aniso.append(b/a)
        elif anisotropy == "kperp":
            aniso.append(b/a)
            aAxis.append(a)
            bAxis.append(b)

    # take the anisotropic parameter average across all of the ellipses
    anisoStd    = np.std(aniso)

    # return both the princple axis the centers and the ansio parametee
    return aniso, anisoStd, prinAxisKeep, centerKeep, aAxis, bAxis

def ContourCoords(data,sigma,isobar):
    """
    Produces the contour coordinate for Ellipse fitting.

    INPTS:
    ----------
    data    - the 2D data that the contours will be created on.
    sigma   - the standard deviations of the Gaussin kernel for smoothing.
    isobar  - a single isobar for the contour to go around.

    OUTPUTS:
    ----------
    contours    - the coordinates for the contours

    """

    # Smooth the 2D data, and create the contours on the plot
    contours = measure.find_contours( filters.gaussian(data, sigma), isobar)

    for n, contour in enumerate(contours):
        return (contour[:, 1], contour[:, 0])
