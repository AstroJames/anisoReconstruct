############################################################################################################################################
# Data handling
import os
import yt
import h5py
import numpy as np
import argparse
import pickle

# Math
import imageio                             # reading in image data
import h5py                                # importing in hdf5 files
import skimage                             # import image data
from skimage import measure, filters       # for drawing contours and Gaussian filters
from scipy import fftpack, misc, optimize  # fourier transform
from scipy.integrate import quad
from scipy.stats import kurtosis, skew
from scipy.optimize import curve_fit
from numpy.linalg import eig, inv

# Visualisations
import matplotlib as mpl
import matplotlib.pyplot as plt            # visualisation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc, ticker, colors  # nicer text in matplotlib and custum ticks for cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patheffects as PathEffects
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',r'\usepackage{mathrsfs}']
rc('font', **{'family': 'DejaVu Sans'})
rc('text', usetex=True)
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

############################################################################################################################################

fs = 16
