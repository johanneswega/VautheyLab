import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
plt.style.use('/Users/wega/Documents/matplotlib/test.mplstyle')
import os

def fit(Aex, grad):
    return grad*Aex

def get_abs(fname, cuts):
    # load data
    data = np.loadtxt(fname, skiprows=2, delimiter=',', usecols=(0,1))
    # get wavelength
    wl = data[:,0]
    # absorbance
    A = data[:,1]
    # subtract baseline
    base = np.mean(A[:150])
    A = A - base
    # cut data
    A = A[(wl>cuts[0]) & (wl<cuts[1])]
    wl = wl[(wl>cuts[0]) & (wl<cuts[1])]
    return wl, A

def get_em(fname, cuts):
    # load data
    data = np.loadtxt(fname, skiprows=2, delimiter=',')
    # get wavelength
    wl = data[:,0]
    # get intensity
    I = data[:,1]
    # subtract baseline for ANS
    # cut data
    I = I[(wl>cuts[0]) & (wl<cuts[1])]
    wl = wl[(wl>cuts[0]) & (wl<cuts[1])]    
    # load corr
    corr = np.loadtxt('corr.csv', skiprows=1, delimiter=',')
    # interpolate
    c = np.interp(wl, corr[:,0], corr[:,1])
    # correct 
    I = I*c
    return wl, I
    
