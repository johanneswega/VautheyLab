import numpy as np
from matplotlib import pyplot as plt
np.seterr(divide='ignore')
import os
plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'style.mplstyle'))

def figure_spectum(figsize=(8,5), IR=False):
    fig, ax = plt.subplots(1,1,figsize=figsize)
    if IR==False:
        ax.set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
        ax.invert_xaxis()    
        ax.axhline(y=0, color='k')
        axsec = ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
        axsec.set_xlabel(r'$\lambda / $ nm') 
    else:
        ax.set_xlabel(r'$\tilde{\nu} / \text{cm}^{-1}$')
    ax.axhline(y=0, color='k')
    return fig, ax

def figure_kinetics(figsize=(8,5), experiment='femto'):
    fig, ax = plt.subplots(1,1,figsize=figsize)
    if experiment=='femto':
        ax.set_xlabel(r'$\Delta t$ / ps')
    else:
        ax.set_xlabel(r'$\Delta t$ / ns')
    return fig, ax

def figure_contour(t, wl, dA, scale, figsize=(8,5), experiment='femto',
                   ylim=None, lines=True, nlines=25):
    wn = (1/wl)*10**4
    # make figure
    fig, ax = plt.subplots(1,1,figsize=figsize)
    if lines==True:
        ax.contour(wn, t, dA, levels=nlines, linewidths=0.12, colors='k')
    D = ax.pcolormesh(wn, t, dA, vmin=scale[0], vmax=scale[1], cmap="RdBu_r")
    cbar = plt.colorbar(D, ax=ax)
    cbar.set_label(r'$\Delta{A} / 10^{-3}$')  
    ax.set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
    ax.invert_xaxis()   
    axsec = ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
    axsec.set_xlabel(r'$\lambda / $ nm') 
    if experiment=='femto':
        ax.set_ylabel(r'$\Delta t / \text{ps}$')
    if ylim!=None:
        ax.set_ylim(ylim)
    return fig, ax