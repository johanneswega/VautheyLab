import numpy as np
from matplotlib import pyplot as plt
from VautheyLab.miscellaneous import rainbow, moving_average
from VautheyLab.fit_functions import *
np.seterr(divide='ignore')
import os
plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'style.mplstyle'))

def load_pdat(file):
    data = np.loadtxt(file, skiprows=1, delimiter=',')
    t = data[1:, 0]
    wl = data[0, 1:]
    dA = data[1:, 1:]  
    return t, wl, dA  

class Kinetics:
    def __init__(self, file, cuts=None, norm=False, normat=None, minnorm=False, ylim=None, xlim=None, xticks=True, xscale='linear', devide=None,
                 yticks=True, wavelengths=None, experiment='femto', colors=None, export=False, MA=False, MA_npoints=10, tightlayout=True, 
                 plot_now = True, savefig=False, figsize=None, title=None,  alpha=1, fmt='.', markersize=4, fontsize=10, outside=False, ylabel=None, legend=True):       
        # set relevant parameters
        self.figsize = figsize
        # get files
        self.file = file
        # get wavelength cuts
        self.cuts = cuts    
        # get experiment
        self.experiment = experiment
        # decide whether to normalize or not
        self.norm = norm
        # if you want to norm on GSB
        self.minnorm = minnorm
        # to normalize at a specific wl
        self.normat = normat
        if self.normat!=None and self.normat[1]<100:
            self.normat[1] = (1/self.normat[1])*10**4
        # specify x and y limits for plot
        self.ylim = ylim
        self.xlim = xlim
        # specify xticks
        self.xticks = xticks
        self.yticks = yticks
        self.xscale = xscale
        # specify list of wavelengths
        self.wavelengths = wavelengths
        # if kK convert to nm
        if self.wavelengths[0]<100:
            for i in range(len(self.wavelengths)):
                self.wavelengths[i] = (1/self.wavelengths[i])*10**4
        # specify list of colors
        self.colors = colors
        #if colors==None:
            #self.colors = rainbow(self.wavelengths, r=True)
        # sepcify title of your plot
        self.title = title
        # whether to export data or not
        self.export = export
        # moving average?
        self.MA = MA
        self.MA_npoints = MA_npoints
        # tight layout
        self.tightlayout = tightlayout
        # save figure
        self.savefig = savefig
        self.alpha = alpha
        self.fmt = fmt
        self.markersize = markersize
        self.fontsize=fontsize
        self.outside=outside
        self.ylabel=ylabel
        self.legend=legend
        # if you want to fit kinetics 
        self.devide=devide

        # intiialize figure
        self.fig, self.ax = plt.subplots(1,1)
        # read data
        self.read_data()
        # plot data
        if plot_now==True:
            self.plot()       

    def read_data(self):
        # load data
        if '.npy' in self.file:
            self.t, self.wl, self.dA = np.load(self.file, allow_pickle=True)
        if '.pdat' in self.file:
            self.t, self.wl, self.dA = load_pdat(self.file)
        self.wn = (1/self.wl)*10**4
        # cut data 
        if self.cuts!=None:
            self.dA = self.dA[(self.t>self.cuts[0])&(self.t<self.cuts[1]), :]
            self.t = self.t[(self.t>self.cuts[0])&(self.t<self.cuts[1])]
        # normalize
        if self.minnorm==False:
            if self.norm == True: 
                if self.normat==None:
                    self.dA = self.dA/np.max(self.dA)
                else:
                    self.dA = self.dA/self.dA[np.argmin(np.abs(self.t - self.normat[0])), np.argmin(np.abs(self.wl - self.normat[1]))]
        else:
            if self.norm == True: 
                if self.normat==None:
                    self.dA = -1*self.dA/np.min(self.dA)
                else:
                    self.dA = -1*self.dA/self.dA[np.argmin(np.abs(self.t - self.normat[0])), np.argmin(np.abs(self.wl - self.normat[1]))]
        if self.export==True:
            for i in range(len(self.wavelengths)):
                y = self.dA[:, np.argmin(np.abs(self.wl - self.wavelengths[i]))]
                np.savetxt('kinetics_at_%.3g_nm.txt'%(self.wavelengths[i]), np.column_stack([self.t, y]), delimiter=',', header='time / ps,     TA / mOD')

    def plot(self):
        # got through all delays and plot
        x = self.t
        for i in range(len(self.wavelengths)):
            if self.norm==True:
                y = self.dA[:, np.argmin(np.abs(self.wl - self.wavelengths[i]))]/self.dA[0, np.argmin(np.abs(self.wl - self.wavelengths[i]))]
            else:
                y = self.dA[:, np.argmin(np.abs(self.wl - self.wavelengths[i]))]
            if self.devide!=None:
                y = y/self.devide[i]
            lab = '%.3g nm'%(self.wavelengths[i])        
            # plot
            if self.MA==False:
                self.ax.plot(x, y, self.fmt, color=self.colors[i], alpha=self.alpha, label=lab, markersize=self.markersize)
            else:
                self.ax.plot(x, y, '-', color=self.colors[i], alpha=0.2)
                self.ax.plot(moving_average(x, self.MA_npoints), moving_average(y, self.MA_npoints), '-', color=self.colors[i], label=lab)

    def show(self):
        if self.figsize!=None:
            self.fig.set_size_inches(self.figsize)
        
        if self.experiment=='femto':
            self.ax.set_xlabel(r'$\Delta t / $ ps')
        elif self.experiment=='nano':
            self.ax.set_xlabel(r'$\Delta t / $ ns')        

        self.ax.set_xscale(self.xscale)
        if self.legend==True:
            if self.outside==True:
                self.ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=self.fontsize)
            else:
                self.ax.legend(fontsize=self.fontsize)
        self.ax.axhline(y=0, color='k')
        if self.ylabel==None:
            if self.norm == True:
                self.ax.set_ylabel(r'norm. $\Delta A$')
            else:
                self.ax.set_ylabel(r'$\Delta A / 10^{-3}$')        
        else:
            self.ax.set_ylabel(self.ylabel)       
        if self.yticks == False:
            self.ax.set_yticks([])
        if self.xticks == False:
            self.ax.set_xticks([])
        if self.ylim != None:
            self.ax.set_ylim(self.ylim)
        if self.xlim != None: 
            self.ax.set_xlim(self.xlim)
        if self.title != None:
            self.ax.set_title(self.title)
        if self.tightlayout == True:
            self.fig.tight_layout()
        if self.savefig == True: 
            self.fig.savefig('%s_kinetics.svg'%(self.file[:self.file.find('.')]), transparent=True)
        plt.show()
