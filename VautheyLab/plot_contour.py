import numpy as np
from matplotlib import pyplot as plt
np.seterr(divide='ignore')
import os
plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'style.mplstyle'))

def load_pdat(file):
    data = np.loadtxt(file, skiprows=1, delimiter=',')
    t = data[1:, 0]
    wl = data[0, 1:]
    dA = data[1:, 1:]  
    return t, wl, dA    

class Contour:
    def __init__(self, file, units='wn', inv=True, t_cuts=None, wl_cuts=None, norm=False, normat=None, minnorm=False, scatter=None, ylim=None, xlim=None, xticks=True,
                 yticks=True, experiment='femto', export=False, MA=False, MA_npoints=10, tightlayout=True, cmap='RdBu_r', lines=False, yscale='linear', IR=False, kk=False,
                 savefig=False, figsize=None, title=None, scale=[-10, 10], nlevels=51, white=False, nwhite=2, zeroline=False, secax=True, arcsinh=False, showcbar=True):          
        # set relevant parameters
        self.figsize = figsize
        # get files
        self.file = file
        self.IR = IR
        self.secax = secax
        # get wavelength cuts
        self.t_cuts = t_cuts 
        self.wl_cuts = wl_cuts
        if not self.wl_cuts==None:
            if self.wl_cuts[0]<100:
                self.wl_cuts[0] = (1/self.wl_cuts[0])*10**4   
                self.wl_cuts[1] = (1/self.wl_cuts[1])*10**4  
        # get scatter 
        self.scatter = scatter
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
        # specify whether you want to invert the x-axis when using kK or not
        self.inv = inv
        # sepcify title of your plot
        self.title = title
        # sepecify your units
        self.units = units
        # whether to export data or not
        self.export = export
        # moving average?
        self.MA = MA
        self.MA_npoints = MA_npoints
        # contour specific things
        self.nlevels = nlevels
        self.white = white
        self.scale = scale
        self.nwhite = nwhite
        self.cmap = cmap
        self.lines = lines
        self.yscale = yscale
        self.zeroline = zeroline
        self.tightlayout = tightlayout
        self.arcsinh = arcsinh
        self.showcbar = showcbar
        self.kk = kk
        # save figure
        self.savefig = savefig
        # intiialize figure
        self.fig, self.ax = plt.subplots(1,1)
        # read data
        self.read_data()
        self.plot()

    def read_data(self):
        # load data
        if '.npy' in self.file:
            self.t, self.wl, self.dA = np.load(self.file, allow_pickle=True)
            self.wn = (1/self.wl)*10**4
        if '.pdat' in self.file:
            if self.IR==False:
                self.t, self.wl, self.dA = load_pdat(self.file)
                self.wn = (1/self.wl)*10**4
            else:
                self.t, self.wn, self.dA = load_pdat(self.file)
                self.wl = (1/self.wn)*10**4
        # cut data 
        if self.wl_cuts!=None:
            if self.IR==False:
                self.dA = self.dA[:, (self.wl>self.wl_cuts[0])&(self.wl<self.wl_cuts[1])]
                self.wn = self.wn[(self.wl>self.wl_cuts[0])&(self.wl<self.wl_cuts[1])]
                self.wl = self.wl[(self.wl>self.wl_cuts[0])&(self.wl<self.wl_cuts[1])]
            else:
                self.dA = self.dA[:, (self.wn>self.wl_cuts[0])&(self.wn<self.wl_cuts[1])]
                self.wl = self.wl[(self.wn>self.wl_cuts[0])&(self.wn<self.wl_cuts[1])]
                self.wn = self.wn[(self.wn>self.wl_cuts[0])&(self.wn<self.wl_cuts[1])]  
                if self.kk == True:              
                    self.wn = self.wn/1000
        # remove scatter
        if self.scatter!=None:
             self.dA[:, (self.wl >= self.scatter[0]) & (self.wl <= self.scatter[1])] = np.nan
        # normalize
        if self.minnorm==False:
            if self.norm == True: 
                if self.normat==None:
                    self.dA = self.dA/np.nanmax(self.dA)
                else:
                    self.dA = self.dA/self.dA[np.argmin(np.abs(self.t - self.normat[0])), np.argmin(np.abs(self.wl - self.normat[1]))]
        else:
            if self.norm == True: 
                if self.normat==None:
                    self.dA = -1*self.dA/np.nanmin(self.dA)
                else:
                    self.dA = -1*self.dA/self.dA[np.argmin(np.abs(self.t - self.normat[0])), np.argmin(np.abs(self.wl - self.normat[1]))]
        if self.export==True:
            for i in range(len(self.delays)):
                y = self.dA[np.argmin(np.abs(self.t - self.delays[i])), :]
                np.savetxt('spectrum_at_%.3g_ps.txt'%(self.delays[i]), np.column_stack([self.wl, self.wn, y]), delimiter=',',
                           header='wavelength / nm,     wavenumber / kK,    TA / mOD')
        if self.arcsinh==True:
            self.dA = np.arcsinh(self.dA)
                
    def plot(self):
        levels = np.linspace(self.scale[0], self.scale[1], self.nlevels)
        if self.cmap=='Reds':
            colors = plt.cm.Reds(np.linspace(0, 1, self.nlevels))   
        if self.cmap=='RdBu_r':
            colors = plt.cm.RdBu_r(np.linspace(0, 1, self.nlevels))  
        if self.cmap=='seismic':
            colors = plt.cm.seismic(np.linspace(0, 1, self.nlevels))  
        if self.white==True:
            colors[int(self.nlevels/2)] = 0
            colors[int(self.nlevels/2)-1] = 0
            colors[int(self.nlevels/2)-2] = 0
        if self.units=='wl':
            D = self.ax.contourf(self.wl, self.t, self.dA, levels=levels, colors=colors)
            if self.lines==True:
                self.ax.contour(self.wl, self.t, self.dA, levels=levels, colors='k', linestyles='-', linewidths=0.5, alpha=0.2)
        else:
            D = self.ax.contourf(self.wn, self.t, self.dA, levels=levels, colors=colors)
            if self.lines==True:
                self.ax.contour(self.wn, self.t, self.dA, levels=levels, colors='k', linestyles='-', linewidths=0.5, alpha=0.2)
        self.cbar = plt.colorbar(D, ax=self.ax)

    def show(self):
        if self.figsize!=None:
            self.fig.set_size_inches(self.figsize)
        if not self.norm==True:
            self.cbar.set_label(r'$\Delta{A} \, \times 10^{3}$')  
        else:
            self.cbar.set_label(r'norm. $\Delta{A}$')
        self.cbar.set_ticks([self.scale[0], self.scale[1]])

        if self.units=='wn':
            if self.IR==False:
                self.ax.set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
                if self.inv==True:
                    self.ax.invert_xaxis()    
                if self.secax==True:
                    axsec = self.ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
                    axsec.set_xlabel(r'$\lambda / $ nm')
            else:
                if self.kk == False:
                    self.ax.set_xlabel(r'$\tilde{\nu} / \text{cm}^{-1}$')
                else:
                    self.ax.set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
                if self.secax==True:
                    axsec = self.ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
                    axsec.set_xlabel(r'$\lambda / $ µm')

        if self.units=='wl':
            self.ax.set_xlabel(r'$\lambda / $ nm')

        if self.experiment=='femto':
            self.ax.set_ylabel(r'$\Delta t$ / ps')
        else:
            self.ax.set_ylabel(r'$\Delta t$ / ns')

        if self.zeroline==True:
            self.ax.axhline(y=0, color='k')

        if self.showcbar==False:
            self.cbar.remove()

        self.ax.set_yscale(self.yscale)
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
            self.fig.savefig('%s_contour.pdf'%(self.file[:self.file.find('.')]), transparent=True)
        plt.show()