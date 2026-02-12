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

class Compare_Contours:
    def __init__(self, files, units='wn', inv=True, t_cuts=None, wl_cuts=None, norm=False, normat=None, minnorm=False, scatter=None, ylim=None, xlim=None, xticks=False,
                 yticks=True, experiment='femto', export=False, MA=[False], MA_npoints=10, tightlayout=True, cmap='RdBu_r', lines=False, yscale='linear', arcsinh=False,
                 savefig=False, figsize=None, titles=None, scale=[-10, 10], nlevels=51, white=False, nwhite=2, zeroline=False, colorbar=False, IR=False,
                 devide=None):          
        # set relevant parameters
        self.figsize = figsize
        # get files
        self.files = files
        self.arcsinh = arcsinh
        self.devide = devide
        self.IR = IR
        # get cuts
        self.t_cuts = t_cuts
        self.wl_cuts = wl_cuts
        if self.wl_cuts!=None:
            if len(self.wl_cuts)==1 and len(self.files)>1:
                # use same range if only one cuts specified
                for i in range(len(self.files) - 1):
                    self.wl_cuts.append(self.wl_cuts[0])  
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
        self.titles = titles
        self.savefig = savefig
        self.colorbar = colorbar
        self.units = units
        # whether to export data or not
        self.export = export
        # moving average?
        self.MA = MA
        if len(self.MA)==1 and len(self.files)>1:
            for i in range(len(self.files)):
                self.MA.append(self.MA[0])
        self.MA_npoints = MA_npoints
        # tight layout
        self.tightlayout = tightlayout
        # intiialize figure
        self.fig, self.ax = plt.subplots(1, len(self.files), sharey=True)
        # read data
        self.read_data()
        self.plot()

    def read_data(self):
        # initialize list of absorbances
        self.dA = []
        self.wl = []
        self.wn = []
        self.t = []
        # go through all files
        for i in range(len(self.files)):
            # load data
            if '.npy' in self.files[i]:
                t, wl, dA = np.load(self.files[i], allow_pickle=True)
                wn = (1/wl)*10**4
            if '.pdat' in self.files[i]:
                if self.IR==False:
                    t, wl, dA = load_pdat(self.files[i])
                    wn = (1/wl)*10**4
                else:
                    t, wn, dA = load_pdat(self.files[i])
                    wl = (1/wn)*10**4
                
            # cut data 
            if self.t_cuts!=None:
                dA = dA[(t>self.t_cuts[0])&(t<self.t_cuts[1]), :]
                t = t[(t>self.t_cuts[0])&(t<self.t_cuts[1])]
            if self.wl_cuts!=None:
                if self.IR==False:
                    dA = dA[:, (wl>self.wl_cuts[0])&(wl<self.wl_cuts[1])]
                    wn = wn[(wl>self.wl_cuts[0])&(wl<self.wl_cuts[1])]
                    wl = wl[(wl>self.wl_cuts[0])&(wl<self.wl_cuts[1])]
                else:
                    dA = dA[:, (wn>self.wl_cuts[0])&(wn<self.wl_cuts[1])]
                    wl = wl[(wn>self.wl_cuts[0])&(wn<self.wl_cuts[1])]
                    wn = wn[(wn>self.wl_cuts[0])&(wn<self.wl_cuts[1])]                    
            # remove scatter
            if self.scatter!=None:
                dA[:, (wl >= self.scatter[i][0]) & (wl <= self.scatter[i][1])] = np.nan
            # normalize
            if self.minnorm==False:
                if self.norm == True: 
                    if self.normat==None:
                        dA = dA/np.nanmax(dA[0,:])
                    else:
                        dA = dA/dA[np.argmin(np.abs(t - self.normat[0])), np.argmin(np.abs(wl - self.normat[1]))]
            else:
                if self.norm == True: 
                    if self.normat==None:
                        dA = -1*dA/np.nanmin(dA[0,:])
                    else:
                        dA = -1*dA/dA[np.argmin(np.abs(t - self.normat[0])), np.argmin(np.abs(wl - self.normat[1]))]
            # append
            if self.arcsinh==True:
                dA = np.arcsinh(dA)
                self.scale = np.arcsinh(self.scale)
            if self.devide!=None:
                dA = dA/self.devide[i]
            self.dA.append(dA)
            self.t.append(t)
            self.wl.append(wl)
            self.wn.append(wn)

    def plot(self):
        # set levels and contour design
        levels = np.linspace(self.scale[0], self.scale[1], self.nlevels)
        if self.cmap=='Reds':
            colors = plt.cm.Reds(np.linspace(0, 1, self.nlevels))   
        if self.cmap=='RdBu_r':
            colors = plt.cm.RdBu_r(np.linspace(0, 1, self.nlevels))  
        if self.cmap=='seismic':
            colors = plt.cm.seismic(np.linspace(0, 1, self.nlevels))  
        if self.white==True:
            colors[int(self.nlevels/2)-1] = 0
            for i in range(2):
                colors[int(self.nlevels/2)-1+i] = 0
                colors[int(self.nlevels/2)-1-i] = 0

        # go through all files
        for i in range(len(self.files)):
            if self.units=='wl':
                D = self.ax[i].contourf(self.wl[i], self.t[i], self.dA[i], levels=levels, colors=colors)
                if self.lines==True:
                    self.ax[i].contour(self.wl[i], self.t[i], self.dA[i], levels=levels, colors='k', linestyles='-', linewidths=0.5, alpha=0.1)
            else:
                D = self.ax[i].contourf(self.wn[i], self.t[i], self.dA[i], levels=levels, colors=colors)
                if self.lines==True:
                    self.ax[i].contour(self.wn[i], self.t[i], self.dA[i], levels=levels, colors='k', linestyles='-', linewidths=0.5, alpha=0.1)
        if self.colorbar==True:
            self.cbar = self.fig.colorbar(D, ax=self.ax[len(self.files)-1])

    def show(self):
        if self.figsize!=None:
            self.fig.set_size_inches(self.figsize)
        # colorbar
        if self.colorbar==True:
            if not self.norm==True:
                self.cbar.set_label(r'$\Delta{A} / 10^{-3}$')  
                if self.arcsinh==True:
                    self.cbar.set_label(r'arcsinh($\Delta{A} / 10^{-3}$)')  
            else:
                self.cbar.set_label(r'norm. $\Delta{A}$')
                if self.arcsinh==True:
                    self.cbar.set_label(r'norm. arcsinh($\Delta{A}$)')
            self.cbar.set_ticks([self.scale[0], self.scale[1]])
            self.fig.savefig('colorbar.svg', transparent=True)

        # go through files
        for i in range(len(self.files)):
            if self.zeroline==True:
                self.ax[i].axhline(y=0, color='k')

            if self.units=='wn':
                if self.IR==False:
                    self.ax[i].set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
                    if self.inv==True:
                        self.ax[i].invert_xaxis()    
                    axsec = self.ax[i].secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
                    axsec.set_xlabel(r'$\lambda / $ nm') 
                else:
                    self.ax[i].set_xlabel(r'$\tilde{\nu} / \text{cm}^{-1}$')
            if self.units=='wl':
                self.ax[i].set_xlabel(r'$\lambda / $ nm')
            if self.xticks != False:
                self.ax[i].set_xticks(self.xticks)
            if self.xlim != None: 
                self.ax[i].set_xlim(self.xlim)
            if type(self.titles) != None:
                self.ax[i].set_title(self.titles[i])

        if self.experiment=='femto':
            self.ax[0].set_ylabel(r'$\Delta t$ / ps')
        else:
            self.ax[0].set_ylabel(r'$\Delta t$ / ns')

        self.ax[0].set_yscale(self.yscale)
        if self.yticks == False:
            self.ax[0].set_yticks([])
        if self.ylim != None:
            self.ax[0].set_ylim(self.ylim)
        if self.tightlayout == True:
            self.fig.tight_layout()
        if self.savefig == True: 
            self.fig.savefig('contour_comp.pdf', transparent=True)
        plt.show()