import numpy as np
from matplotlib import pyplot as plt
from VautheyLab.miscellaneous import rainbow, moving_average
np.seterr(divide='ignore')
import os
plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'style.mplstyle'))

def load_pdat(file):
    data = np.loadtxt(file, skiprows=1, delimiter=',')
    t = data[1:, 0]
    wl = data[0, 1:]
    dA = data[1:, 1:]  
    return t, wl, dA  

class Overview:

    def __init__(self, file, units='wn', inv=True, cuts=None, norm=False, normat=None, minnorm=False, scatter=None, ylim=None, xlim=None, xticks=True,
                 yticks=True, delays=[None], experiment='femto', colors=None, export=False, MA=False, MA_npoints=10, tightlayout=True, steady_state=None,
                 savefig=False, figsize=None, title=None, plot_now=True, legend=True, ylabel=None, IR=False):       
        # set relevant parameters
        self.figsize = figsize
        # get files
        self.file = file
        # get wavelength cuts
        self.cuts = cuts    
        self.ylabel = ylabel
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
        self.steady_state = steady_state
        # specify whether you want to invert the x-axis when using kK or not
        self.inv = inv
        # specify list of delays
        if delays[0]==None:
            if self.experiment == 'femto':
                self.delays = [-1, 0.2, 0.5, 1, 2, 5, 10, 50, 100, 250, 500, 1000, 1200, 1800]
            elif self.experiment == 'nano':
                self.delays = [2, 5, 10, 20, 50, 80, 100, 250, 200, 500, 1e+3, 5e+3, 20e+3, 50e+3, 100e+3, 200e+3, 500e+3]
        else:
            self.delays = delays
        # specify list of colors
        if colors!=None:
            self.colors = colors
        else:
            self.colors = rainbow(self.delays, r=True)
        # sepcify title of your plot
        self.title = title
        # wheter TRIR
        self.IR = IR
        # sepecify your units
        self.units = units
        # whether to export data or not
        self.export = export
        # moving average?
        self.MA = MA
        self.MA_npoints = MA_npoints
        # tight layout
        self.tightlayout = tightlayout
        # save figure
        self.savefig = savefig
        self.legend = legend
        # intiialize figure
        self.fig, self.ax = plt.subplots(1,1)
        # read data
        self.read_data()
        # plot data immedeatly if you want to
        if plot_now==True:
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
        if self.cuts!=None:
            if self.IR==False:
                self.dA = self.dA[:, (self.wl>self.cuts[0])&(self.wl<self.cuts[1])]
                self.wn = self.wn[(self.wl>self.cuts[0])&(self.wl<self.cuts[1])]
                self.wl = self.wl[(self.wl>self.cuts[0])&(self.wl<self.cuts[1])]
            else:
                self.dA = self.dA[:, (self.wn>self.cuts[0])&(self.wn<self.cuts[1])]
                self.wl = self.wl[(self.wn>self.cuts[0])&(self.wn<self.cuts[1])]
                self.wn = self.wn[(self.wn>self.cuts[0])&(self.wn<self.cuts[1])]                
        # remove scatter
        if self.scatter!=None:
             self.dA[:, (self.wl >= self.scatter[0]) & (self.wl <= self.scatter[1])] = np.nan
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
            for i in range(len(self.delays)):
                y = self.dA[np.argmin(np.abs(self.t - self.delays[i])), :]
                np.savetxt('spectrum_at_%.3g_ps.txt'%(self.delays[i]), np.column_stack([self.wl, self.wn, y]), delimiter=',',
                           header='wavelength / nm,     wavenumber / kK,    TA / mOD')

    def plot(self):
        # choose unit
        if self.units == 'wn':
            x = self.wn
        else: 
            x = self.wl 
        # got through all delays and plot
        for i in range(len(self.delays)):
            y = self.dA[np.argmin(np.abs(self.t - self.delays[i])), :]
            # make label
            if self.experiment == 'femto':
                if np.abs(self.delays[i])<1:
                    lab = r'%.3g fs'%(self.delays[i]*1000)
                elif 1<=np.abs(self.delays[i])<1000:
                    lab = r'%.3g ps'%(self.delays[i])
                else:
                    lab = r'%.3g ns'%(self.delays[i]/1000)
            if self.experiment == 'nano':
                if np.abs(self.delays[i])<1:
                    lab = r'%.3g ps'%(self.delays[i]*1000)
                elif 1<=np.abs(self.delays[i])<1000:
                    lab = r'%.3g ns'%(self.delays[i])
                else:
                    lab = r'%.3g µs'%(self.delays[i]/1000)            
            # plot
            if self.MA==False:
                self.ax.plot(x, y, '-', color=self.colors[i], label=lab)
            else:
                self.ax.plot(x, y, '-', color=self.colors[i], alpha=0.2)
                self.ax.plot(moving_average(x, self.MA_npoints), moving_average(y, self.MA_npoints), '-', color=self.colors[i], label=lab)

        # plot steady state spectra
        if self.steady_state!=None:
            for i in range(len(self.steady_state)):
                # read data
                print(self.steady_state[i][0])
                data = np.loadtxt(self.steady_state[i][0], delimiter=',', skiprows=1)
                s_wl = data[:,0]
                s_wn = data[:,1]
                s = data[:,2]
                # multiply intensity by lambda^4 if emission
                if 'em' in self.steady_state[i][0]:
                    s = s*s_wl**4
                    print(self.steady_state[i][0], ' has been multiplied by lambda^4 to convert to SE.')
                if self.IR==False:
                    s = s[(s_wl>self.steady_state[i][1][0]) & (s_wl<self.steady_state[i][1][1])]
                    s_wn = s_wn[(s_wl>self.steady_state[i][1][0]) & (s_wl<self.steady_state[i][1][1])]
                    s_wl = s_wl[(s_wl>self.steady_state[i][1][0]) & (s_wl<self.steady_state[i][1][1])]
                else:
                    s = s[(s_wn>self.steady_state[i][1][0]) & (s_wn<self.steady_state[i][1][1])]
                    s_wl = s_wl[(s_wn>self.steady_state[i][1][0]) & (s_wn<self.steady_state[i][1][1])]
                    s_wn = s_wn[(s_wn>self.steady_state[i][1][0]) & (s_wn<self.steady_state[i][1][1])]                    
                s = self.steady_state[i][2]*s/np.nanmax(s)
                if self.units == 'wl':
                    self.ax.fill_between(s_wl, 0, s, color=self.steady_state[i][3], label=self.steady_state[i][4], alpha=0.1)
                else:
                   self.ax.fill_between(s_wn, 0, s, color=self.steady_state[i][3], label=self.steady_state[i][4], alpha=0.1) 

    def show(self):
        if self.figsize!=None:
            self.fig.set_size_inches(self.figsize)
        
        if self.units=='wn':
            if self.IR==False:
                self.ax.set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
                axsec = self.ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
                axsec.set_xlabel(r'$\lambda / $ nm') 
                self.ax.invert_xaxis()
            else:
                self.ax.set_xlabel(r'$\tilde{\nu} / \text{cm}^{-1}$')

        if self.units=='wl':
            self.ax.set_xlabel(r'$\lambda / $ nm')

        if self.legend==True:
            self.ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
        self.ax.axhline(y=0, color='k')
        if self.norm == True:
            self.ax.set_ylabel(r'norm. $\Delta A$')
        else:
            self.ax.set_ylabel(r'$\Delta A / 10^{-3}$')                
        if self.yticks == False:
            self.ax.set_yticks([])
        if self.xticks == False:
            self.ax.set_xticks([])
        if self.ylim != None:
            self.ax.set_ylim(self.ylim)
        if self.xlim != None: 
            self.ax.set_xlim(self.xlim)
        if self.ylabel != None:
            self.ax.set_ylabel(self.ylabel)
        if self.title != None:
            self.ax.set_title(self.title)
        if self.tightlayout == True:
            self.fig.tight_layout()
        if self.savefig == True: 
            self.fig.savefig('%s.svg'%(self.file[:self.file.find('.')]), transparent=True)
        plt.show()

        