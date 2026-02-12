import numpy as np
from matplotlib import pyplot as plt
from VautheyLab.miscellaneous import rainbow, find_index, moving_average
np.seterr(divide='ignore')
import matplotlib.patches as patches
import os
plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'style.mplstyle'))

def load_pdat(file):
    data = np.loadtxt(file, skiprows=1, delimiter=',')
    t = data[1:, 0]
    wl = data[0, 1:]
    dA = data[1:, 1:]  
    return t, wl, dA    

class Compare_Spectra:
    def __init__(self, files, cuts=None, norm=False, normat=None, minnorm=False, ylim=None, xlim=None, xticks=True, units='wn', inv=True,
                 yticks=True, delays=None, experiment='femto', colors=None, export=False, MA=False, MA_npoints=10, tightlayout=True, legend=True,
                 savefig=False, figsize=None, title=None, labels=[''], plot_now=True, alpha=1, fmt='-', markersize=5, fontsize=10, outside=False,
                 zeroline=True, scatter=[None], steady_state=None, ylabel='', scatterbar=False, devide=None, IR=False):       
        # set relevant parameters
        self.figsize = figsize
        # get files
        self.files = files
        self.IR = IR
        # get wavelength cuts
        self.cuts = cuts    
        # get experiment
        self.experiment = experiment
        # decide whether to normalize or not
        # devide spectra by factor 
        self.devide = devide
        self.norm = norm
        # follows the tuple structure (file, cuts, scale, color, label)
        self.steady_state = steady_state
        # get scatter 
        self.scatter = scatter
        if len(self.files)>1 and not self.scatter[0]==None:
            for i in range(len(self.files) - 1):
                self.scatter.append(self.scatter[0])
        # if you want to norm on GSB
        self.minnorm = minnorm
        # to normalize at a specific wl
        self.normat = normat
        if self.normat!=None and self.normat[1]<100:
            self.normat[1] = (1/self.normat[1])*10**4
        # specify cuts
        self.cuts = cuts
        if self.cuts!=None:
            if len(self.cuts)==1 and len(self.files)>1:
                # use same range if only one cuts specified
                for i in range(len(self.files) - 1):
                    self.cuts.append(self.cuts[0])
        # specify x and y limits for plot
        self.ylim = ylim
        self.xlim = xlim
        # specify xticks
        self.xticks = xticks
        self.yticks = yticks
        self.units = units
        self.ylabel = ylabel
        # specify list of delays
        self.delays = delays
        # take same if only one specified
        if len(self.delays)==1 and len(self.files)>1:
            # use same range if only one cuts specified
            for i in range(len(self.files) - 1):
                self.delays.append(self.delays[0])        
        # sepcify title of your plot
        self.title = title
        # treat labels
        self.labels=labels
        if self.labels==[''] and len(self.files)>1:
            for i in range(len(self.files) - 1):
                self.labels.append('')
        # whether to export data or not
        self.export = export
        # specifie colors
        self.colors = colors
        # moving average?
        if MA!=False:
            self.MA = MA
            self.MA_npoints = MA_npoints
        else:
            self.MA = [False for _ in range(len(self.files))]      
        # tight layout
        self.tightlayout = tightlayout
        # save figure
        self.savefig = savefig
        self.alpha = alpha
        self.fmt = fmt
        self.markersize = markersize
        self.fontsize=fontsize
        self.outside=outside
        self.zeroline=zeroline
        self.inv = inv
        self.legend = legend
        self.scatterbar = scatterbar
        # intiialize figure
        self.fig, self.ax = plt.subplots(1,1)
        # read data
        self.read_data()
        # plot data
        if plot_now==True:
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
            if self.cuts!=None:
                if self.IR==False:
                    dA = dA[:, (wl>self.cuts[0])&(wl<self.cuts[1])]
                    wn = wn[(wl>self.cuts[0])&(wl<self.cuts[1])]
                    wl = wl[(wl>self.cuts[0])&(wl<self.cuts[1])]
                else:
                    dA = dA[:, (wn>self.cuts[0])&(wn<self.cuts[1])]
                    wl = wl[(wn>self.cuts[0])&(wn<self.cuts[1])]
                    wn = wn[(wn>self.cuts[0])&(wn<self.cuts[1])]      
            # remove scatter
            if self.scatter[0]!=None:
                dA[:, (wl >= self.scatter[i][0]) & (wl <= self.scatter[i][1])] = np.nan
            # normalize
            if self.minnorm==False:
                if self.norm == True: 
                    if self.normat==None:
                        if self.MA[i]==True:
                            dA = dA/np.nanmax(moving_average(dA[np.argmin(np.abs(t - self.delays[i])), :], self.MA_npoints[i]))
                        else:
                            dA = dA/np.nanmax(dA[np.argmin(np.abs(t - self.delays[i])), :])
                    else:
                        dA = dA/dA[np.argmin(np.abs(t - self.normat[0])), np.argmin(np.abs(wl - self.normat[1]))]
            else:
                if self.norm == True: 
                    if self.normat==None:
                        dA = -1*dA/np.nanmin(dA[np.argmin(np.abs(t - self.delays[i])), :])
                    else:
                        dA = -1*dA/dA[np.argmin(np.abs(t - self.normat[0])), np.argmin(np.abs(wl - self.normat[1]))]
            # devide if wanted
            if self.devide!=None:
                dA = dA/self.devide[i]
            # append
            self.dA.append(dA)
            self.t.append(t)
            self.wl.append(wl)
            self.wn.append(wn)

    def plot(self):
        # go through all files
        for i in range(len(self.files)):
            if self.units=='wn':
                x = self.wn[i]
            else:
                x = self.wl[i]
            y = self.dA[i][np.argmin(np.abs(self.t[i] - self.delays[i])), :]  
            # plot
            if self.MA[i]==False:
                self.ax.plot(x, y, self.fmt, color=self.colors[i], alpha=self.alpha, label=self.labels[i], markersize=self.markersize)
            else:
                self.ax.plot(x, y, '-', color=self.colors[i], alpha=0.2)
                self.ax.plot(moving_average(x, self.MA_npoints[i]), moving_average(y, self.MA_npoints[i]), '-', color=self.colors[i], label=self.labels[i])
        # plot steady state spectra
        if self.steady_state!=None:
            for i in range(len(self.steady_state)):
                # read data
                data = np.loadtxt(self.steady_state[i][0], delimiter=',', skiprows=1)
                s_wl = data[:,0]
                s_wn = data[:,1]
                s = data[:,2]
                # multiply intensity by lambda^4 if emission
                if 'em' in self.steady_state[i][0]:
                    s = s*s_wl**4
                    print(self.steady_state[i][0], ' has been multiplied by lambda^4 to convert to SE.')
                s = s[(s_wl>self.steady_state[i][1][0]) & (s_wl<self.steady_state[i][1][1])]
                s_wn = s_wn[(s_wl>self.steady_state[i][1][0]) & (s_wl<self.steady_state[i][1][1])]
                s_wl = s_wl[(s_wl>self.steady_state[i][1][0]) & (s_wl<self.steady_state[i][1][1])]
                s = self.steady_state[i][2]*s/np.max(s)
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
                if self.inv==True:
                    self.ax.invert_xaxis()    
                axsec = self.ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
                axsec.set_xlabel(r'$\lambda / $ nm') 
            else:
                self.ax.set_xlabel(r'$\tilde{\nu} / \text{cm}^{-1}$')

        if self.units=='wl':
            self.ax.set_xlabel(r'$\lambda / $ nm')

        if self.legend==True:
            if self.outside==True:
                self.ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=self.fontsize)
            else:
                self.ax.legend(fontsize=self.fontsize)

        if len(set(self.delays)) == 1:
            if self.experiment == 'femto':
                if np.abs(self.delays[0])<1:
                    self.ax.set_title(r'$\Delta t = %.3g\,\text{fs}$'%(self.delays[0]*1000))
                elif 1<=np.abs(self.delays[0])<1000:
                    self.ax.set_title(r'$\Delta t = %.3g\,\text{ps}$'%(self.delays[0]))
                else:
                    self.ax.set_title(r'$\Delta t = %.3g\,\text{ns}$'%(self.delays[0]/1000))
            else:
                if self.delays[0]<1:
                    self.ax.set_title(r'$\Delta t = %.3g\,\text{ps}$'%(self.delays[0]*1000))
                if self.delays[0]>=1000:
                    self.ax.set_title(r'$\Delta t = %.3g\,\text{µs}$'%(self.delays[0]/1000))
                else:
                    self.ax.set_title(r'$\Delta t = %.3g\,\text{ns}$'%self.delays[0])                

        if not self.scatterbar==False:
            rect = [(1/self.scatter[0][1])*10**4, (1/self.scatter[0][0])*10**4]
            self.ax.add_patch(patches.Rectangle((rect[0], self.ax.get_ylim()[0]), rect[1]-rect[0], self.ax.get_ylim()[1] - self.ax.get_ylim()[0], facecolor='white'))

        if self.zeroline==True:
            self.ax.axhline(y=0, color='k')
        if self.norm == True:
            self.ax.set_ylabel(r'norm. $\Delta A$')
        else:
            self.ax.set_ylabel(r'$\Delta A / 10^{-3}$')      
        if self.ylabel!='':
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
            self.fig.savefig('%s.svg'%(self.files[0][:self.files[0].find('.')]), transparent=True)
        plt.show()    
         