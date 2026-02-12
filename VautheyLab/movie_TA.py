import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, writers
from VautheyLab.miscellaneous import moving_average
import os
plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'style.mplstyle'))
np.seterr(divide='ignore')

def load_pdat(file):
    data = np.loadtxt(file, skiprows=1, delimiter=',')
    t = data[1:, 0]
    wl = data[0, 1:]
    dA = data[1:, 1:]  
    return t, wl, dA    

class Movie:
    # initalize class
    def __init__(self, files, units='wn', inv=True, t_cuts=None, wl_cuts=None, norm=False, normat=None, minnorm=False, scatter=None, normall=False,
                 ylim=None, xlim=None, xticks=True, yticks=True, experiment='femto', colors=['r'], export=False, MA=[False], labels=[''],
                 MA_npoints=[10], tightlayout=True, figsize=(8,5), title=None, before=True, time=5, movname='movie.mp4', ylabel=None, outside=False,
                 steady_state=None, devide=None, IR=False, vlines=None, secax=True, merged=False):     
        # set relevant parameters
        self.figsize = figsize
        # get files
        self.files = files
        self.IR = IR
        # get cuts
        self.t_cuts = t_cuts
        self.wl_cuts = wl_cuts
        # get scatter 
        self.scatter = scatter
        self.secax = secax
        # get experiment
        self.experiment = experiment
        # decide whether to normalize or not
        self.norm = norm
        # if you want to norm on GSB
        self.minnorm = minnorm
        self.vlines = vlines
        # to normalize at a specific wl
        self.normat = normat
        # specify x and y limits for plot
        self.ylim = ylim
        self.xlim = xlim
        # specify xticks
        self.devide = devide
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
        self.colors = colors
        self.movname = movname
        self.ylabel = ylabel
        self.outside = outside
        self.normall = normall
        self.merged = merged
        # follows the tuple structure (file, wl_cuts, scale, color, label)
        self.steady_state = steady_state
        # treat labels
        self.labels=labels
        if self.labels==[''] and len(self.files)>1:
            for i in range(len(self.files) - 1):
                self.labels.append('')
        # moving average?
        self.MA = MA
        if len(self.MA)==1 and len(self.files)>1:
            for i in range(len(self.files)):
                self.MA.append(self.MA[0])
        self.MA_npoints = MA_npoints
        # tight layout
        self.tightlayout = tightlayout
        self.before = before
        self.time = time
        # intiialize figure
        self.fig, self.ax = plt.subplots(1,1,figsize=self.figsize)
        # read data
        self.read_data()

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
                if len(self.scatter[i])==2:
                    dA[:, (wl >= self.scatter[i][0]) & (wl <= self.scatter[i][1])] = np.nan
                if len(self.scatter[i])==4:
                    dA[:, (wl >= self.scatter[i][0]) & (wl <= self.scatter[i][1])] = np.nan
                    dA[:, (wl >= self.scatter[i][2]) & (wl <= self.scatter[i][3])] = np.nan
            # normalize
            if self.minnorm==False:
                if self.norm == True: 
                    if self.normat==None:
                        # normalize at earliest delay
                        dA = dA/np.nanmax(dA[0,:])
                    else:
                        if self.MA[i]==False:
                            dA = dA/dA[np.argmin(np.abs(t - self.normat[0])), np.argmin(np.abs(wl - self.normat[1]))]
                        else:
                            # normalize on mov av spectrum
                            wl_av = moving_average(wl, self.MA_npoints[i])
                            dA_av = moving_average(dA[np.argmin(np.abs(t - self.normat[0])),:], self.MA_npoints[i])
                            dA = dA/dA_av[np.argmin(np.abs(wl_av - self.normat[1]))]
            else:
                if self.norm == True: 
                    if self.normat==None:
                        dA = -1*dA/np.nanmin(dA[0,:])
                    else:
                        if self.MA[i]==False:
                            dA = dA/dA[np.argmin(np.abs(t - self.normat[0])), np.argmin(np.abs(wl - self.normat[1]))]
                        else:
                            # normalize on mov av spectrum
                            wl_av = moving_average(wl, self.MA_npoints[i])
                            dA_av = moving_average(dA[np.argmin(np.abs(t - self.normat[0])),:], self.MA_npoints[i])
                            dA = -1*dA/dA_av[np.argmin(np.abs(wl_av - self.normat[1]))]
            # append
            # devide if wanted
            if self.devide!=None:
                dA = dA/self.devide[i]
            # if nano cut size and only take every 5th spectrum
            if self.experiment=='nano' or self.merged==True:
                dA = dA[0::5,:]
                t = t[0::5]
            self.dA.append(dA)
            self.t.append(t)
            self.wl.append(wl)
            self.wn.append(wn)

    def animate(self, i):
        self.ax.clear()
        print("%i / %i"%(i, len(self.t[0])))
        # go through all files and plot
        for j in range(len(self.files)):
            if self.units == 'wl':
                x = self.wl[j]
            else:
                x = self.wn[j]
            if self.before==True:
                before = self.dA[j][self.t[j]<self.t[j][i],:]
                for k in range(len(before)):
                    self.ax.plot(x, before[k], '-', color=self.colors[j], alpha=0.01)
            y = self.dA[j][i, :]
            if self.MA[j]==False:
                if self.normall==False:
                    self.ax.plot(x, y, '-', color=self.colors[j], label=self.labels[j])
                else:
                    self.ax.plot(x, y/np.nanmax(y), '-', color=self.colors[j], label=self.labels[j])
            else:
                if self.normall==True:
                    norm = np.nanmax(moving_average(y, self.MA_npoints[j]))
                    self.ax.plot(x, y/norm, '-', color=self.colors[j], alpha=0.2)
                    self.ax.plot(moving_average(x, self.MA_npoints[j]), moving_average(y, self.MA_npoints[j])/norm, '-', color=self.colors[j], label=self.labels[j])
                else:
                    self.ax.plot(x, y, '-', color=self.colors[j], alpha=0.2)
                    self.ax.plot(moving_average(x, self.MA_npoints[j]), moving_average(y, self.MA_npoints[j]), '-', color=self.colors[j], label=self.labels[j])                    
            if self.experiment=='nano':
                if self.t[j][i]<1:
                    self.ax.set_title(r'$\Delta t = %.3g\,\text{ps}$'%(self.t[j][i]*1000))
                if self.t[j][i]>=1000:
                    self.ax.set_title(r'$\Delta t = %.3g\,\text{µs}$'%(self.t[j][i]/1000))
                else:
                    self.ax.set_title(r'$\Delta t = %.3g\,\text{ns}$'%self.t[j][i])
            else:
                if np.abs(self.t[j][i])<1:
                    self.ax.set_title(r'$\Delta t = %.3g\,\text{fs}$'%(self.t[j][i]*1000))
                elif 1<=np.abs(self.t[j][i])<1000:
                    self.ax.set_title(r'$\Delta t = %.3g\,\text{ps}$'%(self.t[j][i]))
                elif self.t[j][i]>=1000 and self.t[j][i]<1e6:
                    self.ax.set_title(r'$\Delta t = %.3g\,\text{ns}$'%(self.t[j][i]/1000))
                else:
                    self.ax.set_title(r'$\Delta t = %.3g\,\text{µs}$'%(self.t[j][i]/10**6))

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

        # stylistic stuff
        if self.units=='wn':
            if self.IR==False:
                self.ax.set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
                if self.secax==True:
                    axsec = self.ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
                    axsec.set_xlabel(r'$\lambda / $ nm') 
                    #axsec.set_xticks([500, 600, 700, 800, 900, 1200])
                self.ax.invert_xaxis()
            else:
                self.ax.set_xlabel(r'$\tilde{\nu} / \text{cm}^{-1}$')

        if self.units=='wl':
            self.ax.set_xlabel(r'$\lambda / $ nm')
            self.ax.set_xscale('function', functions=(lambda x: 1/x, lambda x: 1/x))

        if self.vlines!=None:
            for v in range(len(self.vlines)):
                self.ax.axvline(x=self.vlines[v], linestyle='--', color='k')

        if self.outside==True:
            self.ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
        else:
            self.ax.legend(loc='upper right')
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
        if self.tightlayout == True:
            self.fig.tight_layout()

    def render(self):
        anim = FuncAnimation(self.fig, func=self.animate, frames=len(self.t[0]), interval=1)
        Writer = writers['ffmpeg']
        writer = Writer(fps=round(len(self.t[0])/self.time), metadata={'artist': 'Me'}, bitrate=2500)
        anim.save(self.movname, writer)
        plt.show()