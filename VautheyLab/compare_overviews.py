from VautheyLab.standard import *
from VautheyLab.miscellaneous import *

def load_pdat(file):
    data = np.loadtxt(file, skiprows=1, delimiter=',')
    t = data[1:, 0]
    wl = data[0, 1:]
    dA = data[1:, 1:]  
    return t, wl, dA    

class Compare_Overviews:
    def __init__(self, files, delays=[None], legend=True, scatcol=None,
                 units='wn', inv=True, colors=None, steady_state=None, 
                 t_cuts=None, wl_cuts=None, horizontal=True, ylabel=False,
                 norm=False, normat=None, minnorm=False, 
                 scatter=None, ylim=None, xlim=None, xticks=True,
                 yticks=True, experiment='femto', export=False, MA=[False], 
                 MA_npoints=[10], tightlayout=True, yscale='linear', 
                 savefig=False, figsize=None, titles=None, zeroline=False,
                 devide=None):          
        # set relevant parameters
        self.figsize = figsize
        # get files
        self.files = files
        self.devide = devide
        self.legend = legend
        self.ylabel = ylabel
        self.scatcol = scatcol
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
        self.steady_state = steady_state
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
        # delays 
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
        self.yscale = yscale
        self.zeroline = zeroline
        self.tightlayout = tightlayout
        self.titles = titles
        self.savefig = savefig
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
        if horizontal==True:
            self.fig, self.ax = plt.subplots(1, len(self.files), sharey=True)
        else:
            self.fig, self.ax = plt.subplots(len(self.files), 1, sharex=True)
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
            if '.pdat' in self.files[i]:
                t, wl, dA = load_pdat(self.files[i])
            wn = (1/wl)*10**4
            # cut data 
            if self.t_cuts!=None:
                dA = dA[(t>self.t_cuts[0])&(t<self.t_cuts[1]), :]
                t = t[(t>self.t_cuts[0])&(t<self.t_cuts[1])]
            if self.wl_cuts!=None:
                dA = dA[:, (wl>self.wl_cuts[0])&(wl<self.wl_cuts[1])]
                wn = wn[(wl>self.wl_cuts[0])&(wl<self.wl_cuts[1])]
                wl = wl[(wl>self.wl_cuts[0])&(wl<self.wl_cuts[1])]
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
            if self.devide!=None:
                dA = dA/self.devide[i]
            self.dA.append(dA)
            self.t.append(t)
            self.wl.append(wl)
            self.wn.append(wn)

    def plot(self):
        for j in range(len(self.files)):
            # plot steady state spectra
            if self.steady_state!=None:
                for k in range(len(self.steady_state[j])):
                    # read data
                    data = np.loadtxt(self.steady_state[j][k][0], delimiter=',', skiprows=1)
                    s_wl = data[:,0]
                    s_wn = data[:,1]
                    s = data[:,2]
                    # multiply intensity by lambda^4 if emission
                    if 'em' in self.steady_state[j][k][0]:
                        s = s*s_wl**4
                        print(self.steady_state[j][k][0], ' has been multiplied by lambda^4 to convert to SE.')
                    s = s[(s_wl>self.steady_state[j][k][1][0]) & (s_wl<self.steady_state[j][k][1][1])]
                    s_wn = s_wn[(s_wl>self.steady_state[j][k][1][0]) & (s_wl<self.steady_state[j][k][1][1])]
                    s_wl = s_wl[(s_wl>self.steady_state[j][k][1][0]) & (s_wl<self.steady_state[j][k][1][1])]
                    s = self.steady_state[j][k][2]*s/np.nanmax(s)
                    if self.units == 'wl':
                        self.ax[j].fill_between(s_wl, 0, s, color=self.steady_state[j][k][3], label=self.steady_state[j][k][4], alpha=0.1)
                    else:
                        self.ax[j].fill_between(s_wn, 0, s, color=self.steady_state[j][k][3], label=self.steady_state[j][k][4], alpha=0.1) 

            if self.scatcol!=None:
                self.ax[j].axvspan(nm_to_kk(self.scatter[j][1]), nm_to_kk(self.scatter[j][0]), color=self.scatcol[j], alpha=0.08)

            # choose unit
            if self.units == 'wn':
                x = self.wn[j]
            else: 
                x = self.wl[j]
            # got through all delays and plot
            for i in range(len(self.delays)):
                y = self.dA[j][np.argmin(np.abs(self.t[j] - self.delays[i])), :]
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
                if self.MA[j]==False:
                    self.ax[j].plot(x, y, '-', color=self.colors[i], label=lab)
                else:
                    self.ax[j].plot(x, y, '-', color=self.colors[i], alpha=0.2)
                    self.ax[j].plot(moving_average(x, self.MA_npoints[j]), moving_average(y, self.MA_npoints[j]), '-', color=self.colors[i], label=lab)
 
    def show(self):
        if self.figsize!=None:
            self.fig.set_size_inches(self.figsize)

        for i in range(len(self.files)):
            if self.units=='wn':
                self.ax[i].set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
                if self.inv==True:
                    self.ax[i].invert_xaxis()    
                axsec = self.ax[i].secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
                axsec.set_xlabel(r'$\lambda / $ nm') 

            if self.units=='wl':
                self.ax[i].set_xlabel(r'$\lambda / $ nm')

            if self.legend==True:
                self.ax[i].legend(loc='center left')
            self.ax[i].axhline(y=0, color='k')
            if self.norm == True:
                self.ax[0].set_ylabel(r'norm. $\Delta A$')
            else:
                self.ax[0].set_ylabel(r'$\Delta A / 10^{-3}$')                
            if self.yticks == False:
                self.ax[0].set_yticks([])
            if self.xticks == False:
                self.ax[0].set_xticks([])
            if self.ylim != None:
                self.ax[0].set_ylim(self.ylim)
            if self.xlim != None: 
                self.ax[i].set_xlim(self.xlim)
            if self.titles != None:
                self.ax[i].set_title(self.titles[i])
            if self.ylabel != False:
                self.ax[0].set_ylabel(self.ylabel)
            if self.tightlayout == True:
                self.fig.tight_layout()
            if self.savefig == True: 
                self.fig.savefig('%s.svg'%(self.files[i][:self.files[i].find('.')]), transparent=True)
        plt.show()