from VautheyLab.standard import *
from VautheyLab.miscellaneous import rainbow, moving_average
from VautheyLab.fit_functions import *

def load_pdat(file):
    data = np.loadtxt(file, skiprows=1, delimiter=',')
    t = data[1:, 0]
    wl = data[0, 1:]
    dA = data[1:, 1:]  
    return t, wl, dA    

class Compare_Kinetics:
    def __init__(self, files, cuts=None, norm=False, normat=None, minnorm=False, ylim=None, xlim=None, xticks=True, xscale='linear',
                 yticks=True, wavelengths=None, experiment='femto', colors=None, export=False, MA=False, MA_npoints=10, tightlayout=True,
                 savefig=False, figsize=None, title=None, labels=[''], plot_now=True, alpha=1, fmt='.', markersize=5, fontsize=10, outside=False,
                 devide=None, shift=None,
                 zeroline=False, fit=None, ylabel='', IR=False):       
        # set relevant parameters
        self.figsize = figsize
        self.IR = IR
        # get files
        self.files = files
        # get wavelength cuts
        self.cuts = cuts    
        # get experiment
        self.experiment = experiment
        # devide traces by factor 
        self.devide = devide
        # decide whether to normalize or not
        self.norm = norm
        self.shift = shift
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
        self.xscale = xscale
        self.ylabel = ylabel
        # specify list of wavelengths
        self.wavelengths = wavelengths
        # if kK convert to nm
        if self.IR==False:
            if self.wavelengths[0]<100:
                for i in range(len(self.wavelengths)):
                    self.wavelengths[i] = (1/self.wavelengths[i])*10**4
        # take same if only one specified
        if len(self.wavelengths)==1 and len(self.files)>1:
            # use same range if only one cuts specified
            for i in range(len(self.files) - 1):
                self.wavelengths.append(self.wavelengths[0])        
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
        if colors==None:
            self.colors = rainbow(self.files)
        else:
            self.colors = colors
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
        self.zeroline=zeroline
        # if you want to fit --> follows tuple structure (model, p0)
        self.fit = fit
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
                dA = dA[(t>self.cuts[i][0])&(t<self.cuts[i][1]), :]
                t = t[(t>self.cuts[i][0])&(t<self.cuts[i][1])]
            # normalize
            if self.minnorm==False:
                if self.norm == True: 
                    if self.normat==None:
                        if self.IR==False:
                            dA = dA/np.max(dA[:, np.argmin(np.abs(wl - self.wavelengths[i]))])
                        else:
                            dA = dA/np.max(dA[:, np.argmin(np.abs(wn - self.wavelengths[i]))]) 
                    else:
                        if self.IR==False:
                            dA = dA/dA[np.argmin(np.abs(t - self.normat[0])), np.argmin(np.abs(wl - self.normat[1]))]
                        else:
                            dA = dA/dA[np.argmin(np.abs(t - self.normat[0])), np.argmin(np.abs(wn - self.normat[1]))]
            else:
                if self.norm == True: 
                    if self.normat==None:
                        dA = -1*dA/np.min(dA)
                    else:
                        dA = -1*dA/dA[np.argmin(np.abs(t - self.normat[0])), np.argmin(np.abs(wl - self.normat[1]))]
            # devide if wanted
            if self.devide!=None:
                dA = dA/self.devide[i]
            if self.shift!=None:
                dA = dA + self.shift[i]
            # append
            self.dA.append(dA)
            self.t.append(t)
            self.wl.append(wl)
            self.wn.append(wn)

    def plot(self):
        # go through all files
        for i in range(len(self.files)):
            x = self.t[i]
            if self.IR==False:
                y = self.dA[i][:, np.argmin(np.abs(self.wl[i] - self.wavelengths[i]))]  
            else:
                y = self.dA[i][:, np.argmin(np.abs(self.wn[i] - self.wavelengths[i]))]  
            # plot
            if self.MA==False:
                self.ax.plot(x, y, self.fmt, color=self.colors[i], alpha=self.alpha, label=self.labels[i], markersize=self.markersize)
            else:
                self.ax.plot(x, y, '-', color=self.colors[i], alpha=0.2)
                self.ax.plot(moving_average(x, self.MA_npoints), moving_average(y, self.MA_npoints), '-', color=self.colors[i], label=self.labels[i])

    def do_fit(self):
        # go through all files
        for i in range(len(self.files)):
            x = self.t[i]
            xfine = np.linspace(x[0], x[-1], 1000)
            if self.IR==False:
                y = self.dA[i][:, np.argmin(np.abs(self.wl[i] - self.wavelengths[i]))] 
            else:
                y = self.dA[i][:, np.argmin(np.abs(self.wn[i] - self.wavelengths[i]))] 
            p0 = self.fit[i][1] 

            # simple monoexponential fit
            if self.fit[i][0]=='mono':
                popt, pcov = curve_fit(mono_exp, x, y, p0=p0)
                tau = popt[1]
                err_tau = pcov[1,1]**(0.5)
                if self.experiment=='femto':
                    lab = r'mono exp. fit with $\tau = (%.4g \pm %.1g)$ ps'%(tau, err_tau)
                else:
                    lab = r'mono exp. fit with $\tau = (%.4g \pm %.1g)$ ns'%(tau, err_tau)
                 # plot the fit
                self.ax.plot(xfine, mono_exp(xfine, *popt), '-', color=self.colors[i], label=lab, linewidth=2.5)

            # monoexponential with background
            if self.fit[i][0]=='mono_bg':
                popt, pcov = curve_fit(mono_exp_with_bg, x, y, p0=p0)
                tau = popt[1]
                err_tau = pcov[1,1]**(0.5)
                if self.experiment=='femto':
                    lab = r'mono exp. fit with $\tau = (%.4g \pm %.1g)$ ps'%(tau, err_tau)
                else:
                    lab = r'mono exp. fit with $\tau = (%.4g \pm %.1g)$ ns'%(tau, err_tau)
                 # plot the fit
                self.ax.plot(xfine, mono_exp_with_bg(xfine, *popt), '-', color=self.colors[i], label=lab, linewidth=2.5)

            # simple bi-exponential fit
            if self.fit[i][0]=='bi':
                popt, pcov = curve_fit(bi_exp, x, y, p0=p0)
                tau1 = popt[1]
                err_tau1 = pcov[1,1]**(0.5)
                tau2 = popt[3]
                err_tau2 = pcov[3,3]**(0.5)
                f1 = (popt[1]*np.abs(popt[0]))/(popt[1]*np.abs(popt[0]) + popt[3]*np.abs(popt[2]))
                f2 = (popt[3]*np.abs(popt[2]))/(popt[1]*np.abs(popt[0]) + popt[3]*np.abs(popt[2]))
                if popt[0]<0:
                    print('tau1 is a rise!')
                if popt[2]<0:
                    print('tau2 is a rise!')
                if self.experiment=='femto':
                    lab = r'bi exp. fit with $\tau_1 = (%.4g \pm %.1g)$ ps (%.2g), $\tau_2 = (%.4g \pm %.1g)$ ps (%.2g)'%(tau1, err_tau1, f1, tau2, err_tau2, f2)   
                else:
                    lab = r'bi exp. fit with $\tau_1 = (%.4g \pm %.1g)$ ns (%.2g), $\tau_2 = (%.4g \pm %.1g)$ ns (%.2g)'%(tau1, err_tau1, f1, tau2, err_tau2, f2)          
                # plot the fit
                self.ax.plot(xfine, bi_exp(xfine, *popt), '-', color=self.colors[i], label=lab, linewidth=2.5)

    def show(self):
        if self.figsize!=None:
            self.fig.set_size_inches(self.figsize)
        
        if self.experiment=='femto':
            self.ax.set_xlabel(r'$\Delta t / $ ps')
        elif self.experiment=='nano':
            self.ax.set_xlabel(r'$\Delta t / $ ns')        

        if self.zeroline==True:
            self.ax.axhline(y=0, color='k')

        self.ax.set_xscale(self.xscale)
        if self.outside==True:
            self.ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=self.fontsize)
        else:
            self.ax.legend(fontsize=self.fontsize)
        self.ax.axhline(y=0, color='k')
        if self.norm == True:
            self.ax.set_ylabel(r'norm. $\Delta A$')
            if len(set(self.wavelengths))==1:
                if self.IR==False:
                    self.ax.set_ylabel(r'norm. $\Delta A$ @ %.3g nm'%(self.wavelengths[0]))
                else:
                   self.ax.set_ylabel(r'norm. $\Delta A$ @ %i cm$^{-1}$'%(round(self.wavelengths[0]))) 
        else:
            self.ax.set_ylabel(r'$\Delta A / 10^{-3}$ @ %.3g nm'%(self.wavelengths[0]))   
        if not self.ylabel=='':
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
