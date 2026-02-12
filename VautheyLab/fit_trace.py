import numpy as np
import matplotlib.pyplot as plt
import os
plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'style.mplstyle'))
from scipy.optimize import curve_fit
from scipy.signal import convolve as conv
from scipy.signal import fftconvolve
import statsmodels.api as sm
from VautheyLab.fit_functions import gaussian
np.seterr(divide='ignore')

class fitClass:
    def __init__(self):
        pass

    def mono_exp(self, t, tau, A, shift):
        IRF = self.irf
        IRF_int = np.interp(t, t + shift, IRF)
        # area normalize IRF to get rid off scaling artefacts
        #IRF_int = IRF_int / np.trapz(IRF_int, t)
        fit = A*conv(A*np.exp(-t/tau), IRF_int)/np.max(conv(A*np.exp(-t/tau), IRF_int))
        return fit[:len(t)] + self.bg
    
    def mono_exp_tail(self, t, tau, A):
        return A*np.exp(-t/tau) + self.bg

    def bi_exp(self, t, tau1, A1, tau2, A2, shift):
        IRF = self.irf
        IRF_int = np.interp(t, t + shift, IRF)
        fit = A1*conv(np.exp(-t/tau1), IRF_int)/np.max(conv(np.exp(-t/tau1), IRF_int)) + A2*conv(np.exp(-t/tau2), IRF_int)/np.max(conv(np.exp(-t/tau2), IRF_int))
        return fit[:len(t)] + self.bg
    
    def tri_exp(self, t, tau1, A1, tau2, A2, tau3, A3, shift):
        IRF = self.irf
        IRF_int = np.interp(t, t + shift, IRF)
        fit = A1*conv(np.exp(-t/tau1), IRF_int)/np.max(conv(np.exp(-t/tau1), IRF_int)) + A2*conv(np.exp(-t/tau2), IRF_int)/np.max(conv(np.exp(-t/tau2), IRF_int)) + A3*conv(np.exp(-t/tau3), IRF_int)/np.max(conv(np.exp(-t/tau3), IRF_int))
        return fit[:len(t)] + self.bg
    
class Fit_trace:
    def __init__(self, file, irf_file, cleanup=5, model='mono', p0=None, tail=None, notime=False, irfnotime=False, xlim=None,
                 figsize=(3.6*1.7, 3*1.7), tightlayout=True, savefig=False, yscale='log', title=None, Ulf=None, t_cuts=None, restrict=None):
        self.file = file
        self.irf_file = irf_file
        self.cleanup = cleanup
        self.model = model
        self.p0 = p0
        self.notime = notime
        self.irfnotime = irfnotime
        self.figsize = figsize
        self.tightlayout = tightlayout
        self.savefig = savefig
        self.yscale = yscale
        self.title = title
        self.tail = tail
        self.xlim = xlim

        if Ulf==None:
            self.t, self.counts = self.read_data(self.file, self.notime)
            self.t_IRF, self.counts_IRF = self.read_data(self.irf_file, self.irfnotime)
        else:
            self.t, self.counts = self.read_data_ulf(self.file)
            self.t_IRF, self.counts_IRF = self.read_data_ulf(self.irf_file)            

        # redefine time zero to be on the peak of the IRF
        self.t = self.t - self.t[self.counts_IRF==np.max(self.counts_IRF)][0]

        # set zero counts to one 
        #self.counts[self.counts==0] = 1
        self.counts += 1
        #self.t = self.t[self.counts!=0]
        #self.counts_IRF = self.counts_IRF[self.counts!=0]
        #self.counts = self.counts[self.counts!=0]

        # get background 
        self.bg = np.mean(self.counts[self.t<-3])
        print(f'background: {self.bg}')
        self.reject = np.std(self.counts[self.t<-3])
        self.restrict = restrict

        # get normalized counting errors
        self.sigma = np.sqrt(self.counts)

        # restrict
        if restrict!=None:
            # find start time
            t_start = self.t[0]
            i = 0
            while self.counts[i]<self.bg+restrict*self.reject:
                t_start = self.t[i]
                i += 1
            self.counts_IRF_for_fit = self.counts_IRF[self.t>t_start]
            self.counts_for_fit = self.counts[self.t>t_start]
            self.t_for_fit = self.t[self.t>t_start]
            self.sigma_for_fit = self.sigma[self.t>t_start]

            # find end time
            t_end = 0
            i = 0
            tcut = self.t[self.counts==np.max(self.counts)]
            while (self.counts[self.t>tcut][i]>self.bg+restrict*self.reject):
                t_end = self.t[self.t>tcut][i]
                i += 1
            self.counts_IRF_for_fit = self.counts_IRF[self.t<t_end]
            self.counts_for_fit = self.counts[self.t<t_end]
            self.t_for_fit = self.t[self.t<t_end]    
            self.sigma_for_fit = self.sigma[self.t<t_end] 

        # cut data if needed
        if t_cuts!=None:
            self.counts = self.counts[(self.t>t_cuts[0]) & (self.t<t_cuts[1])]
            self.counts_IRF = self.counts_IRF[(self.t>t_cuts[0]) & (self.t<t_cuts[1])]
            self.sigma = self.sigma[(self.t>t_cuts[0]) & (self.t<t_cuts[1])]
            self.t = self.t[(self.t>t_cuts[0]) & (self.t<t_cuts[1])]

        # put IRF on peak value of counts
        self.counts_IRF = self.counts_IRF*(np.max(self.counts)/np.max(self.counts_IRF))

        # clean up IRF
        self.counts_IRF[self.counts_IRF<self.cleanup] = 0

        # instanciate fitting class
        self.inst = fitClass()
        if self.restrict==None:
            self.inst.irf = self.counts_IRF 
        else:
            self.inst.irf = self.counts_IRF_for_fit
        self.inst.bg = self.bg

        # fit to model
        self.fit()

        # export fit
        self.export()

        # plot
        self.fig, self.ax = plt.subplots(ncols=1,nrows=3, sharex=True,
                        gridspec_kw={'height_ratios':[1,1,3]})

    def read_data(self, fname, notime):
        if notime==False:
            fh = open(fname, 'r')
            counts = []
            t = []
            for line in fh: 
                if line[0]!=' ' and line[0]!='*':
                    line = line.split()
                    if line[1]=='0':
                        line[1] = line[1].replace('0','1')
                    t.append(float(line[0]))
                    counts.append(float(line[1]))
        else:
            fh = open(fname, 'r')
            counts = []
            t = []
            for line in fh: 
                if line[0]!=' ' and line[0]!='*':
                    line = line.replace('0','1')
                    counts.append(float(line))
            # for normal length = corresponds to 50 ns in total 
            if len(counts)==4096:
                dt = 50/(len(counts) -1)
                t = [i*dt for i in range(len(counts))]
        return np.array(t), np.array(counts)
    
    def read_data_ulf(self, fname):
        data = np.loadtxt(fname, skiprows=10)
        I = data
        # ns per channel 
        t = np.array([i*0.0040 for i in range(len(I))])  
        return t, I

    def fit_IRF(self):
        fig, ax = plt.subplots(1,1, figsize=(5, 3.5))
        p, pcov = curve_fit(gaussian, self.t, self.counts_IRF, p0=[1e3, 0, 0.5])
        FWHM = p[-1]*1000
        errFWHM = (pcov[-1, -1])**(0.5)*1000
        ax.plot(self.t, self.counts_IRF, 'ok', alpha=0.3)
        tfine = np.linspace(-1.2, 1.2, 1000)
        ax.plot(tfine, gaussian(tfine, *p), '-r')
        ax.set_title(r'FHWHM = $(%.3g \pm %.3g)\,\text{ps}$'%(FWHM, errFWHM))
        ax.set_ylabel('counts')
        ax.set_xlabel(r'$\Delta t$ / ns')
        ax.set_xlim([-1.2, 1.2])
        fig.tight_layout()

    def fit(self):
        if self.model=='mono':
            self.fitted_model = self.inst.mono_exp
            if self.restrict==None:
                self.popt, self.pcov = curve_fit(self.inst.mono_exp, self.t, self.counts, p0=self.p0, sigma=self.sigma, absolute_sigma=True,
                                                bounds=((0, 0, -3), (np.inf, np.inf, 3)))
            else:
                self.popt, self.pcov = curve_fit(self.inst.mono_exp, self.t_for_fit, self.counts_for_fit, p0=self.p0, sigma=self.sigma_for_fit, absolute_sigma=True,
                                                bounds=((0, 0, -3), (np.inf, np.inf, 3)))        
                        
        if self.model=='mono_tail':
            self.fitted_model = self.inst.mono_exp_tail
            if self.restrict==None:
                x = self.t
                y = self.counts
                s = self.sigma
            else:
                x = self.t_for_fit
                y = self.counts_for_fit
                s = self.sigma_for_fit                
            yf = y[x>self.tail]
            sf = s[x>self.tail]
            xf = x[x>self.tail]
            self.popt, self.pcov = curve_fit(self.inst.mono_exp_tail, xf, yf, p0=self.p0, sigma=sf)  

        if self.model=='bi':
            self.fitted_model = self.inst.bi_exp
            if self.restrict==None:
                self.popt, self.pcov  = curve_fit(self.inst.bi_exp, self.t, self.counts, p0=self.p0, sigma=self.sigma, absolute_sigma=True,
                                                bounds=((0,-np.inf,0,-np.inf,-np.inf),(np.inf,np.inf,np.inf,np.inf,np.inf)))
            else:
                self.popt, self.pcov  = curve_fit(self.inst.bi_exp, self.t_for_fit, self.counts_for_fit, p0=self.p0, sigma=self.sigma_for_fit, absolute_sigma=True,
                                                bounds=((0,-np.inf,0,-np.inf,-np.inf),(np.inf,np.inf,np.inf,np.inf,np.inf)))                
            
        if self.model=='tri':
            self.fitted_model = self.inst.tri_exp
            if self.restrict==None:
                self.popt, self.pcov  = curve_fit(self.inst.tri_exp, self.t, self.counts, p0=self.p0, sigma=self.sigma,
                                                bounds=((0,-np.inf,0,-np.inf, 0,-np.inf, -np.inf),(np.inf,np.inf,np.inf,np.inf,np.inf, np.inf,np.inf)))  
            else:
                self.popt, self.pcov  = curve_fit(self.inst.tri_exp, self.t_for_fit, self.counts_for_fit, p0=self.p0, sigma=self.sigma_for_fit,
                                                bounds=((0,-np.inf,0,-np.inf, 0,-np.inf, -np.inf),(np.inf,np.inf,np.inf,np.inf,np.inf, np.inf,np.inf)))                  
        if self.restrict!=None:
            self.fitted = self.fitted_model(self.t_for_fit, *self.popt)
        else:
            self.fitted = self.fitted_model(self.t, *self.popt)
        return self.popt, self.pcov    

    def export(self):
        # model params 
        if self.model == 'mono' or self.model == 'mono_tail':
            header = 'model: ' + self.model + '\n'
            header += r'lifetime: $\tau = (%.4g \pm %.1g)\,\text{ns}$'%(self.popt[0], self.pcov[0,0]**0.5) + '\n' 
            header += '\n'
        if self.model == 'bi':
            header = 'model: ' + self.model + '\n'
            # percentage of species contributing to the decay
            per1 = np.abs(self.popt[1]) / (np.abs(self.popt[1]) + np.abs(self.popt[3]) )
            per2 = np.abs(self.popt[3]) / (np.abs(self.popt[1]) + np.abs(self.popt[3]) )
            header += r'''$\tau_1 = (%.3g \pm %.1g)\,\text{ns} \,\,\,\,\,\, \alpha_1 = %.2g$ \n
                    $\tau_2 = (%.3g \pm %.1g)\,\text{ns} \,\,\,\,\,\, \alpha_2 = %.2g$'''%(self.popt[0], self.pcov[0,0]**0.5, per1, self.popt[2], self.pcov[2,2]**0.5, per2)
        header += '\n'
        header += 'time / ns, counts, fit'
        if self.restrict!=None:
            np.savetxt(self.file[:self.file.find('.')] + '_export.txt', np.column_stack([self.t_for_fit, self.counts_for_fit, self.fitted]),
                    header=header, delimiter=',')
        else:
            np.savetxt(self.file[:self.file.find('.')] + '_export.txt', np.column_stack([self.t, self.counts, self.fitted]),
                    header=header, delimiter=',')            

    def show(self):
        self.fig.set_size_inches(self.figsize)

        # plot irf
        self.ax[2].plot(self.t, self.counts_IRF, '-k', alpha=0.4, linewidth=0.5, label=r'\text{IRF}')
        # plot counts
        self.ax[2].plot(self.t, self.counts, '.b', alpha=0.4, linewidth=1, label=r'\text{data}')

        # if restrict data to reject*sigma of bg show limits
        if self.restrict!=None:
            self.ax[2].axhline(y=self.bg + (self.reject * self.restrict), linestyle='--', color='k')
            self.ax[2].axhline(y=self.bg, linestyle='-', color='k')
            self.ax[2].plot(self.t_for_fit, self.counts_for_fit, '.g', alpha=0.4, linewidth=1, label=r'\text{data}')
            # look if irf fit or tail fit
            x = self.t_for_fit
            y = self.counts_for_fit
            s = self.sigma_for_fit
        else:
            # take all data
            x = self.t
            y = self.counts
            s = self.sigma

        #print('estimated colorshift of IRF: %.3g ns'%(x[y==np.max(y)] - x[self.counts_IRF==np.max(self.counts_IRF)]))
        #print('fitted color-shift of IRF: %.3g ns'%(self.popt[-1]))

        if '_tail' in self.model:
            y = y[x>self.tail]
            s = s[x>self.tail]
            x = x[x>self.tail]   
            self.ax[2].axvline(x=self.tail, linestyle='--', color='k', alpha=0.3)

        # plot fit
        self.ax[2].plot(x, self.fitted_model(x, *self.popt), '-r', linewidth=1)
        self.res = (y - self.fitted_model(x, *self.popt))/(s)
        # start taking residuals from maximum of peak + 0.2 ns
        self.res = self.res[x>x[y==np.max(y)]+0.2]
        x = x[x>x[y==np.max(y)]+0.2]
        self.chi2 = np.sum(self.res**2)/(len(x)-len(self.popt))
        print(f'chi2: {self.chi2}')
        self.ax[1].plot(x, self.res, '-r', linewidth=0.4, label=r'$\chi_{v}^2 = %.3g$'%(self.chi2))
        self.ax[1].axhline(y=0, color='k')
        self.ax[1].set_ylabel(r'$\text{residuals}$')
        self.ax[1].legend()
        # print results to title
        self.print_fit_params()
        # plot autocorrelation
        self.ax[0].plot(x, sm.tsa.acf(self.res, nlags=len(x)), '-g', linewidth=0.4)
        self.ax[0].axhline(y=0, color='k')
        self.ax[0].set_ylim([-0.15, 0.15])
        self.ax[0].set_ylabel(r'$\text{autocorr.}$')
        self.ax[2].set_xlabel(r'$\Delta t$ / ns')
        self.ax[2].set_ylabel(r'counts / a.u.')
        self.ax[2].set_yscale(self.yscale)
        if self.xlim !=None:
            self.ax[2].set_xlim(self.xlim)
        if self.title!=None:
            self.ax[0].set_title(self.title)
        if self.tightlayout == True:
            self.fig.tight_layout()
        if self.savefig == True: 
            self.fig.savefig('%s_fit.svg'%(self.file[:self.file.find('.')]), transparent=True)
        plt.show()

    def print_fit_params(self):
        # print results
        if self.model=='mono' or self.model=='mono_tail':
            self.ax[0].set_title(r'$\tau = (%.4g \pm %.1g)$ ns'%(self.popt[0], self.pcov[0,0]**0.5), size=12)
        if self.model=='bi':
            # fractional contribution (fi) of each decay time to the steadystate intensity --> see Lakowitch
            f1 = (self.popt[0]*np.abs(self.popt[1]))/(self.popt[0]*np.abs(self.popt[1]) + self.popt[2]*np.abs(self.popt[3]))
            f2 = (self.popt[2]*np.abs(self.popt[3]))/(self.popt[0]*np.abs(self.popt[1]) + self.popt[2]*np.abs(self.popt[3]))
            # percentage of species contributing to the decay
            per1 = np.abs(self.popt[1]) / (np.abs(self.popt[1]) + np.abs(self.popt[3]) )
            per2 = np.abs(self.popt[3]) / (np.abs(self.popt[1]) + np.abs(self.popt[3]) )
            self.ax[0].set_title(r'''$\tau_1 = (%.3g \pm %.1g)\,\text{ns} \,\,\,\,\,\, \alpha_1 = %.2g (%.2g) \,\,\,\,\,\,  f_1 = %.2g$ \\
                            $\tau_2 = (%.3g \pm %.1g)\,\text{ns} \,\,\,\,\,\, \alpha_2 = %.2g (%.2g) \,\,\,\,\,\,  f_2 = %.2g$'''%(self.popt[0], self.pcov[0,0]**0.5, self.popt[1], per1, f1, self.popt[2], self.pcov[2,2]**0.5, self.popt[3], per2, f2), size=12)
        if self.model=='tri': 
            f1 = (self.popt[0]*np.abs(self.popt[1]))/(self.popt[0]*np.abs(self.popt[1]) + self.popt[2]*np.abs(self.popt[3]) + self.popt[4]*np.abs(self.popt[5]))
            f2 = (self.popt[2]*np.abs(self.popt[3]))/(self.popt[0]*np.abs(self.popt[1]) + self.popt[2]*np.abs(self.popt[3]) + self.popt[4]*np.abs(self.popt[5]))
            f3 = (self.popt[4]*np.abs(self.popt[5]))/(self.popt[0]*np.abs(self.popt[1]) + self.popt[2]*np.abs(self.popt[3]) + self.popt[4]*np.abs(self.popt[5]))
            self.ax[0].set_title(r'''$\tau_1 = (%.3g \pm %.1g)\,\text{ns} \,\,\,\,\,\, \alpha_1 = %.2g \,\,\,\,\,\,  f_1 = %.2g$ \\
                            $\tau_2 = (%.3g \pm %.1g)\,\text{ns} \,\,\,\,\,\, \alpha_2 = %.2g \,\,\,\,\,\,  f_2 = %.2g$ \\
                            $\tau_3 = (%.3g \pm %.1g)\,\text{ns} \,\,\,\,\,\, \alpha_3 = %.2g \,\,\,\,\,\,  f_3 = %.2g$'''%(self.popt[0], self.pcov[0,0]**0.5, self.popt[1], f1, self.popt[2], self.pcov[2,2]**0.5, self.popt[3], f2,  self.popt[4], self.pcov[4,4]**0.5, self.popt[5], f3), size=12)
