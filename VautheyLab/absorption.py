from VautheyLab.miscellaneous import moving_average, solvs, find_index, nm_to_kk
from VautheyLab.standard import *
from VautheyLab.fit_functions import lorentzian

# function for reading theoretical IR spectrum
def get_peaks(fname, UV=False):
    x = []
    y = []
    with open(fname, "r") as file: 
        for line in file:
            # only read lines starting with #
            if line.startswith('#'):
                # split line
                line = line.split()
                # if second entry can be converted to a number add it
                try:
                    x.append(float(line[1]))
                    if UV==False:
                        y.append(float(line[2]))
                    else:
                        y.append(float(line[3]))
                except ValueError:
                    continue
    return np.array(x), np.array(y)

# Class to plot absorption spectra
class Absorption:
    # initialize class
    def __init__(self, files, figsize=None, zeroline=True, devide=None, which=None,
                 cuts=None, norm=None, ylim=None, xlim=None, xticks=True, ylabel=None,
                 yticks=True, fill=False, inv=False, TDM=False, colors=['r'], alpha=[],
                 labels=[''], title='', units='wn', conc=None, pathlegnths=None, plot_now=True,
                 eps=None, export=False, MA=False, MA_npoints=10, tightlayout=True, savefig=False, IR=False,
                 linestyle=None, baseline=None, normat=None, baselineat=None, waterfall=False, outside=False, fontsize=10, legend=True, 
                 secaxticks=None, dipUlf=None, yscale=None):
        
        # set relevant parameters
        self.figsize = figsize
        # get files
        self.files = files
        # get wavelength cuts
        self.cuts = cuts
        if self.cuts!=None:
            if len(self.cuts)==1 and len(self.files)>1:
                # use same range if only one cuts specified
                for i in range(len(self.files) - 1):
                    self.cuts.append(self.cuts[0])
        # decide whether to normalize or not
        self.norm = norm
        # to normalize at a specific wl
        self.normat = normat
        # specify x and y limits for plot
        self.ylim = ylim
        self.xlim = xlim
        # specify xticks
        self.xticks = xticks
        self.yticks = yticks
        self.zeroline = zeroline
        self.devide = devide
        # decide whether to have spectra filled
        self.fill = fill
        # specify whether you want to invert the x-axis when using kK or not
        self.inv = inv
        # transition dipole moment representation
        self.TDM = TDM
        # specify list of colors
        self.colors = colors
        # specify list of labels
        self.labels = labels
        # sepcify title of your plot
        self.title = title
        # sepecify your units
        self.units = units
        self.outside = outside
        self.legend = legend
        self.fontsize = fontsize
        self.ylabel = ylabel
        self.IR = IR
        # specify which files to plot
        self.which = which 
        if self.which==None:
            self.which = [i for i in range(len(self.files))]

        # alpha
        self.alpha = alpha
        if self.alpha == []:
            for i in range(len(self.files)):
                self.alpha.append(1)

        # list of concentrations to calculate eps
        self.conc = conc
        self.plength = pathlegnths

        # list of epsilons at certains wavelengths
        self.eps = eps

        # whether to export data or not
        self.export = export

        # moving average?
        self.MA = MA
        self.MA_npoints = MA_npoints

        # baseline file
        self.baseline = baseline
        self.baselineat = baselineat

        # tight layout
        self.tightlayout = tightlayout
        self.savefig = savefig

        # for waterfall plot
        self.waterfall = waterfall
        self.secaxticks = secaxticks

        # to correct for dip at 800 nm for cary5000
        self.dipUlf = dipUlf
        self.yscale = yscale

        # line style 
        if linestyle==None:
            self.linestyle = ['-' for _ in range(len(self.files))]
        else:
            self.linestyle = linestyle

        # if no labels specified
        if self.labels==[''] and len(self.files)>1:
            self.labels = ['' for _ in range(len(self.files))]

        # if only one baseline use same for all files
        if self.baseline!=None:
            if len(self.baseline)==1 and len(self.files)>1:
                self.baseline = [baseline[0] for _ in range(len(self.files))]

        # intiialize figure
        self.fig, self.ax = plt.subplots(1,1)
        # read data
        self.read_data()
        # plot data
        if plot_now==True:
            self.plot()

    # read in data
    def read_data(self):
        # initialize list of absorbances
        self.A = []
        self.wl = []
        self.wn = []
        # go through all files
        for i in range(len(self.files)):
            if self.IR==False:
                # read file
                if not '.txt' in self.files[i]:
                    data = np.loadtxt('%s'%self.files[i], delimiter=',', usecols=[0,1], skiprows=2)
                else:
                    data = np.loadtxt('%s'%self.files[i], delimiter=',', usecols=[0,2], skiprows=1)
                
                # get wavelength
                if data[:,0][0]<100:
                    self.wl.append((1/data[:,0])*10**4)
                    self.wn.append((1/self.wl[i])*10**4)
                else:
                    self.wl.append(data[:,0])
                    self.wn.append((1/data[:,0])*10**4)
                # get absorbance
                self.A.append(data[:,1]) 
            else:
                if not '.txt' in self.files[i]:
                    data = np.loadtxt('%s'%self.files[i], delimiter=',')
                    self.wn.append(data[:,0])
                    self.wl.append((1/data[:,0])*10**(4))
                    self.A.append(data[:,1])
                else:
                    data = np.loadtxt('%s'%self.files[i], delimiter=',')
                    self.wl.append(data[:,0])
                    self.wn.append(data[:,1])
                    self.A.append(data[:,2])     

            # correct for dip by applying offset to NIR part
            if self.dipUlf!= None:
                a = self.A[i][self.wl[i]>=800][-1]
                b = self.A[i][self.wl[i]<800][0]
                dip = b/a
                #print(dip, a, b, self.wl[i][self.wl[i]>=800][-1], self.wl[i][self.wl[i]<800][0])
                self.A[i] = np.concatenate([self.A[i][self.wl[i]>=800]*dip, self.A[i][self.wl[i]<800]])               
            
            # correct baseline 
            if self.baseline!=None:
                if self.baseline[i]!=None:
                    bdata = np.loadtxt(self.baseline[i], delimiter=',', usecols=[0,1], skiprows=2)
                    self.A[i] = self.A[i] - bdata[:,1]
            if self.baselineat!=None and self.baseline==None:
                self.A[i] = self.A[i] - self.A[i][np.argmin(np.abs(self.wl[i] - self.baselineat))]


            # cut data 
            if self.cuts!=None:
                if self.IR==False:
                    self.A[i] = self.A[i][(self.wl[i]>self.cuts[i][0])&(self.wl[i]<self.cuts[i][1])]
                    self.wl[i] = self.wl[i][(self.wl[i]>self.cuts[i][0])&(self.wl[i]<self.cuts[i][1])]
                    self.wn[i] = (1/self.wl[i])*10**4
                else:
                    self.A[i] = self.A[i][(self.wn[i]>self.cuts[i][0])&(self.wn[i]<self.cuts[i][1])]
                    self.wn[i] = self.wn[i][(self.wn[i]>self.cuts[i][0])&(self.wn[i]<self.cuts[i][1])]      
                    self.wl[i] = (1/self.wn[i])*10**(4)   

            # normalize
            if self.norm == True: 
                if self.normat==None:
                    self.A[i] = self.A[i]/np.nanmax(self.A[i])
                else:
                    self.A[i] = self.A[i]/self.A[i][np.argmin(np.abs(self.wl[i] - self.normat[i]))]
            # devide if wanted
            if self.devide!=None:
                self.A[i] = self.A[i]/self.devide[i]
            # export if wanted
            if self.export==True:
                if self.IR==False:
                    if self.norm == True: 
                        head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  norm. absorbance'
                    else:
                        if self.conc == None:
                            head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  absorbance'
                            np.savetxt('%s.txt'%(self.files[i][:self.files[i].find('.')]), 
                                        np.column_stack([self.wl[i], self.wn[i], self.A[i]]),
                                        header=head, delimiter=',')
                        else:
                            head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  extinction coefficient / M-1 s-1'
                            np.savetxt('%s.txt'%(self.files[i][:self.files[i].find('.')]), 
                                        np.column_stack([self.wl[i], self.wn[i], self.A[i]/(self.conc[i] * self.plength[i])]),
                                        header=head, delimiter=',')                            
                else:
                    if self.norm == True: 
                        head = r'wavenlength / µm,    wavenumber / cm-1,  norm. absorbance'
                    else:
                        if self.conc == None:
                            head = r'wavenlength / µm,    wavenumber / cm-1,  absorbance'
                            np.savetxt('%s.txt'%(self.files[i][:self.files[i].find('.')]), 
                                        np.column_stack([self.wl[i], self.wn[i], self.A[i]]),
                                        header=head, delimiter=',')
                        else:
                            head = r'wavenlength / µm,    wavenumber / cm-1,  extinction coefficient / M-1 s-1'
                            np.savetxt('%s.txt'%(self.files[i][:self.files[i].find('.')]), 
                                        np.column_stack([self.wl[i], self.wn[i], self.A[i]/(self.conc[i] * self.plength[i])]),
                                        header=head, delimiter=',')   

    # find maximum
    def find_max(self, file_index):
        Amax = np.max(self.A[file_index])
        wlmax = self.wl[file_index][self.A[file_index] == np.max(self.A[file_index])]
        wnmax = self.wn[file_index][self.A[file_index] == np.max(self.A[file_index])]
        print("%s --> Amax = %.3g at wl = %.3g nm / wn = %.3g kK"%(self.files[file_index], Amax, wlmax, wnmax))
        if self.conc!=None:
            eps = Amax / (self.plength[file_index] * self.conc[file_index])
            print(r"%s --> \varepsilon_{\text{%.3g nm}} = \SI{%.3g}{\text{M}^{-1}\,\text{cm}^{-1}} at wl = %.3g nm / wn = %.3g kK"%(self.files[file_index], wlmax, eps, wlmax, wnmax))
        return Amax, wlmax, wnmax


    def epsilon_regression(self, M, m, Vstock, Vcuv, Vadd, l, wl):
        # M - molar mass in g/mol
        # m - weighted masses in mg 
        # Vstock - volume used to dissolve weighted masses in mL 
        # Vcuv - volume [mL] of solution in which Vadd [uL] of stock is added
        # l - pathlength [cm]
        # wl - wavelength at which to analyze absorbance
        
        # get mass in g 
        m = m/1000
        # get volumia in L 
        Vstock = Vstock/1000
        Vadd = Vadd/10**6
        Vcuv = Vcuv/1000

        # calculate stock concentrations 
        cstock = m / (Vstock * M)
        # error for mass = 0.01 mg
        err_m = 0.01e-3
        # error of Vstock = 0.05 mL
        err_Vstock = 0.05e-3
        dcdm = 1 / (Vstock * M)
        dcdV = -m / (Vstock**2 * M) 
        cstockerr = np.sqrt( (dcdm * err_m)**2 + (dcdV * err_Vstock)**2)

        # print stock concentrations with error
        print('\n Stock concentrations: ')
        for i in range(len(cstock)):
            print(f'\n c = ({cstock[i]*1000 : .3f} +- {cstockerr[i]*1000 : .3f}) mM')

        # calculate diluted concentration 
        Vfull = Vcuv + Vadd
        err_Vcuv = (Vcuv/Vstock)*err_Vstock
        err_Vadd = 1e-6
        err_Vfull = np.sqrt(err_Vcuv**2 + err_Vadd**2)
        c = (cstock * Vadd) / Vfull 
        dcdc = Vadd/Vfull 
        dcdVa = cstock/Vfull 
        dcdVc = -(cstock * Vadd) / (Vfull **2)
        cerr = np.sqrt( (dcdc * cstockerr)**2 + (dcdVa * err_Vadd)**2 + (dcdVc * err_Vfull)**2 )

        # print concentration in cuvette with error
        print('\n Cuvette concentrations: ')
        for i in range(len(cstock)):
            print(f'\n c = ({c[i]*10**6 : .3f} +- {cerr[i]*10**6 : .3f}) µM')

        # do fit 
        fig, ax = plt.subplots(1,1, figsize=(5, 3.5))
        A = []
        for i in range(len(self.files)):
            Aat = self.A[i][find_index(self.wl[i], wl)]
            A.append(Aat)
        ax.errorbar(x=A, y=c*10**6, yerr=cerr*10**6, ecolor='k', capsize=3, fmt='or')
        # fit to lamber beer law 
        p, pcov = curve_fit(lambda A, eps: A/(eps*l), A, c, sigma=cerr, absolute_sigma=True)
        eps = p[0]
        err_eps = pcov[0,0]**(0.5)
        Afine = np.linspace(A[0]-0.1, A[-1]+0.1, 100)
        ax.plot(Afine, Afine/(eps*10**-6 * l), '--k')
        s = f"{eps:.1e}"  
        power = int(s.split('e')[1])
        ax.set_title(r'$\varepsilon_{%i \text{nm}} = (%.2g \pm %.2g) \times 10^%i\,\text{M}^{-1}\,\text{cm}^{-1}$'%(wl, eps*10**-power, err_eps*10**-power, power))
        ax.set_xlabel(r'$A_{%i \text{nm}}$'%(int(round(wl))))
        self.ax.axvline(x = wl, color='k', linestyle='--', alpha=0.3)
        ax.set_ylabel(r'$c$ / µM')
        fig.tight_layout()
        fig.savefig('eps_regression.svg', transparent=True)
        print("")
        print(r"$\varepsilon_{%i \text{nm}} = (%.2g \pm %.2g) \times 10^%i\,\text{M}^{-1}\,\text{cm}^{-1}$"%(wl, eps*10**-power, err_eps*10**-power, power))
        print(f"relative error = {err_eps*100/eps : .2f} percent")
        print("")
 

    # function to calculate concentration of solution
    # need to procvide eps = [ (eps_at_wl, wl) ] at intialization
    def get_concentration(self, file_index):
        wl_look = self.eps[file_index][1]
        conc = self.A[file_index][[np.argmin(np.abs(self.wl[file_index] - wl_look))]] / (self.eps[file_index][0] * self.plength[file_index])
        if 10**(-3) <= conc < 1:
            print("%s -> c = %.3g mM"%(self.files[file_index], conc*10**3))
        elif 10**(-6) <= conc < 10**(-3):
            print("%s -> c = %.3g µM"%(self.files[file_index], conc*10**6))
        else:
            print("%s -> c = %.3g M"%(self.files[file_index], conc))
        return conc
    
    # function to calculate oscillator strength and krad
    # need to provide: limits in nm over which to calculate the band integral
    # refractive index n, main transition freq. in kK
    def calc_oscillator_strength(self, limits, n, nu0, file_index):
        # plot integral region
        x = self.wn[file_index][ (self.wl[file_index] > limits[0]) & (self.wl[file_index] < limits[1]) ]
        y = self.A[file_index][ (self.wl[file_index] > limits[0]) & (self.wl[file_index] < limits[1]) ] / (self.plength[file_index] * self.conc[file_index])
        self.ax.fill_between(x, 0, y, color='r', alpha=0.1)
        
        # flip arrays if not in acending order 
        if x[-1]<x[0]:
            x = np.flip(x)
            y = np.flip(y)

        # integrate 
        # convert to cm-1
        x = x*1000
        f = 4.3e-9 * np.trapz(y, x=x)
        # convert nu0 to from kK to s-1
        nu0 = sc.c*100 * nu0*1000
        print('%.3g'%nu0)
        krad = 7.4e-22 * ((n**2 + 2)/3) * nu0**2 * n * f
        tau_rad = (1/krad)*10**9
        # calculate TDM
        fL = 3*n**2/(2*n**2 + 1)
        TDM = 9.584e-2 * (n/fL**2)**(0.5) * (np.trapz(y/x, x=x))**(0.5)
        print('The oscillator strenth calculated over the region is: f = %.3g'%(f))
        print(r'The TDM is: $\mu_{\text{TDM}} = %.3g\,\text{D}$'%(TDM))
        print(r'Radiative rate constant $k_{\text{rad}} = \SI{%.3g}{\per \second}$'%(krad))
        print(r'Radiative lifetime $\tau_\text{rad} = %.3g\,\text{ns}'%(tau_rad))

    # function to plot calculated spectrum
    def plot_calculated(self, fname, color, label, fwhm=10, shift=0, scale=1, UV=False, norm=True, velox=False):
        # load peak information
        if velox==False:
            wnp, Ip = get_peaks(fname, UV=UV)
            if UV==True:
                # convert to cm–1
                wnp = 10**4/wnp
        else:
            data = np.loadtxt(fname, delimiter=',', skiprows=1)
            wnp = data[:,0]
            Ip = data[:,1]

        wnp = scale*wnp

        # take only those peaks that lie within the range of cuts of exp. 
        if UV==False:
            Ip = Ip[(wnp>self.cuts[0][0]) & (wnp<self.cuts[0][1])]
            wnp = wnp[(wnp>self.cuts[0][0]) & (wnp<self.cuts[0][1])]
        
        # make frequency axis
        if UV==False:
            wn = np.linspace(self.cuts[0][0], self.cuts[0][1], 5000)
        else:
            wn = np.linspace(nm_to_kk(self.cuts[0][1]), nm_to_kk(self.cuts[0][0]), 5000)
        I = np.zeros(len(wn))

        # construct spectrum
        for i in range(len(Ip)):
            I += lorentzian(wn, Ip[i], wnp[i]+shift, fwhm)

        if norm==True:
            norm = np.max(I)
        else:
            norm = 1
        label += r'/ scaling = %.3g'%(scale)
        self.ax.plot(wn, I/norm, '-', color=color, label=label)

        for i in range(len(Ip)):
            self.ax.plot([wnp[i]+shift, wnp[i]+shift], [0, Ip[i]/norm - 0.2*Ip[i]/norm], '-', color=color, alpha=0.5)

    # function to plot absorption solvatochromism
    def solvchrom(self, solv, lim=None, save=False):
        fig1, ax1 = plt.subplots(1,1,figsize=(5, 3.5))
        fig2, ax2 = plt.subplots(1,1,figsize=(5, 3.5))
        fig3, ax3 = plt.subplots(1,1,figsize=(5, 3.5))
        for i in range(len(solv)):
            n = solvs[solv[i]][0]
            er = solvs[solv[i]][1]
            fe = (2*(er - 1)/(2*er+1))
            fn = (2*(n**2 - 1)/(2*n**2+1))
            df = fe - fn
            ax1.plot(df, self.wn[i][self.A[i]==np.max(self.A[i])], 'o', color=self.colors[i])
            ax2.plot(fe, self.wn[i][self.A[i]==np.max(self.A[i])], 'o', color=self.colors[i])
            ax3.plot(fn, self.wn[i][self.A[i]==np.max(self.A[i])], 'o', color=self.colors[i])
        ax1.set_ylabel(r'$\Tilde{\nu}_{\text{max, abs}} / 10^3$ cm$^{-1}$')
        ax2.set_ylabel(r'$\Tilde{\nu}_{\text{max, abs}} / 10^3$ cm$^{-1}$')
        ax3.set_ylabel(r'$\Tilde{\nu}_{\text{max, abs}} / 10^3$ cm$^{-1}$')
        ax1.set_xlabel(r'$\Delta f = f(\varepsilon_r) - f(n^2)$')
        ax2.set_xlabel(r'$f(\varepsilon_r)$')
        ax3.set_xlabel(r'$f(n^2)$')
        if lim!=None:
            ax1.set_ylim(lim)
            ax2.set_ylim(lim)
            ax3.set_ylim(lim)
        fig1.tight_layout()
        fig2.tight_layout()
        fig3.tight_layout()
        if save!=False:
            fig1.savefig('df.svg', transparent=True)
            fig2.savefig('f_epsr.svg', transparent=True)
            fig3.savefig('f_n2.svg', transparent=True)

    # function to plot difference spectrum
    def plot_diff(self, file1, file2, scale=False, export=False):
        wn = self.wn[0]
        A1 = self.A[file1]
        A2 = self.A[file2]
        if scale!=False:
            dA = A1*scale[0] - A2*scale[1]
            self.ax.plot(wn, A1*scale[0] - A2*scale[1], '--k', label=r'Difference $\times$ %i'%(round(scale[1])))
        else:
            dA = A1 - A2
            self.ax.plot(wn, A1 - A2, '--k', label=r'Difference')
        if export==True:
            head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  absorbance'
            np.savetxt('difference.txt', 
                        np.column_stack([(1/wn)*10**4, wn, dA]),
                        header=head, delimiter=',')

    # find abs at certian wavelength / wavenumber
    def find(self, wl_look, file_index):
        self.wl_look = wl_look
        if self.wl_look < 100:
            self.wn_look = self.wl_look
            self.wl_look = (1/self.wn_look)*10**4
        else:
            self.wn_look = (1/self.wl_look)*10**4
        A_look = self.A[file_index][np.argmin(np.abs(self.wl[file_index] - self.wl_look))]
        if self.conc!=None:
            A_look = A_look / (self.plength[file_index] * self.conc[file_index])
            print(r"%s --> \varepsilon_{\text{%.3g nm}} = \SI{%.3g}{\text{M}^{-1}\,\text{cm}^{-1}} at wl = %.3g nm / wn = %.3g kK"%(self.files[file_index],self.wl_look, A_look, self.wl_look, self.wn_look))
        else:
            print("%s --> A = %.3g at wl = %.3g nm / wn = %.3g kK"%(self.files[file_index], A_look, self.wl_look, self.wn_look))
        return A_look
    
    def plot(self):
        # for waterfall plot
        shift = 0
        # go thourgh all spectra and plot them
        for i in self.which:
            if self.units=='wl':
                x = self.wl[i]
            else:
                x = self.wn[i]

            # transition dipole moment representation
            if self.TDM==True:
                y = self.A[i]/x
            # calculate extinction coefficient if c not None
            elif self.conc!=None:
                y = self.A[i]/(self.conc[i] * self.plength[i])
            else:
                y = self.A[i]

            # shift if waterfall
            if self.waterfall!=None:
                y += shift*self.waterfall
                shift+=1

            # plot
            if self.MA==False:
                self.ax.plot(x, y, color=self.colors[i], label=self.labels[i], linestyle=self.linestyle[i], alpha=self.alpha[i])
                if self.fill == True:
                    self.ax.fill_between(x, shift, y, color=self.colors[i], alpha=self.alpha[i]/4)  
            else:
                self.ax.plot(x, y, '-', color=self.colors[i], alpha=0.2)
                self.ax.plot(moving_average(x, self.MA_npoints), moving_average(y, self.MA_npoints), '-', color=self.colors[i], label=self.labels[i])
                if self.fill == True:
                    self.ax.fill_between(moving_average(x, self.MA_npoints), 0, moving_average(y, self.MA_npoints), color=self.colors[i], alpha=0.05)                  
    
    def show(self):
        if self.figsize!=None:
            self.fig.set_size_inches(self.figsize)
        
        if self.units=='wn':
            if self.IR==False:
                self.ax.invert_xaxis()
                self.ax.set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
                axsec = self.ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
                axsec.set_xlabel(r'$\lambda / $ nm') 
                if self.secaxticks!=None:
                    axsec.set_xticks(self.secaxticks)
            else:
                self.ax.set_xlabel(r'$\tilde{\nu} / \text{cm}^{-1}$')
                #axsec = self.ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
                #axsec.set_xlabel(r'$\lambda / $ µm')                 

        if self.inv==True:
            self.ax.invert_xaxis()    
        if self.units=='wl':
            if self.IR==False:
                self.ax.set_xlabel(r'$\lambda / $ nm')
            else:
                self.ax.set_xlabel(r'$\lambda / $ µm')

        if self.legend==True:
            if self.outside==True:
                self.ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=self.fontsize)
            else:
                self.ax.legend(fontsize=self.fontsize)

        if self.yscale!=None:
            self.ax.set_yscale(self.yscale)
            
        if self.zeroline==True:
            self.ax.axhline(y=0, color='k')
        if self.norm == True:
            self.ax.set_ylim([-0.1, 1.2])
            if not self.TDM == True:
                self.ax.set_ylabel('norm. $A$')
            else:
                self.ax.set_ylabel(r'$\varepsilon(\tilde{\nu})/\tilde{\nu}$')                
        elif self.TDM == True:
            self.ax.set_ylabel(r'$\varepsilon(\tilde{\nu})/\tilde{\nu}$')
        elif self.conc!=None:
            self.ax.set_ylabel(r'$\varepsilon / \text{M}^{-1} \, \text{cm}^{-1}$')
        else:
            self.ax.set_ylabel('$A$')
        if self.ylabel!=None:
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
        if self.savefig != False: 
            self.fig.savefig('%s.svg'%(self.savefig), transparent=True)
        plt.show()