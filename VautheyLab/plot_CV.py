from VautheyLab.standard import *

def cotrell(t, k, tprime, m, b):
    return k/((t - tprime)**0.5)  + m*t + b

class CV:
    def __init__(self, files, colors, labels, fc_cal=None, figsize=None, which=None, ref_el=None, savefig=None, conv='IUPAC',
                 column=None, devide=None, waterfall=None, E0=None, cuts_E=None, ylim=None, xlim=None, yticks=None, 
                 xlabel=None, ylabel=None, xticks=None, legtextcolor=None):
        self.files = files
        self.labels = labels
        self.conv = conv
        self.colors = colors 
        self.fc_cal = fc_cal
        self.xlim = xlim
        self.which = which 
        self.E0 = E0
        self.cuts_E = cuts_E
        self.savefig = savefig
        self.figsize = figsize
        self.column = column
        self.ylim = ylim
        self.yticks = yticks
        self.devide = devide
        self.waterfall = waterfall
        self.ylabel = ylabel
        self.xlabel = xlabel
        self.xticks = xticks
        self.legtextcolor = legtextcolor
        if self.waterfall == None:
            self.waterfall = [0 for _ in self.files]
        if self.column==None:
            self.column = [1 for _ in self.files]
        if self.which==None:
            self.which = [i for i in range(len(self.files))]
        self.ref_el = ref_el

        self.read_data()
        self.fig, self.ax = plt.subplots(1,1)
        self.plot()

    def read_data(self):
        self.E = []
        self.I = []
        for i in range(len(self.files)):
            data = np.genfromtxt(self.files[i], delimiter=',', skip_header=6, 
                                 skip_footer=1, encoding='utf-16', filling_values=np.nan) 
            E = data[:,self.column[i]-1]
            I = data[:,self.column[i]]
            # shift potential by ferrocene reference
            E = E - self.fc_cal[i]
            if self.fc_cal[i]!=0:
                if self.ref_el=='SCE':
                    E += 0.38
            # go to different reference electrode if needed
            if self.cuts_E!=None:
                I = I[(E>=self.cuts_E[i][0])&(E<=self.cuts_E[i][1])]
                E = E[(E>=self.cuts_E[i][0])&(E<=self.cuts_E[i][1])]
            self.E.append(E)
            if self.conv=='burger':
                I = I*(-1)
            if self.devide!=None:
                I = I / self.devide[i]
            self.I.append(I)

    def peak_separation(self, scan_rates, E_range):
        # file index and range
        dE = []
        for i in range(len(scan_rates)):
            E = self.E[i]
            I = self.I[i]
            # cut in desried voltage range
            I = I[(E>=E_range[0])&(E<=E_range[1])]
            E = E[(E>=E_range[0])&(E<=E_range[1])]      
            Epa = E[np.argmax(I)]
            Epc = E[np.argmin(I)]  
            dE.append(np.abs(Epa - Epc)*1000)
        fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
        ax.plot(scan_rates, dE, 'ob')
        ax.set_ylabel(r'$\Delta E_p$ / mV')
        ax.set_xlabel(r'$v / \text{V} \, \text{s}^{-1}$')
        ax.set_ylim([35, 155])
        fig.tight_layout()
        fig.savefig('peak_sep.svg', transparent=True)

    def cotrell_fit(self, scan_rates, E_range, plot=True):
        ratio = []
        # make big overview figure
        if plot==True:
            figbig, axbig = plt.subplots(len(scan_rates), 1, figsize=(5, len(scan_rates)*3))
        # go through all scan rates
        for i in range(len(scan_rates)):
            v = scan_rates[i]
            E = self.E[i]
            I = self.I[i]
            # cut in desried voltage range
            I = I[(E>=E_range[0])&(E<=E_range[1])]
            E = E[(E>=E_range[0])&(E<=E_range[1])]    
            # get scan direction
            Estart = E[0]
            if E[1]<E[5]:
                scan_dir = 'ox'
                # get switching potential
                Elam = np.max(E)
            # reductive scan
            else:
                scan_dir = 'red'
                # get switching potential
                Elam = np.min(E) 
            # calculate switching time
            tlam = np.abs(Elam-Estart)/v
            # potentialstep
            dE = np.abs(E[1]-E[0])
            # calculate time axis
            t = np.array([i*dE/v for i in range(len(E))])
            # find forward peak
            if scan_dir=='ox':
                tpf = t[I==np.max(I)]
                Ipf0 = I[I==np.max(I)]
            else:
                tpf = t[I==np.min(I)]
                Ipf0 = I[I==np.min(I)]    
            # fit linear baseline on the first 20 points
            m, b = np.polyfit(t[:20], I[:20], 1)
            # make fine axis to extrapolate baseline
            lin_for = np.linspace(0, tpf, 1000)
            # calculate baseline
            I_for = m*lin_for + b
            # find time at peak potential on baseline
            linpf = I_for[np.argmin(np.abs(lin_for - tpf))]
            # calculate actual peak current
            Ipf = Ipf0 - linpf

            # specify fitting range for Cotrell-fit
            # as +0.05 V above peak and -0.05 V before switching potential
            lim = 0.05/v
            tfit = t[(t>tpf+lim)&(t<tlam-lim)]
            Ifit = I[(t>tpf+lim)&(t<tlam-lim)]

            # find peak values for back peak
            if scan_dir=='ox':
                tpb = t[I==np.min(I)]
                Ipb0 = I[I==np.min(I)]
            else:
                tpb = t[I==np.max(I)]
                Ipb0 = I[I==np.max(I)]    

            # do curve fit to obtain cotrell fit
            p, pcov = curve_fit(lambda t, k, tprime: cotrell(t, k, tprime, m, b), tfit, Ifit, p0=[Ipf0, tpf])
            # make axis to extrapolate cotrell fit
            t_fit_long = np.linspace(tfit[0], tpb, 1000)
            k = p[0]
            tprime = p[1]
            # calculate fit current
            Ifit_long = []
            for n in range(len(t_fit_long)):
                if t_fit_long[n]<=tlam:
                    Ifit_long.append(k/(np.sqrt(t_fit_long[n] - tprime)) + m*t_fit_long[n] + b)
                else:
                    Ifit_long.append(k/(np.sqrt(t_fit_long[n] - tprime)) - m*(t_fit_long[n] - tlam) - b)
            Ifit_long = np.array(Ifit_long)
            # calculate back current
            cotpb = Ifit_long[np.argmin(np.abs(t_fit_long - tpb))]
            Ipb = Ipb0 - cotpb
            ratio.append(np.abs(Ipb/Ipf))

            # plot if wanted
            if plot==True:
                axbig[i].plot(t, I, '-b')

                axbig[i].plot(lin_for, m*lin_for + b, '--k')
                axbig[i].plot([tpf, tpf], [linpf, Ipf0], 'o--g', label=r'$i_{pf}$ = %.4g µA'%(Ipf))
                axbig[i].plot(tfit, Ifit, '.k')

                axbig[i].plot([tpb, tpb], [cotpb, Ipb0], 'o--r', label=r'$i_{pr}$ = %.4g µA'%(Ipb))
                axbig[i].plot(t_fit_long, Ifit_long, '--k')
                axbig[i].axvline(x=tlam, color='k', linestyle='--', alpha=0.2, label=r'$t_{\lambda}$')

                if v>=1:
                    axbig[i].set_title(r'$v = %.3g$ V/s   $|i_{pr}/i_{pf}|$ = %.3g'%(v, np.abs(Ipb/Ipf)))
                else:
                    axbig[i].set_title(r'$v = %.3g$ mV/s   $|i_{pr}/i_{pf}|$ = %.3g'%(v*1000, np.abs(Ipb/Ipf)))
                axbig[i].plot(lin_for, m*lin_for + b, '--k')
                axbig[i].set_xlabel(r'$t$ / s')
                axbig[i].set_ylabel(r'$i$ / µA')
                axbig[i].legend(loc='upper left', fontsize=10)
        if plot==True:
            figbig.tight_layout()
            figbig.savefig('fits.pdf', transparent=True)

        # make figure for peak current ratio 
        fig, ax = plt.subplots(1,1,figsize=(5,3.5))
        ax.axhline(y=1, color='k', linestyle='--', alpha=0.2)
        ax.plot(scan_rates, ratio, 'og')
        ax.set_ylabel(r'$|i_{pr}/i_{pf}|$')
        ax.set_xlabel(r'$v$ / V$\cdot$s$^{-1}$')
        ax.set_ylim([-0.1, 2.1])
        fig.tight_layout()
        fig.savefig('peak_current_ratio.svg', transparent=True)

    def randles_sevcik(self, scan_rates, scan_direction, E_range, conc=None):
        Ip = []
        # file index and range
        for i in range(len(scan_rates)):
            E = self.E[i]
            I = self.I[i]
            # cut in desried voltage range
            I = I[(E>=E_range[0])&(E<=E_range[1])]
            E = E[(E>=E_range[0])&(E<=E_range[1])]
            if scan_direction=='ox':
                Ip.append(np.max(I))
            else:
                Ip.append(-1*np.min(I))
        fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
        p, pcov = np.polyfit(np.sqrt(scan_rates), Ip, 1, cov=True)
        ax.plot(np.sqrt(scan_rates), p[0]*np.sqrt(scan_rates) + p[1], '--k')
        for i in range(len(scan_rates)):
            ax.plot(np.sqrt(scan_rates[i]), Ip[i], 'o', color=self.colors[i])
        if scan_direction=='red':
            ax.set_ylabel(r'$|i_{pc}|$ / µA')
        else:
            ax.set_ylabel(r'$|i_{pa}|$ / µA')
        # conc in mM
        if conc!=None:
            # convert slope to A s^{-1/2}
            slope = p[0]*10**-6
            err_slope = (pcov[0,0]**(0.5))*10**-6
            # WE electrode glassy carbon d = 3 mm = 0.3 cm | CH Instruments 
            d = 0.3
            # electrode area in cm2
            A = np.pi*(d/2)**2
            # number of e– transfered
            n = 1
            # constant in the equation C mol^-1 V^(-0.5)
            cont = 2.69e5
            # we need to convert the concentration from mM to mol/cm^3
            conc = conc*10**(-6)
            # calculate diffusion coefficient
            D = (slope/(cont * n**(1.5) * A * conc))**2
            Derr = 2*slope*err_slope/(cont * n**(1.5) * A * conc)**2
            print(r'D = (%.3g \pm %.3g) \times \, 10^{-5} \text{cm}^2 \cdot \text{s}^{-1}'%(D*10**5, Derr*10**5))
            ax.set_title(r'$D = (%.3g \pm %.3g) \times \, 10^{-5} \text{cm}^2 \cdot \text{s}^{-1}$'%(D*10**5, Derr*10**5))
        ax.set_xlabel(r'$v^{1/2} / \text{V}^{1/2} \, \text{s}^{-1/2}$')
        fig.tight_layout()
        fig.savefig('randles_sevcik.svg', transparent=True)

    def EC_pot_shift(self, scan_rates, scan_direction, E_range):
        Ep = []
        # file index and range
        for i in range(len(scan_rates)):
            E = self.E[i]
            I = self.I[i]
            # cut in desried voltage range
            I = I[(E>=E_range[0])&(E<=E_range[1])]
            E = E[(E>=E_range[0])&(E<=E_range[1])]
            if scan_direction=='ox':
                Ep.append(E[I==np.max(I)]*1000)
            else:
                Ep.append(E[I==np.min(I)]*1000)
        Ep = np.array(Ep)
        # make figure
        fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
        for i in range(len(scan_rates)):
            ax.plot(np.log10(scan_rates[i]), Ep[i], 'o', color=self.colors[i])
        p, cov = np.polyfit(np.log10(scan_rates), Ep, 1, cov=True)
        ax.plot(np.log10(scan_rates), p[0]*np.log10(scan_rates) + p[1], '--k')
        ax.set_title(r'$\partial E_\text{p} / \partial \log v = (%.3g \pm %.3g)\,\text{mV}$'%(p[0], cov[0,0]**(0.5)))
        ax.set_xlabel(r'$\log(v / \text{V}^{-1} \cdot \text{s})$')
        ax.set_ylabel(r'$E_\text{p}$ / mV')
        fig.tight_layout()
        fig.savefig('EC_peak_shift.svg', transparent=True)

    def export(self):
        for i in range(len(self.files)):
            if self.fc_cal[i]==0:
                if self.ref_el=='SCE':
                    header_x = 'E / V vs. SCE'
                else:
                    header_x = 'E / V vs. Ag/AgCl'
            else:
                header_x = 'E / V vs. Fc+/Fc'
            if self.devide==None:
                header_y = 'i / uA'
            else:
                header_y = 'norm. i'
            np.savetxt('exported/' + self.files[i][:self.files[i].find('.')] + '.txt', np.column_stack([self.E[i], self.I[i]]), delimiter=',',
                       header=header_x + ' , ' + header_y)

    # function to find peak and potentially also E0
    def find_peak(self, index, E_range, plot=False):
        # file index and range
        E = self.E[index]
        I = self.I[index]
        # cut in desried voltage range
        I = I[(E>=E_range[0])&(E<=E_range[1])]
        E = E[(E>=E_range[0])&(E<=E_range[1])]
        ipa = np.max(I)
        ipc = np.min(I)
        print('ipc = %.5g µA'%(ipc))
        Epa = E[np.argmax(I)]
        Epc = E[np.argmin(I)]
        E0 = (np.array(Epa) + np.array(Epc))/2
        print('')
        print(f'E0 = {E0 : .4f} V')
        if plot==True:
            self.ax.plot(E, I, '-g')
            self.ax.plot(Epa, ipa, 'ok')
            self.ax.plot(Epc, ipc, 'ok')
            self.ax.axvline(x=E0, color='k', linestyle='--')
            self.ax.set_title(r'$E_{1/2} = %.4g\,\text{V}$'%E0)
        return E0, ipc

    def plot(self):
        for i in self.which:
            self.ax.plot(self.E[i], self.I[i] + self.waterfall[i], '-', color=self.colors[i], label=self.labels[i])
            if self.E0!=None:
                self.ax.plot([self.E0[i], self.E0[i]], [-10, 10], '--', alpha=0.2, color=self.colors[i])

    def show(self):
        if self.figsize!=None:
            self.fig.set_size_inches(self.figsize)

        self.ax.set_ylabel(r'$i$ / µA')     
        if self.ref_el=='SCE':
            self.ax.set_xlabel(r'$E$ vs. SCE / V')
        if self.ref_el==None:
            self.ax.set_xlabel(r'$E$ vs. Fc$^+$/Fc / V')
        if self.fc_cal[0]==0:
            self.ax.set_xlabel(r'$E$ vs. Ag/AgCl / V')   
        if self.conv=='burger':
            self.ax.invert_xaxis()
        if self.devide!=None:
            self.ax.set_ylabel(r'norm. $i$')
            self.ax.set_yticks([])
        if self.ylim!=None:
            self.ax.set_ylim(self.ylim)
        if self.xlim!=None:
            self.ax.set_xlim(self.xlim)
        if self.yticks!=None:
            self.ax.set_yticks(self.yticks)
        if self.xticks!=None:
            self.ax.set_xticks(self.xticks)
            plt.xticks(fontsize=12)
        if self.xlabel!=None:
            self.ax.set_xlabel(self.xlabel, size=14)
        if self.ylabel!=None:
            self.ax.set_ylabel(self.ylabel, size=14)
        self.fig.tight_layout()
        self.ax.legend()
        if self.legtextcolor!=None:
            for text, color in zip(self.ax.legend().get_texts(), self.colors):
                text.set_color(color)
        if self.savefig!=None: 
            self.fig.savefig(self.savefig, transparent=True)
        plt.show()