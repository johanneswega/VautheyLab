from VautheyLab.standard import *
from VautheyLab.miscellaneous import *

class Ex_Em_Map:
    # initialize class
    def __init__(self, S_file, R_file, corr=True, figsize=(8,5), time=10,
                 MA=True, MA_npoints=5, outside=False, ylim=None, xlim=None, 
                 yticks=False, xticks=True, inv=False, movname_em='movie_emission.mp4',
                 units='wn', abs_spec=None, yscale='linear', movname_ex='excitation_movie.mp4',
                 em_spec=None, cuts_em=None, cuts_ex=None, before=True, Rbg=None, Sbg=None):
        self.S_file = S_file
        self.R_file = R_file
        self.corr = corr
        self.figsize = figsize
        self.time = time
        self.MA = MA
        self.MA_npoints = MA_npoints
        self.outside = outside
        self.ylim = ylim
        self.xlim = xlim
        self.yticks = yticks
        self.xticks = xticks
        self.inv = inv
        self.movname_em = movname_em
        self.units = units
        self.abs_spec = abs_spec
        self.yscale = yscale
        self.movname_ex = movname_ex
        self.em_spec = em_spec
        self.cuts_em = cuts_em
        self.cuts_ex = cuts_ex
        self.before = before
        self.Rbg = Rbg
        self.Sbg = Sbg
        self.get_data()
        self.fig, self.ax = plt.subplots(1,1,figsize=self.figsize)

    # get data
    def get_data(self):
        # load data
        data = np.loadtxt(self.S_file, delimiter=',', dtype=str)
        # get excitation wavelengths
        self.wl_ex = np.array(data[0,1:].astype(float))
        self.wn_ex = (1/self.wl_ex)*10**4
        # get emission wavelengths
        self.wl_em = np.array(data[2:,0].astype(float))
        self.wn_em = (1/self.wl_em)*10**4
        # get emission intensity
        S = data[2:, 1:].astype(float)
        
        # same for R 
        data = np.loadtxt(self.R_file, delimiter=',', dtype=str)
        R = data[2:, 1:].astype(float)      

        # subtract background
        if self.Rbg!=None:
            print('subtract background')
            data = np.loadtxt(self.Rbg, delimiter=',', dtype=str)
            Rb = data[2:, 1:].astype(float)
            data = np.loadtxt(self.Sbg, delimiter=',', dtype=str)
            Sb = data[2:, 1:].astype(float)      
            self.I_wl_bg = (Sb/Rb)
            #S = S - Sb

        # calculate intensity
        self.I_wl = S/R

        self.I_wn = np.zeros((len(self.wl_em), len(self.wl_ex)))
        if self.Rbg!=None:
            self.I_wn_bg = np.zeros((len(self.wl_em), len(self.wl_ex)))

        # correct data if needed
        if self.corr==True:
            corr = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'FMax_lamp_20151217.txt'))
            c = np.interp(self.wl_em, corr[:,0], corr[:,1])
            for i in range(len(self.wl_ex)):
                self.I_wl[:,i] = self.I_wl[:,i] - np.min(self.I_wl[:,i])
                self.I_wl[:,i] = c*self.I_wl[:,i]
                if self.Rbg!=None:
                    self.I_wl_bg[:,i] = self.I_wl_bg[:,i] - np.min(self.I_wl_bg[:,i])
                    self.I_wl_bg[:,i] = c*self.I_wl_bg[:,i]                    

        # calculate intensity for wn spacing
        for i in range(len(self.wl_ex)):
            self.I_wn[:,i] = self.I_wl[:,i]*self.wl_em**2  
            if self.Rbg!=None:
                self.I_wn_bg[:,i] = self.I_wl_bg[:,i]*self.wl_em**2 

        # cut data accordingly
        if self.cuts_em!=None:
            self.I_wn = self.I_wn[(self.wl_em>self.cuts_em[0]) & (self.wl_em<self.cuts_em[1]), :] 
            self.I_wl = self.I_wl[(self.wl_em>self.cuts_em[0]) & (self.wl_em<self.cuts_em[1]), :] 
            self.wn_em = self.wn_em[(self.wl_em>self.cuts_em[0]) & (self.wl_em<self.cuts_em[1])]
            self.wl_em = self.wl_em[(self.wl_em>self.cuts_em[0]) & (self.wl_em<self.cuts_em[1])]
        if self.cuts_ex!=None:
            self.I_wn = self.I_wn[:, (self.wl_ex>self.cuts_ex[0]) & (self.wl_ex<self.cuts_ex[1])] 
            self.I_wl = self.I_wl[:, (self.wl_ex>self.cuts_ex[0]) & (self.wl_ex<self.cuts_ex[1])] 
            self.wn_ex = self.wn_ex[(self.wl_ex>self.cuts_ex[0]) & (self.wl_ex<self.cuts_ex[1])] 
            self.wl_ex = self.wl_ex[(self.wl_ex>self.cuts_ex[0]) & (self.wl_ex<self.cuts_ex[1])]     

        if self.Rbg!=None:
            if self.cuts_em!=None:
                self.I_wn_bg = self.I_wn_bg[(self.wl_em>self.cuts_em[0]) & (self.wl_em<self.cuts_em[1]), :] 
                self.I_wl_bg = self.I_wl_bg[(self.wl_em>self.cuts_em[0]) & (self.wl_em<self.cuts_em[1]), :] 
                self.wn_em = self.wn_em[(self.wl_em>self.cuts_em[0]) & (self.wl_em<self.cuts_em[1])]
                self.wl_em = self.wl_em[(self.wl_em>self.cuts_em[0]) & (self.wl_em<self.cuts_em[1])]
            if self.cuts_ex!=None:
                self.I_wn_bg = self.I_wn_bg[:, (self.wl_ex>self.cuts_ex[0]) & (self.wl_ex<self.cuts_ex[1])] 
                self.I_wl_bg = self.I_wl_bg[:, (self.wl_ex>self.cuts_ex[0]) & (self.wl_ex<self.cuts_ex[1])] 
                self.wn_ex = self.wn_ex[(self.wl_ex>self.cuts_ex[0]) & (self.wl_ex<self.cuts_ex[1])] 
                self.wl_ex = self.wl_ex[(self.wl_ex>self.cuts_ex[0]) & (self.wl_ex<self.cuts_ex[1])]    

    def export_emission(self, wl_exp):
        head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  Intensity for wavelength / a.u., Intensity for wavenumber / a.u.'
        lim1 = wl_exp + 15
        lim2 = 2*wl_exp - 15
        np.savetxt('em_%.3g_nm.txt'%(wl_exp), 
                    np.column_stack([self.wl_em[(self.wl_em>lim1) & (self.wl_em<lim2)], 
                                     self.wn_em[(self.wl_em>lim1) & (self.wl_em<lim2)], 
                                     self.I_wl[(self.wl_em>lim1) & (self.wl_em<lim2), find_index(self.wl_ex, wl_exp)], 
                                     self.I_wn[(self.wl_em>lim1) & (self.wl_em<lim2), find_index(self.wl_ex, wl_exp)]]),
                    header=head, delimiter=',')      
        
    def export_excitation(self, wl_exp):
        fname = 'ex_at_%.3g_nm_em.txt'%(wl_exp)  
        head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  Intensity for wavelength / a.u., Intensity for wavenumber / a.u.'  
        lim1 = wl_exp/2 + 15
        lim2 = wl_exp - 15
        np.savetxt(fname, 
                    np.column_stack([self.wl_ex[(self.wl_ex>lim1) & (self.wl_ex<lim2)],
                                        self.wn_ex[(self.wl_ex>lim1) & (self.wl_ex<lim2)],
                                        self.I_wl[find_index(self.wl_em, wl_exp), (self.wl_ex>lim1) & (self.wl_ex<lim2)],
                                        self.I_wn[find_index(self.wl_em, wl_exp), (self.wl_ex>lim1) & (self.wl_ex<lim2)]]),
                    header=head, delimiter=',')  
        
    def export_all(self):
        if not 'emission_data' in os.listdir():
            os.mkdir('emission_data')
        if not 'excitation_data' in os.listdir():
            os.mkdir('excitation_data')

        for i in range(len(self.wl_ex)):
            fname = 'emission_data/em_at_%.3g_nm_ex.txt'%(self.wl_ex[i])
            head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  Intensity for wavelength / a.u., Intensity for wavenumber / a.u.'
            lim1 = self.wl_ex[i] + 15
            lim2 = 2*self.wl_ex[i]- 15
            np.savetxt(fname, 
                        np.column_stack([self.wl_em[(self.wl_em>lim1) & (self.wl_em<lim2)], 
                                        self.wn_em[(self.wl_em>lim1) & (self.wl_em<lim2)], 
                                        self.I_wl[(self.wl_em>lim1) & (self.wl_em<lim2), find_index(self.wl_ex, self.wl_ex[i])], 
                                        self.I_wn[(self.wl_em>lim1) & (self.wl_em<lim2), find_index(self.wl_ex, self.wl_ex[i])]]),
                        header=head, delimiter=',')  

        for i in range(len(self.wl_em)):
            fname = 'excitation_data/ex_at_%.3g_nm_em.txt'%(self.wl_em[i])  
            head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  Intensity for wavelength / a.u., Intensity for wavenumber / a.u.'  
            lim1 = self.wl_em[i]/2 + 15
            lim2 = self.wl_em[i] - 15
            np.savetxt(fname, 
                        np.column_stack([self.wl_ex[(self.wl_ex>lim1) & (self.wl_ex<lim2)],
                                         self.wn_ex[(self.wl_ex>lim1) & (self.wl_ex<lim2)],
                                         self.I_wl[find_index(self.wl_em, self.wl_em[i]), (self.wl_ex>lim1) & (self.wl_ex<lim2)],
                                         self.I_wn[find_index(self.wl_em, self.wl_em[i]), (self.wl_ex>lim1) & (self.wl_ex<lim2)]]),
                        header=head, delimiter=',')  
            
    def plot_emission_slices(self, slices, bg=False, showbg=False, scale=None, export=True):
        self.ax.clear()
        cols = rainbow(slices)

        for i in range(len(slices)):
            # make sure to get rid of rayleigh region
            lim1 = self.wl_ex[find_index(self.wl_ex, slices[i])] + 20
            lim2 = 2*self.wl_ex[find_index(self.wl_ex, slices[i])] - 20

            # get units
            if self.units=='wn':
                x = self.wn_em[(self.wl_em>lim1) & (self.wl_em<lim2)]
                y = self.I_wn[(self.wl_em>lim1) & (self.wl_em<lim2), find_index(self.wl_ex, slices[i])] 
                if bg==True:
                    ybg = self.I_wn_bg[(self.wl_em>lim1) & (self.wl_em<lim2), find_index(self.wl_ex, slices[i])] 
                self.ax.axvline(x=(1/slices[i])*10**4, color=cols[i], linestyle='--', alpha=0.7)
            else:
                x = self.wl_em[(self.wl_em>lim1) & (self.wl_em<lim2)]
                y = self.I_wl[(self.wl_em>lim1) & (self.wl_em<lim2), find_index(self.wl_ex, slices[i])] 
                if bg==True:
                    ybg = self.I_wl[(self.wl_em>lim1) & (self.wl_em<lim2), find_index(self.wl_ex, slices[i])] 
                self.ax.axvline(x=slices[i], color=cols[i], linestyle='--', alpha=0.7)
            
            # plot background
            if bg==True:
                sub = y/np.max(y) - scale[i]*ybg/np.max(ybg)
                if showbg==True:
                    self.ax.plot(x, scale[i]*ybg/np.max(ybg), '-k')
                    self.ax.plot(x, sub/np.max(sub), '-r')
                else:
                    y = sub
                
            # plot slice
            if self.MA==True:
                norm = np.max(moving_average(y, self.MA_npoints))
                y = y/norm
                self.ax.plot(x, y, '-', alpha=0.2, color=cols[i])      
                self.ax.plot(moving_average(x, self.MA_npoints), moving_average(y, self.MA_npoints)/norm, 
                            '-', color=cols[i], label=r'$\lambda_{\text{ex}} = %i$ nm'%(round(slices[i]))) 
            else:
                y = y/np.max(y)
                self.ax.plot(x, y, '-', color=cols[i], label=r'$\lambda_{\text{ex}} = %i$ nm'%(round(slices[i])))

            # export
            if export==True:
                head = 'wavenumber / 10^3 cm-1, norm. intensity'
                np.savetxt('em_%i_nm.txt'%(slices[i]), np.column_stack([x, y]), header=head, delimiter=',')    

        # plot absorption spectrum
        if self.abs_spec!=None:
            data = np.loadtxt(self.abs_spec[0], delimiter=',', skiprows=1)
            s_wl = data[:,0]
            s_wn = data[:,1]
            s = data[:,2]
            s = s[(s_wl>self.abs_spec[1][0]) & (s_wl<self.abs_spec[1][1])]
            s_wn = s_wn[(s_wl>self.abs_spec[1][0]) & (s_wl<self.abs_spec[1][1])]
            s_wl = s_wl[(s_wl>self.abs_spec[1][0]) & (s_wl<self.abs_spec[1][1])]
            s = self.abs_spec[2]*s/np.max(s)
            if self.units == 'wl':
                self.ax.fill_between(s_wl, 0, s, color=self.abs_spec[3], label=self.abs_spec[4], alpha=0.1)
            else:
                self.ax.fill_between(s_wn, 0, s, color=self.abs_spec[3], label=self.abs_spec[4], alpha=0.1) 

        # stylistic stuff
        self.ax.axhline(y=0, color='k')

        if self.units=='wn':
            self.ax.invert_xaxis() 
            self.ax.set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
            axsec = self.ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
            axsec.set_xlabel(r'$\lambda / $ nm') 

        if self.inv==True:
            self.ax.invert_xaxis()    

        if self.units=='wl':
            self.ax.set_xlabel(r'$\lambda / $ nm')

        #if self.outside==True:
            #self.ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
        #else:
            #self.ax.legend(loc='upper right')

        self.ax.axhline(y=0, color='k')
        self.ax.set_ylabel(r'norm. Emission Intensity')           
        self.ax.set_yscale(self.yscale) 
        if self.yticks == False:
            self.ax.set_yticks([])
        if self.xticks == False:
            self.ax.set_xticks([])
        if self.ylim != None:
            self.ax.set_ylim(self.ylim)
        if self.xlim != None:
            if self.units=='wl':
                self.ax.set_xlim(self.xlim)
            else:
                self.ax.set_xlim([nm_to_kk(self.xlim[0]), nm_to_kk(self.xlim[1])])
        self.fig.tight_layout()
        self.fig.savefig('emission_slices.svg', transparent=True)
        plt.show()

    def emission_movie(self, i):
        self.ax.clear()
        cols = rainbow(self.wl_ex)
        print("%i / %i"%(i, len(self.wl_ex)))
        # make sure to get rid of rayleigh region
        lim1 = self.wl_ex[i] + 20
        lim2 = 2*self.wl_ex[i] - 20
        # get units
        if self.units=='wn':
            x = self.wn_em[(self.wl_em>lim1) & (self.wl_em<lim2)]
            y = self.I_wn[(self.wl_em>lim1) & (self.wl_em<lim2), find_index(self.wl_ex, self.wl_ex[i])] 
            self.ax.axvline(x=(1/self.wl_ex[i])*10**4, color=cols[i])
        else:
            x = self.wl_em[(self.wl_em>lim1) & (self.wl_em<lim2)]
            y = self.I_wl[(self.wl_em>lim1) & (self.wl_em<lim2), find_index(self.wl_ex, self.wl_ex[i])] 
            self.ax.axvline(x=self.wl_ex[i], color=cols[i])

        # plot previous emission spectra
        if self.before==True:
            for j in range(len(self.wl_ex)):
                lim1 = self.wl_ex[j] + 20
                lim2 = 2*self.wl_ex[j] - 20
                if self.wl_ex[j]<self.wl_ex[i]:
                    if self.units=='wl':
                        before_x = self.wl_em[(self.wl_em>lim1) & (self.wl_em<lim2)]
                        before = self.I_wl[(self.wl_em>lim1) & (self.wl_em<lim2), find_index(self.wl_ex, self.wl_ex[j])] 
                    else:
                        before_x = self.wn_em[(self.wl_em>lim1) & (self.wl_em<lim2)]
                        before = self.I_wn[(self.wl_em>lim1) & (self.wl_em<lim2), find_index(self.wl_ex, self.wl_ex[j])] 
                    self.ax.plot(before_x, before/np.max(moving_average(before, 5)), '-k', alpha=0.01)

        # plot current emission sepctrum
        if self.MA==True:
            norm = np.max(moving_average(y, self.MA_npoints))
            self.ax.plot(x, y/norm, '-', alpha=0.2, color=cols[i])      
            self.ax.plot(moving_average(x, self.MA_npoints), moving_average(y, self.MA_npoints)/norm, 
                         '-', color=cols[i]) 
        else:
           self.ax.plot(x, y/np.max(y), '-', color=cols[i])  

        # plot absorption spectrum
        if self.abs_spec!=None:
            data = np.loadtxt(self.abs_spec[0], delimiter=',', skiprows=1)
            s_wl = data[:,0]
            s_wn = data[:,1]
            s = data[:,2]
            s = s[(s_wl>self.abs_spec[1][0]) & (s_wl<self.abs_spec[1][1])]
            s_wn = s_wn[(s_wl>self.abs_spec[1][0]) & (s_wl<self.abs_spec[1][1])]
            s_wl = s_wl[(s_wl>self.abs_spec[1][0]) & (s_wl<self.abs_spec[1][1])]
            s = self.abs_spec[2]*s/np.max(s)
            if self.units == 'wl':
                self.ax.fill_between(s_wl, 0, s, color=self.abs_spec[3], label=self.abs_spec[4], alpha=0.1)
            else:
                self.ax.fill_between(s_wn, 0, s, color=self.abs_spec[3], label=self.abs_spec[4], alpha=0.1) 

        # stylistic stuff
        self.ax.set_title(r'$\lambda_{\text{ex}} = %.3g$ nm'%(self.wl_ex[i]))
        self.ax.axhline(y=0, color='k')

        if self.units=='wn':
            self.ax.invert_xaxis() 
            self.ax.set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
            axsec = self.ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
            axsec.set_xlabel(r'$\lambda / $ nm') 

        if self.inv==True:
            self.ax.invert_xaxis()    

        if self.units=='wl':
            self.ax.set_xlabel(r'$\lambda / $ nm')

        if self.outside==True:
            self.ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
        else:
            self.ax.legend(loc='upper right')
        self.ax.axhline(y=0, color='k')
        self.ax.set_ylabel(r'norm. Emission Intensity')            
        if self.yticks == False:
            self.ax.set_yticks([])
        if self.xticks == False:
            self.ax.set_xticks([])
        if self.ylim != None:
            self.ax.set_ylim(self.ylim)
        if self.xlim != None:
            if self.units=='wl':
                self.ax.set_xlim(self.xlim)
            else:
                self.ax.set_xlim([nm_to_kk(self.xlim[0]), nm_to_kk(self.xlim[1])])
        self.ax.set_yscale(self.yscale)
        self.fig.tight_layout()

    def plot_excitation_slices(self, slices, bg=False, scale=None, showbg=False, export=True):
        self.ax.clear()
        cols = rainbow(slices)

        for i in range(len(slices)):
            # make sure to get rid of rayleigh region
            lim1 = self.wl_em[find_index(self.wl_em, slices[i])]/2 + 15
            lim2 = self.wl_em[find_index(self.wl_em, slices[i])] - 15
            # get units
            if self.units=='wn':
                x = self.wn_ex[(self.wl_ex>lim1) & (self.wl_ex<lim2)]
                y = self.I_wn[find_index(self.wl_em, slices[i]), (self.wl_ex>lim1) & (self.wl_ex<lim2)] 
                if bg==True:
                    ybg = self.I_wn_bg[find_index(self.wl_em, slices[i]), (self.wl_ex>lim1) & (self.wl_ex<lim2)] 
                self.ax.axvline(x=(1/slices[i])*10**4, color=cols[i], linestyle='--', alpha=0.7)
            else:
                x = self.wl_ex[(self.wl_ex>lim1) & (self.wl_ex<lim2)]
                y = self.I_wl[find_index(self.wl_em, slices[i]), (self.wl_ex>lim1) & (self.wl_ex<lim2)] 
                if bg==True:
                    ybg = self.I_wl_bg[find_index(self.wl_em, slices[i]), (self.wl_ex>lim1) & (self.wl_ex<lim2)] 
                self.ax.axvline(x=slices[i], color=cols[i], linestyle='--', alpha=0.7)

            # plot background
            if bg==True:
                sub = y/np.max(y) - scale[i]*ybg/np.max(ybg)
                if showbg==True:
                    self.ax.plot(x, scale[i]*ybg/np.max(ybg), '-k')
                    self.ax.plot(x, sub/np.max(sub), '-r')
                else:
                    y = sub
                
            # plot current excitation sepctrum
            if self.MA==True:
                norm = np.max(moving_average(y, self.MA_npoints))
                y = y/norm
                self.ax.plot(x, y, '-', alpha=0.2, color=cols[i])      
                self.ax.plot(moving_average(x, self.MA_npoints), moving_average(y, self.MA_npoints)/norm, 
                            '-', color=cols[i]) 
            else:
                y = y/np.max(y)
                self.ax.plot(x, y, '-', color=cols[i])  

            # export
            if export==True:
                head = 'wavenumber / 10^3 cm-1, intensity'
                np.savetxt('ex_at_%i_nm_em.txt'%(slices[i]), np.column_stack([x, y]), header=head, delimiter=',')        

        # plot absorption spectrum
        if self.abs_spec!=None:
            data = np.loadtxt(self.abs_spec[0], delimiter=',', skiprows=1)
            s_wl = data[:,0]
            s_wn = data[:,1]
            s = data[:,2]
            s = s[(s_wl>self.abs_spec[1][0]) & (s_wl<self.abs_spec[1][1])]
            s_wn = s_wn[(s_wl>self.abs_spec[1][0]) & (s_wl<self.abs_spec[1][1])]
            s_wl = s_wl[(s_wl>self.abs_spec[1][0]) & (s_wl<self.abs_spec[1][1])]
            s = self.abs_spec[2]*s/np.max(s)
            if self.units == 'wl':
                self.ax.fill_between(s_wl, 0, s, color=self.abs_spec[3], label=self.abs_spec[4], alpha=0.1)
            else:
                self.ax.fill_between(s_wn, 0, s, color=self.abs_spec[3], label=self.abs_spec[4], alpha=0.1) 

        # plot em spectrum
        if self.em_spec!=None:
            data = np.loadtxt(self.em_spec[0], delimiter=',', skiprows=1)
            s_wl = data[:,0]
            s_wn = data[:,1]
            s = data[:,2]
            s_e = data[:,3]
            s = s[(s_wl>self.em_spec[1][0]) & (s_wl<self.em_spec[1][1])]
            s_e = s_e[(s_wl>self.em_spec[1][0]) & (s_wl<self.em_spec[1][1])]
            s_wn = s_wn[(s_wl>self.em_spec[1][0]) & (s_wl<self.em_spec[1][1])]
            s_wl = s_wl[(s_wl>self.em_spec[1][0]) & (s_wl<self.em_spec[1][1])]
            s = self.em_spec[2]*s/np.max(s)
            s_e = self.em_spec[2]*s_e/np.max(s_e)
            if self.units == 'wl':
                self.ax.fill_between(s_wl, 0, s, color=self.em_spec[3], label=self.em_spec[4], alpha=0.1)
            else:
                self.ax.fill_between(s_wn, 0, s_e, color=self.em_spec[3], label=self.em_spec[4], alpha=0.1) 

        # stylistic stuff
        self.ax.axhline(y=0, color='k')

        if self.units=='wn':
            self.ax.invert_xaxis() 
            self.ax.set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
            axsec = self.ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
            axsec.set_xlabel(r'$\lambda / $ nm') 

        if self.inv==True:
            self.ax.invert_xaxis()    

        if self.units=='wl':
            self.ax.set_xlabel(r'$\lambda / $ nm')

        #if self.outside==True:
            #self.ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
        #else:
            #self.ax.legend(loc='upper right')

        self.ax.axhline(y=0, color='k')
        self.ax.set_ylabel(r'norm. Emission Intensity')           
        self.ax.set_yscale(self.yscale) 
        if self.yticks == False:
            self.ax.set_yticks([])
        if self.xticks == False:
            self.ax.set_xticks([])
        if self.ylim != None:
            self.ax.set_ylim(self.ylim)
        if self.xlim != None:
            if self.units=='wl':
                self.ax.set_xlim(self.xlim)
            else:
                self.ax.set_xlim([nm_to_kk(self.xlim[0]), nm_to_kk(self.xlim[1])])
        self.fig.tight_layout()
        self.fig.savefig('emission_slices.svg', transparent=True)
        plt.show()

    def excitation_movie(self, i):
        self.ax.clear()
        cols = rainbow(self.wl_em)
        print("%i / %i"%(i, len(self.wl_em)))
        # make sure to get rid of rayleigh region
        lim1 = self.wl_em[i]/2 + 15
        lim2 = self.wl_em[i] - 15
        # get units
        if self.units=='wn':
            x = self.wn_ex[(self.wl_ex>lim1) & (self.wl_ex<lim2)]
            y = self.I_wn[find_index(self.wl_em, self.wl_em[i]), (self.wl_ex>lim1) & (self.wl_ex<lim2)] 
            self.ax.axvline(x=(1/self.wl_em[i])*10**4, color=cols[i])
        else:
            x = self.wl_ex[(self.wl_ex>lim1) & (self.wl_ex<lim2)]
            y = self.I_wl[find_index(self.wl_em, self.wl_em[i]), (self.wl_ex>lim1) & (self.wl_ex<lim2)] 
            self.ax.axvline(x=self.wl_em[i], color=cols[i])

        # plot previous excitation spectra
        if self.before==True:
            for j in range(len(self.wl_em)):
                lim1 = self.wl_em[j]/2 + 15
                lim2 = self.wl_em[j] - 15
                if self.wl_em[j]<self.wl_em[i]:
                    if self.units=='wl':
                        before_x = self.wl_ex[(self.wl_ex>lim1) & (self.wl_ex<lim2)]
                        before = self.I_wl[find_index(self.wl_em, self.wl_em[j]), (self.wl_ex>lim1) & (self.wl_ex<lim2)] 
                    else:
                        before_x = self.wn_ex[(self.wl_ex>lim1) & (self.wl_ex<lim2)]
                        before = self.I_wn[find_index(self.wl_em, self.wl_em[j]), (self.wl_ex>lim1) & (self.wl_ex<lim2)] 
                    self.ax.plot(before_x, before/np.max(moving_average(before, 5)), '-k', alpha=0.05)

        # plot current excitation sepctrum
        if self.MA==True:
            norm = np.max(moving_average(y, self.MA_npoints))
            self.ax.plot(x, y/norm, '-', alpha=0.2, color=cols[i])      
            self.ax.plot(moving_average(x, self.MA_npoints), moving_average(y, self.MA_npoints)/norm, 
                         '-', color=cols[i]) 
        else:
           self.ax.plot(x, y/np.max(y), '-', color=cols[i])  

        # plot absorption spectrum
        if self.abs_spec!=None:
            data = np.loadtxt(self.abs_spec[0], delimiter=',', skiprows=1)
            s_wl = data[:,0]
            s_wn = data[:,1]
            s = data[:,2]
            s = s[(s_wl>self.abs_spec[1][0]) & (s_wl<self.abs_spec[1][1])]
            s_wn = s_wn[(s_wl>self.abs_spec[1][0]) & (s_wl<self.abs_spec[1][1])]
            s_wl = s_wl[(s_wl>self.abs_spec[1][0]) & (s_wl<self.abs_spec[1][1])]
            s = self.abs_spec[2]*s/np.max(s)
            if self.units == 'wl':
                self.ax.fill_between(s_wl, 0, s, color=self.abs_spec[3], label=self.abs_spec[4], alpha=0.1)
            else:
                self.ax.fill_between(s_wn, 0, s, color=self.abs_spec[3], label=self.abs_spec[4], alpha=0.1) 

        # plot em spectrum
        if self.em_spec!=None:
            data = np.loadtxt(self.em_spec[0], delimiter=',', skiprows=1)
            s_wl = data[:,0]
            s_wn = data[:,1]
            s = data[:,2]
            s_e = data[:,3]
            s = s[(s_wl>self.em_spec[1][0]) & (s_wl<self.em_spec[1][1])]
            s_e = s_e[(s_wl>self.em_spec[1][0]) & (s_wl<self.em_spec[1][1])]
            s_wn = s_wn[(s_wl>self.em_spec[1][0]) & (s_wl<self.em_spec[1][1])]
            s_wl = s_wl[(s_wl>self.em_spec[1][0]) & (s_wl<self.em_spec[1][1])]
            s = self.em_spec[2]*s/np.max(s)
            s_e = self.em_spec[2]*s_e/np.max(s_e)
            if self.units == 'wl':
                self.ax.fill_between(s_wl, 0, s, color=self.em_spec[3], label=self.em_spec[4], alpha=0.1)
            else:
                self.ax.fill_between(s_wn, 0, s_e, color=self.em_spec[3], label=self.em_spec[4], alpha=0.1) 

        # stylistic stuff
        self.ax.set_title(r'$\lambda_{\text{em}} = %.3g$ nm'%(self.wl_em[i]))
        self.ax.axhline(y=0, color='k')

        if self.units=='wn':
            self.ax.invert_xaxis() 
            self.ax.set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
            axsec = self.ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
            axsec.set_xlabel(r'$\lambda / $ nm') 

        if self.inv==True:
            self.ax.invert_xaxis()    

        if self.units=='wl':
            self.ax.set_xlabel(r'$\lambda / $ nm')

        if self.outside==True:
            self.ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
        else:
            self.ax.legend(loc='upper right')
        self.ax.axhline(y=0, color='k')
        self.ax.set_ylabel(r'norm. Emission Intensity')            
        if self.yticks == False:
            self.ax.set_yticks([])
        if self.xticks == False:
            self.ax.set_xticks([])
        if self.ylim != None:
            self.ax.set_ylim(self.ylim)
        if self.xlim != None:
            if self.units=='wl':
                self.ax.set_xlim(self.xlim)
            else:
                self.ax.set_xlim([nm_to_kk(self.xlim[0]), nm_to_kk(self.xlim[1])])
        self.ax.set_yscale(self.yscale)
        self.fig.tight_layout()

    def render_em(self):
        anim = FuncAnimation(self.fig, func=self.emission_movie, frames=len(self.wl_ex), interval=1)
        Writer = writers['ffmpeg']
        writer = Writer(fps=round(len(self.wl_ex)/self.time), metadata={'artist': 'Me'}, bitrate=2500)
        anim.save(self.movname_em, writer)
        plt.show()

    def render_ex(self):
        anim = FuncAnimation(self.fig, func=self.excitation_movie, frames=len(self.wl_em), interval=1)
        Writer = writers['ffmpeg']
        writer = Writer(fps=round(len(self.wl_em)/self.time), metadata={'artist': 'Me'}, bitrate=2500)
        anim.save(self.movname_ex, writer)
        plt.show()