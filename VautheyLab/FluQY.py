from VautheyLab.standard import *
from VautheyLab.miscellaneous import find_index, solvs

class FluQY:
    # initialize class
    def __init__(self, files_abs, files_em, ex_wl, solv, phi_R, baselineat=None, OD_filter=None, slit=None, corr=False):
        # files_abs = [sample, ref]
        self.files_abs = files_abs
        self.files_em = files_em
        self.ex_wl = ex_wl
        self.baselineat = baselineat
        self.OD_filter = OD_filter
        # excitation mono slit in nm
        self.slit = slit
        self.corr = corr

        # get refractive indices
        n_S, n_R = solvs[solv[0]][0], solvs[solv[1]][0]
        print(f'\n refractive indices at 25°C of \n nS ({solv[0]}) = {n_S} \n nR ({solv[1]}) = {n_R} \n')

        # get abs at ex_wl
        self.get_Aex()

        # get integrated emission intensities
        self.get_Int()

        # calculate fluoresence quantum yield
        fS = 1 - 10**(-self.AexS)
        fR = 1 - 10**(-self.AexR)
        phiQflu = phi_R * (fR/fS) * (self.int_S / self.int_R) * (n_S/n_R)**2
        print(f'\n Fluorescence Quantum Yield = {phiQflu*100} percent')

        ### error estimation ### 
        # take error in Aex as 0.002
        errAex = 0.002
        # take error in intensities as relative error of 2 % 
        errint_S = 0.02*self.int_S
        errint_R = 0.02*self.int_R
        # take error in refractive indicies as 0.0002
        errn = 0.0002
        # take error in reference quantum yield as 2 % 
        errphiR = 0.02

        # partial derivatives
        A_S, A_R, I_S, I_R = self.AexS, self.AexR, self.int_S, self.int_R
        dphiR = 10**(-A_R + self.AexS)*I_S*n_S**2*(10**A_R - 1)/(I_R*n_R**2*(10**A_S - 1))
        dI_S = 10**(-A_R + A_S)*n_S**2*phi_R*(10**A_R - 1)/(I_R*n_R**2*(10**A_S - 1))
        dI_R = -10**(-A_R + A_S)*I_S*n_S**2*phi_R*(10**A_R - 1)/(I_R**2*n_R**2*(10**A_S - 1))
        dA_R = 10**(-A_R + A_S)*I_S*n_S**2*phi_R*np.log(10)/(I_R*n_R**2*(10**A_S - 1))
        dA_S = -10**(-A_R + A_S)*I_S*n_S**2*phi_R*(10**A_R - 1)*np.log(10)/(I_R*n_R**2*(10**A_S - 1)**2)
        dn_S = 2*10**(-A_R + A_S)*I_S*n_S*phi_R*(10**A_R - 1)/(I_R*n_R**2*(10**A_S - 1))
        dn_R = -2*10**(-A_R + A_S)*I_S*n_S**2*phi_R*(10**A_R - 1)/(I_R*n_R**3*(10**A_S - 1))

        # calculate error 
        errQflu = np.sqrt( (dphiR * errphiR)**2 + (dI_S * errint_S)**2 + (dI_R * errint_R)**2 + (dA_R * errAex)**2 + (dA_S * errAex)**2 + (dn_S * errn)**2 + (dn_R * errn)**2 )
        print(f'\n Fluorescence Quantum Yield = ({phiQflu*100 : .3f} +- {errQflu*100 : .3f}) percent \n ({100*errQflu/phiQflu : .3f} perc. rel error)')

    # function to get absorption and excitation wavelength
    def get_Aex(self):
        # for sample
        data = np.loadtxt(self.files_abs[0], delimiter=',', usecols=[0,1], skiprows=2)
        wl_S, A_S = data[:,0], data[:,1]
        if self.baselineat[0]!=None:
            A_S = A_S - A_S[find_index(wl_S, self.baselineat[0])]
        indS = find_index(wl_S, self.ex_wl[0])
        # use the method proposed by Demas et. al to average abs at ex. wavelength taking into account the slitwidth of the excitation mono
        if self.slit==None:
            self.AexS = A_S[indS]
        else:
            self.AexS = (A_S[find_index(wl_S, self.ex_wl[0])] + A_S[find_index(wl_S, self.ex_wl[0] + self.slit/2)] + A_S[find_index(wl_S, self.ex_wl[0] - self.slit/2)])/3
            print(f'AS = {self.AexS} (Demas) vs. {A_S[indS]} (value at wl_ex)')

        # for reference
        data = np.loadtxt(self.files_abs[1], delimiter=',', usecols=[0,1], skiprows=2)
        wl_R, A_R = data[:,0], data[:,1]
        if self.baselineat[1]!=None:
            A_R = A_R - A_R[find_index(wl_R, self.baselineat[1])]
        indR = find_index(wl_R, self.ex_wl[1])
        if self.slit==None:
            self.AexR = A_R[indS]
        else:
            self.AexR = (A_R[find_index(wl_R, self.ex_wl[1])] + A_R[find_index(wl_R, self.ex_wl[1] + self.slit/2)] + A_R[find_index(wl_R, self.ex_wl[1] - self.slit/2)])/3
            print(f'AR = {self.AexR} (Demas) vs. {A_R[indR]} (value at wl_ex)')

        # plot absorption spectra 
        fig, ax = plt.subplots(1, 2, sharey=True, sharex=True)
        ax[0].plot(wl_S, A_S, '-r')
        ax[0].plot(wl_S[indS], A_S[indS], 'ok', markersize=5)
        ax[1].plot(wl_R, A_R, '-b')
        ax[1].plot(wl_R[indR], A_R[indR], 'ok', markersize=5)
        ax[0].set_xlabel(r'$\lambda / $ nm')
        ax[1].set_xlabel(r'$\lambda / $ nm')
        ax[0].set_ylabel('$A$')
        ax[0].set_title(r'sample, $A_{%i nm} = %.3g $'%(self.ex_wl[0], self.AexS))
        ax[1].set_title(r'reference, $A_{%i nm} = %.3g $'%(self.ex_wl[1], self.AexR))
        fig.tight_layout()
        plt.show()

    # function to get integrated emission intensities
    def get_Int(self):
        # for sample
        data = np.loadtxt(self.files_em[0], delimiter=',', usecols=[0,1], skiprows=2)
        wl_S, I_S = data[:,0], data[:,1]

        # for reference
        data = np.loadtxt(self.files_em[1], delimiter=',', usecols=[0,1], skiprows=2)
        wl_R, I_R = data[:,0], data[:,1]

        if self.corr!=False: 
            corr = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'FMax_lamp_20151217.txt'))
            c = np.interp(wl_S, corr[:,0], corr[:,1])
            I_S = c*I_S
            c = np.interp(wl_R, corr[:,0], corr[:,1])
            I_R = c*I_R

        # get scaling factor if OD_filter was used to reduce intensity of reference emission
        if self.OD_filter!=None:
            data = np.loadtxt(self.OD_filter[0], delimiter=',', usecols=[0,1], skiprows=2)
            wl, I_noOD = data[:,0], data[:,1]
            data = np.loadtxt(self.OD_filter[1], delimiter=',', usecols=[0,1], skiprows=2)
            wl, I_OD = data[:,0], data[:,1]

            # calculate integrals
            int_noOD = np.trapz(I_noOD, x=wl)
            int_OD = np.trapz(I_OD, x=wl)

            # get scaling factor
            factor = int_noOD/int_OD
            # apply scaling factor 
            I_R = factor*I_R

            fig, ax = plt.subplots(1,1)
            ax.plot(wl, I_noOD, '-b', label='without OD filter')
            ax.plot(wl, I_OD, '--b', label='with OD filter')
            ax.set_title('scaling factor = %.3g'%(factor))
            ax.set_xlabel(r'$\lambda / $ nm')
            ax.set_ylabel('counts')
            ax.legend()
            fig.tight_layout()
            plt.show()
        
        # get integrated intensites
        self.int_S = np.trapz(I_S, x=wl_S)
        self.int_R = np.trapz(I_R, x=wl_R)

        fig, ax = plt.subplots(1, 2, sharey=True, sharex=True)
        ax[0].plot(wl_S, I_S, '-r')
        ax[1].plot(wl_R, I_R, '-b')
        ax[0].set_xlabel(r'$\lambda / $ nm')
        ax[1].set_xlabel(r'$\lambda / $ nm')
        ax[0].set_ylabel('counts')
        ax[0].set_title(r'sample, $I_{S} = %.3g $'%(self.int_S))
        ax[1].set_title(r'reference, $I_{R} = %.3g $'%(self.int_R))
        fig.tight_layout()
        plt.show()

