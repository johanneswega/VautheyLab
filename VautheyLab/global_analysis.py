from VautheyLab.standard import *
from scipy.special import erf
import matplotlib.patches as patches
import re
import numpy as np

# funciton to load pdat
def load_pdat(file):
    data = np.loadtxt(file, skiprows=1, delimiter=',')
    t = data[1:, 0]
    wl = data[0, 1:]
    dA = data[1:, 1:]  
    return t, wl, dA    

# class to perform global kinetic analysis
class Global_Analysis:
    # initialize class
    def __init__(self, file, p0, t_cuts=None, wl_cuts=None, scatter=None, experiment='femto', scale=[-10, 10],
                 model='sequential', IRF=False, wavelengths=np.arange(350, 750, 20), IR=False, nlevels=51, K=None, C0=None):
        # file to be loaded 
        self.file = file
        # initial values for the fit 
        self.p0 = p0
        # experiment type
        self.experiment = experiment
        # scale for plots
        self.scale = scale
        # time cuts
        self.t_cuts = t_cuts
        # wavelength cuts
        self.wl_cuts = wl_cuts
        # scatter as (x nm, x nm) to exclude
        self.scatter = scatter
        # kinetic model to use 
        self.model = model
        # whether to do an analytical gaussian IRF fit
        self.IRF = IRF
        # wavelengths to plot fitted kinetics at 
        self.wavelengths = wavelengths
        # set True for TRIR data
        self.IR = IR
        # number of contour levels
        self.nlevels = nlevels
        # target matrix for target kinetic analysis
        self.K = K
        # initial conditions for target model
        self.C0 = C0

        # read in data
        self.read_data()

    # method to read in data
    def read_data(self):
        # load data 
        if '.npy' in self.file:
            self.t, self.wl, self.dA = np.load(self.file, allow_pickle=True)
        if '.pdat' in self.file:
            self.t, self.wl, self.dA = load_pdat(self.file)

        # cut the data 
        if self.t_cuts!=None:
            self.dA = self.dA[(self.t>self.t_cuts[0])&(self.t<self.t_cuts[1]), :]
            self.t = self.t[(self.t>self.t_cuts[0])&(self.t<self.t_cuts[1])]
        if self.wl_cuts!=None:
            self.dA = self.dA[:, (self.wl>self.wl_cuts[0])&(self.wl<self.wl_cuts[1])]
            self.wl = self.wl[(self.wl>self.wl_cuts[0])&(self.wl<self.wl_cuts[1])]
        
        # exclude scatter
        if self.scatter!=None:
            lindex = []
            for i in range(len(self.wl)):
                if self.wl[i]>=self.scatter[1] or self.wl[i]<=self.scatter[0]:
                    continue
                else:
                    lindex.append(i)
            self.wl = np.delete(self.wl, lindex)
            self.dA = np.delete(self.dA, lindex, axis=1)

        # calculate wavenumber in kK
        if self.IR==False:
            self.wn = (1/self.wl)*10**4
        else:
            self.wn = self.wl

    # method to plot 2D map
    def plot_2D_map(self, ax, t, wn, dA, scale, white=False):
        if abs(scale[0])<=abs(scale[1]):
            levels = np.linspace(-1*np.abs(self.scale[1]), np.abs(self.scale[1]), self.nlevels)
        else:
            levels = np.linspace(-1*np.abs(self.scale[0]), np.abs(self.scale[0]),self.nlevels)
        colors = plt.cm.RdBu_r(np.linspace(0,1, self.nlevels))
        if white==True:
            res_scale = scale[1]/10
            levels = np.linspace(-1*res_scale , res_scale , self.nlevels)
            colors = plt.cm.seismic(np.linspace(0,1, self.nlevels))
            colors[int(self.nlevels/2)-1] = 0
            for i in range(2):
                colors[int(self.nlevels/2)-1+i] = 0
                colors[int(self.nlevels/2)-1-i] = 0
        D = ax.contourf(wn, t, dA, levels=levels, colors=colors)
        cbar = plt.colorbar(D, ax=ax)
        if white==False:
            if abs(scale[0])<=abs(scale[1]):
                cbar.set_ticks([-scale[1], 0, scale[1]])
            else:
                cbar.set_ticks([-scale[0], 0, scale[0]])
        if white==True: 
            cbar.set_ticks([-1*res_scale , 0, res_scale])
            cbar.set_label(r'$\Delta{A} / 10^{-4}$') 
        else:
            cbar.set_label(r'$\Delta{A} / 10^{-3}$')
        if self.experiment=='nano':
            ax.set_ylabel(r'$\Delta t / \text{ns}$')
            if np.min(t)>0:
                ax.set_yscale('log')
            else:
                ax.set_yscale('symlog')
        else:
            ax.set_ylabel(r'$\Delta t / \text{ps}$')
            if np.min(t)>0:
                ax.set_yscale('log')
            else:
                ax.set_yscale('symlog')
        if self.IR==False:
            ax.invert_xaxis()
            ax.set_xlabel(r'$\tilde{\nu} / 10^{3} \, \text{cm}^{-1}$')
            ax2 = ax.secondary_xaxis("top", functions=(lambda x: (1/x)*10**+4 ,lambda x: (1/x)*10**-4))
            ax2.set_xlabel(r'$\lambda / \text{nm}$')
        else:
            ax.set_xlabel(r'$\tilde{\nu} / \, \text{cm}^{-1}$')
        # if scatter not None draw a white rectangle over scatter region
        if self.scatter!=None:
            rect = [(1/self.scatter[1])*10**4, (1/self.scatter[0])*10**4]
            ax.add_patch(patches.Rectangle((rect[0], ax.get_ylim()[0]), rect[1]-rect[0], ax.get_ylim()[1] - ax.get_ylim()[0], facecolor='white'))

    # method to calculate and plot singular values
    def SVD(self):
        S = np.linalg.svd(self.dA)
        sigma = S[1][:25]
        N = np.linspace(1,25,25)
        fig, ax = plt.subplots(1,1,figsize=(5,3.5))
        ax.plot(N, np.log10(sigma), '.b')
        ax.set_ylabel(r'$\log{(\sigma_i)}$')
        ax.set_xlabel(r'$i$')
        ax.set_title("Singular Values")
        fig.tight_layout()
        plt.show()

    # methods to get bounds of the fit
    def get_bounds(self, p0, IRF):
        b = []
        if IRF==True: 
            # lower bound for tau
            b0 = [0 for i in range(len(p0)-2)]
            # upper bound for tau
            b1 = [np.inf for i in range(len(p0)-2)]
            # add bounds for t0 and fwhm
            b0.append(-10)
            b0.append(-10)
            b1.append(10)
            b1.append(10)
        else:
            b0 = [0 for i in range(len(p0))]
            b1 = [np.inf for i in range(len(p0))]
        b.append(b0)
        b.append(b1)
        return b

    # helper function to build target matrix from string
    def _build_target_K(self, k):
        Kt = np.zeros(self.K.shape, dtype=float)

        for i in range(self.K.shape[0]):
            for j in range(self.K.shape[1]):

                entry = self.K[i, j].replace(" ", "")

                if entry == "0":
                    continue

                tokens = re.findall(r'[+-]?k\d+', entry)

                val = 0.0
                for tok in tokens:
                    sign = -1 if tok.startswith('-') else 1
                    idx = int(tok.lstrip('+-')[1:])
                    val += sign * k[idx]

                Kt[i, j] = val

        return Kt

    # method to calculate concentration matrix based on fit
    def calculate_C(self, p):
        # get parameters
        if self.IRF==False:
            k = 1/p
        else:
            k = 1/p[:-2]
            t0 = p[-2]
            fwhm = p[-1]
            sigma = fwhm/(2*(2*np.log(2))**0.5)

        # get K-matrix for parallel model (DADS)
        if self.model=='parallel':
            K = np.diag(-k)
            # initial values
            C0 = [1 for _ in range(len(K[0]))]

        # get K-matrix for sequential model (EADS)
        if self.model=='sequential':
            K = np.zeros((len(k), len(k)))
            K[0,0] = -1*k[0]
            for i in range(1, len(K)):
                K[i,i-1] = k[i-1]
                K[i,i] = -1*k[i]
            # initial values
            C0 = [1 if i == 0 else 0 for i in range(len(K[0]))]

        if self.model=='target':
            # build the target matrix from the string matrix 
            K = self._build_target_K(k)
            print(K)
            C0 = self.C0

        # Calculate eigenvalues and eigenvectors
        eigenvalues, eigenvectors = np.linalg.eig(K)

        # Create the diagonal matrix of exponential eigenvalues
        # expm_Lambda = np.diag(np.exp(eigenvalues))

        # Matrices of eigenvectors
        U = eigenvectors

        # Calculate the inverse of U
        U_inv = np.linalg.inv(U)

        # Calculate matrix exponential
        # expm_K = np.dot(U, np.dot(expm_Lambda, U_inv))

        # calculate alpha
        alpha = np.dot(U_inv, C0)

        # calculate theta
        if self.IRF==False:
            theta = np.zeros((len(eigenvalues), len(self.t)))
            for i in range(len(eigenvalues)):
                theta[i,:] = alpha[i]*np.exp(eigenvalues[i]*self.t)
        else:
            psi = np.zeros((len(eigenvalues), len(self.t)))
            theta = np.zeros((len(eigenvalues), len(self.t)))
            for i in range(len(eigenvalues)):
                psi[i,:] = 0.5*np.exp(eigenvalues[i]*(self.t - t0 + eigenvalues[i]*sigma**2/2))*(1 + erf((t-t0+eigenvalues[i]*sigma**2)/(sigma*2**0.5)))
            for i in range(len(eigenvalues)):
                theta[i,:] = alpha[i]*psi[i,:]

        # calculate C
        C = np.transpose(np.dot(U, theta))
        return C

    # method to perform least squares fit
    def least_squares_fit(self, p):
        C = self.calculate_C(p)
        self.Sim = np.dot(C, np.dot(np.linalg.pinv(C), self.dA))
        return (self.Sim - self.dA).flatten()
    
    # method to plot the results of the fit
    def plot_fit(self, dA, l, t, p):
        # Calculate C matrix using found parameters
        C = self.calculate_C(p)
        # Calculate simulated spectrum
        Sim = np.dot(C, np.dot(np.linalg.pinv(C), dA))

        # plot simulated spectra and residuals
        fig,ax = plt.subplots(ncols=3,nrows=1, figsize=(13, 4), sharey=True)

        # experimental
        self.plot_2D_map(ax[0], t, l, dA, self.scale, False)
        ax[0].set_title("Experimental Data")

        # simulated
        self.plot_2D_map(ax[1], t, l, Sim, self.scale, False)
        ax[1].set_title("Simulated")
        ax[1].set_ylabel('')

        # residuals 
        self.plot_2D_map(ax[2], t, l, (Sim-dA), self.scale, True)
        residuals = (Sim-dA)**2
        np.save('Res.npy', np.array([t, l, Sim-dA], dtype=object))
        chi2 = np.sum((residuals.flatten()))/(len((Sim-dA).flatten()) - len(p))
        ax[2].set_title(r"Residuals ($\chi^2 = %.3g$)"%(chi2))
        ax[2].set_ylabel('')

        fig.tight_layout()
        fig.savefig('fit.pdf', transparent=True)
        plt.show()
        return chi2
    
    # method to print errors
    def print_errors(self, res, p, IRF):
        print("")
        dic = {0:'A', 1:'B', 2:'C', 3:'D', 4:'E', 5:'F', 6:'G', 7:'H', 8:'I'}
        # Calculate the covariance matrix
        covariance_matrix = np.linalg.inv(np.dot(res.jac.T, res.jac))
        # Calculate the parameter errors
        p_err = np.sqrt(np.diagonal(covariance_matrix))

        # print errors:
        label = ''
        if IRF==True: 
            for i in range(len(p)-2):
                if self.experiment=='nano':
                    if p[i]>=1000:
                        print('tau_%s = (%.3g +- %.3g) µs'%(dic[i], p[i]/1000, p_err[i]/1000))
                        label += r'$\tau_%s$ = (%.3g $\pm$ %.3g) µs \\'%(dic[i], p[i]/1000, p_err[i]/1000)
                    else:
                        print('tau_%s = (%.3g +- %.3g) ns'%(dic[i], p[i], p_err[i]))
                        label += r'$\tau_%s$ = (%.3g $\pm$ %.3g) ns \\'%(dic[i], p[i], p_err[i])
                else: 
                    if p[i]>=1000:
                        print('tau_%s = (%.3g +- %.3g) ns'%(dic[i], p[i]/1000, p_err[i]/1000))
                        label += r'$\tau_%s$ = (%.3g $\pm$ %.3g) ns \\'%(dic[i], p[i]/1000, p_err[i]/1000)
                    else:
                        print('tau_%s = (%.3g +- %.3g) ps'%(dic[i], p[i], p_err[i]))
                        label += r'$\tau_%s$ = (%.3g $\pm$ %.3g) ps \\'%(dic[i], p[i], p_err[i])
            if self.experiment=='nano':
                print('t0 = %.3g ns'%(p[-2]))
                print('FWHM = %.3g ns'%(p[-1]))
            else:
                print('t0 = %.3g ps'%(p[-2]))
                print('FWHM = %.3g ps'%(p[-1]))
            # save lifetimes + error
            np.savetxt('lifetimes.txt', np.column_stack((p[:-2], p_err[:-2])), delimiter=',')
        else:
            for i in range(len(p)):
                if self.experiment=='nano':
                    if p[i]>=1000:
                        print('tau_%s = (%.3g +- %.3g) µs'%(dic[i], p[i]/1000, p_err[i]/1000))
                        label += r'$\tau_%s$ = (%.3g $\pm$ %.3g) µs \\'%(dic[i], p[i]/1000, p_err[i]/1000)
                    else:
                        print('tau_%s = (%.3g +- %.3g) ns'%(dic[i], p[i], p_err[i]))
                        label += r'$\tau_%s$ = (%.3g $\pm$ %.3g) ns \\'%(dic[i], p[i], p_err[i])
                else: 
                    if p[i]>=1000:
                        print('tau_%s = (%.3g +- %.3g) ns'%(dic[i], p[i]/1000, p_err[i]/1000))
                        label += r'$\tau_%s$ = (%.3g $\pm$ %.3g) ns \\'%(dic[i], p[i]/1000, p_err[i]/1000)
                    else:
                        print('tau_%s = (%.3g +- %.3g) ps'%(dic[i], p[i], p_err[i]))  
                        label += r'$\tau_%s$ = (%.3g $\pm$ %.3g) ps \\'%(dic[i], p[i], p_err[i])
            # save lifetimes + error
            np.savetxt('lifetimes.txt', np.column_stack((p, p_err)), delimiter=',')
        print("")
        return label
    
    # method to plot fitted kinetics
    def plot_kinetics(self, dA, l, t, p):
        # Calculate C-matrix using found parameters
        C = self.calculate_C(p)
        # Calculate simulated dA 
        Sim = np.dot(C, np.dot(np.linalg.pinv(C), dA))
        self.wavelengths = np.array(self.wavelengths)
        if self.IR==False:
            self.wavelengths = (1/self.wavelengths)*10**4
        # make rainbow colormap for all chosen kinectics
        col = plt.cm.rainbow_r(np.linspace(0, 1, len(self.wavelengths)))
        # generate figure
        fig, ax = plt.subplots(2, 1, figsize=(8, 6), gridspec_kw={'height_ratios':[1,3]}, sharex=True)
        # go through chosen wavenumbers and plot raw data, fit as well as the residuals
        for i in range(len(self.wavelengths)):
            # find wavenumber index closest to the wavenumber you picked 
            index = min(range(len(l)), key=lambda j: abs(l[j]-self.wavelengths[i]))
            # plot raw data
            ax[1].plot(t, dA[:, index], 'o', alpha=0.2, markersize=5, color=col[i])
            # plot simulated data
            if self.IR==False:
                ax[1].plot(t, Sim[:, index], '-', color=col[i], linewidth=2, label=r'%.3g nm'%((1/l[index])*10**4))
            else:
                ax[1].plot(t, Sim[:, index], '-', color=col[i], linewidth=2, label=r'%i cm$^{-1}$'%(round(l[index])))
            # plot residuals
            ax[0].plot(t, Sim[:, index]-dA[:, index], '-', color=col[i])
        # set x-lables
        if self.experiment=='nano':
            ax[1].set_xlabel(r'$\Delta t / \text{ns}$')
            ax[1].set_xscale('log')
        else:
            ax[1].set_xlabel(r'$\Delta t / \text{ps}$')
            if np.min(t)<0.1:
                ax[1].set_xscale('symlog')
                ax[1].set_xticks([1, 2, 10, 100, 1000])
                ax[1].set_xticklabels(['1', '2', '10', '100', '1000'])
            else:
                ax[1].set_xscale('log')
        # stylistic stuff
        ax[0].axhline(y=0, color='k')
        ax[0].set_ylim([-2, 2])
        ax[0].set_ylabel(r'residuals')
        ax[1].set_ylabel(r'$\Delta{A} / 10^{-3}$')
        ax[1].axhline(y=0, color='k')
        ax[1].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)
        fig.tight_layout()  
        fig.savefig('kinetics.svg', transparent=True)
        plt.show()

    # metod to plot EADS / SADS / DADS 
    def plot_EAS(self, dA, l, t, p, model, labelt):
        fig,ax = plt.subplots(ncols=2,nrows=1, figsize=(10, 4)) 
        dic = {0:'A', 1:'B', 2:'C', 3:'D', 4:'E', 5:'F', 6:'G', 7:'H', 8:'I'}
        col = ['r', 'b', 'g', 'orange', 'purple', 'k']
        C = self.calculate_C(p)
        EAS = np.dot(np.linalg.pinv(C), dA)

        # plot results
        # time evolution
        for i in range(len(C[0,:])):
            ax[0].plot(t, C[:,i], '-', color=col[i], label='%s'%(dic[i]))
        ax[0].legend()
        ax[0].set_ylabel('rel. concentration')
        if self.experiment=='nano':
            ax[0].set_xlabel(r'$\Delta t / \text{ns}$')
            if np.min(t)>0:
                ax[0].set_xscale('log')
            else:
                ax[0].set_xscale('symlog')
        else:
            ax[0].set_xlabel(r'$\Delta t / \text{ps}$')
            if np.min(t)>0:
                ax[0].set_xscale('log')
            else:
                ax[0].set_xscale('symlog')
        if model=='sequential':
            ax[1].set_ylabel('EADS / mOD')
        if model=='parallel':
            ax[1].set_ylabel('DADS / mOD')
        if model=='target':
            ax[1].set_ylabel('SADS / mOD')
        ax[0].set_title(labelt, size=11)
        ax[0].set_ylim([-0.1,1.1])

        # spectra
        for i in range(len(C[0,:])):
            ax[1].plot(l, EAS[i,:], '-', color=col[i], label='%s'%(dic[i]), zorder=1)
            np.savetxt('%s.txt'%(dic[i]), np.column_stack([(1/l)*10**4, l, EAS[i,:]]), delimiter=',',
                        header='wavenumber / 10^3 cm-1, wavelength / nm, EADS / mOD')
        ax[1].legend()
        ax[1].axhline(y=0, color='k')
        if self.IR==False:
            ax[1].invert_xaxis()
            ax[1].set_xlabel(r'$\tilde{\nu} / 10^{3} \, \text{cm}^{-1}$')
            ax2 = ax[1].secondary_xaxis("top", functions=(lambda x: (1/x)*10**+4,lambda x: (1/x)*10**-4))
            ax2.set_xlabel(r'$\lambda / \text{nm}$')
        else:
            ax[1].set_xlabel(r'$\tilde{\nu} / \text{cm}^{-1}$')
        fig.tight_layout()
        if self.scatter!=None:
            rect = [(1/self.scatter[1])*10**4, (1/self.scatter[0])*10**4]
            ax[1].add_patch(patches.Rectangle((rect[0], ax[1].get_ylim()[0]), rect[1]-rect[0], ax[1].get_ylim()[1]-ax[1].get_ylim()[0], facecolor='white'))
        fig.savefig('global.svg', transparent=True)
        plt.show()

    # method to start fit
    def fit(self):
        # plot raw data 
        fig, ax = plt.subplots(1,1)
        self.plot_2D_map(ax, self.t, self.wn, self.dA, self.scale)
        fig.tight_layout()
        plt.show()

        # calculate singular values
        self.SVD()

        # estimate fit bounds
        self.b = self.get_bounds(self.p0, self.IRF)

        # perform fit
        res = least_squares(self.least_squares_fit, self.p0, bounds=self.b)
        p = res.x

        # plot fit
        self.plot_fit(self.dA, self.wn, self.t, p)

        # plot kinetics
        self.plot_kinetics(self.dA, self.wn, self.t, p)
        label = self.print_errors(res, p, self.IRF)
        self.plot_EAS(self.dA, self.wn, self.t, p, self.model, label)

