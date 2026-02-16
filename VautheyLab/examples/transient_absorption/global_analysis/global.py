import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from VautheyLab.standard import *

# file name
fname = 'dA.npy'

# experiment
experiment = 'femto'

# define scale
scale = np.array([-25, 25])

# specify number of contours
nlevels = 50

# specify cuts 
# time in pss
delay = [5, 1800]
# wavelength in nm / or cm-1 for TRIR
wl = [350, 750]
# specify scatter region to exclude
#scatter = None
scatter = [520, 560]

# specify model
model = 'sequential'

# convolution fit
IRF = False

# TRIR data set
IR = False

# initial guess
p0 = [10, 100]

# choose wl to plot kinetics
wlkin = np.arange(350, 750, 20)
#wlkin = np.arange(2100, 2250, 10)


def cut_data(t, l, dA):
    # cut the data 
    dA = dA[(t>delay[0])&(t<delay[1]), :]
    dA = dA[:, (l>wl[0])&(l<wl[1])]
    t = t[(t>delay[0])&(t<delay[1])]
    l = l[(l>wl[0])&(l<wl[1])]
    if scatter!=None:
        # exclude scatter
        lindex = []
        for i in range(len(l)):
            if l[i]>=scatter[1] or l[i]<=scatter[0]:
                continue
            else:
                lindex.append(i)
        l = np.delete(l, lindex)
        dA = np.delete(dA, lindex, axis=1)
    return t, l, dA

def calculate_C(t, p, model, IRF=True):
    # get parameters
    if IRF==False:
        k = 1/p
    else:
        k = 1/p[:-2]
        t0 = p[-2]
        fwhm = p[-1]
        sigma = fwhm/(2*(2*np.log(2))**0.5)

    # get K-matrix for parallel model (DADS)
    if model=='parallel':
        K = np.diag(-k)
        # initial values
        C0 = [1 for i in range(len(K[0]))]

    # get K-matrix for sequential model (EADS)
    if model=='sequential':
        K = np.zeros((len(k), len(k)))
        K[0,0] = -1*k[0]
        for i in range(1, len(K)):
            K[i,i-1] = k[i-1]
            K[i,i] = -1*k[i]
        # initial values
        C0 = [1 if i == 0 else 0 for i in range(len(K[0]))]

    if model=='target':
        # let's try it with the branching model
        # A -->k0 0
        # A -->k1 B -->k3 0
        # A -->k2 C -->k4 0
        # specify your own K-matrix
        K = np.array([[-1*(k[0] + k[1] + k[2]), 0, 0],
                    [k[1], -1*k[3], 0],
                    [k[2], 0, -1*k[4]]])
    
        # initial values
        #C0 = [1 if i == 0 else 0 for i in range(len(K[0]))]
        C0 = [1, 1, 0]

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
    if IRF==False:
        theta = np.zeros((len(eigenvalues), len(t)))
        for i in range(len(eigenvalues)):
            theta[i,:] = alpha[i]*np.exp(eigenvalues[i]*t)
    else:
        psi = np.zeros((len(eigenvalues), len(t)))
        theta = np.zeros((len(eigenvalues), len(t)))
        for i in range(len(eigenvalues)):
            psi[i,:] = 0.5*np.exp(eigenvalues[i]*(t - t0 + eigenvalues[i]*sigma**2/2))*(1 + erf((t-t0+eigenvalues[i]*sigma**2)/(sigma*2**0.5)))
        for i in range(len(eigenvalues)):
            theta[i,:] = alpha[i]*psi[i,:]

    # calculate C
    C = np.transpose(np.dot(U, theta))
    return C

def plot_2D_map(ax, t, l, dA, scale, white=False):
    if abs(scale[0])<=abs(scale[1]):
        levels = np.linspace(-1*np.abs(scale[1]), np.abs(scale[1]), nlevels)
    else:
        levels = np.linspace(-1*np.abs(scale[0]), np.abs(scale[0]), nlevels)
    colors = plt.cm.RdBu_r(np.linspace(0,1, nlevels))
    if white==True:
        res_scale = scale[1]/10
        levels = np.linspace(-1*res_scale , res_scale , nlevels)
        colors = plt.cm.seismic(np.linspace(0,1, nlevels))
        colors[int(nlevels/2)-1] = 0
        for i in range(2):
            colors[int(nlevels/2)-1+i] = 0
            colors[int(nlevels/2)-1-i] = 0
    D = ax.contourf(l, t, dA, levels=levels, colors=colors)
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
    if experiment=='nano':
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
    if IR==False:
        ax.invert_xaxis()
        ax.set_xlabel(r'$\tilde{\nu} / 10^{3} \, \text{cm}^{-1}$')
        ax2 = ax.secondary_xaxis("top", functions=(lambda x: (1/x)*10**+4 ,lambda x: (1/x)*10**-4))
        ax2.set_xlabel(r'$\lambda / \text{nm}$')
    else:
        ax.set_xlabel(r'$\tilde{\nu} / \, \text{cm}^{-1}$')
    # if scatter not None draw a white rectangle over scatter region
    if scatter!=None:
        rect = [(1/scatter[1])*10**4, (1/scatter[0])*10**4]
        ax.add_patch(patches.Rectangle((rect[0], ax.get_ylim()[0]), rect[1]-rect[0], ax.get_ylim()[1] - ax.get_ylim()[0], facecolor='white'))

def SVD(dA):
    S = np.linalg.svd(dA)
    sigma = S[1][:25]
    N = np.linspace(1,25,25)
    fig,ax = plt.subplots(ncols=1,nrows=1)
    ax.plot(N, np.log10(sigma), '.b')
    ax.set_ylabel(r'$\log{(\sigma_i)}$')
    ax.set_xlabel(r'$i$')
    ax.set_title("Singular Values")
    fig.tight_layout()
    plt.show()

def get_bounds(p0, IRF):
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

def least_squares_fit(p):
    C = calculate_C(t, p, model, IRF)
    Sim = np.dot(C, np.dot(np.linalg.pinv(C), dA))
    return (Sim - dA).flatten()

def plot_fit(dA, l, t, p, model):
    # Calculate C matrix using found parameters
    C = calculate_C(t, p, model, IRF)
    # Calculate simulated spectrum
    Sim = np.dot(C, np.dot(np.linalg.pinv(C), dA))

    # plot simulated spectra and residuals
    fig,ax = plt.subplots(ncols=3,nrows=1, figsize=(13, 4), sharey=True)

    # experimental
    plot_2D_map(ax[0], t, l, dA, scale, False)
    ax[0].set_title("Experimental Data")

    # simulated
    plot_2D_map(ax[1], t, l, Sim, scale, False)
    ax[1].set_title("Simulated")
    ax[1].set_ylabel('')

    # residuals 
    plot_2D_map(ax[2], t, l, (Sim-dA), scale, True)
    residuals = (Sim-dA)**2
    chi2 = np.sum((residuals.flatten()))/(len((Sim-dA).flatten()) - len(p))
    ax[2].set_title(r"Residuals ($\chi^2 = %.3g$)"%(chi2))
    ax[2].set_ylabel('')

    fig.tight_layout()
    plt.show()
    return chi2

def plot_kinetics(dA, l, t, p, model):
    # Calculate C-matrix using found parameters
    C = calculate_C(t, p, model, IRF)
    # Calculate simulated dA 
    Sim = np.dot(C, np.dot(np.linalg.pinv(C), dA))
    # get 15 equally spaced wavenumbers through dA
    if IR==False:
        wl = (1/wlkin)*10**4
    else:
        wl = wlkin
    # make rainbow colormap for all chosen kinectics
    col = plt.cm.rainbow_r(np.linspace(0, 1, len(wl)))
    # generate figure
    fig, ax = plt.subplots(2, 1, figsize=(8, 6), gridspec_kw={'height_ratios':[1,3]}, sharex=True)
    # go through chosen wavenumbers and plot raw data, fit as well as the residuals
    for i in range(len(wl)):
        # find wavenumber index closest to the wavenumber you picked 
        index = min(range(len(l)), key=lambda j: abs(l[j]-wl[i]))
        # plot raw data
        ax[1].plot(t, dA[:, index], '-', alpha=0.5, color=col[i])
        # plot simulated data
        if IR==False:
            ax[1].plot(t, Sim[:, index], '-', color=col[i], linewidth=2, label=r'%.3g nm'%((1/l[index])*10**4))
        else:
            ax[1].plot(t, Sim[:, index], '-', color=col[i], linewidth=2, label=r'%i cm$^{-1}$'%(round(l[index])))
        # plot residuals
        ax[0].plot(t, Sim[:, index]-dA[:, index], '-', color=col[i])
    # set x-lables
    if experiment=='nano':
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

def plot_EAS(dA, l, t, p, model, labelt):
    fig,ax = plt.subplots(ncols=2,nrows=1, figsize=(10, 4)) 
    dic = {0:'A', 1:'B', 2:'C', 3:'D', 4:'E', 5:'F', 6:'G', 7:'H', 8:'I'}
    col = ['r', 'b', 'g', 'orange', 'purple', 'k']
    C = calculate_C(t, p, model, IRF)
    EAS = np.dot(np.linalg.pinv(C), dA)

    # plot results
    # time evolution
    for i in range(len(C[0,:])):
        ax[0].plot(t, C[:,i], '-', color=col[i], label='%s'%(dic[i]))
    ax[0].legend()
    ax[0].set_ylabel('rel. concentration')
    if experiment=='nano':
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
        np.savetxt('%s.txt'%(dic[i]), np.column_stack([(1/l)*10**4, EAS[i,:]]), delimiter=',')
    ax[1].legend()
    ax[1].axhline(y=0, color='k')
    if IR==False:
        ax[1].invert_xaxis()
        ax[1].set_xlabel(r'$\tilde{\nu} / 10^{3} \, \text{cm}^{-1}$')
        ax2 = ax[1].secondary_xaxis("top", functions=(lambda x: (1/x)*10**+4,lambda x: (1/x)*10**-4))
        ax2.set_xlabel(r'$\lambda / \text{nm}$')
    else:
        ax[1].set_xlabel(r'$\tilde{\nu} / \text{cm}^{-1}$')
    fig.tight_layout()
    if scatter!=None:
        rect = [(1/scatter[1])*10**4, (1/scatter[0])*10**4]
        ax[1].add_patch(patches.Rectangle((rect[0], ax[1].get_ylim()[0]), rect[1]-rect[0], ax[1].get_ylim()[1]-ax[1].get_ylim()[0], facecolor='white'))
    fig.savefig('global.svg', transparent=True)
    plt.show()

def print_errors(res, p, IRF):
    print("")
    # Calculate the covariance matrix
    covariance_matrix = np.linalg.inv(np.dot(res.jac.T, res.jac))
    # Calculate the parameter errors
    p_err = np.sqrt(np.diagonal(covariance_matrix))

    # print errors:
    label = ''
    if IRF==True: 
        for i in range(len(p)-2):
            if experiment=='nano':
                if p[i]>=1000:
                    print('tau_%i = (%.3g +- %.3g) µs'%(i, p[i]/1000, p_err[i]/1000))
                    label += r'$\tau_%i$ = (%.3g $\pm$ %.3g) µs \\'%(i, p[i]/1000, p_err[i]/1000)
                else:
                    print('tau_%i = (%.3g +- %.3g) ns'%(i, p[i], p_err[i]))
                    label += r'$\tau_%i$ = (%.3g $\pm$ %.3g) ns \\'%(i, p[i], p_err[i])
            else: 
                if p[i]>=1000:
                    print('tau_%i = (%.3g +- %.3g) ns'%(i, p[i]/1000, p_err[i]/1000))
                    label += r'$\tau_%i$ = (%.3g $\pm$ %.3g) ns \\'%(i, p[i]/1000, p_err[i]/1000)
                else:
                    print('tau_%i = (%.3g +- %.3g) ps'%(i, p[i], p_err[i]))
                    label += r'$\tau_%i$ = (%.3g $\pm$ %.3g) ps \\'%(i, p[i], p_err[i])
        if experiment=='nano':
            print('t0 = %.3g ns'%(p[-2]))
            print('FWHM = %.3g ns'%(p[-1]))
        else:
            print('t0 = %.3g ps'%(p[-2]))
            print('FWHM = %.3g ps'%(p[-1]))
        # save lifetimes + error
        np.savetxt('lifetimes.txt', np.column_stack((p[:-2], p_err[:-2])), delimiter=',')
    else:
        for i in range(len(p)):
            if experiment=='nano':
                if p[i]>=1000:
                    print('tau_%i = (%.3g +- %.3g) µs'%(i, p[i]/1000, p_err[i]/1000))
                    label += r'$\tau_%i$ = (%.3g $\pm$ %.3g) µs \\'%(i, p[i]/1000, p_err[i]/1000)
                else:
                    print('tau_%i = (%.3g +- %.3g) ns'%(i, p[i], p_err[i]))
                    label += r'$\tau_%i$ = (%.3g $\pm$ %.3g) ns \\'%(i, p[i], p_err[i])
            else: 
                if p[i]>=1000:
                    print('tau_%i = (%.3g +- %.3g) ns'%(i, p[i]/1000, p_err[i]/1000))
                    label += r'$\tau_%i$ = (%.3g $\pm$ %.3g) ns \\'%(i, p[i]/1000, p_err[i]/1000)
                else:
                    print('tau_%i = (%.3g +- %.3g) ps'%(i, p[i], p_err[i]))  
                    label += r'$\tau_%i$ = (%.3g $\pm$ %.3g) ps \\'%(i, p[i], p_err[i])
        # save lifetimes + error
        np.savetxt('lifetimes.txt', np.column_stack((p, p_err)), delimiter=',')
    print("")
    return label

def load_pdat(file):
    data = np.loadtxt(file, skiprows=1, delimiter=',')
    t = data[1:, 0]
    wl = data[0, 1:]
    dA = data[1:, 1:]  
    return t, wl, dA    

# load data 
if '.npy' in fname:
    t, l, dA = np.load(fname, allow_pickle=True)
if '.pdat' in fname:
    t, l, dA = load_pdat(fname)

# cut data 
t, l, dA = cut_data(t, l, dA)

# convert to 10^3 cm-1
if IR==False:
    l = (1/l)*10**4

# bounds for parametersc
b = get_bounds(p0, IRF)

# plot raw data 
fig, ax = plt.subplots(1,1)
plot_2D_map(ax, t, l, dA, scale)
fig.tight_layout()
plt.show()

# calculate Singular Values
SVD(dA)

# do the fit
res = least_squares(least_squares_fit, p0, bounds=b)
p = res.x

# plot the results
chi2 = plot_fit(dA, l, t, p, model)
plot_kinetics(dA, l, t, p, model)
label = print_errors(res, p, IRF)
plot_EAS(dA, l, t, p, model, label)
