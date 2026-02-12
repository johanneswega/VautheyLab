from helper import *
from VautheyLab.miscellaneous import *

# sample name
sample = 'NP'
ref = 'R6G'

# excitaiton wavelength
ex_wl = [470]
phi = []
phi_err = []

for i in range(len(ex_wl)):
    # initalize figure
    fig, ax = plt.subplots(2, 3, figsize=(15,6))

    # initialize empty Aex and integrated fluorescene
    Aex_S = []
    Iint_S = []
    Aex_R = []
    Iint_R = []

    # load absorption spectra
    files_S = [i for i in np.sort(os.listdir('abs/%s'%sample)) if not '.DS_Store' in i]
    files_R = [i for i in np.sort(os.listdir('abs/%s'%ref)) if not '.DS_Store' in i]
    colors = rainbow(files_S)

    # plot abs spectra for S
    for j in range(len(files_S)):
        wl, A = get_abs('abs/%s/'%sample + files_S[j], [300, 500])
        Aex_S.append(A[find_index(wl, ex_wl[i])])
        ax[0,0].plot(wl, A, '-', color=colors[j])

    # plot abs spectra for R
    for j in range(len(files_R)):
        wl, A = get_abs('abs/%s/'%ref +files_R[j], [300, 500])
        Aex_R.append(A[find_index(wl, ex_wl[i])])
        ax[1,0].plot(wl, A, '-', color=colors[j])    

    # load emission spectra
    files_S = [i for i in np.sort(os.listdir('Flu/%s'%sample)) if not '.DS_Store' in i]
    files_R = [i for i in np.sort(os.listdir('Flu/%s'%ref)) if not '.DS_Store' in i]

    # plot em spectra for S
    for j in range(len(files_S)):
        wl, I = get_em('Flu/%s/'%sample + files_S[j], [300, 800])
        Iint_S.append(np.trapz(I, wl))
        ax[0,1].plot(wl, I, '-', color=colors[j])

    # plot em spectra for R
    for j in range(len(files_R)):
        wl, I = get_em('Flu/%s/'%ref + files_R[j], [300, 800])
        Iint_R.append(np.trapz(I, wl))
        ax[1,1].plot(wl, I, '-', color=colors[j])
    print(Iint_R)
    print(Aex_R)
    # plot lines
    ax[0,2].plot(Aex_S, Iint_S, 'or')
    ax[1,2].plot(Aex_R, Iint_R, 'ob')

    # fit lines and calculate QY
    xfine = np.linspace(0, 0.13, 1000)
    p, pcov = curve_fit(fit, Aex_S, Iint_S)
    mS = p[0]
    err_mS = pcov[0,0]**0.5
    ax[0,2].plot(xfine, fit(xfine, *p), '--r')

    xfine = np.linspace(0, 0.13, 1000)
    p, pcov = curve_fit(fit, Aex_R, Iint_R, p0=[1e7])
    mR = p[0]
    err_mR = pcov[0,0]**0.5
    ax[1,2].plot(xfine, fit(xfine, *p), '--b')

    print(mR, mS)

    phi_R = 0.93
    err_phi_R = 0.02
    phi_S = phi_R * (mS/mR) 
    err_phi_S = phi_S * ( (err_mS/mS)**2 + (err_mR/mR)**2 + (err_phi_R/phi_R)**2)**(0.5)

    #for l in range(len(Aex_S)):
        #phi = (Iint_S[l] / Iint_R[l])*((1 - 10**(-Aex_R[l]))/(1 - 10**(-Aex_S[l])))
        #print(phi)

    phi.append(phi_S*100)
    phi_err.append( err_phi_S*100)
    print(phi_S, phi_R)

    # make labels
    ax[0,1].set_title('Sample = %s'%sample)
    ax[1,1].set_title('Refrence = %s'%ref)
    for k in [0,1]:
        ax[k,0].set_xlabel(r'$\lambda$ / nm')
        ax[k,0].set_xlabel(r'$\lambda$ / nm')
        ax[k,1].set_ylabel(r'counts / a.u.')
        ax[k,1].set_ylabel(r'counts / a.u.')
        ax[k,0].set_ylabel(r'$A$')
        ax[k,0].set_ylabel(r'$A$')
        ax[k,2].set_ylabel(r'$I_{int}$ / a.u.')
        ax[k,2].set_xlabel(r'$A_{ex}$')
    fig.tight_layout()
fig.savefig('QY.svg')
plt.show()

fig, ax = plt.subplots(1,1)
ax.errorbar(ex_wl, phi, yerr=phi_err, fmt='ob', capsize=3)
fig.tight_layout()
plt.show()
    
