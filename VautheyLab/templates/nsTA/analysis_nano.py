import numpy as np
from scipy.optimize import least_squares
import os
import shutil
from PIL import Image
from PyPDF2 import PdfMerger    
import matplotlib.pyplot as plt
np.seterr(divide='ignore')
plt.style.use('/Users/wega/Documents/matplotlib/test.mplstyle')

def pixel_to_lambda(Ho_TA_file, WL_file, Ho_ref_file):
    print('### Starting Analysis ###')
    Ho_TA = np.loadtxt(Ho_TA_file, skiprows=1)
    Ho_ref = np.loadtxt(Ho_ref_file, skiprows=1, delimiter=',', usecols=[0,1])
    WL = np.loadtxt(WL_file, skiprows=1)
    fig, ax = plt.subplots(ncols=1,nrows=1, figsize=(8,5))
    # calculate Absorbance 
    pixel = Ho_TA[:,0] 
    A_TA = -np.log10(np.abs(Ho_TA[:,1])/np.abs(WL[:,1]))
    # subtract baseline and normalize
    A_TA = (A_TA - A_TA[0])/(np.max(A_TA))
    # reference spectrum
    wl_ref = Ho_ref[:,0]
    A_ref = Ho_ref[:,1]
    A_ref = (A_ref[wl_ref>300] - A_ref[wl_ref>300][-1]) /np.max(A_ref[wl_ref>300])
    wl_ref = wl_ref[wl_ref>300]
    # least squares fit
    # limit pixel range for fit
    px_lim = [50, 250]
    res = least_squares(conversion, x0=[0.7, 360],
                        args=(pixel[(pixel<px_lim[1])&(pixel>px_lim[0])], 
                              A_TA[(pixel<px_lim[1])&(pixel>px_lim[0])], wl_ref, A_ref))
    p = res.x
    scale = p[0]
    shift = p[1]
    ax.plot(scale*pixel+shift, A_TA, '-b', linewidth=1.5, label='TA setup')
    ax.plot(wl_ref, A_ref, '-r', linewidth=1.5, label='Ref.')
    ax.set_xlabel(r'$\lambda$ / nm')
    ax.set_ylabel(r'norm. Abs.')
    ax.set_xticks(np.linspace(300,800,6))
    ax.set_yticks([])
    ax.legend()
    ax2 = ax.secondary_xaxis("top", functions=(lambda x: (x-shift)/scale, lambda x: x*scale + shift))
    ax2.set_xlabel('pixel')
    ax2.set_xticks(np.linspace(0,500,6))
    ax.set_title('scale = %.3g, shift = %.3g'%(scale,shift))
    fig.tight_layout()
    if os.path.exists('Results')==True:
        shutil.rmtree('Results')
    os.mkdir('Results')
    fig.savefig('Results/01_Holm.png')
    print('[#......] 1/7', end="\r")
    plt.show()
    return scale*pixel+shift
    
def conversion(p, pixel, A_TA, wl_ref, A_ref):
    wl = p[0]*pixel+p[1]
    A_int = np.interp(wl_ref, wl, A_TA)
    return A_ref - A_int

def get_individual_dA(TA_file):
    fh = open(TA_file)
    # initialize array for data
    data = []
    a = []
    b = []
    # list for nsamples in each delay bin
    nsamples = []
    for line in fh:
        # split line
        line = line.split()
        # if line not empty, i.e. not blank
        if len(line)!=0:
            # reject header
            if line[0]!='%':
                    # append time delay
                    a.append(float(line[0]))
                    # append TA signal
                    a.append(float(line[2]))
                    # append samples in bin
                    b.append(float(line[-1]))
        # if line blank --> next delay
        if len(line)==0:
            # add a-list to data 
            data.append(a)
            # add n samples to b list
            nsamples.append(b[0])
            # reset place holder lists
            a = []
            b = []
    fh.close()

    # the data array is now organized as:
    # [delay0, TA@px0, delay0, TA@px1, delay0, ...,delay0, TA@px524]
    # [delay1, TA@px0, delay1, TA@px1, delay1, ...,delay1, TA@px524]
    # and so on until the last delay

    # make empty list for delays
    t = []
    for i in range(len(data)):
        if len(data[i])!=0:
            t.append(data[i][0]*10**9)

    # number of pixel of CCD
    npx = 524
    # initialize empty dA-array
    dA = np.zeros((len(t), npx))
    # add values into dA-array
    for i in range(len(t)):
        spec = np.array(data[i])
        # only take the TA@px values 
        # i.e. reject every second entry or where data==delay_i
        spec = spec[spec!=data[i][0]]
        dA[i,:] = spec

    return np.array(t), dA, np.array(nsamples)

def plot_2D_map(t, l, dA, ax, lab):
    l = (1/l)*10**4
    nlevels = 51
    levels = np.linspace(scale[0], scale[1], nlevels)
    colors = plt.cm.RdBu_r(np.linspace(0,1, nlevels))
    #ax.contour(l, t, dA, levels=levels, linewidths=0.05, colors='k')
    D = ax.contourf(l, t, dA, levels=levels, colors=colors)
    cbar = plt.colorbar(D, ax=ax)
    cbar.set_label(lab)  
    cbar.set_ticks([scale[0], 0, scale[1]])
    ax.set_xlabel(r'$\tilde{\nu} / 10^{3} \,$ cm$^{-1}$')
    ax.set_ylabel(r'$\Delta t / \text{ns}$')
    ax.set_yscale('log')
    ax.invert_xaxis()
    ax2 = ax.secondary_xaxis("top", functions=(lambda x: (1/x)*10**+4, lambda x: (1/x)*10**-4))
    ax2.set_xlabel(r'$\lambda / \text{nm}$')
    ax.set_yscale('symlog')
    ax.set_ylim([-3, 5*10**4])
    ax.axhline(y=0, color='k')

def plot_background_spectra(ax, l, dA, ex_wl):
    l = (1/l)*10**4
    ax.plot(l, dA, '-')
    ax.set_title('Background')
    ax.set_ylabel(r'$\Delta{A} / 10^{-3}$')  
    ax.set_xlabel(r'$\tilde{\nu} / 10^{3} \,$ cm$^{-1}$')
    ax.invert_xaxis()
    ax.axvline(x=(1/ex_wl)*10**4, color='k')
    ax2 = ax.secondary_xaxis("top", functions=(lambda x: (1/x)*10**+4, lambda x: (1/x)*10**-4))
    ax2.set_xlabel(r'$\lambda / \text{nm}$')
    
def plot_overview(ax, t, l, dA, delay):
    l = (1/l)*10**4
    c = plt.cm.rainbow_r(np.linspace(0, 1, len(delay)))
    for i in range(len(delay)):
        index = min(range(len(t)), key=lambda j: abs(t[j]-delay[i]))
        if delay[i]>=1000:
            lab='%.3g µs'%(delay[i]/1000)
        else:
            lab='%.3g ns'%(delay[i])
        ax.plot(l, dA[index,:], '-', linewidth=1, color=c[i], label=lab)
    ax.invert_xaxis()
    ax.set_xlim([29, 13])
    ax.set_xticks(np.linspace(28,14,8))
    # figure out limit
    ax.set_ylim(scale)
    ax.axhline(y=0, color='k', linewidth=1)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
    ax.set_ylabel(r'$\Delta{A} / 10^{-3}$')  
    ax.set_xlabel(r'$\tilde{\nu} / 10^{3} \,$ cm$^{-1}$')
    ax2 = ax.secondary_xaxis("top", functions=(lambda x: (1/x)*10**+4, lambda x: (1/x)*10**-4))
    ax2.set_xlabel(r'$\lambda$ / nm')

def load_TA_data(ex_wl, TA_file, Ho_TA_file, WL_file, Ho_ref_file):
    # do px-to-lambda conversion and get wavelength axis
    wl = pixel_to_lambda(Ho_TA_file, WL_file, Ho_ref_file) 

    # get TA data
    t, dA, nsamples = get_individual_dA(TA_file)

    # plot counting statistics
    fig,ax = plt.subplots(ncols=1,nrows=1, figsize=(8,5))
    ax.plot(t, nsamples, '.r')
    ax.set_xlabel(r'$\Delta t$ / ns')
    ax.set_ylabel(r'samples')
    ax.set_yscale('log')
    ax.set_xscale('symlog')
    ax.set_title('Counting Statistics')
    fig.savefig("Results/02_counting.png")
    print('[##.....] 2/7', end="\r")
    fig.tight_layout()
    plt.show()

    # apply shift
    t = t + offset

    # plot on entire timescale
    fig,ax = plt.subplots(ncols=1,nrows=1, figsize=(8,5))
    plot_2D_map(t, wl, dA, ax, r'$\Delta{A}$ / mOD')
    ax.set_title(r'raw spectrum')
    fig.tight_layout()
    fig.savefig("Results/03_raw_spec.png")
    print('[###....] 3/7', end="\r")
    plt.show()

    # background 
    bg = np.mean(dA[t<-100,:], axis=0)
    fig,ax = plt.subplots(ncols=1,nrows=1, figsize=(8,5))
    plot_background_spectra(ax, wl, bg, ex_wl)
    fig.tight_layout()
    fig.savefig("Results/04_background.png")
    print('[####...] 4/7', end="\r")  
    plt.show()

    # subtract background
    for j in range(len(t)):
        dA[j,:] = dA[j,:] - bg

    # plot corrected spec
    fig, ax = plt.subplots(ncols=1,nrows=1, figsize=(8,5))
    plot_2D_map(t, wl, dA, ax, r'$\Delta{A}$ / mOD')
    ax.set_title(r'background corrected')
    fig.tight_layout()
    fig.savefig("Results/05_spec.png")
    print('[#####..] 5/7', end="\r")  
    plt.show()

    # plot early times
    fig, ax = plt.subplots(ncols=1,nrows=1, figsize=(8,5))
    plot_2D_map(t, wl, dA, ax, r'$\Delta{A}$ / mOD')
    ax.set_title(r'time zero set correctly?')
    ax.set_ylim([-1, 1])
    ax.axhline(y=0, color='k')
    ax.set_yscale('linear')
    fig.tight_layout()
    fig.savefig("Results/06_spec.png")
    print('[#####..] 6/7', end="\r")  
    plt.show()

    # plot overview spectrum
    fig,ax = plt.subplots(ncols=1, nrows=1, figsize=(8,5))
    delay = np.array([2, 5, 10, 20, 50, 80, 100, 200, 250, 500, 1e+3, 5e+3, 20e+3, 50e+3, 100e+3, 200e+3, 500e+3])
    plot_overview(ax, t, wl, dA, delay)
    fig.tight_layout()
    fig.savefig("Results/07_overview.png")
    print('[#######] 7/7', end="\r")  
    plt.show()
    
    # generate results file
    if os.path.exists('Results.pdf')==True:
            os.remove('Results.pdf')
    names = os.listdir('Results')
    names = np.sort(names)
    for i in range(len(names)):
        image = Image.open('Results/%s'%names[i])
        im = image.convert('RGB')
        im.save('Results/%s.pdf'%(names[i][:-4]))
    merger = PdfMerger()
    for i in range(len(names)):
        merger.append('Results/%s.pdf'%(names[i][:-4]))
    merger.write('Results.pdf') 
    merger.close()

    # save data
    np.save('dA.npy', np.array([t, wl, dA], dtype=object))
    # convert to txt
    # Assuming you have numpy arrays x, y, and z
    t, l, dA = np.load('dA.npy', allow_pickle=True)
    # Save array a to a text file
    np.savetxt('wavelength.txt', l, header='wavelength / nm')
    # Save array b to a text file
    np.savetxt('time.txt', t, header = 'time / ns')
    # Save array c to a text file
    np.savetxt('TA.txt', dA, header='Transient Absorption / mOD')

# define scale for 2D TA plots in mOD
scale = [-70,70] 
# define excitation wavelength in nm
ex_wl = 355
# offset for kinetics
offset = -1.3

# get data files
files = [f for f in os.listdir() if f.endswith('.dat')]
for i in range(len(files)):
    if 'SAMPLE' in files[i]:
        sample = files[i]
    if 'WL' in files[i]:
        WL = files[i]
    if 'HOLMIUM' in files[i]:
        HO = files[i]

load_TA_data(ex_wl,
             sample,
             HO,
             WL,
             'HOLM_ref.csv')

