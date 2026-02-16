from VautheyLab.plot_overview import Overview
from VautheyLab.plot_kinetics import Kinetics
from VautheyLab.compare_kinetics import Compare_Kinetics
from VautheyLab.compare_contours import Compare_Contours
from VautheyLab.compare_spectra import Compare_Spectra
from VautheyLab.compare_overviews import Compare_Overviews
from VautheyLab.global_analysis import Global_Analysis
from VautheyLab.plot_2DIR import twoDIR
from VautheyLab.plot_contour import Contour
from VautheyLab.movie_TA import Movie
from VautheyLab.miscellaneous import find_index
import numpy as np

def save(fname, t, wl, dA):
    np.save(fname, np.array([t, wl, dA], dtype=object))

def load(fname, t_cuts=None, wl_cuts=None):
    t, wl, dA = np.load(fname, allow_pickle=True)

    if wl_cuts!=None:
        dA = dA[:,(wl>wl_cuts[0])&(wl<wl_cuts[1])]
        wl = wl[(wl>wl_cuts[0])&(wl<wl_cuts[1])]

    if t_cuts!=None:
        dA = dA[(t>t_cuts[0])&(t<t_cuts[1]), :]
        t = t[(t>t_cuts[0])&(t<t_cuts[1])]
    
    return t, wl, dA

def interpolate_TA(file1, file2):
    # file for time axis
    t1, wl1, dA1 = load(file1)
    # file to be interpolated
    t2, wl2, dA2 = load(file2)
    # make empty dA frame
    dAint = np.zeros((len(t1), len(wl1)))   
    # go through wavelengths and interpolate
    for i in range(len(wl1)):
        # find index 
        index = find_index(wl2, wl1[i])
        dAint[:, i] = np.interp(t1, t2, dA2[:,index])
    # save interpolated file
    save('%s_int.npy'%(file2[:file2.find('.')]), t1, wl1, dAint)

# function to interpolate UVvis TA and TRIR on the same time grid
def interpolate_TA_TRIR(file1, file2):
    # load TRIR
    tIR, wnIR, dAIR = load_pdat(file1)
    # load vis
    tvis, wlvis, dAvis = load(file2)
    # for each wavelength we want to interplote the vis dA on the tIR axis
    dAint = np.zeros((len(tIR), len(wlvis)))
    for i in range(len(wlvis)):
        # get trace
        trace = dAvis[:, i]
        # interpolate
        int_trace = np.interp(tIR, tvis, trace)
        # add to array
        dAint[:, i] = int_trace
    save('%s_int.npy'%(file2[:file2.find('.')]), tIR, wlvis, dAint)

# function to interpolate two TRIR spectra on the same time grid
def interpolate_TRIR(file1, file2):
    # load TRIR
    tIR, wnIR, dAIR = load_pdat(file1)
    # load vis
    t, wn, dA = load_pdat(file2)
    # for each wavelength we want to interplote the vis dA on the tIR axis
    dAint = np.zeros((len(tIR), len(wn)))
    for i in range(len(wn)):
        # get trace
        trace = dA[:, i]
        # interpolate
        int_trace = np.interp(tIR, t, trace)
        # add to array
        dAint[:, i] = int_trace
    save_pdat('%s_int.pdat'%(file2[:file2.find('.')]), tIR, wn, dAint)

def save_pdat(name, t, wl, dA):
    pdat = np.zeros((len(t)+1, len(wl)+1))
    pdat[1:,0] = t
    pdat[0, 1:] = wl
    pdat[1:, 1:] = dA
    np.savetxt(name, pdat, header='ps*nm', delimiter=',')    

def load_pdat(file):
    data = np.loadtxt(file, skiprows=1, delimiter=',')
    t = data[1:, 0]
    wl = data[0, 1:]
    dA = data[1:, 1:]  
    return t, wl, dA    

def average_pdats(files):
    dAs = []
    for i in range(len(files)):
        t, wl, dA = load_pdat(files[i])
        dAs.append(dA)
    dAs = np.array(dAs)
    dAmean = np.mean(dAs, axis=0)
    save_pdat('mean.pdat', t, wl, dAmean)

def npy_to_pdat(file):
    t, wl, dA = load(file)
    pdat = np.zeros((len(t)+1, len(wl)+1))
    pdat[1:,0] = t
    pdat[0, 1:] = wl
    pdat[1:, 1:] = dA
    np.savetxt('%s.pdat'%(file[:file.find('.')]), pdat, header='ps*nm', delimiter=',')

def pdat_to_npy(file):
    t, wl, dA = load_pdat(file)  
    np.save('%s.npy'%(file[:file.find('.')]), np.array([t, wl, dA], dtype=object))

def txt_to_npy(t_file, wl_file, dA_file):
    t = np.loadtxt(t_file, skiprows=1)
    wl = np.loadtxt(wl_file, skiprows=1)
    dA = np.loadtxt(dA_file, skiprows=1)
    save('dA_from_txt.npy', t, wl, dA)

def convert_to_txt(fname, experiment='femto', t_cuts=None, wl_cuts=None):
    if '.npy' in fname:
        t, l, dA = load(fname, t_cuts=t_cuts, wl_cuts=wl_cuts)
        # Save array a to a text file
        np.savetxt('wavelength.txt', l, header='wavelength / nm')

    elif '.pdat' in fname:
        t, l, dA = load_pdat(fname)
        # Save array a to a text file
        np.savetxt('wavenumber.txt', l, header='wavenumber / cm-1')

    if experiment=='nano':
        headert = 'time / ns'
    else:
        headert = 'time / ps'

    # Save array b to a text file
    np.savetxt('time.txt', t, header = headert)

    # Save array c to a text file
    np.savetxt('TA.txt', dA, header='Transient Absorption / mOD')        

def export(fname, spectra=None, kinetics=None, experiment='femto', t_cuts=None, wl_cuts=None):
    # load data 
    t, wl, dA = load(fname, t_cuts=t_cuts, wl_cuts=wl_cuts)
    wn = (1/wl)*10**4

    # go through all spectra and export 
    if not spectra==None:
        for i in range(len(spectra)):
            if experiment=='femto':
                if spectra[i] > 1000:
                    label = '%s_spectrum_at_%3g_ns.txt'%(fname[:fname.find('.')], t[find_index(t, spectra[i])]/1000)
                elif np.abs(spectra[i]) < 1:
                    label = '%s_spectrum_at_%3g_fs.txt'%(fname[:fname.find('.')], t[find_index(t, spectra[i])]*1000)
                else:
                    label = '%s_spectrum_at_%3g_ps.txt'%(fname[:fname.find('.')], t[find_index(t, spectra[i])])
            else:
                if spectra[i] > 1000:
                    label = '%s_spectrum_at_%3g_us.txt'%(fname[:fname.find('.')], t[find_index(t, spectra[i])]/1000)
                elif np.abs(spectra[i]) < 1:
                    label = '%s_spectrum_at_%3g_ps.txt'%(fname[:fname.find('.')], t[find_index(t, spectra[i])]*1000)  
                else:
                    label = '%s_spectrum_at_%3g_ns.txt'%(fname[:fname.find('.')], t[find_index(t, spectra[i])])           
                
            np.savetxt(label, np.column_stack([wl, wn, dA[find_index(t, spectra[i]), :]]),
                        header='wavelength / nm,     wavenumber / kK,        TA / mOD', delimiter=',')
        
    # go through all kinetics and export
    if not kinetics==None:
        for i in range(len(kinetics)):
            label = '%s_kinetics_at_%3g_nm.txt'%(fname[:fname.find('.')], wl[find_index(wl, kinetics[i])])
            if experiment=='femto':
                header = 'time / ps,      TA / mOD'
            else:
                header = 'time / ns,      TA / mOD'
            np.savetxt(label, np.column_stack([t, dA[:, find_index(wl, kinetics[i])]]),
            header=header, delimiter=',')

                